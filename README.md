# Hierarchical Bloom filter Index (HBI)

A compact sliding-window index for **exact pattern matching over data streams**, where queries are ad hoc: patterns are not known in advance and can arrive at any time.

HBI keeps the most recent `w` symbols of a stream and answers "where does pattern `P` occur in the current window?" exactly. Every occurrence is reported, and nothing else. It is the implementation accompanying the master's thesis [*Pattern Matching on Big Data Streams*](https://research.tue.nl/nl/studentTheses/pattern-matching-on-big-data-streams/) (Eindhoven University of Technology, 2025; [PDF](https://research.tue.nl/files/415744451/Papadimitriou_D.pdf)).

## How it works

- **Hierarchy of Bloom filters.** The window is split into a binary hierarchy of intervals. Level `L` splits each tree into `2^L` intervals, and each level has one Bloom filter holding keys `<L, interval, symbol>`, meaning "this symbol occurs somewhere in this interval". An arriving symbol adds one key per level, so the update cost is `O(log w)`.
- **Coarse-to-fine search with exact verification.** A query walks the hierarchy depth-first and tests the pattern's symbols against each interval's filter. A negative answer prunes the whole interval. Surviving intervals are checked symbol by symbol against the buffered window, so Bloom filter false positives never reach the answer.
- **Sliding window by whole trees.** The window is covered by a few implicit trees. A tree is dropped in `O(1)` as soon as none of its positions is inside the window, so no per-item deletion is needed.
- **Search strategies.** All strategies return the same answers and differ only in how much work they do:
  - **Block Search (BS):** starts at the root.
  - **Confidence-Primed BS:** starts at the level where the rarest pattern symbol stops being almost certainly present.
  - **Confidence-Gated BS:** tests each symbol only from its own level downward.
  - **Cost-Function BS:** picks the start level that minimizes the expected number of filter probes.

  The last three use a small per-tree Count Sketch (width 2048, depth 8) to estimate symbol frequencies.
- **n-grams.** Indexing overlapping n-grams instead of single characters makes each probe far more selective on low-entropy data such as natural language. In exchange, patterns must be at least `n` characters long.
- **Provisioning policies.** The reactive and predictive policies drop the top levels that offer little pruning. They use a constant-time sampler and a low-percentile frequency estimate.

HBI uses far less memory than suffix-tree based indexes such as SSWSI or an online suffix tree, and inserts take well under a microsecond per symbol. Query speed depends on symbol selectivity. On very low-entropy streams, such as a repeated de Bruijn sequence, pruning disappears and suffix-based indexes answer faster.

## Using the index

`PMIndex.HierarchicalBloomFilterIndex` is the public entry point:

```java
import PMIndex.HierarchicalBloomFilterIndex;
import PMIndex.HierarchicalBloomFilterIndex.SearchStrategy;

HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
        .windowLength(1 << 20)                  // search the last 1,048,576 characters
        .nGram(4)                               // index 4-grams; patterns need >= 4 characters
        .searchStrategy(SearchStrategy.CONFIDENCE_PRIMED)
        .build();

index.append(chunkOfText);                      // or index.append(char) per character
long[] hits = index.search("needle");           // absolute offsets, ascending
boolean seen = index.contains("haystack");
long oldest = index.windowStart();              // offset of the first character still in the window
```

Builder options. Only `windowLength` is required.

| Option | Default | Meaning |
| --- | --- | --- |
| `windowLength` | required | Number of most recent characters that queries search. |
| `nGram` | 1 | Characters per indexed symbol. Larger values prune better on text. |
| `treeLength` | largest power of two ≤ window | Symbols per tree (power of two). Smaller trees mean finer expiration. |
| `expectedDistinctSymbols` | `treeLength` | Upper bound on distinct n-grams per tree. It only sizes the Bloom filters, and underestimating costs speed, not correctness. |
| `fpRate` | 0.2 | Target false positive rate of each Bloom filter level. |
| `searchStrategy` | `BLOCK_SEARCH` | `BLOCK_SEARCH`, `CONFIDENCE_PRIMED`, `CONFIDENCE_GATED`, or `COST_FUNCTION`. |
| `strides` | false | Probe only every n-th n-gram of a pattern. Faster for large `n`. |
| `provisioningPolicy` | `NONE` | `NONE`, `REACTIVE`, or `PREDICTIVE`. |
| `confidence` | 0.99 | Confidence threshold used by the confidence-based strategies and the policies. |

Results are exact. N-grams of up to four characters are stored losslessly. Longer n-grams are hashed, and every reported match is confirmed against the buffered characters. Instances are not thread-safe.

## Building and testing

Requirements: Java 21+ and Maven.

```bash
mvn test          # differential tests against a brute-force scan, plus regression tests
mvn package       # shaded jar: target/Hierarchical-Bloom-filter-Index-1.0-SNAPSHOT.jar
```

The tests in `src/test/java` run every strategy, n-gram size, stride mode and provisioning policy on streams longer than the window, and compare each answer with a brute-force scan.

## Reproducing the experiments

The shaded jar runs `HBIDatasetBenchmarkMulti`. It loads each dataset window under `data/<window>`, runs the query workloads under `queries/<window>`, and compares HBI with the SSWSI and online suffix tree baselines:

```bash
java -jar target/Hierarchical-Bloom-filter-Index-1.0-SNAPSHOT.jar \
     --window w21 --ngrams 1,2,4,8 --fp 0.2 --algorithm bs --type all --runs 5 --warmup 1
```

Useful flags:
- `--algorithm bs|cp|cgbs|cf`
- `--policy NONE|REACTIVE|PREDICTIVE`
- `--strides true|false`
- `--fp-list 0.1,0.2,0.3`
- `--tree-power 20`
- `--mode segments`, for one token per line, as in CAIDA-style records
- `--run-suffix false`, to skip the baselines
- `--verify true`, to check that every index reports the same occurrences for every query (all of them are exact, so any difference is a bug)
- `--sswsi-mode amortized|legacy`, `--sswsi-builder dc3|fc` and `--sswsi-dict fks|fastutil`, which SSWSI runs (also in `MemoryBenchmark`). The default, `amortized` with `dc3` and `fks`, is the amortized structure of Bille et al. (CPM 2023, Sec. 3.3) with the paper's FKS dictionaries; it needs a power-of-two window. `fc` is the paper's Farach-Colton construction: it builds the same trees, 1.1–1.4 times slower. `legacy` is the thesis-era structure, with fixes for deep trees, n-gram remapping and position overflow, and takes no builder or dictionary. The setting is in the CSV name and in its `sswsi_*` columns
- `--out-dir <dir>`, where the CSV goes (by default the working directory)

The window length comes from `--window-length`, then `--window-power`, then the number in the window name (`w19` gives 2^19). Insert and query times are measured in nanoseconds, after the window and the query files have been read into memory. `avgQueryMs` (`avg_query_ms_*` in the CSV) is the time for a whole query file of 500 queries, not for one query: queries per second are `500 * 1000 / avgQueryMs`.

### Query workloads

Every benchmark query set (`queries/w19|w20|w21`, `queries/wzipf_*`, `queries/adversarial/*` and `queries/mawi*`) is generated by `datagenerators.GenerateQueries`, 500 queries per file:
- `uniform`: substrings at random positions of the window;
- `missing`: the same substrings, half of them changed in 30% of their symbols, drawn from the window's own symbol frequencies, until they no longer occur, so they are realistic near misses;
- `rare`: substrings that contain one of the window's rarest symbols.

Window `i` of a family with windows of 2^p symbols uses seed `100 * p + i`, for example:

```bash
java -Xmx6g -cp "$CP" datagenerators.GenerateQueries data/w21/1/1_W21.txt queries/w21/1 chars 2101 10,20,40,80,160,320,450,640
```

The thesis used earlier query files, which are in the git history. The `unique_substrings_*` files and `queries/pg2701` are separate sets for the start-level estimator experiments.

### Network trace datasets (MAWI)

The packet-trace windows are built from a MAWI trace (`data/traces/202501011400.pcap.gz`, the 15-minute sample of 2025-01-01, 136.9M packets) and are not in git. Each IPv4/IPv6 packet becomes one record `src,dst,proto`, the 3-tuple symbol the thesis uses for CAIDA. Regenerate them with:

```bash
java -cp "$CP" datagenerators.PcapTokenizer data/traces/202501011400.pcap.gz data/traces/202501011400.tokens.txt
java -cp "$CP" datagenerators.SampleWindows data/traces/202501011400.tokens.txt data mawi 1 5 19 20 21
for p in 19 20 21; do for i in 1 2 3 4 5; do
  java -Xmx4g -cp "$CP" datagenerators.GenerateQueries data/mawi$p/$i/${i}_W$p.txt queries/mawi$p/$i segments $((p*100+i))
done; done
```

This gives `data/mawi19|mawi20|mawi21/{1..5}` (windows of 2^19, 2^20 and 2^21 records at random, non-overlapping positions; `SOURCE.txt` lists them) and the matching `queries/`. Run them with `--window mawi21 --mode segments`. A 2^21 window holds 670K–770K distinct records, so pass `--alphabet 800000` to size the Bloom filters. On this data `rare` behaves almost like `uniform`, because a third of all records occur only once.

Other entry points:

| Class | Purpose |
| --- | --- |
| `MemoryBenchmark` | Memory footprint, measured with JOL. |
| `TrackIndexingSpeed` | Insertion speed. |
| `ConfidenceExperiment` | Accuracy of the start-level estimators. |
| `SamplingSketchExperiment` | The percentile sampler. |
| `GenerateDatasets`, `GenerateAdversarialDatasets` | Synthetic and adversarial streams. |

## Repository layout

| Path | Contents |
| --- | --- |
| `src/main/java/PMIndex` | The index (`HBI`, `HierarchicalBloomFilterIndex`, configuration, stats) and the suffix-based baselines. |
| `src/main/java/tree` | Implicit tree layout, per-level filter directory, raw symbol buffer, and the suffix tree and array builders used by the baselines. |
| `src/main/java/search` | Block Search and its variants, candidate verification, and start-level planning. |
| `src/main/java/membership` | Bloom filter with double hashing. |
| `src/main/java/estimators` | Count Sketch, the percentile sampler (HOPS), and the start-level cost functions. |
| `src/main/java/utilities` | Benchmark orchestration, dataset readers, reporting, and hashing helpers. |
| `data/`, `queries/` | Dataset windows (Gutenberg text, Zipf, uniform, adversarial) and query workloads per pattern length. |
| `perf_data/`, `scripts/`, `*.sh` | Hardware counter measurements and the scripts that summarize them. |

## License

MIT, see [LICENSE](LICENSE).
