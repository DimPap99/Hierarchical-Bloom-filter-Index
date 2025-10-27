import PMIndex.*;
import search.*;
import estimators.*;
import membership.BloomFilter;
import membership.Membership;
import utilities.*;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.function.Supplier;

/**
 * HBIDatasetBenchmark
 *
 * This benchmark compares HBI versus SuffixTreeIndex on the same dataset and query set.
 *
 * We now support two streaming modes, exactly like in ConfidenceExperiment:
 *
 *   mode = "chars"
 *     Treat the dataset file as one continuous character stream.
 *     Build n-grams of characters. Queries are substrings, exactly like Experiment.run().
 *
 *   mode = "segments"
 *     Treat the dataset file as one token per line (for example one packet per line).
 *     Each token is a logical "symbol". We slide a RingBuffer<String> of length NGRAMS
 *     over those symbols to build patterns. We index those patterns in HBI or SuffixTreeIndex.
 *     Queries are sequences of tokens separated by spaces.
 *
 * You can pick the mode via --mode chars or --mode segments at runtime.
 *
 * Everything else in this benchmark (timing, result comparison, averages) is preserved.
 */
public class HBIDatasetBenchmark {

    /** Default input paths and parameters. Change these as you like. */
    private static String DATA_FILE =
            "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/zipf_21_1.txt";

    private static String QUERY_FILE =
            "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/zipf21_1/unique_substrings_zipf21_1_10.txt";

    private static final int WINDOW_LEN   = 1 << 21;
    private static final int TREE_LEN     = 1 << 21;
    private static int ALPHABET           = 74;
    private static final double FP_RATE   = 0.005;
    private static final int RUNS         = 1;
    private static final boolean USE_STRIDES = true;
    private static int NGRAMS             = 4;

    /**
     * Utility method from your original code.
     * Compare match results across two indexes for each query.
     */
    public static void compared(ArrayList<ArrayList<Integer>> arr1,
                                ArrayList<ArrayList<Integer>> arr2) {

        Objects.requireNonNull(arr1, "First index results must not be null.");
        Objects.requireNonNull(arr2, "Second index results must not be null.");

        boolean resultsMatch = arr1.size() == arr2.size();
        if (!resultsMatch) {
            System.out.printf(
                    "Query count mismatch: first=%d second=%d%n",
                    arr1.size(),
                    arr2.size());
        }

        int comparisons = Math.min(arr1.size(), arr2.size());
        for (int i = 0; i < comparisons; i++) {
            List<Integer> row1 = arr1.get(i);
            List<Integer> row2 = arr2.get(i);

            List<Integer> a = (row1 == null) ? new ArrayList<>() : new ArrayList<>(row1);
            List<Integer> b = (row2 == null) ? new ArrayList<>() : new ArrayList<>(row2);
            Collections.sort(a);
            Collections.sort(b);

            if (!a.equals(b)) {
                resultsMatch = false;
                ArrayList<Integer> onlyInFirst = new ArrayList<>();
                ArrayList<Integer> onlyInSecond = new ArrayList<>();
                int i1 = 0, i2 = 0;
                while (i1 < a.size() || i2 < b.size()) {
                    if (i1 >= a.size()) { onlyInSecond.add(b.get(i2++)); continue; }
                    if (i2 >= b.size()) { onlyInFirst.add(a.get(i1++)); continue; }
                    int va = a.get(i1), vb = b.get(i2);
                    if (va == vb) { i1++; i2++; }
                    else if (va < vb) { onlyInFirst.add(va); i1++; }
                    else { onlyInSecond.add(vb); i2++; }
                }
                System.out.printf(
                        "Mismatch at query %d: onlyInFirst=%s onlyInSecond=%s%n",
                        i,
                        onlyInFirst,
                        onlyInSecond);
            }
        }

        if (resultsMatch) {
            System.out.println("Indexes returned identical match results for all compared queries.");
        } else {
            System.out.println("Indexes produced differing results.");
        }
    }

    /**
     * Helper for "segments" mode only.
     * We read queries as lines and trim blanks. Same logic as ConfidenceExperiment.loadQueries.
     */
    private static List<String> loadQueries(String queriesFile) throws IOException {
        List<String> lines = Files.readAllLines(Path.of(queriesFile), StandardCharsets.UTF_8);
        ArrayList<String> out = new ArrayList<>(lines.size());
        for (String l : lines) {
            if (l == null) continue;
            String s = l.strip();
            if (!s.isEmpty()) out.add(s);
        }
        return out;
    }

    /**
     * Stream and index a dataset in "segments" mode, using the given IPMIndexing implementation.
     *
     * This is analogous to Experiment.run(...) for character mode,
     * and analogous to runStreamingSegments(...) from ConfidenceExperiment,
     * but returns ExperimentRunResult with matchRes filled so we can diff
     * HBI vs SuffixTreeIndex.
     *
     * Semantics:
     *   - Treat each nonempty line from DATA_FILE as one token (String).
     *   - Maintain a RingBuffer<String> of length ngram.
     *   - For every filled buffer snapshot, call index.insert(snapshotString).
     *   - Every WINDOW_LEN tokens, run all queries.
     *
     * Query semantics in segments mode:
     *   - Each line in QUERY_FILE is "tok0 tok1 tok2 ...".
     *   - We split by spaces, build new Pattern(splitArray, ngram),
     *     and call index.report(pat).
     *
     * We measure:
     *   totalInsertTimeMs
     *   totalRunTimeMs (query time sum)
     *   avgQueryLen in n-grams
     * We also build matchRes, parallel to what Experiment.run(...) returns.
     */
    private static ExperimentRunResult runSegmentsMode(
            String datasetPath,
            String queriesPath,
            IPMIndexing index,
            int ngram,
            int windowLen,
            boolean verbose
    ) throws IOException {

        // load queries up front
        List<String> queries = loadQueries(queriesPath);

        // compute average query length in n-grams
        int sumQueryLen = 0;
        for (String q : queries) {
            if (q == null || q.isEmpty()) continue;
            String[] splitArray = q.split(" ");
            Pattern p = new Pattern(splitArray, ngram);
            sumQueryLen += p.nGramToLong.length;
        }
        double avgQueryLen = queries.isEmpty() ? 0.0 : (sumQueryLen * 1.0 / queries.size());

        SegmentReader seg = new SegmentReader(
                datasetPath,
                queriesPath,
                '\n',
                false,  // do not include delimiter in token
                true    // skip empty segments
        );
        SegmentReaderAdapter adapter = new SegmentReaderAdapter(seg);

        RingBuffer<String> window = new RingBuffer<>(ngram);

        long totalInsertNanos   = 0L;
        long insertionEvents    = 0L;
        long tokensSeen         = 0L;

        long nextQueryAt        = Math.max(1, windowLen);
        long sumQueryLoadMs     = 0L;
        long queryLoads         = 0L;

        ArrayList<ArrayList<Integer>> allQueryMatches = new ArrayList<>();

        while (adapter.hasNext()) {
            String tok = adapter.next();
            window.append(tok);
            if (!window.isFilled()) {
                continue;
            }

            String s = window.snapshot().toString();

            long t0 = System.nanoTime();
            index.insert(s);
            long t1 = System.nanoTime();

            totalInsertNanos += (t1 - t0);
            insertionEvents++;
            tokensSeen++;

            if (tokensSeen >= nextQueryAt - ngram + 1) {
                long qMs = runAllQueriesSegments(
                        index,
                        queries,
                        ngram,
                        verbose,
                        allQueryMatches
                );
                sumQueryLoadMs += qMs;
                queryLoads++;
                nextQueryAt += windowLen;
            }
        }

        double totalInsertMs = totalInsertNanos / 1_000_000.0;
        double totalQueryMs  = sumQueryLoadMs;
        double avgQueryLoadMs =
                (queryLoads == 0) ? 0.0 : (totalQueryMs / queryLoads);

        return new ExperimentRunResult(
                totalQueryMs,
                totalInsertMs,
                /* patternResults = */ new ArrayList<>(),
                avgQueryLen,
                /* avgInsertMsPerSymbol = */ (insertionEvents == 0)
                ? 0.0
                : (totalInsertMs / insertionEvents),
                avgQueryLoadMs,
                allQueryMatches
        );
    }

    /**
     * Helper for runSegmentsMode.
     * Execute ALL queries (segments mode interpretation) against index, collect matches,
     * and measure wall clock time in milliseconds.
     *
     * allQueryMatches grows one entry per query. Each entry is the match result
     * list returned by index.report(pat).
     */
    private static long runAllQueriesSegments(
            IPMIndexing index,
            List<String> queries,
            int ngram,
            boolean verbose,
            ArrayList<ArrayList<Integer>> allQueryMatches
    ) {
        long start = System.currentTimeMillis();

        for (String q : queries) {
            if (q == null || q.isEmpty()) continue;

            String[] splitArray = q.split(" ");
            Pattern pat = new Pattern(splitArray, ngram);

            ArrayList<Integer> report = index.report(pat);
            allQueryMatches.add(report);

            if (verbose) {
                if (report.size() < 20 && q.length() < 200) {
                    System.out.println(q + ":" + report);
                }
            }
        }

        return System.currentTimeMillis() - start;
    }

    /**
     * Build a fresh HBI with the configured parameters.
     * This is unchanged structurally from your code, except for:
     *  - we keep NGRAMS as the last constructor parameter so it respects your chosen n-gram size
     */
    private static HBI newHbi(double conf) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(TREE_LEN);
        Supplier<Membership> memFactory = () -> new BloomFilter();
        Supplier<PruningPlan> prFactory = () -> new MostFreqPruning(conf);
        Verifier v = new VerifierLinearLeafProbe();

        return new HBI(
                new BlockSearch(),
                WINDOW_LEN,
                FP_RATE,
                ALPHABET,
                TREE_LEN,
                estFactory,
                memFactory,
                prFactory,
                v,
                /* cost function */ null,
                conf,
                NGRAMS
        );
    }

    public static void main(String[] args) throws IOException {

        // default mode is "chars" to preserve original behavior
        String mode = "chars";

        // light weight argument parsing for overrides
        for (int i = 0; i < args.length - 1; i += 2) {
            switch (args[i]) {
                case "--mode" -> mode = args[i + 1]; // "chars" or "segments"
                case "--data" -> DATA_FILE = args[i + 1];
                case "--queries" -> QUERY_FILE = args[i + 1];
            }
        }

        System.out.println("Benchmark starting…");
        System.out.println("Mode: " + mode);

        System.out.println("N-gram: " + NGRAMS);
        System.out.println("Window Size: " + WINDOW_LEN);
        System.out.println("Tree Length: " + TREE_LEN);

        // expand alphabet based on NGRAMS exactly as before
        ALPHABET = (int) Math.pow(ALPHABET, NGRAMS);
        ALPHABET = Math.min(ALPHABET, TREE_LEN);

        System.out.println("Alphabet: " + ALPHABET);
        System.out.println();

        double hbiTotalMs = 0.0;
        double hbiTotalMsInsert = 0.0;
        double suffixTotalMs = 0.0;
        double suffixTotalMsInsert = 0.0;

        double lpShareSum = 0.0;
        double avgQueryTimeSum = 0.0;
        double avgLpTimeSum = 0.0;
        int statsSamples = 0;

        // ------------------
        // Warm-up iteration
        // ------------------
        {
            HBI hbi = newHbi(0.999);
            hbi.strides = USE_STRIDES;
            hbi.stats().setCollecting(false);
            hbi.stats().setExperimentMode(false);

            ExperimentRunResult warmHbi;
            if ("segments".equalsIgnoreCase(mode)) {
                warmHbi = runSegmentsMode(
                        DATA_FILE, QUERY_FILE,
                        hbi,
                        NGRAMS,
                        WINDOW_LEN,
                        /* verbose */ false
                );
            } else {
                warmHbi = Experiment.run(
                        DATA_FILE, QUERY_FILE,
                        hbi,
                        NGRAMS,
                        /* verbose */ false,
                        /* queryResults */ false
                );
            }
            ArrayList<ArrayList<Integer>> warmHbiMatches = warmHbi.matchRes();

            IPMIndexing suffix = new StreamingSlidingWindowIndex(WINDOW_LEN);//SuffixTreeIndex(ALPHABET, 0.0001, WINDOW_LEN);

            ExperimentRunResult warmSuffix;
            if ("segments".equalsIgnoreCase(mode)) {
                warmSuffix = runSegmentsMode(
                        DATA_FILE, QUERY_FILE,
                        suffix,
                        /* force ngram=1? or use NGRAMS */ NGRAMS,
                        WINDOW_LEN,
                        false
                );
            } else {
                warmSuffix = Experiment.run(
                        DATA_FILE, QUERY_FILE,
                        suffix,
                        NGRAMS,
                        false,
                        false
                );
            }
            ArrayList<ArrayList<Integer>> warmSuffixMatches = warmSuffix.matchRes();

            compared(warmHbiMatches, warmSuffixMatches);
            int b = 2;
        }

        // -----------
        // Timed runs
        // -----------
        for (int i = 0; i < RUNS; i++) {

            // HBI pass
            HBI hbi = newHbi(0.99);
            hbi.strides = USE_STRIDES;
            HbiStats stats = hbi.stats();
            stats.setCollecting(true);
            stats.setExperimentMode(false);

            ExperimentRunResult hbiRes;
            if ("segments".equalsIgnoreCase(mode)) {
                hbiRes = runSegmentsMode(
                        DATA_FILE, QUERY_FILE,
                        hbi,
                        NGRAMS,
                        WINDOW_LEN,
                        false
                );
            } else {
                hbiRes = Experiment.run(
                        DATA_FILE, QUERY_FILE,
                        hbi,
                        NGRAMS,
                        false,
                        false
                );
            }

            hbiTotalMs       += hbiRes.totalRunTimeMs();
            hbiTotalMsInsert += hbiRes.totalInsertTimeMs();

            if (stats.totalQueryCount() > 0) {
                lpShareSum      += stats.lpShareOfQuery();
                avgQueryTimeSum += stats.averageQueryTimeMillis();
                avgLpTimeSum    += stats.averageLpTimeMillis();
                statsSamples++;
            }

            ArrayList<ArrayList<Integer>> hbiMatches = hbiRes.matchRes();

            // Suffix tree pass
            IPMIndexing suffix = new StreamingSlidingWindowIndex(WINDOW_LEN); //new SuffixTreeIndex(ALPHABET, 0.0001, WINDOW_LEN);

            ExperimentRunResult suffixRes;
            if ("segments".equalsIgnoreCase(mode)) {
                suffixRes = runSegmentsMode(
                        DATA_FILE, QUERY_FILE,
                        suffix,
                        NGRAMS,
                        WINDOW_LEN,
                        false
                );
            } else {
                suffixRes = Experiment.run(
                        DATA_FILE, QUERY_FILE,
                        suffix,
                        NGRAMS,
                        false,
                        false
                );
            }

            suffixTotalMs       += suffixRes.totalRunTimeMs();
            suffixTotalMsInsert += suffixRes.totalInsertTimeMs();

            ArrayList<ArrayList<Integer>> suffixMatches = suffixRes.matchRes();

            // direct correctness comparison for this run
            compared(hbiMatches, suffixMatches);
        }

        // --------------------------
        // Print aggregate statistics
        // --------------------------
        if (RUNS > 0) {
            System.out.printf(
                    Locale.ROOT,
                    "HBI avg (ms): %.3f%n",
                    hbiTotalMs / RUNS);
            System.out.printf(
                    Locale.ROOT,
                    "HBI Insert avg (ms): %.3f%n",
                    hbiTotalMsInsert / RUNS);
            System.out.printf(
                    Locale.ROOT,
                    "HBI Insert avg per symbol (ms): %.4f%n",
                    (hbiTotalMsInsert / RUNS) / WINDOW_LEN
            );

            if (statsSamples > 0) {
                System.out.printf(
                        Locale.ROOT,
                        "HBI avg query time per pattern (ms): %.3f%n",
                        (avgQueryTimeSum / statsSamples));
                System.out.printf(
                        Locale.ROOT,
                        "HBI avg LP computation time per pattern (ms): %.3f%n",
                        (avgLpTimeSum / statsSamples));
                System.out.printf(
                        Locale.ROOT,
                        "HBI LP time share of query (%%): %.2f%n",
                        (lpShareSum / statsSamples) * 100.0
                );
            }

            System.out.printf(
                    Locale.ROOT,
                    "SuffixTreeIndex avg (ms): %.3f%n",
                    suffixTotalMs / RUNS);
            System.out.printf(
                    Locale.ROOT,
                    "SuffixTreeIndex Insert avg (ms): %.3f%n",
                    suffixTotalMsInsert / RUNS);
            System.out.printf(
                    Locale.ROOT,
                    "SuffixTreeIndex Insert avg (ms) per char: %.3f%n",
                    (suffixTotalMsInsert / RUNS)/WINDOW_LEN);        }
    }
}
