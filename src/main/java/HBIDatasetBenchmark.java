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
    private static final String DEFAULT_DATA_FILE =
            "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/w20/1/1_W20.txt";

    private static final String DEFAULT_QUERY_FILE =
            "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/w20/1/40.uniform.txt";

    private static final int WINDOW_LEN       = 1 << 20;
    private static final int TREE_LEN         = 1 << 20;
    private static final int ALPHABET_BASE    = 500;
    private static final double DEFAULT_FP_RATE = 0.25;
    private static final int DEFAULT_RUNS     = 20;
    private static final double Confidence = 0.99;
    private static final boolean USE_STRIDES  = true;
    private static int NGRAMS                 = 8;

    private record BenchmarkOptions(String mode,
                                     String dataFile,
                                     String queryFile,
                                     double fpRate,
                                     boolean runSuffix,
                                     boolean skipQueries,
                                     int warmupRuns,
                                     int runs) {

        static BenchmarkOptions parse(String[] args) {
            String mode = "chars";
            String dataFile = DEFAULT_DATA_FILE;
            String queryFile = DEFAULT_QUERY_FILE;
            double fpRate = DEFAULT_FP_RATE;
            boolean runSuffix = false;
            boolean skipQueries = false;
            int warmupRuns = 2;
            int runs = DEFAULT_RUNS;

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (!arg.startsWith("--")) {
                    continue;
                }

                String key;
                String value;
                int eq = arg.indexOf('=');
                if (eq >= 0) {
                    key = arg.substring(2, eq);
                    value = arg.substring(eq + 1);
                } else {
                    key = arg.substring(2);
                    if (i + 1 >= args.length) {
                        throw new IllegalArgumentException("Missing value for option --" + key);
                    }
                    value = args[++i];
                }

                switch (key) {
                    case "mode" -> mode = value;
                    case "data" -> dataFile = value;
                    case "queries" -> queryFile = value;
                    case "fp" -> fpRate = Double.parseDouble(value);
                    case "run-suffix" -> runSuffix = Boolean.parseBoolean(value);
                    case "skip-queries" -> skipQueries = Boolean.parseBoolean(value);
                    case "warmup" -> warmupRuns = Integer.parseInt(value);
                    case "runs" -> runs = Integer.parseInt(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            if (runs < 0) {
                throw new IllegalArgumentException("--runs must be non-negative");
            }
            if (warmupRuns < 0) {
                throw new IllegalArgumentException("--warmup must be non-negative");
            }
            if (fpRate <= 0.0 || fpRate >= 1.0) {
                throw new IllegalArgumentException("--fp must be between 0 and 1 (exclusive)");
            }

            return new BenchmarkOptions(mode, dataFile, queryFile, fpRate, runSuffix, skipQueries, warmupRuns, runs);
        }
    }

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


    private static ExperimentRunResult runSegmentsMode(
            String datasetPath,
            String queriesPath,
            IPMIndexing index,
            int ngram,
            int windowLen,
            boolean verbose,
            boolean runQueries
    ) throws IOException {

        List<String> queries;
        if (runQueries) {
            queries = loadQueries(queriesPath);
        } else {
            queries = Collections.emptyList();
        }

        int sumQueryLen = 0;
        if (!queries.isEmpty()) {
            for (String q : queries) {
                if (q == null || q.isEmpty()) continue;
                String[] splitArray = q.split(" ");
                Pattern p = new Pattern(splitArray, ngram);
                sumQueryLen += p.nGramToLong.length;
            }
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
        boolean shouldRunQueries = runQueries && !queries.isEmpty();

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

            if (shouldRunQueries && tokensSeen >= nextQueryAt - ngram + 1) {
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
    private static HBI newHbi(double conf, int alphabet, double fpRate) {
        //ε=0.05, δ=7.5e-4 → w=2048, d=8 → ~64 K
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(TREE_LEN);
        Supplier<Membership> memFactory = () -> new BloomFilter();
        Supplier<PruningPlan> prFactory = () -> new MostFreqPruning(conf, fpRate);
        Verifier v = new VerifierLinearLeafProbe();

        return new HBI(
                new BlockSearch(),
                WINDOW_LEN,
                fpRate,
                alphabet,
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
        BenchmarkOptions options = BenchmarkOptions.parse(args);

        String mode = options.mode();
        String dataFile = options.dataFile();
        String queryFile = options.queryFile();
        double fpRate = options.fpRate();
        boolean runSuffix = options.runSuffix();
        boolean skipQueries = options.skipQueries();
        boolean runQueries = !skipQueries;
        int warmupRuns = options.warmupRuns();
        int runs = options.runs();

        System.out.println("Benchmark starting…");
        System.out.println("Mode: " + mode);

        System.out.println("N-gram: " + NGRAMS);
        System.out.println("Window Size: " + WINDOW_LEN);
        System.out.println("Tree Length: " + TREE_LEN);
        System.out.printf(Locale.ROOT, "False Positive Rate: %.4f%n", fpRate);
        System.out.println("Warmup runs: " + warmupRuns);
        System.out.println("Timed runs: " + runs);
        System.out.println("Suffix enabled: " + runSuffix);
        System.out.println("Run queries: " + runQueries);

        // expand alphabet based on NGRAMS exactly as before
        int alphabet = (int) Math.pow(ALPHABET_BASE, NGRAMS);
        alphabet = Math.min(alphabet, TREE_LEN);

        System.out.println("Alphabet: " + alphabet);
        System.out.println();

        // ------------------
        // Warm-up iterations
        // ------------------
        for (int i = 0; i < warmupRuns; i++) {
            HBI hbi = newHbi(Confidence, alphabet, fpRate);
            hbi.strides = USE_STRIDES;
            HbiStats stats = hbi.stats();
            stats.setCollecting(false);
            stats.setExperimentMode(false);

            if ("segments".equalsIgnoreCase(mode)) {
                runSegmentsMode(
                        dataFile,
                        queryFile,
                        hbi,
                        NGRAMS,
                        WINDOW_LEN,
                        false,
                        runQueries
                );
            } else {
                Experiment.run(
                        dataFile,
                        queryFile,
                        hbi,
                        NGRAMS,
                        false,
                        false,
                        runQueries
                );
            }

            if (runSuffix) {
                IPMIndexing suffix = new StreamingSlidingWindowIndex(WINDOW_LEN, alphabet);
                if ("segments".equalsIgnoreCase(mode)) {
                    runSegmentsMode(
                            dataFile,
                            queryFile,
                            suffix,
                            1,
                            WINDOW_LEN,
                            false,
                            runQueries
                    );
                } else {
                    Experiment.run(
                            dataFile,
                            queryFile,
                            suffix,
                            1,
                            false,
                            false,
                            runQueries
                    );
                }
            }
        }

        double hbiTotalMs = 0.0;
        double hbiTotalMsInsert = 0.0;
        double suffixTotalMs = 0.0;
        double suffixTotalMsInsert = 0.0;

        double lpShareSum = 0.0;
        double avgLpLevelSum = 0.0; // average chosen LP level per run
        double avgQueryTimeSum = 0.0;
        double avgLpTimeSum = 0.0;
        int statsSamples = 0;
        int suffixRuns = 0;

        // -----------
        // Timed runs
        // -----------
        for (int i = 0; i < runs; i++) {

            // HBI pass
            HBI hbi = newHbi(Confidence, alphabet, fpRate);
            hbi.strides = USE_STRIDES;
            HbiStats stats = hbi.stats();
            stats.setCollecting(true);
            stats.setExperimentMode(false);

            ExperimentRunResult hbiRes;
            if ("segments".equalsIgnoreCase(mode)) {
                hbiRes = runSegmentsMode(
                        dataFile,
                        queryFile,
                        hbi,
                        NGRAMS,
                        WINDOW_LEN,
                        false,
                        runQueries
                );
            } else {
                hbiRes = Experiment.run(
                        dataFile,
                        queryFile,
                        hbi,
                        NGRAMS,
                        false,
                        false,
                        runQueries
                );
            }

            hbiTotalMs += hbiRes.totalRunTimeMs();
            hbiTotalMsInsert += hbiRes.totalInsertTimeMs();

            if (stats.totalQueryCount() > 0) {
                lpShareSum += stats.lpShareOfQuery();
                avgQueryTimeSum += stats.averageQueryTimeMillis();
                avgLpTimeSum += stats.averageLpTimeMillis();
                // compute average LP level chosen this run
                double runAvgLp = 0.0;
                List<Integer> lps = stats.lpLevels();
                if (!lps.isEmpty()) {
                    long sum = 0L;
                    for (int lp : lps) sum += lp;
                    runAvgLp = (double) sum / lps.size();
                }
                avgLpLevelSum += runAvgLp;
                statsSamples++;
            }

            ArrayList<ArrayList<Integer>> hbiMatches = hbiRes.matchRes();

            if (runSuffix) {
                IPMIndexing suffix = new StreamingSlidingWindowIndex(WINDOW_LEN, alphabet);

                ExperimentRunResult suffixRes;
                if ("segments".equalsIgnoreCase(mode)) {
                    suffixRes = runSegmentsMode(
                            dataFile,
                            queryFile,
                            suffix,
                            1,
                            WINDOW_LEN,
                            false,
                            runQueries
                    );
                } else {
                    suffixRes = Experiment.run(
                            dataFile,
                            queryFile,
                            suffix,
                            1,
                            false,
                            false,
                            runQueries
                    );
                }

                suffixTotalMs += suffixRes.totalRunTimeMs();
                suffixTotalMsInsert += suffixRes.totalInsertTimeMs();
                suffixRuns++;

                ArrayList<ArrayList<Integer>> suffixMatches = suffixRes.matchRes();
                if (runQueries) {
                    compared(hbiMatches, suffixMatches);
                }
            }
        }

        // --------------------------
        // Print aggregate statistics
        // --------------------------
        if (runs > 0) {
            System.out.printf(
                    Locale.ROOT,
                    "HBI avg (ms): %.3f%n",
                    hbiTotalMs / runs);
            System.out.printf(
                    Locale.ROOT,
                    "HBI Insert avg (ms): %.3f%n",
                    hbiTotalMsInsert / runs);
            System.out.printf(
                    Locale.ROOT,
                    "HBI Insert avg per symbol (ms): %.4f%n",
                    (hbiTotalMsInsert / runs) / WINDOW_LEN
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
                System.out.printf(
                        Locale.ROOT,
                        "HBI avg chosen LP level: %.3f%n",
                        (avgLpLevelSum / statsSamples)
                );
            }
        }

        if (suffixRuns > 0) {
            System.out.printf(
                    Locale.ROOT,
                    "SuffixTreeIndex avg (ms): %.3f%n",
                    suffixTotalMs / suffixRuns);
            System.out.printf(
                    Locale.ROOT,
                    "SuffixTreeIndex Insert avg (ms): %.3f%n",
                    suffixTotalMsInsert / suffixRuns);
            System.out.printf(
                    Locale.ROOT,
                    "SuffixTreeIndex Insert avg (ms) per char: %.3f%n",
                    (suffixTotalMsInsert / suffixRuns) / WINDOW_LEN
            );
        } else if (!runSuffix) {
            System.out.println("Suffix index was disabled; no suffix statistics collected.");
        }
    }
}
