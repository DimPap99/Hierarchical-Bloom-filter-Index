import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import estimators.*;
import membership.BloomFilter;
import membership.Membership;
import membership.MockMembership;
import search.*;
import utilities.*;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * ConfidenceExperiment runs multiple end to end workloads and collects accuracy and cost model diagnostics.
 *
 * After the changes in this version, we ALSO support two streaming modes:
 *
 *   mode="chars"
 *     Treat the dataset as a continuous character stream and queries as raw substrings.
 *     This is the original Experiment.run(...) pipeline.
 *
 *   mode="segments"
 *     Treat the dataset as a token-per-line stream. Each line (for example each packet record)
 *     is one symbol. We index sliding n-grams of those symbols. Queries are sequences of those
 *     symbols separated by spaces. This matches ProcessStream.
 *
 * You select the mode with --mode chars or --mode segments. Default is chars, so if you do not pass
 * --mode you get your original behavior unmodified.
 *
 * Everything else in this class (accuracy stats, timing summaries, CSV export) remains intact.
 *
 * The rest of this header comment is the original explanation of timing metrics:
 *
 * We measure per-run:
 *    averageQueryTimeMillis   average total query time per pattern
 *    averageLpTimeMillis      average pruning level estimation time per pattern
 *    lpShareOfQuery           fraction of query time spent only on pruning level estimation
 */
public class ConfidenceExperiment {

    /** Default file locations for convenience. */
    private static final String DATA_FILE    = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/pg2701.txt";
    private static final String QUERIES_FILE = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/pg2701/10.txt";

    private static final int TextSize = 20;
    private static final int WINDOW_LEN   = 1 << TextSize;
    private static final int TREE_LEN     = 1 << TextSize;
    private static final double FP_RATE   = 0.0;
    private static final boolean LINEBREAK_DATASET = false; // not directly changed here, kept for compatibility

    // How many times we rerun the whole workload to smooth variability
    private static final int RUNS     = (int) (TextSize - Math.ceil(Math.log(10)/Math.log(2)));

    // n-grams for this experiment
    private static int NGRAMS = 1;

    // size of base alphabet for one symbol. We raise it to NGRAMS in main
    private static int ALPHABET = 89;

    /** ---------------------------------------
     *  Data structures for accuracy reporting
     *  ---------------------------------------
     */
    private record PatternRow(
            int runIdx,
            int lp, int cfLp,
            int actualProbes,
            double estProbes,
            double relError,
            int patLen
    ) {}

    private static final class RunStats {
        int    runIdx;
        int    patterns;             // how many patterns with at least one probe
        double sumActual;
        double sumEstimated;
        double sumAbsError;
        double sumAbsRelError;
        double sumSqError;
        int    lpMatches;
        double minRelError = Double.POSITIVE_INFINITY;
        double maxRelError = 0.0;
        PatternRow minRow  = null;
        PatternRow maxRow  = null;

        int bloomProbes;
        int leafProbes;
        int overEstimations = 0;

        double overestimation(){return (double)overEstimations/patterns;}
        double mape() { return (patterns == 0) ? 0.0 : (sumAbsRelError / patterns); }
        double overallRelError() {
            return (sumActual <= 0.0) ? 0.0 : Math.abs(1.0 - (sumEstimated / sumActual));
        }
        double rmse() {
            return (patterns == 0) ? 0.0 : Math.sqrt(sumSqError / patterns);
        }
        double lpMatchRate() {
            return (patterns == 0) ? 0.0 : (lpMatches * 1.0 / patterns);
        }
    }

    /**
     * NEW helper:
     * Read queries file into memory as strings, trimming blank lines.
     * We use this in "segments" mode, and it mirrors how Experiment.run reads queries for "chars" mode.
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
     * NEW helper:
     * In "segments" mode, run the streaming experiment where each symbol is a whole line / packet record.
     *
     * This is conceptually ProcessStream.run(...) but modified to:
     *   1. collect PatternResult objects (index.getLatestStats()) so summarizeRun(...) still works
     *   2. compute avgQueryLen the same way ConfidenceExperiment expects it
     *
     * Important details for "segments" mode:
     *   - We treat the dataset as a sequence of tokens (Strings), one per segment.
     *   - We slide a RingBuffer<String> of length NGRAMS over those tokens.
     *   - Each time the buffer is full, we call index.insert(currentNgramAsString).
     *   - Every WINDOW_LEN tokens we run ALL queries.
     *   - A query line from the queries file is assumed to be "tok0 tok1 tok2 ...".
     *     We split by spaces, and build new Pattern(splitArray, NGRAMS).
     */
    private static ExperimentRunResult runStreamingSegments(
            String datasetPath,
            String queriesPath,
            IPMIndexing index,
            int ngram,
            int windowLen,
            boolean verbose,
            boolean collectResults
    ) throws IOException {

        // Load queries once
        List<String> queries = loadQueries(queriesPath);

        // average query length in n-grams
        int sumQueryLen = 0;
        for (String q : queries) {
            if (q == null || q.isEmpty()) continue;
            String[] splitArray = q.split(" ");
            Pattern p = new Pattern(splitArray, ngram);
            sumQueryLen += p.nGramToLong.length;
        }
        double avgQueryLen = queries.isEmpty() ? 0.0 : (sumQueryLen * 1.0 / queries.size());

        // SegmentReader gives you tokens (Strings), not characters.
        SegmentReader seg = new SegmentReader(
                datasetPath,
                queriesPath,
                '\n',    // delimiter for segments
                false,   // do not include delimiter in returned token
                true     // skip empty segments
        );
        SegmentReaderAdapter adapter = new SegmentReaderAdapter(seg);

        RingBuffer<String> window = new RingBuffer<>(ngram);

        long totalInsertNanos = 0L;
        long insertionEvents  = 0L;  // how many n-grams we actually inserted
        long tokensSeen       = 0L;  // how many segment tokens consumed

        long nextQueryAt      = Math.max(1, windowLen); // schedule query workload every windowLen tokens
        long sumQueryLoadMs   = 0L;
        long queryLoads       = 0L;

        ArrayList<PatternResult> patternResultsList = new ArrayList<>();

        while (adapter.hasNext()) {
            String tok = adapter.next();
            window.append(tok);

            if (!window.isFilled()) {
                continue;
            }

            // snapshot().toString() should give the canonical string form we index.
            // In your ProcessStream code you were doing exactly this:
            //   String s = window.snapshot().toString();
            String s = window.snapshot().toString();

            long t0 = System.nanoTime();
            index.insert(s);
            long t1 = System.nanoTime();

            totalInsertNanos += (t1 - t0);
            insertionEvents++;
            tokensSeen++;

            if (tokensSeen >= nextQueryAt) {
                long qMs = runQueryLoadSegments(index, queries, ngram, verbose, collectResults, patternResultsList);
                queryLoads++;
                sumQueryLoadMs += qMs;
                nextQueryAt += windowLen;
            }
        }

        // If the stream was shorter than the window threshold, make sure we still run the queries once
        if (tokensSeen > 0 && queryLoads == 0) {
            long qMs = runQueryLoadSegments(index, queries, ngram, verbose, collectResults, patternResultsList);
            queryLoads++;
            sumQueryLoadMs += qMs;
        }

        double totalInsertMs = totalInsertNanos / 1_000_000.0;
        double avgInsertMsPerSymbol =
                insertionEvents == 0 ? 0.0 : (totalInsertMs / insertionEvents);

        double totalQueryMs = sumQueryLoadMs;
        double avgQueryLoadMs =
                queryLoads == 0 ? 0.0 : (totalQueryMs / queryLoads);

        return new ExperimentRunResult(
                totalQueryMs,
                totalInsertMs,
                patternResultsList,
                avgQueryLen,
                avgInsertMsPerSymbol,
                avgQueryLoadMs,
                null // match results omitted here, but could be added later if you wish
        );
    }

    /**
     * NEW helper:
     * Issue the full query workload for "segments" mode.
     *
     * For each query line:
     *   - split the line on spaces to get an array of tokens
     *   - build Pattern(splitArray, ngram)
     *   - ask the index to report()
     *   - if collectResults is true, read index.getLatestStats() and append to outList
     *
     * This mirrors runQueryLoad(...) in ProcessStream, but now we are explicitly
     * accumulating PatternResult objects so that summarizeRun(...) still works.
     */
    private static long runQueryLoadSegments(
            IPMIndexing index,
            List<String> queries,
            int ngram,
            boolean verbose,
            boolean collectResults,
            ArrayList<PatternResult> outList
    ) {

        long start = System.currentTimeMillis();

        for (String q : queries) {
            if (q == null || q.isEmpty()) continue;

            // segments mode: query is space separated tokens
            String[] splitArray = q.split(" ");
            Pattern pat = new Pattern(splitArray, ngram);

            ArrayList<Integer> report = index.report(pat);

            if (verbose) {
                if (report.size() < 20 && q.length() < 200) {
                    System.out.println(q + ":" + report);
                }
            }

            if (collectResults) {
                PatternResult rr = index.getLatestStats();
                if (rr != null) {
                    outList.add(rr);
                }
            }
        }

        return System.currentTimeMillis() - start;
    }

    /**
     * Helper from your original code to manage cost model accuracy accounting.
     * I am leaving summarizeRun(...) and PatternAccuracy exactly the same.
     */

    public static void main(String[] args) throws IOException {
        System.out.println("Starting experiment…");

        // Default mode is character streaming, which is your current behavior
        String mode = "chars";

        // Allow minimal override via command line: --mode chars|segments
        for (int i = 0; i < args.length - 1; i += 2) {
            if ("--mode".equalsIgnoreCase(args[i])) {
                mode = args[i + 1];
            }
        }

        // Expand ALPHABET to n-gram space (same as your current code does with ALPHABET, NGRAMS)
        ALPHABET = (int) Math.pow(ALPHABET, NGRAMS);

        System.out.printf(
                "Window=%d, Tree=%d, σ=%d, FP=%.3g, n-gram=%d, runs=%d, mode=%s%n",
                WINDOW_LEN, TREE_LEN, ALPHABET, FP_RATE, NGRAMS, (RUNS + 1), mode);

        // Per-run summaries for runs_summary.csv
        List<List<?>> runRows = new ArrayList<>();
        runRows.add(List.of(
                "run",
                "patterns",
                "sumActualProbes",
                "sumEstimatedProbes",
                "overallRelError",
                "MAPE",
                "RMSE",
                "LpMatchRate",
                "avgQueryMs",
                "avgLpMs",
                "lpShare",
                "queryTimeMsTotal",
                "insertTimeMsTotal",
                "avgQueryLen"));

        // Optional per-pattern CSV rows
        List<List<?>> patternRows = new ArrayList<>();
        patternRows.add(List.of(
                "run",
                "Lp",
                "cfLp",
                "actualProbes",
                "estProbes",
                "relError",
                "patternLen"));

        // Cross-run aggregates
        double aggOverallRelErr = 0.0;
        double aggMAPE          = 0.0;
        double aggRMSE          = 0.0;
        double aggAvgMinCostLpMillis = 0.0;

        double aggAvgQueryMs = 0.0;
        double aggAvgLpMs    = 0.0;
        double aggLpShare    = 0.0;

        Map<String, PatternAccuracy> patternAccuracy = new HashMap<>();
        double avgOverallOverestimation = 0.0;
        long start = System.currentTimeMillis();

        for (int run = 0; run <= RUNS; run++) {

            // ---------------------------------
            // Build a fresh HBI index for this run
            // ---------------------------------
            HBI hbi = newHbi(0.99);
            hbi.isMarkov = false;
            hbi.strides = false;
            hbi.setLpOverride(run);
            hbi.stats().setCollecting(true);
            hbi.stats().setExperimentMode(true);

            // ---------------------------------
            // Run either chars-mode pipeline or segments-mode pipeline
            // ---------------------------------
            ExperimentRunResult result;
            if ("segments".equalsIgnoreCase(mode)) {
                result = runStreamingSegments(
                        DATA_FILE,
                        QUERIES_FILE,
                        hbi,
                        NGRAMS,
                        WINDOW_LEN,
                        false,   // verbose
                        true     // collect PatternResult objects
                );
            } else {
                // chars mode (original behavior)
                result = Experiment.run(
                        DATA_FILE,
                        QUERIES_FILE,
                        hbi,
                        NGRAMS,
                        false,   // verbose
                        true     // collect PatternResult objects
                );
            }

            // ---------------------------------
            // Summarize model accuracy for this run
            // ---------------------------------
            RunStats stats = summarizeRun(
                    run,
                    result,
                    patternRows,
                    true,
                    patternAccuracy);

            // Pull timing metrics from HBI
            double avgQueryMs = hbi.stats().averageQueryTimeMillis();
            double avgLpMs    = hbi.stats().averageLpTimeMillis();
            double lpShare    = hbi.stats().lpShareOfQuery(); // fraction in [0,1]

            // Print concise per-run summary
            System.out.printf(
                    Locale.ROOT,
                    "Run %2d: patterns=%4d  sumActual=%.0f  sumEst=%.1f  overallRelErr=%.4f  MAPE=%.4f  RMSE=%.2f  LpMatch=%.2f%n",
                    run,
                    stats.patterns,
                    stats.sumActual,
                    stats.sumEstimated,
                    stats.overallRelError(),
                    stats.mape(),
                    stats.rmse(),
                    stats.lpMatchRate());

            System.out.println(
                    "Leaf probes: " + stats.leafProbes + " Bloom Probes: " + stats.bloomProbes);

            double avgMinCostLpMillis = hbi.stats().averageMinCostLpTimeMillis();
            System.out.printf(
                    Locale.ROOT,
                    "Avg minCostLp time = %.3f ms%n",
                    avgMinCostLpMillis);

            // Also print timing share
            System.out.printf(
                    Locale.ROOT,
                    "Avg query time per pattern = %.3f ms, Lp estimation time = %.3f ms (%.2f%% of query time)%n",
                    avgQueryMs,
                    avgLpMs,
                    lpShare * 100.0);

            // Add this run to runs_summary.csv
            runRows.add(List.of(
                    run,
                    stats.patterns,
                    (long) stats.sumActual,
                    stats.sumEstimated,
                    stats.overallRelError(),
                    stats.mape(),
                    stats.rmse(),
                    stats.lpMatchRate(),
                    avgQueryMs,
                    avgLpMs,
                    lpShare,
                    result.totalRunTimeMs(),
                    result.totalInsertTimeMs(),
                    result.avgQuerySize()
            ));

            // Accumulate cross-run aggregates
            aggOverallRelErr       += stats.overallRelError();
            aggMAPE                += stats.mape();
            aggRMSE                += stats.rmse();
            avgOverallOverestimation += stats.overestimation();
            aggAvgMinCostLpMillis  += avgMinCostLpMillis;

            aggAvgQueryMs          += avgQueryMs;
            aggAvgLpMs             += avgLpMs;
            aggLpShare             += lpShare;
        }

        int runCount = RUNS + 1;

        // Final cross-run summary
        double avgOverallRelErr = aggOverallRelErr / runCount;
        double avgMAPE          = aggMAPE / runCount;
        double avgRMSE          = aggRMSE / runCount;
        double avgMinCostLpMillis = aggAvgMinCostLpMillis / runCount;

        double finalAvgQueryMs = aggAvgQueryMs / runCount;
        double finalAvgLpMs    = aggAvgLpMs / runCount;
        double finalAvgLpShare = aggLpShare  / runCount;

        System.out.println("\n=== Cross-run summary ===");
        System.out.printf(Locale.ROOT, "Avg overallRelError = %.4f%n", avgOverallRelErr);
        System.out.printf(Locale.ROOT, "Avg MAPE            = %.4f%n", avgMAPE);
        System.out.printf(Locale.ROOT, "Avg RMSE (probes)   = %.2f%n", avgRMSE);
        System.out.printf(Locale.ROOT, "Avg minCostLp time  = %.3f ms%n", avgMinCostLpMillis);

        System.out.printf(
                Locale.ROOT,
                "Avg query time per pattern (across runs) = %.3f ms%n",
                finalAvgQueryMs);
        System.out.printf(
                Locale.ROOT,
                "Avg Lp estimation time per pattern (across runs) = %.3f ms%n",
                finalAvgLpMs);
        System.out.printf(
                Locale.ROOT,
                "Mean fraction of query time spent in Lp estimation (across runs) = %.2f%%%n",
                finalAvgLpShare * 100.0);

        // Probe accuracy diagnostics across runs
        Map<String, PatternAccuracy> patternAccuracyMap = patternAccuracy;
        int predictedMatches = 0;
        int predictedStrictOne = 0;
        int predictedWithinOneInclusive = 0;
        long cfVsArbitraryGain = 0L;
        int cfVsArbitraryCount = 0;
        long arbitraryVsRootGain = 0L;
        int arbitraryVsRootCount = 0;
        int mispredictedTotal = 0;
        int mispredictedComparable = 0;
        int mispredictedPredBetter = 0;
        int mispredictedPredWorse = 0;

        for (PatternAccuracy accuracy : patternAccuracyMap.values()) {
            boolean matches = accuracy.predictedMatchesActual();
            if (matches) {
                predictedMatches++;
            }
            boolean withinStrictOne = accuracy.predictedWithinExactTolerance(1);
            boolean withinZeroOrOne = matches || withinStrictOne;
            if (withinStrictOne) {
                predictedStrictOne++;
            }
            if (withinZeroOrOne) {
                predictedWithinOneInclusive++;
            }

            Integer cfProbes = accuracy.probesForPredictedLp();
            Integer arbitraryProbes = accuracy.probesForArbitraryLp();
            Integer rootProbes = accuracy.probesForLp(0);
            boolean hasComparison = cfProbes != null && arbitraryProbes != null;

            if (hasComparison) {
                cfVsArbitraryGain += (long) arbitraryProbes - cfProbes;
                cfVsArbitraryCount++;
            }
            if (!withinZeroOrOne && accuracy.hasPredictions()) {
                mispredictedTotal++;
                if (hasComparison) {
                    mispredictedComparable++;
                    if (cfProbes < arbitraryProbes) {
                        mispredictedPredBetter++;
                    } else if (cfProbes > arbitraryProbes) {
                        mispredictedPredWorse++;
                    }
                }
            }
            if (arbitraryProbes != null && rootProbes != null) {
                arbitraryVsRootGain += (long) rootProbes - arbitraryProbes;
                arbitraryVsRootCount++;
            }
        }

        int evaluatedPatterns = patternAccuracyMap.size();
        double predictedMatchRate = evaluatedPatterns == 0 ? 0.0 : (predictedMatches * 1.0 / evaluatedPatterns);
        double predictedNearRate  = evaluatedPatterns == 0 ? 0.0 : (predictedStrictOne * 1.0 / evaluatedPatterns);
        double predictedNearInclusiveRate =
                evaluatedPatterns == 0 ? 0.0 : (predictedWithinOneInclusive * 1.0 / evaluatedPatterns);

        double avgCfVsArbitraryGain =
                cfVsArbitraryCount == 0 ? 0.0 : (cfVsArbitraryGain * 1.0 / cfVsArbitraryCount);
        double avgArbitraryVsRootGain =
                arbitraryVsRootCount == 0 ? 0.0 : (arbitraryVsRootGain * 1.0 / arbitraryVsRootCount);

        double mispredBetterRate =
                mispredictedComparable == 0 ? 0.0 : (mispredictedPredBetter * 1.0 / mispredictedComparable);
        double mispredWorseRate  =
                mispredictedComparable == 0 ? 0.0 : (mispredictedPredWorse * 1.0 / mispredictedComparable);

        avgOverallOverestimation /= runCount;
        double u = 1f - avgOverallOverestimation;

        System.out.printf(
                Locale.ROOT,
                "Predicted optimal Lp matches actual best: %d/%d (%.2f%%)%n",
                predictedMatches,
                evaluatedPatterns,
                predictedMatchRate * 100.0);
        System.out.printf(
                Locale.ROOT,
                "Predicted optimal Lp exactly ±1 level: %d/%d (%.2f%%)%n",
                predictedStrictOne,
                evaluatedPatterns,
                predictedNearRate * 100.0);
        System.out.printf(
                Locale.ROOT,
                "Predicted optimal Lp within {0,±1}: %d/%d (%.2f%%)%n",
                predictedWithinOneInclusive,
                evaluatedPatterns,
                predictedNearInclusiveRate * 100.0);
        System.out.printf(
                Locale.ROOT,
                "Avg probe reduction (CF vs arbitrary): %.2f over %d patterns%n",
                avgCfVsArbitraryGain,
                cfVsArbitraryCount);
        System.out.printf(
                Locale.ROOT,
                "Avg probe reduction (arbitrary vs Lp=0): %.2f over %d patterns%n",
                avgArbitraryVsRootGain,
                arbitraryVsRootCount);
        System.out.printf(
                Locale.ROOT,
                "Mispredicted cases (>|1| off optimal): %d  with arbitrary comparison: %d  cf better: %.2f%%  arbitrary better: %.2f%%%n",
                mispredictedTotal,
                mispredictedComparable,
                mispredBetterRate * 100.0,
                mispredWorseRate * 100.0);

        System.out.printf(
                Locale.ROOT,
                "Overall Overestimation: " + avgOverallOverestimation + " Underestimation: " + u);

        long end = System.currentTimeMillis() - start;
        System.out.println("\nTime taken: " + end);

        // Write CSVs
        CsvUtil.writeRows(Path.of("runs_summary.csv"), runRows);
        CsvUtil.writeRows(Path.of("patterns_summary.csv"), patternRows);
    }

    private static RunStats summarizeRun(
            int runIdx,
            ExperimentRunResult res,
            List<List<?>> patternRows,
            boolean verbose,
            Map<String, PatternAccuracy> patternAccuracy) {

        RunStats s = new RunStats();
        s.runIdx = runIdx;

        for (PatternResult pr : res.patternResults()) {
            int    actual = pr.probes();
            double est    = pr.predictedCost();
            int    leafprobes = pr.leafProbes();
            if (actual <= 0) continue; // skip degenerate cases

            double diff   = est - actual;
            double relErr = Math.abs(1.0 - (est / actual));

            String patternKey = pr.p().patternTxt;
            PatternAccuracy accuracy = patternAccuracy.computeIfAbsent(
                    patternKey, k -> new PatternAccuracy());
            accuracy.recordPrediction(pr.cfLp());
            accuracy.recordArbitrary(pr.arbitraryConfLp());
            accuracy.recordObservation(pr.Lp(), actual);

            s.patterns++;
            s.sumActual      += actual;
            s.sumEstimated   += est;
            s.sumAbsError    += Math.abs(diff);
            s.sumAbsRelError += relErr;
            s.sumSqError     += diff * diff;
            s.leafProbes     += leafprobes;
            s.bloomProbes    += actual - leafprobes;
            if (est >= actual) s.overEstimations++;

            if (pr.Lp() == pr.cfLp()) s.lpMatches++;

            // track extremes
            if (relErr < s.minRelError) {
                s.minRelError = relErr;
                s.minRow = new PatternRow(
                        runIdx,
                        pr.Lp(),
                        pr.cfLp(),
                        actual,
                        est,
                        relErr,
                        pr.p().nGramToLong.length);
            }
            if (relErr > s.maxRelError) {
                s.maxRelError = relErr;
                s.maxRow = new PatternRow(
                        runIdx,
                        pr.Lp(),
                        pr.cfLp(),
                        actual,
                        est,
                        relErr,
                        pr.p().nGramToLong.length);
            }

            // optional per-pattern CSV
            patternRows.add(List.of(
                    runIdx,
                    pr.Lp(),
                    pr.cfLp(),
                    actual,
                    est,
                    relErr,
                    pr.p().nGramToLong.length
            ));
        }

        if (s.minRow != null && s.maxRow != null && verbose) {
            System.out.printf(
                    Locale.ROOT,
                    "  ↳ minRelErr=%.4f (Lp=%d cf=%d act=%d est=%.1f)   maxRelErr=%.4f (Lp=%d cf=%d act=%d est=%.1f)%n",
                    s.minRow.relError,
                    s.minRow.lp,
                    s.minRow.cfLp,
                    s.minRow.actualProbes,
                    s.minRow.estProbes,
                    s.maxRow.relError,
                    s.maxRow.lp,
                    s.maxRow.cfLp,
                    s.maxRow.actualProbes,
                    s.maxRow.estProbes);

            double pctOver = s.overestimation();
            double under = 1f - pctOver;
            System.out.println("Pct Overstimate: " + pctOver);
            System.out.println("Pct Underestimated: " + under);
        }

        return s;
    }

    private static final class PatternAccuracy {
        private final Set<Integer> predictedLps = new HashSet<>();
        private final Set<Integer> bestActualLps = new HashSet<>();
        private final Set<Integer> arbitraryLps = new HashSet<>();
        private final Map<Integer, Integer> probesByLp = new HashMap<>();
        private int bestActualProbes = Integer.MAX_VALUE;

        void recordPrediction(int predictedLp) {
            predictedLps.add(predictedLp);
        }

        boolean hasPredictions() {
            return !predictedLps.isEmpty();
        }

        void recordArbitrary(int arbitraryLp) {
            arbitraryLps.add(arbitraryLp);
        }

        void recordObservation(int actualLp, int probes) {
            probesByLp.merge(actualLp, probes, Math::min);
            if (probes < bestActualProbes) {
                bestActualProbes = probes;
                bestActualLps.clear();
                bestActualLps.add(actualLp);
            } else if (probes == bestActualProbes) {
                bestActualLps.add(actualLp);
            }
        }

        boolean predictedMatchesActual() {
            if (bestActualLps.isEmpty() || predictedLps.isEmpty()) {
                return false;
            }
            for (int predicted : predictedLps) {
                if (bestActualLps.contains(predicted)) {
                    return true;
                }
            }
            return false;
        }

        boolean predictedWithinTolerance(int tolerance) {
            if (bestActualLps.isEmpty() || predictedLps.isEmpty()) {
                return false;
            }
            for (int predicted : predictedLps) {
                for (int actual : bestActualLps) {
                    if (Math.abs(predicted - actual) <= tolerance) {
                        return true;
                    }
                }
            }
            return false;
        }

        boolean predictedWithinExactTolerance(int tolerance) {
            if (bestActualLps.isEmpty() || predictedLps.isEmpty()) {
                return false;
            }
            for (int predicted : predictedLps) {
                for (int actual : bestActualLps) {
                    if (Math.abs(predicted - actual) == tolerance) {
                        return true;
                    }
                }
            }
            return false;
        }

        Integer probesForLp(int lp) {
            return probesByLp.get(lp);
        }

        Integer probesForPredictedLp() {
            return bestProbesForSet(predictedLps);
        }

        Integer probesForArbitraryLp() {
            return bestProbesForSet(arbitraryLps);
        }

        private Integer bestProbesForSet(Set<Integer> lps) {
            Integer best = null;
            for (int lp : lps) {
                Integer probes = probesByLp.get(lp);
                if (probes == null) continue;
                if (best == null || probes < best) {
                    best = probes;
                }
            }
            return best;
        }
    }

    /**
     * Your existing helper that constructs a new HBI.
     * I have intentionally left this unchanged so that behavior (Bloom filter config,
     * pruning plan, etcetera) remains identical to what you are currently evaluating.
     *
     * NOTE:
     * Right now this hard codes ngrams=3 in the constructor call to HBI.
     * I am not changing that because you asked to keep behavior the same.
     * If you truly want to honor NGRAMS here dynamically, replace that 3 with NGRAMS.
     */
    private static HBI newHbi(double conf) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(TREE_LEN);
        Supplier<Membership> memFactory = MockMembership::new;
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
                new CostFunctionMaxProb(),
                conf,
                1
        );
    }
}
