import PMIndex.HbiConfiguration;
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
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.IntStream;

public class ConfidenceExperiment {

    // The workload: mawi21 part 1, uniform queries of length 10, read as segments, unless the
    // DataCatalog.Selection flags (--dataset, --part, --pattern-length, --data, --queries, ...) say otherwise.
    private static Workload WORKLOAD;

    // log2 of the window and tree length. Override with --log2w.
    private static int TextSize = 21;
    private static int WINDOW_LEN   = 1 << TextSize;
    private static int TREE_LEN     = 1 << TextSize;
    private static final double FP_RATE   = 0.0;
    private static final boolean LINEBREAK_DATASET = false; // not directly changed here, kept for compatibility

    // Run r forces Block Search to start at level r, so RUNS is the deepest start level evaluated.
    // It is recomputed in main from the shortest query in the workload.
    private static int RUNS     = (int) (TextSize - Math.ceil(Math.log(10)/Math.log(2)));

    // Cost function under test: maxprob (biased estimator of the thesis) or occurrence.
    private static String COST_FUNCTION = "maxprob";
    // Prefix for the CSV files, so parallel runs do not overwrite each other.
    private static String OUT_PREFIX = "";

    // n-grams for this experiment
    private static int NGRAMS = 1;

    // size of base alphabet for one symbol. We raise it to NGRAMS in main
    private static int ALPHABET = 130000;

    // Data structures for accuracy reporting.
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

    // Helper to manage cost model accuracy accounting.

    public static void main(String[] args) throws IOException {
        System.out.println("Starting experiment...");

        // Command line overrides: the DataCatalog.Selection flags (--dataset, --part, --pattern-length,
        // --query-type, --reader/--mode, --data, --queries, --data-root), and --log2w, --ngram,
        // --alphabet, --cf maxprob|occurrence, --out <csv prefix>
        DataCatalog.Selection selection = new DataCatalog.Selection("mawi21");
        for (int i = 0; i < args.length - 1; i += 2) {
            String value = args[i + 1];
            String key = args[i].toLowerCase(Locale.ROOT);
            if (key.startsWith("--") && selection.accept(key.substring(2), value)) {
                continue;
            }
            switch (key) {
                case "--log2w" -> {
                    TextSize = Integer.parseInt(value);
                    WINDOW_LEN = 1 << TextSize;
                    TREE_LEN = 1 << TextSize;
                }
                case "--ngram" -> NGRAMS = Integer.parseInt(value);
                case "--alphabet" -> ALPHABET = Integer.parseInt(value);
                case "--cf" -> COST_FUNCTION = value.toLowerCase(Locale.ROOT);
                case "--out" -> OUT_PREFIX = value;
                default -> throw new IllegalArgumentException("unknown option " + args[i]);
            }
        }
        WORKLOAD = selection.resolve();
        RUNS = deepestStartLevel(WORKLOAD);

        // Expand ALPHABET to n-gram space (same as your current code does with ALPHABET, NGRAMS)
        ALPHABET = (int) Math.pow(ALPHABET, NGRAMS);

        System.out.printf(
                "Window=%d, Tree=%d, sigma=%d, FP=%.3g, n-gram=%d, runs=%d, workload=%s, cf=%s%n",
                WINDOW_LEN, TREE_LEN, ALPHABET, FP_RATE, NGRAMS, (RUNS + 1), WORKLOAD, COST_FUNCTION);

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
            // Segments stream and re-run the queries after every window; characters query once at the end.
            WorkloadRunner.Options collect = WorkloadRunner.Options.QUIET.collectStats(true);
            ExperimentRunResult result = (WORKLOAD.mode() == Workload.ReaderMode.SEGMENTS)
                    ? WorkloadRunner.stream(WORKLOAD, hbi, NGRAMS, WINDOW_LEN, collect)
                    : WorkloadRunner.run(WORKLOAD, hbi, NGRAMS, collect);

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

        double sumRegret = 0.0;
        int regretCount = 0;
        double sumAbsErrAtPred = 0.0;
        int overAtPred = 0;
        int errCount = 0;
        for (PatternAccuracy accuracy : patternAccuracyMap.values()) {
            Double regret = accuracy.regret();
            if (regret != null) {
                sumRegret += regret;
                regretCount++;
            }
            Double err = accuracy.signedErrorAtPrediction();
            if (err != null) {
                sumAbsErrAtPred += Math.abs(err);
                if (err > 0) overAtPred++;
                errCount++;
            }
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
                "Predicted optimal Lp exactly +/-1 level: %d/%d (%.2f%%)%n",
                predictedStrictOne,
                evaluatedPatterns,
                predictedNearRate * 100.0);
        System.out.printf(
                Locale.ROOT,
                "Predicted optimal Lp within {0,+/-1}: %d/%d (%.2f%%)%n",
                predictedWithinOneInclusive,
                evaluatedPatterns,
                predictedNearInclusiveRate * 100.0);
        System.out.printf(
                Locale.ROOT,
                "Avg regret (probes at predicted Lp / best): %.3f over %d patterns%n",
                regretCount == 0 ? 0.0 : sumRegret / regretCount,
                regretCount);
        System.out.printf(
                Locale.ROOT,
                "Avg |rel. error| at predicted Lp: %.3f, overestimates: %.1f%% over %d patterns%n",
                errCount == 0 ? 0.0 : sumAbsErrAtPred / errCount,
                errCount == 0 ? 0.0 : 100.0 * overAtPred / errCount,
                errCount);
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
        CsvUtil.writeRows(Path.of(OUT_PREFIX + "runs_summary.csv"), runRows);
        CsvUtil.writeRows(Path.of(OUT_PREFIX + "patterns_summary.csv"), patternRows);
    }

    private static RunStats summarizeRun(
            int runIdx,
            ExperimentRunResult res,
            List<List<?>> patternRows,
            boolean verbose,
            Map<String, PatternAccuracy> patternAccuracy) {

        RunStats s = new RunStats();
        s.runIdx = runIdx;
        int i =0;
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
            accuracy.recordEstimate(pr.Lp(), est);

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
//            System.out.println("p: " + i + " err: " + relErr);
            i++;
        }

        if (s.minRow != null && s.maxRow != null && verbose) {
            System.out.printf(
                    Locale.ROOT,
                    "  -> minRelErr=%.4f (Lp=%d cf=%d act=%d est=%.1f)   maxRelErr=%.4f (Lp=%d cf=%d act=%d est=%.1f)%n",
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
        private final Map<Integer, Double> estimateByLp = new HashMap<>();
        private int bestActualProbes = Integer.MAX_VALUE;

        void recordEstimate(int lp, double estimate) {
            estimateByLp.putIfAbsent(lp, estimate);
        }

        // Probes at the predicted start level divided by the probes at the best start level.
        Double regret() {
            Integer predicted = probesForPredictedLp();
            if (predicted == null || bestActualProbes == Integer.MAX_VALUE || bestActualProbes <= 0) {
                return null;
            }
            return predicted / (double) bestActualProbes;
        }

        // Signed relative error of the estimate at the predicted start level.
        Double signedErrorAtPrediction() {
            for (int lp : predictedLps) {
                Integer actual = probesByLp.get(lp);
                Double estimate = estimateByLp.get(lp);
                if (actual != null && estimate != null && actual > 0) {
                    return (estimate - actual) / actual;
                }
            }
            return null;
        }

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


    private static CostFunction newCostFunction() {
        return switch (COST_FUNCTION) {
            case "maxprob" -> new CostFunctionMaxProb();
            case "occurrence" -> new CostFunctionOccurrence();
            default -> throw new IllegalArgumentException("unknown cost function " + COST_FUNCTION);
        };
    }

    // Deepest level whose nodes can still hold the shortest query (Block Search never starts lower).
    private static int deepestStartLevel(Workload workload) {
        int shortest = Integer.MAX_VALUE;
        for (Workload.Query q : workload.parsedQueries()) {
            int len = (q.tokens() != null) ? q.tokens().length : q.text().length();
            shortest = Math.min(shortest, len);
        }
        if (shortest == Integer.MAX_VALUE) {
            shortest = 10;
        }
        return MathUtils.deepestVisitedLevel(TREE_LEN, shortest);
    }

    private static HBI newHbi(double conf) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(TREE_LEN);
        Supplier<Membership> memFactory = MockMembership::new;
        Supplier<PruningPlan> prFactory = () -> new MostFreqPruning(conf, FP_RATE);
        Verifier v = new VerifierLinearLeafProbe();
        // Experiment mode: every query records a PatternResult with the probe counts per start level.
        return new HBI(HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(WINDOW_LEN)
                .fpRate(FP_RATE)
                .alphabetSize(ALPHABET)
                .treeLength(TREE_LEN)
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(v)
                .costFunction(newCostFunction())
                .confidence(conf)
                .nGram(NGRAMS)
                .experimentMode(true)
                .build());
    }
}
