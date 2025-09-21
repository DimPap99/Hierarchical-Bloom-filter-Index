import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import estimators.CostFunctionMaxProb;
import estimators.Estimator;
import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;
import search.*;
import utilities.*;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Supplier;
import java.util.stream.IntStream;

public class ConfidenceExperiment {

    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/pg2701.txt";
    private static final String QUERIES_FILE= "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/pg2701/unique_substrings_pg2701_1000.txt";

    private static final int TextSize = 20;
    private static final int WINDOW_LEN   = 1 << TextSize;
    private static final int TREE_LEN     = 1 << TextSize;
    private static final double FP_RATE   = 0.001;
    private static final boolean LINEBREAK_DATASET = false;

    // Controls how many times we rerun the whole workload to average out JIT etc.
    private static final int RUNS     = (int) (TextSize - Math.ceil(Math.log(1000)/Math.log(2)));


    // N-grams for this experiment (you can parameterize if needed)
    private static int NGRAMS = 2;

    private static int ALPHABET = 89;

    /** Per-pattern row for (optional) CSV dump. */
    private record PatternRow(
            int runIdx,
            int lp, int cfLp,
            int actualProbes,
            double estProbes,
            double relError,   // |1 - est/act|
            int patLen
    ) {}

    /** Aggregates for one run. */
    private static final class RunStats {
        int    runIdx;
        int    patterns;             // counted (non-zero actual probes)
        double sumActual;
        double sumEstimated;
        double sumAbsError;          // Σ |est - actual|
        double sumAbsRelError;       // Σ |1 - est/actual|
        double sumSqError;           // Σ (est - actual)^2
        int    lpMatches;            // count of patterns with Lp == cfLp
        // extremes
        double minRelError = Double.POSITIVE_INFINITY;
        double maxRelError = 0.0;
        PatternRow minRow  = null;
        PatternRow maxRow  = null;

        int bloomProbes;
        int leafProbes;

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
    //356578
    //258699
    //53372
    public static void main(String[] args) throws IOException {
        System.out.println("Starting experiment…");
        ALPHABET = (int) Math.pow(ALPHABET, NGRAMS);

        System.out.printf("Window=%d, Tree=%d, σ=%d, FP=%.3g, n-gram=%d, runs=%d%n",
                WINDOW_LEN, TREE_LEN, ALPHABET, FP_RATE, NGRAMS, RUNS);

        // Per-run summaries
        List<List<?>> runRows = new ArrayList<>();
        runRows.add(List.of("run",
                "patterns",
                "sumActualProbes",
                "sumEstimatedProbes",
                "overallRelError",
                "MAPE",
                "RMSE",
                "LpMatchRate",
                "queryTimeMs",
                "insertTimeMs",
                "avgQueryLen"));

        // Optional per-pattern CSV (comment out if too big)
        List<List<?>> patternRows = new ArrayList<>();
        patternRows.add(List.of("run", "Lp", "cfLp", "actualProbes", "estProbes", "relError", "patternLen"));

        // Cross-run aggregates (over the per-run summaries)
        double aggOverallRelErr = 0.0;
        double aggMAPE          = 0.0;
        double aggRMSE          = 0.0;

        Map<String, PatternAccuracy> patternAccuracy = new HashMap<>();

        for (int run = 0; run <= RUNS; run++) {
            // --- Build a fresh HBI and stream data ---
            HBI hbi = newHbi(0.99);
            hbi.setLpOverride(run);
            hbi.resetAlphabetMap(ALPHABET);
            hbi.stats().setCollecting(true);

            ExperimentRunResult result = Experiment.run(DATA_FILE, QUERIES_FILE, hbi, NGRAMS, false, true);

            // --- Summarize this run ---
            RunStats stats = summarizeRun(run, result, patternRows, true, patternAccuracy);

            // Print concise per-run summary
            System.out.printf(Locale.ROOT,
                    "Run %2d: patterns=%4d  sumActual=%.0f  sumEst=%.1f  overallRelErr=%.4f  MAPE=%.4f  RMSE=%.2f  LpMatch=%.2f%n",
                    run, stats.patterns, stats.sumActual, stats.sumEstimated,
                    stats.overallRelError(), stats.mape(), stats.rmse(), stats.lpMatchRate());
            System.out.println("Leaf probes: " + stats.leafProbes + " Bloom Probes: " + stats.bloomProbes);


            // Add run row for CSV
            runRows.add(List.of(
                    run,
                    stats.patterns,
                    (long) stats.sumActual,
                    stats.sumEstimated,
                    stats.overallRelError(),
                    stats.mape(),
                    stats.rmse(),
                    stats.lpMatchRate(),
                    result.totalRunTimeMs(),
                    result.totalInsertTimeMs(),
                    result.avgQuerySize()
            ));

            // Accumulate for cross-run
            aggOverallRelErr += stats.overallRelError();
            aggMAPE          += stats.mape();
            aggRMSE          += stats.rmse();
        }

        // --- Final cross-run summary ---
        double avgOverallRelErr = aggOverallRelErr / RUNS;
        double avgMAPE          = aggMAPE / RUNS;
        double avgRMSE          = aggRMSE / RUNS;

        System.out.println("\n=== Cross-run summary ===");
        System.out.printf(Locale.ROOT, "Avg overallRelError = %.4f%n", avgOverallRelErr);
        System.out.printf(Locale.ROOT, "Avg MAPE            = %.4f%n", avgMAPE);
        System.out.printf(Locale.ROOT, "Avg RMSE (probes)   = %.2f%n", avgRMSE);

        // --- Predicted vs actual optimal level statistics ---
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
        for (PatternAccuracy accuracy : patternAccuracy.values()) {
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
        int evaluatedPatterns = patternAccuracy.size();
        double predictedMatchRate = evaluatedPatterns == 0 ? 0.0 : (predictedMatches * 1.0 / evaluatedPatterns);
        double predictedNearRate  = evaluatedPatterns == 0 ? 0.0 : (predictedStrictOne * 1.0 / evaluatedPatterns);
        double predictedNearInclusiveRate = evaluatedPatterns == 0 ? 0.0 : (predictedWithinOneInclusive * 1.0 / evaluatedPatterns);
        double avgCfVsArbitraryGain = cfVsArbitraryCount == 0 ? 0.0 : (cfVsArbitraryGain * 1.0 / cfVsArbitraryCount);
        double avgArbitraryVsRootGain = arbitraryVsRootCount == 0 ? 0.0 : (arbitraryVsRootGain * 1.0 / arbitraryVsRootCount);
        double mispredBetterRate = mispredictedComparable == 0 ? 0.0 : (mispredictedPredBetter * 1.0 / mispredictedComparable);
        double mispredWorseRate = mispredictedComparable == 0 ? 0.0 : (mispredictedPredWorse * 1.0 / mispredictedComparable);
        System.out.printf(Locale.ROOT,
                "Predicted optimal Lp matches actual best: %d/%d (%.2f%%)%n",
                predictedMatches,
                evaluatedPatterns,
                predictedMatchRate * 100.0);
        System.out.printf(Locale.ROOT,
                "Predicted optimal Lp exactly ±1 level: %d/%d (%.2f%%)%n",
                predictedStrictOne,
                evaluatedPatterns,
                predictedNearRate * 100.0);
        System.out.printf(Locale.ROOT,
                "Predicted optimal Lp within {0,±1}: %d/%d (%.2f%%)%n",
                predictedWithinOneInclusive,
                evaluatedPatterns,
                predictedNearInclusiveRate * 100.0);
        System.out.printf(Locale.ROOT,
                "Avg probe reduction (CF vs arbitrary): %.2f over %d patterns%n",
                avgCfVsArbitraryGain,
                cfVsArbitraryCount);
        System.out.printf(Locale.ROOT,
                "Avg probe reduction (arbitrary vs Lp=0): %.2f over %d patterns%n",
                avgArbitraryVsRootGain,
                arbitraryVsRootCount);
        System.out.printf(Locale.ROOT,
                "Mispredicted cases (>|1| off optimal): %d  with arbitrary comparison: %d  cf better: %.2f%%  arbitrary better: %.2f%%%n",
                mispredictedTotal,
                mispredictedComparable,
                mispredBetterRate * 100.0,
                mispredWorseRate * 100.0);

        // --- Write CSVs ---
        CsvUtil.writeRows(Path.of("runs_summary.csv"), runRows);
        CsvUtil.writeRows(Path.of("patterns_summary.csv"), patternRows);
    }

    /** Computes all per-run stats, optionally filling per-pattern rows. */
    private static RunStats summarizeRun(int runIdx,
                                         ExperimentRunResult res,
                                         List<List<?>> patternRows,
                                         boolean verbose,
                                         Map<String, PatternAccuracy> patternAccuracy) {
        RunStats s = new RunStats();
        s.runIdx = runIdx;

        for (PatternResult pr : res.patternResults()) {
            int    actual = pr.probes();
            double est    = pr.predictedCost();   // already in "probes" units
            int leafprobes = pr.leafProbes();
            if (actual <= 0) continue;            // avoid divide-by-zero; skip empty probe cases

            double diff   = est - actual;
            double relErr = Math.abs(1.0 - (est / actual));

            String patternKey = pr.p().patternTxt;
            PatternAccuracy accuracy = patternAccuracy.computeIfAbsent(patternKey, k -> new PatternAccuracy());
            accuracy.recordPrediction(pr.cfLp());
            accuracy.recordArbitrary(pr.arbitraryConfLp());
            accuracy.recordObservation(pr.Lp(), actual);

            s.patterns++;
            s.sumActual     += actual;
            s.sumEstimated  += est;
            s.sumAbsError   += Math.abs(diff);
            s.sumAbsRelError+= relErr;
            s.sumSqError    += diff * diff;
            s.leafProbes += leafprobes;
            s.bloomProbes += actual - leafprobes;
            if (pr.Lp() == pr.cfLp()) s.lpMatches++;

            // track extremes
            if (relErr < s.minRelError) {
                s.minRelError = relErr;
                s.minRow = new PatternRow(runIdx, pr.Lp(), pr.cfLp(), actual, est, relErr, pr.p().nGramToInt.length);
            }
            if (relErr > s.maxRelError) {
                s.maxRelError = relErr;
                s.maxRow = new PatternRow(runIdx, pr.Lp(), pr.cfLp(), actual, est, relErr, pr.p().nGramToInt.length);
            }

            // optional detailed CSV
            patternRows.add(List.of(
                    runIdx, pr.Lp(), pr.cfLp(), actual, est, relErr, pr.p().nGramToInt.length
            ));
        }

        // (Optional) print extremes for debug
        if (s.minRow != null && s.maxRow != null && verbose) {
            System.out.printf(Locale.ROOT,
                    "  ↳ minRelErr=%.4f (Lp=%d cf=%d act=%d est=%.1f)   maxRelErr=%.4f (Lp=%d cf=%d act=%d est=%.1f)%n",
                    s.minRow.relError, s.minRow.lp, s.minRow.cfLp, s.minRow.actualProbes, s.minRow.estProbes,
                    s.maxRow.relError, s.maxRow.lp, s.maxRow.cfLp, s.maxRow.actualProbes, s.maxRow.estProbes);
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

    // Helper: fresh HBI wired to suppliers each time
    private static HBI newHbi(double conf) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(TREE_LEN);
        Supplier<Membership> memFactory = BloomFilter::new;
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
                conf
        );
    }
}
