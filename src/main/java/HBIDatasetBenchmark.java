import PMIndex.HBI;
import PMIndex.HbiStats;
import PMIndex.IPMIndexing;
import PMIndex.SuffixTreeIndex;
import estimators.*;
import membership.BloomFilter;
import membership.Membership;
import search.*;
import utilities.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.Supplier;
import java.util.stream.IntStream;
import org.openjdk.jol.info.ClassLayout;
import org.openjdk.jol.info.GraphLayout;

import javax.xml.stream.FactoryConfigurationError;

public class HBIDatasetBenchmark {

    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/wzipf20_e1/1/1_Wzipf20_e1.txt";


    private static final int WINDOW_LEN   = 1 << 20;//1 << 21;
    private static final int TREE_LEN     = 1 << 20;
    private static int ALPHABET     = 32;
    private static final double FP_RATE   = 0.001;
    private static final int RUNS         = 2;        // set to 0 for a dry run
    private static final boolean USE_STRIDES = true;
    private static int NGRAMS = 1;
    private static String QUERY_FILE = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/wzipf20_e1/1/10.missing.txt";
    private static int NUMQUERIES = 135;

    public static void compared(ArrayList<ArrayList<Integer>> arr1, ArrayList<ArrayList<Integer>> arr2) {
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
    public static void main(String[] args) throws IOException {

        List<String> queryFiles = new ArrayList<>();

        List<List<?>> rows =  new ArrayList<>();
        List<?> header = List.of("nGram", "runAvgMs", "insertAvgMs", "avgQueryLength", "index");
        rows.add(header);

        ExperimentRunResult runResult;
        double suffixTotalMs = 0;
        double suffixTotalMsInsert = 0;


        System.out.println("Running queries for " + QUERY_FILE);
        double hbiTotalMs = 0;
        suffixTotalMs = 0;
        double hbiTotalMsInsert = 0;
        suffixTotalMsInsert = 0;
        double avgLp = 0;
        double lpShareSum = 0;
        double avgQueryTimeSum = 0;
        double avgLpTimeSum = 0;
        int statsSamples = 0;
        System.out.println("N-gram: " + NGRAMS);
        System.out.println("Window Size: " + WINDOW_LEN);
        System.out.println("Tree Length: " + TREE_LEN);
        ALPHABET = (int) Math.pow(ALPHABET, NGRAMS);
        ALPHABET = Math.min(ALPHABET, TREE_LEN);
        System.out.println("Alphabet: " + ALPHABET);
        System.out.println("\n");
        double avgAlpha = 0;
        /* JIT warm-up so HotSpot reaches steady state */

        for (int i = 0; i < 1; i++) {
            HBI hbi = newHbi(0.999);
            hbi.strides = USE_STRIDES;
            hbi.stats().setCollecting(false);
            hbi.stats().setExperimentMode(false);

            ArrayList<ArrayList<Integer>> warmHbi =
                    Experiment.run(DATA_FILE, QUERY_FILE, hbi, NGRAMS, false, false).matchRes();

            SuffixTreeIndex suffix = new SuffixTreeIndex(ALPHABET, 0.0001, WINDOW_LEN);
            ArrayList<ArrayList<Integer>> warmSuffix =
                    Experiment.run(DATA_FILE, QUERY_FILE, suffix, NGRAMS, false, false).matchRes();
////            suffix.debugDumpForQuery(" the repub", 3, 2048);
////            System.out.println(new MemUtil().jolMemoryReportPartitioned(hbi));
//            suffix.compactForQuerying();
//            System.out.println(suffix.jolMemoryReportPartitioned());
            compared(warmHbi, warmSuffix);
            int b = 2;
        }
                ArrayList<Long> timings;
        for (int i = 0; i < RUNS; i++) {

            HBI hbi = newHbi(0.99);
            hbi.strides = USE_STRIDES;
            HbiStats stats = hbi.stats();
            hbi.stats().setCollecting(false);
            hbi.stats().setExperimentMode(false);
//            stats.setCollecting(true);
            runResult = Experiment.run(DATA_FILE, QUERY_FILE, hbi, NGRAMS, false, false);
            ArrayList<ArrayList<Integer>> hbiMatches = runResult.matchRes();
            hbiTotalMs += runResult.totalRunTimeMs();
            hbiTotalMsInsert += runResult.totalInsertTimeMs();
            if (stats.totalQueryCount() > 0) {
                lpShareSum += stats.lpShareOfQuery();
                avgQueryTimeSum += stats.averageQueryTimeMillis();
                avgLpTimeSum += stats.averageLpTimeMillis();
                statsSamples++;
            }
            IPMIndexing suffix = new SuffixTreeIndex(ALPHABET, 0.0001, WINDOW_LEN);
            runResult = Experiment.run(DATA_FILE, QUERY_FILE, suffix, 1, false, false);
            ArrayList<ArrayList<Integer>> suffixMatches = runResult.matchRes();
            suffixTotalMs += runResult.totalRunTimeMs();
            suffixTotalMsInsert += runResult.totalInsertTimeMs();
//            MemUtil memUtil = new MemUtil();
//            compared(hbiMatches, suffixMatches);

        }

        if (RUNS > 0) {
            System.out.printf("HBI avg (ms): %.3f%n", hbiTotalMs / RUNS);
            System.out.printf("HBI Insert avg (ms): %.3f%n", hbiTotalMsInsert / RUNS);
            System.out.printf("HBI Insert avg per symbol (ms): %.4f%n", (hbiTotalMsInsert / RUNS)/WINDOW_LEN);

            if (statsSamples > 0) {
                System.out.printf("HBI avg query time per pattern (ms): %.3f%n", avgQueryTimeSum / statsSamples);
                System.out.printf("HBI avg LP computation time per pattern (ms): %.3f%n", avgLpTimeSum / statsSamples);
                System.out.printf("HBI LP time share of query (%%): %.2f%n", (lpShareSum / statsSamples) * 100.0);
            }
            System.out.println("Avg LP: " + avgLp);
            System.out.println("Avg Alpha: " + avgAlpha);

            System.out.printf("SuffixTreeIndex avg (ms): %.3f%n", suffixTotalMs / RUNS);
            System.out.printf("SuffixTreeIndex Insert avg (ms): %.3f%n", suffixTotalMsInsert / RUNS);
            System.out.println("\n");

        }
//                    rows.add(List.of(n, hbiTotalMs / RUNS,  hbiTotalMsInsert / RUNS, avgQueryLength, "hbi"));

            }
//                rows.add(List.of(1, ipmTotalMs / RUNS,  ipmTotalMsInsert / RUNS, avgQueryLength, "regex"));





//        CsvUtil.writeRows(Path.of("pg2701_arb.csv"), rows);


    // Helper that builds a fresh HBI wired to suppliers each time
    private static HBI newHbi(double conf) {
        Supplier<Estimator> estFactory =
                () -> new HashMapEstimator(TREE_LEN);

        Supplier<Membership> memFactory =
                () -> new BloomFilter();
        Supplier<PruningPlan> prFactory =
                () -> new MostFreqPruning(conf);

        Verifier v = new VerifierLinearLeafProbe();
        return new HBI(new BlockSearch(),

                WINDOW_LEN,
                FP_RATE,
                ALPHABET,
                TREE_LEN,
                estFactory,
                memFactory,
                prFactory,
                v, new CostFunctionDefaultRoot(), conf, NGRAMS);
    }


}
