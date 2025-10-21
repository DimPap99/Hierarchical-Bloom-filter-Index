import PMIndex.HBI;
import PMIndex.HbiStats;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
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
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.IntStream;
import org.openjdk.jol.info.ClassLayout;
import org.openjdk.jol.info.GraphLayout;

import javax.xml.stream.FactoryConfigurationError;

public class HBIDatasetBenchmark {

    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/w21/1/1_W21.txt";


    private static final int WINDOW_LEN   = 1 << 21;//1 << 21;
    private static final int TREE_LEN     = 1 << 21;
    private static int ALPHABET     = 89;
    private static final double FP_RATE   = 0.001;
    private static final int RUNS         = 3;        // set to 0 for a dry run
    private static final boolean USE_STRIDES = true;
    private static int NGRAMS = 10;
    private static String QUERY_FILE = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/w21/1/10.uniform.txt";
    private static int NUMQUERIES = 135;

    public void compared(ArrayList<ArrayList<Integer>> arr1, ArrayList<ArrayList<Integer>> arr2){


        for(int i =0; i<arr1.size(); i++){
            ArrayList<Integer> row1 = arr1.get(i);
            ArrayList<Integer> row2 = arr2.get(i);


        }
    }
    public static void main(String[] args) throws IOException {

        List<String> queryFiles = new ArrayList<>();

        List<List<?>> rows =  new ArrayList<>();
        List<?> header = List.of("nGram", "runAvgMs", "insertAvgMs", "avgQueryLength", "index");
        rows.add(header);

        ExperimentRunResult runResult;
        double ipmTotalMs = 0;
        double ipmTotalMsInsert = 0;


        System.out.println("Running queries for " + QUERY_FILE);
        double hbiTotalMs = 0;
        ipmTotalMs = 0;
        double hbiTotalMsInsert = 0;
        ipmTotalMsInsert = 0;
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

            Experiment.run(DATA_FILE, QUERY_FILE, hbi, NGRAMS, true, false);

            IPMIndexing ipm = new RegexIndex();
            Experiment.run(DATA_FILE, QUERY_FILE, ipm, 1, true, false);
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
            hbiTotalMs += runResult.totalRunTimeMs();
            hbiTotalMsInsert += runResult.totalInsertTimeMs();
            if (stats.totalQueryCount() > 0) {
                lpShareSum += stats.lpShareOfQuery();
                avgQueryTimeSum += stats.averageQueryTimeMillis();
                avgLpTimeSum += stats.averageLpTimeMillis();
                statsSamples++;
            }
            IPMIndexing ipm = new RegexIndex();
            runResult = Experiment.run(DATA_FILE, QUERY_FILE, ipm, 1, false, false);
            ipmTotalMs += runResult.totalRunTimeMs();
            ipmTotalMsInsert += runResult.totalInsertTimeMs();
            MemUtil memUtil = new MemUtil();


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

            System.out.printf("RegexIndex avg (ms): %.3f%n", ipmTotalMs / RUNS);
            System.out.printf("RegexIndex Insert avg (ms): %.3f%n", ipmTotalMsInsert / RUNS);
            System.out.println("\n");

        }
//                    rows.add(List.of(n, hbiTotalMs / RUNS,  hbiTotalMsInsert / RUNS, avgQueryLength, "hbi"));

            }
//                rows.add(List.of(1, ipmTotalMs / RUNS,  ipmTotalMsInsert / RUNS, avgQueryLength, "regex"));





//        CsvUtil.writeRows(Path.of("pg2701_arb.csv"), rows);


    // Helper that builds a fresh HBI wired to suppliers each time
    private static HBI newHbi(double conf) {
        Supplier<Estimator> estFactory =
                () -> new CSEstimator(TREE_LEN, 5, 16384 );//new HashMapEstimator(TREE_LEN);

        Supplier<Membership> memFactory =
                () -> new BloomFilter();
        Supplier<PruningPlan> prFactory =
                () -> new MultiLevelPruning(conf);

        Verifier v = new VerifierLinearLeafProbe();
        return new HBI(new BlockSearchCharSet(),

                WINDOW_LEN,
                FP_RATE,
                ALPHABET,
                TREE_LEN,
                estFactory,
                memFactory,
                prFactory,
                v, new CostFunctionMaxProb(), conf, NGRAMS);
    }


}
