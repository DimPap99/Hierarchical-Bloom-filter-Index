import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import estimators.CostFunctionMaxProb;
import estimators.Estimator;
import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;
import search.*;
import utilities.AlphabetMapGen;
import utilities.CsvUtil;
import utilities.ExperimentRunResult;
import utilities.RunResult;

import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.IntStream;
import org.openjdk.jol.info.ClassLayout;
import org.openjdk.jol.info.GraphLayout;

public class HBIDatasetBenchmark {

    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/zipf_21_1.txt";
    private static final String QUERIES_FILE= "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/unique_substrings_zipf16_1_15.txt";

    private static final int WINDOW_LEN   = 1 << 21;//1 << 21;
    private static final int TREE_LEN     = 1 << 20;
    private static int ALPHABET     = 75;
    private static final double FP_RATE   = 0.001;
    private static final int RUNS         = 5;        // set to 0 for a dry run
    private static int NGRAMS = 4;
    private static String queryDir = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/";
    private static int NUMQUERIES = 135;
    public static void main(String[] args) throws IOException {

        List<String> queryFiles = new ArrayList<>();
        File dir = new File(queryDir + "zipf21_1");
        File[] directoryListing = dir.listFiles();
        List<List<?>> rows =  new ArrayList<>();
        List<?> header = List.of("nGram", "runAvgMs", "insertAvgMs", "avgQueryLength", "index");
        rows.add(header);

        ExperimentRunResult runResult;
        double ipmTotalMs = 0;
        double ipmTotalMsInsert = 0;
        double avgQueryLength = 0;
        if(directoryListing != null) {

            for (File file : directoryListing) {

                System.out.println("Running queries for " + file.toString());
                for (int n = 1; n <= 3; n++) {
                    double hbiTotalMs = 0;
                    ipmTotalMs = 0;
                    double hbiTotalMsInsert = 0;
                    ipmTotalMsInsert = 0;
                    double avgLp = 0;
                    NGRAMS = n;
                    System.out.println("N-gram: " + NGRAMS);
                    System.out.println("Window Size: " + WINDOW_LEN);
                    System.out.println("Tree Length: " + TREE_LEN);
                    ALPHABET = (int) Math.pow(ALPHABET, NGRAMS);

                    System.out.println("Alphabet: " + ALPHABET);
                    System.out.println("\n");
                    double avgAlpha = 0;
                    /* JIT warm-up so HotSpot reaches steady state */
                    for (int i = 0; i < 2; i++) {
                        HBI hbi = newHbi(0.999);

                        hbi.getStats = true;

                        Experiment.run(DATA_FILE, file.toString(), hbi, NGRAMS, false, false);

                        IPMIndexing ipm = new RegexIndex();
                        Experiment.run(DATA_FILE, file.toString(), ipm, 1, false, false);
                        avgLp = hbi.Lp.stream()
                                .mapToDouble(a -> a)
                                .sum() / hbi.Lp.size();
                        avgAlpha = hbi.alphas.stream()
                                .mapToDouble(a -> a)
                                .sum() / hbi.alphas.size();

                    }

                    ArrayList<Long> timings;
                    for (int i = 0; i < RUNS; i++) {

                        HBI hbi = newHbi(0.99);

                        hbi.getStats = true;
                        runResult = Experiment.run(DATA_FILE, file.toString(), hbi, NGRAMS, false, false);
                        hbiTotalMs += runResult.totalRunTimeMs();
                        hbiTotalMsInsert += runResult.totalInsertTimeMs();
                        IPMIndexing ipm = new RegexIndex();
                        runResult = Experiment.run(DATA_FILE, file.toString(), ipm, 1, false, false);
                        ipmTotalMs += runResult.totalRunTimeMs();
                        ipmTotalMsInsert += runResult.totalInsertTimeMs();
                        avgQueryLength = runResult.avgQuerySize();
                    }

                    if (RUNS > 0) {
                        System.out.printf("HBI avg (ms): %.3f%n", hbiTotalMs / RUNS);
                        System.out.printf("HBI Insert avg (ms): %.3f%n", hbiTotalMsInsert / RUNS);
                        System.out.println("Avg LP: " + avgLp);
                        System.out.println("Avg Alpha: " + avgAlpha);

                        System.out.printf("RegexIndex avg (ms): %.3f%n", ipmTotalMs / RUNS);
                        System.out.printf("RegexIndex Insert avg (ms): %.3f%n", ipmTotalMsInsert / RUNS);
                        System.out.println("\n");

                    }
                    rows.add(List.of(n, hbiTotalMs / RUNS,  hbiTotalMsInsert / RUNS, avgQueryLength, "hbi"));

                }
                rows.add(List.of(1, ipmTotalMs / RUNS,  ipmTotalMsInsert / RUNS, avgQueryLength, "regex"));

            }


        }
        CsvUtil.writeRows(Path.of("zipf21cf_set.csv"), rows);
    }

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
                v, new CostFunctionMaxProb(), conf);
    }


}
