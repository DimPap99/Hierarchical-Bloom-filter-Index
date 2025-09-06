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
import utilities.ExperimentRunResult;
import utilities.RunResult;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.IntStream;

public class HBIDatasetBenchmark {


    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/zipf_16_1.txt";
    private static final String QUERIES_FILE= "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/unique_substrings_zipf16.txt";

    private static final int WINDOW_LEN   = 1 << 16;//1 << 21;
    private static final int TREE_LEN     = 1 << 16;
    private static int ALPHABET     = 75;
    private static final double FP_RATE   = 0.001;
    private static final int RUNS         = 150;        // set to 0 for a dry run
    private static int NGRAMS = 4;

    private static int NUMQUERIES = 135;
    public static void main(String[] args) throws IOException {


        List<Character> letters = IntStream.rangeClosed(48,122)
                .mapToObj(c -> (char)c)
                .toList();

        ArrayList<RunResult> results = new ArrayList<>();

        for (double alpha = 0.1; alpha <= 1 + 0.1; alpha += 0.1+ 1e-9) {
            if(alpha >=1) alpha = 0.99;            double hbiTotalMs = 0;
            double ipmTotalMs = 0;
            double hbiTotalMsInsert = 0;
            double ipmTotalMsInsert = 0;
            double avgLp = 0;
            NGRAMS = 1;
            System.out.println("N-gram: " + NGRAMS);
            System.out.println("Window Size: " + WINDOW_LEN);
            System.out.println("Tree Length: " + TREE_LEN);
            AlphabetMapGen<Character> gen = new AlphabetMapGen<>(NGRAMS, letters);
            ALPHABET = gen.alphabetMap.size();
            System.out.println("Alphabet: " + ALPHABET);
            System.out.println("ALPHA: " + alpha);
            int maxLvl;
            double avgAlpha =0;
            /* JIT warm-up so HotSpot reaches steady state */
            for (int i = 0; i < 5; i++) {
                HBI hbi = newHbi(alpha);
                hbi.alphabetMap = gen.alphabetMap;
                hbi.getStats = true;
                //if(i ==0)hbi.getAvgTimes(new Pattern("wZE2bl[cuO", 1));
                Experiment.run(DATA_FILE, QUERIES_FILE, hbi, NGRAMS, false, false);

                IPMIndexing ipm = new RegexIndex();
                Experiment.run(DATA_FILE, QUERIES_FILE, ipm, 1, false, false);
                avgLp = hbi.Lp.stream()
                        .mapToDouble(a -> a)
                        .sum()/hbi.Lp.size();
                avgAlpha = hbi.alphas.stream()
                        .mapToDouble(a -> a)
                        .sum()/hbi.alphas.size();
            }

            ArrayList<Long> timings;
            ExperimentRunResult runResult;
            int probes = 0;
            for (int i = 0; i < RUNS; i++) {
                HBI hbi = newHbi(alpha);
                hbi.alphabetMap = gen.alphabetMap;
                hbi.getStats = true;
                runResult = Experiment.run(DATA_FILE, QUERIES_FILE, hbi, NGRAMS, false, false);
                hbiTotalMs += runResult.totalRunTimeMs();
                hbiTotalMsInsert += runResult.totalInsertTimeMs();
                IPMIndexing ipm = new RegexIndex();
                runResult = Experiment.run(DATA_FILE, QUERIES_FILE, ipm, 1, false, false);
                ipmTotalMs += runResult.totalRunTimeMs();
                ipmTotalMsInsert += runResult.totalInsertTimeMs();
                probes+= hbi.getAllprobes();
//                System.out.println("Run with confidence " + conf);
//                System.out.printf("HBI avg (ms): %.3f%n", hbiTotalMs / RUNS);
//                avgLp = hbi.Lp.stream()
//                        .mapToDouble(a -> a)
//                        .sum()/hbi.Lp.size();
//                System.out.println("Avg Lp for this run " + avgLp);
//                System.out.println("   ");
                //System.out.printf("HBI Insert avg (ms): %.3f%n", hbiTotalMsInsert / RUNS);

            }

            if (RUNS > 0) {
                System.out.printf("HBI avg (ms): %.3f%n", hbiTotalMs / RUNS);
                System.out.printf("HBI Insert avg (ms): %.3f%n", hbiTotalMsInsert / RUNS);
                System.out.println("Avg LP: " + avgLp);
                System.out.println("Avg Alpha: " + avgAlpha);

                System.out.printf("RegexIndex avg (ms): %.3f%n", ipmTotalMs / RUNS);
                System.out.printf("RegexIndex Insert avg (ms): %.3f%n", ipmTotalMsInsert / RUNS);
                System.out.println("\n\n\n");
                results.add(new RunResult(alpha, hbiTotalMs / RUNS, avgLp, probes/RUNS));
            }
            if(alpha == 0.99)break;
        }

        System.out.println("\n\nAlpha,AvgMs,Lp");
        for(RunResult r: results){
            r.print();
        }

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
