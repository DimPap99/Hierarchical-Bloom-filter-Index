
import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;

import search.*;

import estimators.Estimator;
import estimators.HashMapEstimator;

import membership.BloomFilter;
import membership.Membership;
import utilities.AlphabetMapGen;

import java.io.IOException;
import java.util.List;
import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * Simple driver that benchmarks both the refactored HBI and the legacy
 * RegexIndex on the same data / query workload.
 */
public final class Main {

    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/zipf_text.txt";
    private static final String QUERIES_FILE= "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries1.txt";

    private static final int WINDOW_LEN   = 1 << 17;   // 131 072
    private static final int TREE_LEN     = 1 << 16;   // 65 536
    private static int ALPHABET     = 81;
    private static final double FP_RATE   = 0.001;
    private static final int RUNS         = 100;        // set to 0 for a dry run
    private static final int NGRAMS = 1;
    public static void main(String[] args) throws IOException {

        double hbiTotalMs = 0;
        double ipmTotalMs = 0;
        List<Character> letters = IntStream.rangeClosed(48,122)
                .mapToObj(c -> (char)c)
                .toList();

        AlphabetMapGen<Character> gen = new AlphabetMapGen<>(NGRAMS, letters);
        ALPHABET = gen.alphabetMap.size();
        /* ---------- 3× JIT warm-up so HotSpot reaches steady state ----- */
        for (int i = 0; i < 3; i++) {
            HBI hbi = newHbi();
            hbi.alphabetMap = gen.alphabetMap;
            Experiment.run(DATA_FILE, QUERIES_FILE, hbi, NGRAMS, false);

            IPMIndexing ipm = new RegexIndex();
            Experiment.run(DATA_FILE, QUERIES_FILE, ipm, 1, false);
        }

        /* ---------- actual benchmark ----------------------------------- */
//        for (int i = 0; i < RUNS; i++) {
//            HBI hbi = newHbi();
//            hbi.alphabetMap = gen.alphabetMap;
//
//            hbiTotalMs += Experiment.run(DATA_FILE, QUERIES_FILE, hbi, NGRAMS, false);
//
//            IPMIndexing ipm = new RegexIndex();
//            ipmTotalMs += Experiment.run(DATA_FILE, QUERIES_FILE, ipm, NGRAMS, false);
//        }
//
//        if (RUNS > 0) {
//            System.out.printf("HBI avg (ms): %.3f%n", hbiTotalMs / RUNS);
//            System.out.printf("RegexIndex avg (ms): %.3f%n", ipmTotalMs / RUNS);
//        }
    }

    /* ------------------------------------------------------------------ */
    /*  Helper that builds a *fresh* HBI wired to our suppliers each time */
    /* ------------------------------------------------------------------ */
    private static HBI newHbi() {
        /* Estimator supplier – one instance per ImplicitTree */
        Supplier<Estimator> estFactory =
                () -> new HashMapEstimator(TREE_LEN);

        /* Membership supplier – one instance *per tree level* */
        Supplier<Membership> memFactory =
                () -> new BloomFilter();
        Supplier<PruningPlan> prFactory =
                () -> new MostFreqPruning();
        Verifier v = new VerifierLinearLeafProbe();
        return new HBI(new BlockSearch(),
                WINDOW_LEN,
                FP_RATE,
                ALPHABET,
                TREE_LEN,
                estFactory,
                memFactory,
                prFactory,
                v);
    }
}
