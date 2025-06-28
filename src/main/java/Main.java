
import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;

import search.BlockSearch;

import estimators.Estimator;
import estimators.HashMapEstimator;

import membership.BloomFilter;
import membership.Membership;

import java.io.IOException;
import java.util.function.Supplier;

/**
 * Simple driver that benchmarks both the refactored HBI and the legacy
 * RegexIndex on the same data / query workload.
 */
public final class Main {

    /** Adjust to your file locations. */
    private static final String DATA_FILE   = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/zipf_text.txt";
    private static final String QUERIES_FILE= "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/substrings.txt";

    private static final int WINDOW_LEN   = 1 << 17;   // 131 072
    private static final int TREE_LEN     = 1 << 16;   // 65 536
    private static final int ALPHABET     = 81;
    private static final double FP_RATE   = 0.001;
    private static final int RUNS         = 50;        // set to 0 for a dry run

    public static void main(String[] args) throws IOException {

        double hbiTotalMs = 0;
        double ipmTotalMs = 0;

        /* ---------- 3× JIT warm-up so HotSpot reaches steady state ----- */
        for (int i = 0; i < 3; i++) {
            HBI hbi = newHbi();
            Experiment.run(DATA_FILE, QUERIES_FILE, hbi, /*verbose=*/false);

            IPMIndexing ipm = new RegexIndex();
            Experiment.run(DATA_FILE, QUERIES_FILE, ipm, false);
        }

        /* ---------- actual benchmark ----------------------------------- */
        for (int i = 0; i < RUNS; i++) {
            HBI hbi = newHbi();
            hbiTotalMs += Experiment.run(DATA_FILE, QUERIES_FILE, hbi, false);

            IPMIndexing ipm = new RegexIndex();
            ipmTotalMs += Experiment.run(DATA_FILE, QUERIES_FILE, ipm, false);
        }

        if (RUNS > 0) {
            System.out.printf("HBI avg (ms): %.3f%n", hbiTotalMs / RUNS);
            System.out.printf("RegexIndex avg (ms): %.3f%n", ipmTotalMs / RUNS);
        }
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

        return new HBI(new BlockSearch(),
                WINDOW_LEN,
                FP_RATE,
                ALPHABET,
                TREE_LEN,
                estFactory,
                memFactory);
    }
}
