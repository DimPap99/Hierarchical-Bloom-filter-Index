package utilities;

import tree.ImplicitTree;

import java.util.ArrayList;

public final class  MathUtils {
    private MathUtils() {
        throw new AssertionError("MathUtils must not be instantiated");
    }

    /** h(b) = 1 - (1-p)^block, with block = W / 2^L */
    public static double h_b(int width, int lvl, double prob){
        int blockLen = width >> lvl;                 // W / 2^L
        return 1.0 - Math.pow(1.0 - prob, blockLen);
    }

    /** q = h + beta*(1-h) — Bloom YES prob for one symbol at level L */
    public static double q_yes(double prob, int width, int lvl, double bfFp){
        double h = h_b(width, lvl, prob);
        return h + bfFp * (1.0 - h);
    }

    /**
     * Expected probes per node (tail sum):
     * H = 1 + q1 + q1*q2 + ... + q1*...*q_{r-1}
     * Probe order = as given in probs[].
     */
    public static double expectedProbesPerNode(double[] probs,
                                               double bloomFp,
                                               int width,
                                               int level) {
        final int r = probs.length;
        if (r <= 1) return 1.0;

        double[] qs = new double[r];
        for (int i = 0; i < r; i++) {
            qs[i] = q_yes(probs[i], width, level, bloomFp);
        }

        double total = 1.0;      // P(N>=1) = 1
        double prod  = 1.0;
        for (int i = 0; i < r - 1; i++) {  // add q1, q1*q2, ..., q1*...*q_{r-1}
            prod *= qs[i];
            total += prod;
        }
        //  cap by r. The total expected probes cant exceed the length of the pattern
        return Math.min(total, r);
    }

    public static double fp_rate(double[] probs, int width, int level, double bloomFp){
        double f = 1.0;
        for (double p : probs) {
            f *= q_yes(p, width, level, bloomFp);
        }
        return Math.max(0.0, Math.min(1.0, f));
    }

    public static double fp_rate_from_q(double[] q_yes, int width, int level, double bloomFp){
        double f = 1.0;
        for (double p : q_yes) {
            f *= p;
        }
        return Math.max(0.0, Math.min(1.0, f));
    }


    // Can the children of level L still host a pattern of length r? */
    public static boolean childCanHost(int width, int level, int patternLength){
        return (width >> (level + 1)) >= patternLength;
    }

    // Alpha→Lp mapping for a single p
    public static int pruningLevel(ImplicitTree<?> tree, double conf, double prob){
        double bAlpha = Math.log(1.0 - conf) / Math.log(1.0 - prob);   // b_α
        double log2   = Math.log(tree.baseIntervalSize() / bAlpha) / Math.log(2.0);
        int rawLp     = (int) Math.floor(log2) + 1;
        return Math.max(0, Math.min(rawLp, tree.maxDepth() - 1));
    }
}
