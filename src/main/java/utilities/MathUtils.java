package utilities;

import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.HashSet;

public final class  MathUtils {
    private MathUtils() {
        throw new AssertionError("MathUtils must not be instantiated");
    }

    // h(b) = 1 - (1-p)^block, with block = W / 2^L */
    public static double h_b(int width, int lvl, double prob){
        int blockLen = width >> lvl;                 // W / 2^L
        return 1.0 - Math.pow(1.0 - prob, blockLen);
    }
    public static double h_b(int width, int lvl, double prob, int count){
        int blockLen = width >> lvl;                 // W / 2^L
        blockLen -= count;
        return 1.0 - Math.pow(1.0 - prob, blockLen);
    }


    // q = h + beta*(1-h) — Bloom YES prob for one symbol at level L */
    public static double q_yes(double prob, int width, int lvl, double bfFp){
        double h = h_b(width, lvl, prob);
        return h + bfFp * (1.0 - h);
    }

    public static double q_child_given_parent_yes(double prob,
                                                  int width,
                                                  int parentLevel,
                                                  double beta) {
        // Parent at L, child at L+1
        int childLevel = parentLevel + 1;
        double hParent = h_b(width, parentLevel, prob);      // 1 - (1-p)^(W/2^L)
        double hChild  = h_b(width, childLevel,  prob);      // 1 - (1-p)^(W/2^(L+1))
        double qParent = q_yes(prob, width, parentLevel, beta); // β + (1-β)*hParent
        // Conditional child YES given "parent YES"
        return beta + (1.0 - beta) * (hChild / qParent);
    }

    /**
     * Expected probes per node (tail sum):
     * H = 1 + q1 + q1*q2 + ... + q1*...*q_{r-1}
     * Probe order = as given in probs[].
     */
    public static double expectedProbesPerNode(double[] probs,
                                               long[] keySeq,
                                               double bloomFp,
                                               int width,
                                               int level) {
        final int bL  = width >> level;
        final int ell = Math.min(keySeq.length, bL);
        if (ell <= 0) return 0.0;
        if (ell == 1) return 1.0;

        double total = 1.0;     // P(N >= 1) = 1
        double prod  = 1.0;     // product over q's for first occurrences seen so far
        HashSet<Long> seen = new HashSet<>();

        for (int i = 0; i < ell - 1; i++) {
            long sym = keySeq[i];
            if (seen.add(sym)) {
                // first time we see this key at this node
                prod *= q_yes(probs[i], width, level, bloomFp);
            }
            total += prod;      // duplicates add the current prod again
        }
        return total;
        }



    public static double fp_rate(double[] probs, long[] keySeq, int width, int level, double bloomFp){
        final int bL  = width >> level;
        final int ell = Math.min(keySeq.length, bL);
        double prod = 1.0;
        HashSet<Long> seen = new HashSet<>();
        for (int i = 0; i < ell; i++) {
            long sym = keySeq[i];
            if (seen.add(sym)) {
                prod *= q_yes(probs[i], width, level, bloomFp);
            }
        }
        return prod;
    }


    public static double fp_rate_from_q(double[] q_yes, int width, int level, double bloomFp){
        double f = 1.0;
        for (double p : q_yes) {
            f *= p;
        }
        return Math.max(0.0, Math.min(1.0, f));
    }
    /** Empirical CDF of a value x from an ascending array a: fraction of entries ≤ x. */
    private static double empiricalCdfLE(int[] a, int x) {
        if (a.length == 0) return 1.0;
        int lo = 0, hi = a.length; // search first index > x
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (a[mid] <= x) lo = mid + 1; else hi = mid;
        }
        return lo / (double) a.length;
    }


    // Can the children of level L still host a pattern of length r? */
    public static boolean childCanHost(int width, int level, int patternLength){
        return (width >> (level + 1)) >= patternLength;
    }

    // Alpha→Lp mapping for a single p
    public static int pruningLevel(ImplicitTree<?> tree, double conf, double prob){
        double bAlpha = Math.log(1.0 - conf) / Math.log(1.0 - prob); // > 0
        double val    = Math.log(tree.baseIntervalSize() / bAlpha) / Math.log(2.0);
        int L         = (int) Math.ceil(val);
        return Math.max(0, Math.min(L, tree.maxDepth() - 1)) + 1;
    }

}
