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
        double lp = Math.min(L, tree.maxDepth() - 1);
        if(lp > 0) lp = lp+1;
        return (int) Math.max(0, lp);
    }

    public static double clamp01(double value) {
        if (value <= 0.0) {
            return 0.0;
        }
        if (value >= 1.0) {
            return 1.0;
        }
        return value;
    }

    public static int deepestVisitedLevel(int width, int patternLength) {
        int ratio = width / Math.max(patternLength, 1);
        if (ratio <= 0) {
            return 0;
        }
        return 31 - Integer.numberOfLeadingZeros(ratio);
    }

    public static double[] qUncondAtLevel(double[] probs, int width, int level, double beta) {
        int bL   = width >> level;
        int ell  = Math.min(probs.length, bL);
        double[] q = new double[ell];
        double omb = 1.0 - beta;
        for (int i = 0; i < ell; i++) {
            double h = h_b(width, level, probs[i]);
            q[i] = clamp01(beta + omb * h);
        }
        return q;
    }

    public static double[] qCondChildGivenParent(double[] probs,
                                                  int width,
                                                  int level,
                                                  double betaPrev,
                                                  double betaL) {
        int bL   = width >> level;
        int ell  = Math.min(probs.length, bL);
        double[] q = new double[ell];
        for (int i = 0; i < ell; i++) {
            double p     = probs[i];
            double hPrev = h_b(width, level - 1, p);
            double hL    = h_b(width, level,     p);

            double numer = hL + betaL * (hPrev - hL) + betaL * betaPrev * (1.0 - hPrev);
            double denom = betaPrev + (1.0 - betaPrev) * hPrev;

            double qc = (denom > 0.0) ? (numer / denom) : 1.0;
            q[i] = clamp01(qc);
        }
        return q;
    }

    public static final class HF {
        public final double H;
        public final double F;

        public HF(double H, double F) {
            this.H = H;
            this.F = F;
        }
    }

    public static HF HF_uncond_pos_beta(int width, int level,
                                        long[] keySeq, double[] probs, double betaL) {
        return HF_uncond_pos_beta(width, level, keySeq, probs, betaL, Integer.MAX_VALUE);
    }

    public static HF HF_uncond_pos_beta(int width, int level,
                                        long[] keySeq, double[] probs, double betaL, int ieMaxOrder) {
        final int bL  = width >> level;
        final int ell = Math.min(keySeq.length, bL);
        if (ell <= 0) return new HF(0.0, 1.0);

        HashSet<Long> seen = new HashSet<>(ell * 2);
        ArrayList<Integer> first = new ArrayList<>();
        for (int pos = 0; pos < ell; pos++) if (seen.add(keySeq[pos])) first.add(pos);
        final int M = first.size();
        if (M == 0) return new HF(1.0, 1.0);

        int orderLimit = (ieMaxOrder < 0) ? Integer.MAX_VALUE : ieMaxOrder;

        double[] p = new double[M];
        for (int m = 0; m < M; m++) {
            p[m] = clamp01(probs[first.get(m)]);
        }

        double[] Fm = new double[M];
        for (int m = 1; m <= M; m++) {
            Fm[m - 1] = IE_prefix_collapsed_beta(p, m, bL, betaL, orderLimit);
        }

        double H = 1.0;
        for (int m = 0; m < M; m++) {
            int next = (m + 1 < M) ? first.get(m + 1) : (ell - 1);
            int mult = next - first.get(m);
            H += mult * Fm[m];
        }
        return new HF(H, Fm[M - 1]);
    }

    public static HF HF_cond_from_q_pos_beta(int width, int level,
                                             long[] keySeq, double[] qCond, double betaL) {
        return HF_cond_from_q_pos_beta(width, level, keySeq, qCond, betaL, Integer.MAX_VALUE);
    }

    public static HF HF_cond_from_q_pos_beta(int width, int level,
                                             long[] keySeq, double[] qCond, double betaL, int ieMaxOrder) {
        final int bL  = width >> level;
        final int ell = Math.min(keySeq.length, bL);
        if (ell <= 0) return new HF(0.0, 1.0);

        HashSet<Long> seen = new HashSet<>(ell * 2);
        ArrayList<Integer> first = new ArrayList<>();
        for (int pos = 0; pos < ell; pos++) if (seen.add(keySeq[pos])) first.add(pos);
        final int M = first.size();
        if (M == 0) return new HF(1.0, 1.0);

        final double omb = 1.0 - betaL;
        int orderLimit = (ieMaxOrder < 0) ? Integer.MAX_VALUE : ieMaxOrder;

        double[] pEff = new double[M];
        for (int m = 0; m < M; m++) {
            double q = clamp01(qCond[first.get(m)]);
            double g = (omb > 0.0) ? clamp01((q - betaL) / omb) : 1.0;
            double p = 1.0 - Math.pow(1.0 - g, 1.0 / Math.max(1, bL));
            pEff[m] = clamp01(p);
        }

        double[] Fm = new double[M];
        for (int m = 1; m <= M; m++) {
            Fm[m - 1] = IE_prefix_collapsed_beta(pEff, m, bL, betaL, orderLimit);
        }

        double H = 1.0;
        for (int m = 0; m < M; m++) {
            int next = (m + 1 < M) ? first.get(m + 1) : (ell - 1);
            int mult = next - first.get(m);
            H += mult * Fm[m];
        }
        return new HF(H, Fm[M - 1]);
    }

    private static double IE_prefix_collapsed_beta(double[] pFirstOccur, int m, int bL, double betaL) {
        return IE_prefix_collapsed_beta(pFirstOccur, m, bL, betaL, Integer.MAX_VALUE);
    }

    private static double IE_prefix_collapsed_beta(double[] pFirstOccur,
                                                   int m,
                                                   int bL,
                                                   double betaL,
                                                   int maxOrder) {
        final int M = Math.max(0, Math.min(m, pFirstOccur.length));
        if (M == 0) {
            return 1.0;
        }

        final int t = Math.max(0, Math.min(maxOrder, M));
        final double omb = 1.0 - betaL;

        double F = 1.0; // k = 0 term
        for (int k = 1; k <= t; k++) {
            double sumOverSubsetsK = sumCombPowers_k(pFirstOccur, M, k, 0, 0.0, bL);
            double coeff = (((k & 1) == 0) ? 1.0 : -1.0) * Math.pow(omb, k);
            F += coeff * sumOverSubsetsK;
        }
        return clamp01(F);
    }

    private static double sumCombPowers_k(double[] p,
                                          int M,
                                          int k,
                                          int start,
                                          double sumP,
                                          int bL) {
        if (k == 0) {
            double base = clamp01(1.0 - sumP);
            return Math.pow(base, bL);
        }
        double total = 0.0;
        for (int i = start; i <= M - k; i++) {
            total += sumCombPowers_k(p, M, k - 1, i + 1, sumP + p[i], bL);
        }
        return total;
    }

}
