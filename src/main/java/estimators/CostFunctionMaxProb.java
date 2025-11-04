package estimators;

import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

public class CostFunctionMaxProb extends AbstractCostFunction {

    /** Per-pattern cache of level-local summaries that do not depend on the chosen base level Lp. */
    private static final class LpCostCache {
        final int width;
        final int r;
        final int maxLevel;
        final int Ldesc;
        final int minLevel;
        final int upperLevel;
        final long[] keySeqRef;
        final double[] probsRef;
        final double[] beta;
        final double[] H_uncond;
        final double[] F_uncond;
        final double[] H_cond;
        final double[] F_cond;
        final boolean[] childCanHostFrom;
        final int[][] firstIdx;
        final int[][] segMult;

        LpCostCache(int width,
                    int r,
                    int maxLevel,
                    int Ldesc,
                    int minLevel,
                    int upperLevel,
                    long[] keySeqRef,
                    double[] probsRef,
                    double[] beta,
                    double[] H_uncond,
                    double[] F_uncond,
                    double[] H_cond,
                    double[] F_cond,
                    boolean[] childCanHostFrom,
                    int[][] firstIdx,
                    int[][] segMult) {
            this.width = width;
            this.r = r;
            this.maxLevel = maxLevel;
            this.Ldesc = Ldesc;
            this.minLevel = minLevel;
            this.upperLevel = upperLevel;
            this.keySeqRef = keySeqRef;
            this.probsRef = probsRef;
            this.beta = beta;
            this.H_uncond = H_uncond;
            this.F_uncond = F_uncond;
            this.H_cond = H_cond;
            this.F_cond = F_cond;
            this.childCanHostFrom = childCanHostFrom;
            this.firstIdx = firstIdx;
            this.segMult = segMult;
        }

        boolean covers(int level) {
            return level >= minLevel && level <= upperLevel;
        }

        int offset(int level) {
            return level - minLevel;
        }

        double HUncond(int level) {
            return this.H_uncond[offset(level)];
        }

        double FUncond(int level) {
            return this.F_uncond[offset(level)];
        }

        double HCond(int level) {
            return this.H_cond[offset(level)];
        }

        double FCond(int level) {
            return this.F_cond[offset(level)];
        }

        boolean childCanHostFrom(int level) {
            return this.childCanHostFrom[offset(level)];
        }
    }

    private static final class HF {
        final double H;
        final double F;

        HF(double H, double F) {
            this.H = H;
            this.F = F;
        }
    }

    private LpCostCache lpCache = null;

    public CostFunctionMaxProb() {
    }

    @Override
    public double costAtLevel(ImplicitTree<?> tree,
                              double[] probs,
                              long[] keySeq,
                              int Lp,
                              double bloomFp,
                              int stopLp) {
        if (probs == null || keySeq == null || keySeq.length == 0) {
            return 0.0;
        }

        if (lpCache == null
            || lpCache.width != tree.baseIntervalSize()
            || lpCache.r != keySeq.length
            || lpCache.keySeqRef != keySeq
            || lpCache.probsRef != probs
            || !lpCache.covers(Lp)) {
            buildLpCache(tree, probs, keySeq, Lp);
        }

        if (lpCache == null) {
            return 0.0;
        }

        final int maxLevel = lpCache.maxLevel;
        if (Lp < lpCache.minLevel || Lp > maxLevel || !lpCache.covers(Lp)) {
            return 0.0;
        }

        final int Ldesc = lpCache.Ldesc;

        double nodesAtLp = (double) (1L << Lp);
        double total = lpCache.HUncond(Lp) * nodesAtLp;

        if (!lpCache.childCanHostFrom(Lp) || Lp >= Ldesc) {
            return total;
        }

        int level = Lp + 1;
        if (level > lpCache.upperLevel) {
            return total;
        }

        double nodesAtLevel = 2.0 * nodesAtLp * lpCache.FUncond(Lp);
        if (nodesAtLevel <= 0.0) {
            return total;
        }

        total += lpCache.HCond(level) * nodesAtLevel;

        while (level < Ldesc && lpCache.childCanHostFrom(level)) {
            double parents = nodesAtLevel;
            int next = level + 1;
            nodesAtLevel = 2.0 * parents * lpCache.FCond(level);
            if (nodesAtLevel <= 0.0 || next > lpCache.upperLevel) {
                break;
            }
            total += lpCache.HCond(next) * nodesAtLevel;
            level = next;
        }

        return total;
    }

    private void buildLpCache(ImplicitTree<?> tree, double[] probs, long[] keySeq, int requestedLp) {
        final int width = tree.baseIntervalSize();
        final int maxLevel = Math.max(0, tree.maxDepth() - 1);
        final int r = keySeq.length;

        double minProb = Arrays.stream(probs).min().orElse(0.0);
        int minCandidate = 0;
        int maxCandidate = MathUtils.pruningLevel(tree, 0.05, minProb);

        // Convert to level indices and clamp to valid range
        minCandidate = clampLevel(minCandidate, maxLevel);
        maxCandidate = clampLevel(maxCandidate, maxLevel);

        final int Ldesc = Math.min(MathUtils.deepestVisitedLevel(width, r), maxLevel);
        maxCandidate = Math.min(maxCandidate, Ldesc);
        if (maxCandidate < minCandidate) {
            maxCandidate = minCandidate;
        }

        int desired = clampLevel(requestedLp, Ldesc);
        int minLevel = Math.min(minCandidate, desired);
        int upperLevel = Math.max(Math.max(minCandidate, maxCandidate), desired);
        upperLevel = Math.min(upperLevel, Ldesc);
        if (upperLevel < minLevel) {
            upperLevel = minLevel;
        }

        final double[] beta = new double[upperLevel + 1];
        for (int level = 0; level <= upperLevel; level++) {
            beta[level] = tree.getMembershipFpRate(level);
        }

        final int range = upperLevel - minLevel + 1;
        final double[] H_uncond = new double[range];
        final double[] F_uncond = new double[range];
        final double[] H_cond = new double[range];
        final double[] F_cond = new double[range];
        final boolean[] childCanHostFrom = new boolean[range];
        final int[][] firstIdx = new int[range][];
        final int[][] segMult = new int[range][];

        int startLevel = minLevel;
        for (int level = startLevel; level <= upperLevel; level++) {
            int idx = level - startLevel;
            final int bL = width >> level;
            final int ell = Math.min(r, bL);

            childCanHostFrom[idx] = MathUtils.childCanHost(width, level, r);

            final int[] first = firstOccurrencePositions(keySeq, ell);
            final int[] mult = segmentMultiplicities(ell, first);
            firstIdx[idx] = first;
            segMult[idx] = mult;

            final double[] qUn = MathUtils.qUncondAtLevel(probs, width, level, beta[level]);
            HF hfUn = HF_from_q_preindexed(ell, first, mult, qUn);
            H_uncond[idx] = hfUn.H;
            F_uncond[idx] = hfUn.F;

            if (level >= 1) {
                final double[] qCond = MathUtils.qCondChildGivenParent(probs,
                                                                       width,
                                                                       level,
                                                                       beta[level - 1],
                                                                       beta[level]);
                HF hfCond = HF_from_q_preindexed(ell, first, mult, qCond);
                H_cond[idx] = hfCond.H;
                F_cond[idx] = hfCond.F;
            } else {
                H_cond[idx] = 0.0;
                F_cond[idx] = 1.0;
            }
        }

        this.lpCache = new LpCostCache(width,
                                       r,
                                       maxLevel,
                                       Ldesc,
                                       minLevel,
                                       upperLevel,
                                       keySeq,
                                       probs,
                                       beta,
                                       H_uncond,
                                       F_uncond,
                                       H_cond,
                                       F_cond,
                                       childCanHostFrom,
                                       firstIdx,
                                       segMult);
    }

    private static int clampLevel(int level, int maxLevel) {
        if (level < 0) {
            return 0;
        }
        if (level > maxLevel) {
            return maxLevel;
        }
        return level;
    }

    private static HF HF_from_q_preindexed(int ell,
                                           int[] firstIdx,
                                           int[] segMult,
                                           double[] q) {
        if (ell <= 0) {
            return new HF(0.0, 1.0);
        }
        final int M = firstIdx.length;
        if (M == 0) {
            return new HF(1.0, 1.0);
        }

        double H = 1.0;
        double prod = 1.0;
        for (int m = 0; m < M; m++) {
            int pos = firstIdx[m];
            double qi = MathUtils.clamp01(q[pos]);
            prod *= qi;
            H += (double) segMult[m] * prod;
        }
        return new HF(H, prod);
    }

    private static int[] firstOccurrencePositions(long[] keySeq, int ell) {
        if (ell <= 0) {
            return new int[0];
        }
        HashSet<Long> seen = new HashSet<>(ell * 2);
        ArrayList<Integer> first = new ArrayList<>();
        for (int i = 0; i < ell; i++) {
            if (seen.add(keySeq[i])) {
                first.add(i);
            }
        }
        int[] out = new int[first.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = first.get(i);
        }
        return out;
    }

    private static int[] segmentMultiplicities(int ell, int[] firstIdx) {
        final int M = firstIdx.length;
        int[] mult = new int[M];
        if (ell <= 1 || M == 0) {
            return mult;
        }
        for (int m = 0; m < M; m++) {
            int start = firstIdx[m];
            int next = (m + 1 < M) ? firstIdx[m + 1] : (ell - 1);
            int count = next - start;
            mult[m] = Math.max(0, count);
        }
        return mult;
    }
}
