package estimators;

import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.HashSet;

public class CostFunctionMaxProb extends AbstractCostFunction {

    /** Per-pattern cache of level-local summaries that do not depend on the chosen base level Lp. */
    private static final class LpCostCache {
        final int width;
        final int r;
        final int maxLevel;
        final int Ldesc;
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
            || lpCache.probsRef != probs) {
            buildLpCache(tree, probs, keySeq);
        }

        if (lpCache == null) {
            return 0.0;
        }

        final int maxLevel = lpCache.maxLevel;
        if (Lp < 0 || Lp > maxLevel) {
            return 0.0;
        }

        final int Ldesc = lpCache.Ldesc;

        double nodesAtLp = (double) (1L << Lp);
        double total = lpCache.H_uncond[Lp] * nodesAtLp;

        if (!lpCache.childCanHostFrom[Lp] || Lp >= Ldesc) {
            return total;
        }

        int level = Lp + 1;
        if (level > maxLevel) {
            return total;
        }

        double nodesAtLevel = 2.0 * nodesAtLp * lpCache.F_uncond[Lp];
        if (nodesAtLevel <= 0.0) {
            return total;
        }

        total += lpCache.H_cond[level] * nodesAtLevel;

        while (level < Ldesc && lpCache.childCanHostFrom[level]) {
            double parents = nodesAtLevel;
            int next = level + 1;
            nodesAtLevel = 2.0 * parents * lpCache.F_cond[level];
            if (nodesAtLevel <= 0.0 || next > maxLevel) {
                break;
            }
            total += lpCache.H_cond[next] * nodesAtLevel;
            level = next;
        }

        return total;
    }

    private void buildLpCache(ImplicitTree<?> tree, double[] probs, long[] keySeq) {
        final int width = tree.baseIntervalSize();
        final int maxLevel = Math.max(0, tree.maxDepth() - 1);
        final int r = keySeq.length;
        final int Ldesc = Math.min(MathUtils.deepestVisitedLevel(width, r), maxLevel);

        final double[] beta = new double[maxLevel + 1];
        for (int level = 0; level <= maxLevel; level++) {
            beta[level] = tree.getMembershipFpRate(level);
        }

        final double[] H_uncond = new double[maxLevel + 1];
        final double[] F_uncond = new double[maxLevel + 1];
        final double[] H_cond = new double[maxLevel + 1];
        final double[] F_cond = new double[maxLevel + 1];
        final boolean[] childCanHostFrom = new boolean[maxLevel + 1];
        final int[][] firstIdx = new int[maxLevel + 1][];
        final int[][] segMult = new int[maxLevel + 1][];

        for (int level = 0; level <= maxLevel; level++) {
            final int bL = width >> level;
            final int ell = Math.min(r, bL);

            childCanHostFrom[level] = MathUtils.childCanHost(width, level, r);

            final int[] first = firstOccurrencePositions(keySeq, ell);
            final int[] mult = segmentMultiplicities(ell, first);
            firstIdx[level] = first;
            segMult[level] = mult;

            final double[] qUn = MathUtils.qUncondAtLevel(probs, width, level, beta[level]);
            HF hfUn = HF_from_q_preindexed(ell, first, mult, qUn);
            H_uncond[level] = hfUn.H;
            F_uncond[level] = hfUn.F;

            if (level >= 1) {
                final double[] qCond = MathUtils.qCondChildGivenParent(probs,
                                                                       width,
                                                                       level,
                                                                       beta[level - 1],
                                                                       beta[level]);
                HF hfCond = HF_from_q_preindexed(ell, first, mult, qCond);
                H_cond[level] = hfCond.H;
                F_cond[level] = hfCond.F;
            } else {
                H_cond[0] = 0.0;
                F_cond[0] = 1.0;
            }
        }

        this.lpCache = new LpCostCache(width,
                                       r,
                                       maxLevel,
                                       Ldesc,
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
