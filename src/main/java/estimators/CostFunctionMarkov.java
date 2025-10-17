package estimators;

import PMIndex.NgramModel;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

public class CostFunctionMarkov extends AbstractCostFunction {

    private double[] PI;
    private double[][] T;
    private int SIGMA;

    private boolean loggedChainStats;

    private final Map<Long, Integer> markovSymbolIndex = new HashMap<>();

    private boolean singletonReady;
    private double[] singletonColSum;

    @Override
    public double costAtLevel(ImplicitTree<?> tree,
                              double[] probs,
                              long[] keySeq,
                              int Lp,
                              double bloomFp,
                              int stopLp) {

        ensureMarkovFromModel();
        if (PI == null || T == null) {
            return 0.0;
        }

        int[] keySeqIdx = mapToMarkovIndices(keySeq);
        if (keySeqIdx == null) {
            return 0.0;
        }

        final int width = tree.baseIntervalSize();
        final int r = keySeq.length;
        final int Ldesc = MathUtils.deepestVisitedLevel(width, r);

        final int bLp = width >> Lp;
        final int ellP = Math.min(keySeqIdx.length, bLp);
        if (ellP <= 0) {
            return 0.0;
        }

        final double betaLp = tree.getMembershipFpRate(Lp);

        double[] FmLp = Fm_uncond_markov(width, Lp, keySeqIdx, ellP, betaLp);
        int[] firstLp = firstOccurrencePositions(keySeqIdx, ellP);
        double HLp = H_from_Fm_by_positions(keySeqIdx, ellP, firstLp, FmLp);
        double F_prev_all = (FmLp.length == 0 ? 1.0 : FmLp[FmLp.length - 1]);

        double parentsVisited = 1 << Lp;
        double total = HLp * parentsVisited;

        if (Lp >= Ldesc) {
            return total;
        }

        for (int L = Lp + 1; L <= Ldesc; L++) {
            if (!MathUtils.childCanHost(width, L - 1, r)) {
                break;
            }

            double nodesAtL = 2.0 * parentsVisited * F_prev_all;
            if (nodesAtL <= 0.0) {
                break;
            }

            double betaPrev = tree.getMembershipFpRate(L - 1);
            double betaL = tree.getMembershipFpRate(L);

            double[] qCondL = qCondChildGivenParent_Markov(width, L, keySeqIdx, betaPrev, betaL);
            MathUtils.HF ns = MathUtils.HF_cond_from_q_pos_beta(width, L, keySeq, qCondL, betaL);

            total += ns.H * nodesAtL;

            parentsVisited = nodesAtL;
            F_prev_all = ns.F;
        }

        return total;
    }

    @Override
    public void setModel(NgramModel.Model bigramModel) {
        super.setModel(bigramModel);
        replaceMarkovSymbolMapping((bigramModel != null) ? bigramModel.symbolToIndex() : Collections.emptyMap());
        this.PI = null;
        this.T = null;
        this.singletonColSum = null;
        this.singletonReady = false;
        this.loggedChainStats = false;
    }

    private void ensureMarkovFromModel() {
        if (this.bigramModel == null) {
            return;
        }
        if (this.PI != null && this.T != null) {
            if (!loggedChainStats) {
                logChainStats(-1);
            }
            return;
        }

        this.SIGMA = bigramModel.sigma;
        this.PI = new double[SIGMA];
        for (int v = 0; v < SIGMA; v++) {
            this.PI[v] = bigramModel.pi(v);
        }

        this.T = new double[SIGMA][SIGMA];
        Map<Long, Map<Integer, Long>> M1 = bigramModel.ctxMaps.isEmpty()
                ? Collections.emptyMap()
                : bigramModel.ctxMaps.get(0);

        int fallbackRows = 0;
        for (int u = 0; u < SIGMA; u++) {
            long rowSum = 0L;
            Map<Integer, Long> row = (M1 == null) ? null : M1.get((long) u);
            if (row != null) {
                for (long c : row.values()) {
                    rowSum += c;
                }
            }

            if (rowSum > 0L) {
                if (row != null) {
                    for (Map.Entry<Integer, Long> e : row.entrySet()) {
                        int v = e.getKey();
                        this.T[u][v] = (double) e.getValue() / (double) rowSum;
                    }
                }
            } else {
                fallbackRows++;
                double Z = 0.0;
                for (int v = 0; v < SIGMA; v++) {
                    Z += this.PI[v];
                }
                if (Z <= 0.0) {
                    Z = 1.0;
                }
                for (int v = 0; v < SIGMA; v++) {
                    this.T[u][v] = this.PI[v] / Z;
                }
            }
        }

        logChainStats(fallbackRows);
    }

    private void logChainStats(int fallbackRows) {
        if (loggedChainStats || this.PI == null || this.T == null) {
            return;
        }
        double sumPi = 0.0;
        double maxRowErr = 0.0;
        for (int a = 0; a < SIGMA; a++) {
            sumPi += PI[a];
            double rowsum = 0.0;
            for (int b = 0; b < SIGMA; b++) {
                rowsum += T[a][b];
            }
            maxRowErr = Math.max(maxRowErr, Math.abs(rowsum - 1.0));
        }
        loggedChainStats = true;
    }

    private void replaceMarkovSymbolMapping(Map<Long, Integer> mapping) {
        markovSymbolIndex.clear();
        if (mapping != null) {
            markovSymbolIndex.putAll(mapping);
        }
    }

    private int resolveMarkovSymbol(long symbol) {
        Integer idx = markovSymbolIndex.get(symbol);
        if (idx != null) {
            return idx;
        }
        return -1;
    }

    private int[] mapToMarkovIndices(long[] keySeq) {
        if (markovSymbolIndex.isEmpty()) {
            return null;
        }
        int[] mapped = new int[keySeq.length];
        for (int i = 0; i < keySeq.length; i++) {
            int idx = resolveMarkovSymbol(keySeq[i]);
            if (idx < 0) {
                return null;
            }
            mapped[i] = idx;
        }
        return mapped;
    }

    private double noU_block_markov(int b, int[] U_syms) {
        if (b <= 0) {
            return 0.0;
        }
        boolean[] inU = new boolean[SIGMA];
        for (int s : U_syms) {
            if (s >= 0 && s < SIGMA) {
                inU[s] = true;
            }
        }

        double pi_notU = 0.0;
        double numer = 0.0;
        for (int a = 0; a < SIGMA; a++) {
            if (inU[a]) {
                continue;
            }
            double rowStay = 0.0;
            for (int bsym = 0; bsym < SIGMA; bsym++) {
                if (!inU[bsym]) {
                    rowStay += T[a][bsym];
                }
            }
            pi_notU += PI[a];
            numer += PI[a] * rowStay;
        }
        if (pi_notU <= 0.0) {
            return 0.0;
        }
        double T00 = numer / pi_notU;
        if (T00 < 0.0) {
            T00 = 0.0;
        } else if (T00 > 1.0) {
            T00 = 1.0;
        }
        return pi_notU * Math.pow(T00, Math.max(0, b - 1));
    }

    private double Fm_markov_beta(int width, int level, int[] keySeq, int[] firstIdx, int m, double betaL) {
        final int bL = width >> level;
        if (m <= 0) {
            return 1.0;
        }

        int[] sym = new int[m];
        for (int k = 0; k < m; k++) {
            sym[k] = keySeq[firstIdx[k]];
        }

        final int subsets = 1 << m;
        final double omb = 1.0 - betaL;

        double F = 0.0;
        F += 1.0;

        for (int mask = 1; mask < subsets; mask++) {
            int bits = Integer.bitCount(mask);
            int[] U = new int[bits];
            for (int i = 0, t = 0; i < m; i++) {
                if ((mask & (1 << i)) != 0) {
                    U[t++] = sym[i];
                }
            }

            double PnoU = noU_block_markov(bL, U);
            double coeff = ((bits & 1) == 0 ? 1.0 : -1.0) * Math.pow(omb, bits);
            F += coeff * PnoU;
        }
        return MathUtils.clamp01(F);
    }

    private double[] Fm_uncond_markov(int width, int level, int[] keySeq, int ell, double betaL) {
        int[] firstIdx = firstOccurrencePositions(keySeq, ell);
        int M = firstIdx.length;
        double[] Fm = new double[M];
        for (int m = 1; m <= M; m++) {
            Fm[m - 1] = Fm_markov_beta(width, level, keySeq, firstIdx, m, betaL);
        }
        return Fm;
    }

    private static double H_from_Fm_by_positions(int[] keySeq, int ell, int[] firstIdx, double[] Fm) {
        final int M = firstIdx.length;
        if (ell <= 0) {
            return 0.0;
        }
        if (M == 0) {
            return 1.0;
        }

        double H = 1.0;
        int m = 0;
        int nextK = firstIdx[0];
        for (int t = 0; t <= ell - 2; t++) {
            while (m < M && t >= nextK) {
                m++;
                nextK = (m < M ? firstIdx[m] : ell);
            }
            H += (m == 0) ? 1.0 : Fm[m - 1];
        }
        return H;
    }

    private static int[] firstOccurrencePositions(int[] keySeq, int ell) {
        HashSet<Integer> seen = new HashSet<>(ell * 2);
        ArrayList<Integer> first = new ArrayList<>();
        for (int pos = 0; pos < ell; pos++) {
            if (seen.add(keySeq[pos])) {
                first.add(pos);
            }
        }
        int[] out = new int[first.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = first.get(i);
        }
        return out;
    }

    private void ensureSingleton() {
        if (singletonReady) {
            return;
        }
        ensureMarkovFromModel();
        if (PI == null || T == null) {
            return;
        }
        singletonColSum = new double[SIGMA];
        for (int s = 0; s < SIGMA; s++) {
            double cs = 0.0;
            for (int a = 0; a < SIGMA; a++) {
                cs += PI[a] * T[a][s];
            }
            singletonColSum[s] = cs;
        }
        singletonReady = true;
    }

    private double hSingletonMarkovIdx(int b, int idx) {
        ensureSingleton();
        if (b <= 0 || singletonColSum == null || idx < 0 || idx >= SIGMA) {
            return 0.0;
        }

        double pi_s = PI[idx];
        double pi_not = 1.0 - pi_s;
        if (pi_not <= 0.0) {
            return 1.0;
        }

        double exitToS = singletonColSum[idx] - (pi_s * T[idx][idx]);
        double theta = 1.0 - (exitToS / pi_not);
        if (theta < 0.0) {
            theta = 0.0;
        } else if (theta > 1.0) {
            theta = 1.0;
        }

        double noS = pi_not * Math.pow(theta, Math.max(0, b - 1));
        double h = 1.0 - noS;
        return (h < 0.0) ? 0.0 : Math.min(h, 1.0);
    }

    private double[] qCondChildGivenParent_Markov(int width, int level, int[] keySeqIdx,
                                                  double betaPrev, double betaL) {
        final int bPrev = width >> (level - 1);
        final int bL = width >> level;
        final int ell = Math.min(keySeqIdx.length, bL);
        double[] q = new double[ell];

        for (int i = 0; i < ell; i++) {
            int idx = keySeqIdx[i];
            double hPrev = hSingletonMarkovIdx(bPrev, idx);
            double hL = hSingletonMarkovIdx(bL, idx);

            double numer = hL + betaL * (hPrev - hL) + betaL * betaPrev * (1.0 - hPrev);
            double denom = betaPrev + (1.0 - betaPrev) * hPrev;
            double qc = (denom > 0.0) ? (numer / denom) : 1.0;
            q[i] = MathUtils.clamp01(qc);
        }
        return q;
    }
}
