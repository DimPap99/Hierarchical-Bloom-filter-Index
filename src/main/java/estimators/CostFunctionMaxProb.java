package estimators;

import org.apache.commons.math3.util.Pair;
import search.Pattern;
import tree.ImplicitTree;
import utilities.MathUtils;

public class CostFunctionMaxProb implements CostFunction {

    //cost of constant operations. Typically measured in nano seconds.
    double bloomProbeCost = 97;
    double leafSearchCost = 26;
    public double alpha;
    public CostFunctionMaxProb() {

    }




    @Override
    public double getAlpha() {
        return this.alpha;
    }



    /**
     * Branching-aware Cost–3:
     *   C_total = C_hor(Lp) + C_vert(Lp..stop) + C_leaf
     * where:
     *   C_hor  = bc * H(Lp) * 2^Lp
     *   C_vert = sum_{L=Lp+1..L_stop} bc * H(L) * N_L
     *            with N_{Lp+1} = 2 * (1 + F(Lp)*(2^Lp - 1)),
     *                 N_{L+1}  = 2 * (1 + F(L)*(N_L - 1))
     *   C_leaf = lc * ( r + (followed_final - 1) * (r - 1) )
     *
     * F(L) = ∏ q_i(L)  (all-pass),  H(L) = tail-sum expected probes.
     */


    @Override
    public int minCostLp(ImplicitTree tree, double bfFalsePosRate, double confInit, Pattern p, double bfCost, double leafCost) {   // measured lc (ns)
        if (bfCost > 0) this.bloomProbeCost = bfCost;
        if (leafCost > 0) this.leafSearchCost = leafCost;

        // Per-symbol probabilities for the pattern from the estimator
        final double[] probs = tree.estimator.estimateALl(p);
        final int width      = tree.baseIntervalSize();
        final int maxDepth   = tree.maxDepth();
        final int r          = p.nGramToInt.length;

        // Deepest level we might touch (children still host the full pattern)
        int LStopMax = 0;
        for (int L = 0; L < maxDepth; L++) {
            if (MathUtils.childCanHost(width, L, r)) LStopMax = L;
            else break;
        }

        double bestCost = Double.POSITIVE_INFINITY;
        int    bestLp   = 0;

        // Evaluate cost for each candidate Lp (no α sweep needed)
        for (int Lp = 0; Lp <= LStopMax; Lp++) {
            double cost = costAtLevel(tree, probs, Lp, bfFalsePosRate);
            if (cost < bestCost) {
                bestCost = cost;
                bestLp   = Lp;
            }
        }

//            // Record an α consistent with bestLp (purely for reporting)
//            double Fbest = MathUtils.allPass(probs, width, bestLp, bloomFp);
//            double Fprev = (bestLp > 0) ? MathUtils.allPass(probs, width, bestLp - 1, bloomFp) : 1.0;
//            this.alpha   = Math.min(0.99, Math.max(Fbest + 1e-9, 0.5 * (Fbest + Fprev)));

        return bestLp;
    }


    /* ------------ cost pieces ---------------- */

    private double costAtLevel(ImplicitTree<?> tree,
                               double[] probs,
                               int Lp,
                               double bloomFp) {
        final int width    = tree.baseIntervalSize();
        final int maxDepth = tree.maxDepth();
        final int r        = probs.length;

        // Horizontal at Lp
        double H_lp  = MathUtils.expectedProbesPerNode(probs, bloomFp, width, Lp);
        double C_hor = bloomProbeCost * H_lp * (1 << Lp);

        // Vertical from Lp+1 down (Lp work already paid in C_hor)
        Pair<Double, Double> vert = verticalCostBranching(probs, width, maxDepth, r, Lp, bloomFp);
        double C_vert   = vert.getFirst();
        double followed = vert.getSecond();

        // Leaf (worst-case model, using 'followed' from the deepest level reached)
        double C_leaf = leafSearchCost * (r + (followed - 1.0) * (r - 1.0));

        return C_hor + C_vert + C_leaf;
    }

    /**
     * Returns (verticalCost, followedFinal).
     * followedFinal is the expected "1 + F*(N-1)" at the last processed level
     * and is used in the caller's leaf term (to mirror the provided Python).
     */
    private Pair<Double, Double> verticalCostBranching(double[] probs,
                                                       int width,
                                                       int maxDepth,
                                                       int r,
                                                       int Lp,
                                                       double bloomFp) {

//        double minProb =         Arrays.stream(probs).min().getAsDouble();
;
        double F_lp    = MathUtils.fp_rate(probs, width, Lp, bloomFp);
        if (!MathUtils.childCanHost(width, Lp, r)) {
            // No deeper levels to visit; followed = 1 + F(Lp)*(2^Lp - 1)
//            double F_lp    = MathUtils.fp_rate(probs, width, Lp, bloomFp);
            double followed = 1.0 + F_lp * ((1 << Lp) - 1.0);
            return Pair.create(0.0, followed);
        }

        // Start from children of Lp
        double followed = 1.0 + F_lp * ((1 << Lp) - 1.0);
        double N        = 2.0 * followed;               // nodes processed at Lp+1

        // Compute deepest level we could ever reach under the stop rule
        int L_stop = Lp + 1;
        while (L_stop < maxDepth - 1 && MathUtils.childCanHost(width, L_stop, r)) {
            L_stop++;
        }

        double cost = 0.0;
        for (int L = Lp + 1; L <= L_stop; L++) {
            // Pay per-node probes at this level
            double H_L = MathUtils.expectedProbesPerNode(probs, bloomFp, width, L);
            cost += bloomProbeCost * H_L * N;

            // Update branching if we can go deeper
            if (MathUtils.childCanHost(width, L, r)) {
                double F_L = MathUtils.fp_rate(probs, width, L, bloomFp);
                followed   = 1.0 + F_L * (N - 1.0);
                N          = 2.0 * followed;
            } else {
                break;
            }
        }

        return Pair.create(cost, followed);
    }





}
