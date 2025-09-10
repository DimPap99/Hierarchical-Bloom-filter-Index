package estimators;

import membership.Membership;
import org.apache.commons.math3.util.Pair;
import search.Pattern;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.Arrays;

import static utilities.MathUtils.q_yes;

public class CostFunctionMaxProb implements CostFunction {

    //cost of constant operations. Typically measured in nano seconds.
    double bloomProbeCost = 97;
    double leafSearchCost = 26;

    double predictedBloomProbeCost = 0;

    public double alpha;
    public CostFunctionMaxProb() {

    }




    @Override
    public double getAlpha() {
        return this.alpha;
    }

    @Override
    public double getEstimatedProbeCost() {
        return this.predictedBloomProbeCost;
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
     * F(L) = Π q_i(L)  (all-pass),  H(L) = tail-sum expected probes.
     */


    @Override
    public int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost) {   // measured lc (ns)
        if (bfCost > 0) this.bloomProbeCost = bfCost;
        if (leafCost > 0) this.leafSearchCost = leafCost;

        // Per-symbol probabilities for the pattern from the estimator
        final double[] probs = tree.estimator.estimateALl(p);
        final int width      = tree.baseIntervalSize();
        final int maxDepth   = tree.maxDepth();
        final int r          = p.nGramToInt.length;
        double minProb = Arrays.stream(probs).min().getAsDouble();


        //Get the maximum Lp level (Deepest Lp) that can actually FIT the pattern. After that level
        //Our pattern matching algorithm switches to leaf probing so there is no reason to search further
        //e.g. if pattern length is 17 then that level is 5, since 5^2 = 32 --> fits 17. (4^2 = 16 doesnt fit it)
        //11
        int maxLpLevel = MathUtils.pruningLevel(tree, 0.05, minProb);

        //Thats the Lp level that is closer to the root
        int minLp =  MathUtils.pruningLevel(tree, 0.99, minProb);


        double bestCost = Double.POSITIVE_INFINITY;
        int    bestLp   = 0;

        // Evaluate cost for each candidate Lp (no α sweep needed)
        for (int Lp = minLp; Lp <= maxLpLevel; Lp++) {

            double cost = costAtLevel(tree, probs, p.nGramToInt, Lp, tree.getMembershipFpRate(Lp), maxDepth);
            if (cost < bestCost) {
                bestCost = cost;
                bestLp   = Lp;
                this.predictedBloomProbeCost = cost/this.bloomProbeCost;
            }
        }

//            // Record an α consistent with bestLp (purely for reporting)
//            double Fbest = MathUtils.allPass(probs, width, bestLp, bloomFp);
//            double Fprev = (bestLp > 0) ? MathUtils.allPass(probs, width, bestLp - 1, bloomFp) : 1.0;
//            this.alpha   = Math.min(0.99, Math.max(Fbest + 1e-9, 0.5 * (Fbest + Fprev)));

        return bestLp;
    }


//    public double costAtLevel(ImplicitTree<?> tree,
//                              double[] probs,
//                              int[] keySeq,
//                              int Lp,
//                              double _unusedBloomFp,  // keep signature
//                              int _unusedStopLp) {
//
//        final int width    = tree.baseIntervalSize();
//        final int maxLevel = tree.maxDepth() - 1;   // levels are 0..maxLevel
//        final int r        = probs.length;
//
//        // --- Horizontal (Lp) stays as you have it ---
//        double betaLp       = tree.getMembershipFpRate(Lp);
//        double H_lp         = MathUtils.expectedProbesPerNode(probs, keySeq, betaLp, width, Lp);
//        double parentsVisited = 1 << Lp;            // nodes we probe at Lp
//        double C_hor        = bloomProbeCost * H_lp * parentsVisited;
//
//        // We'll carry the conditional pass prob for the current level.
//        // At Lp there is no conditioning from above: qCond(Lp) = qUncond(Lp).
//        double qCondPrev;
//        {
//            double qUncondLp = MathUtils.fp_rate(probs, keySeq, width, Lp, betaLp);
//            qCondPrev = qUncondLp;
//        }
//
//        double C_vert = 0.0;
//
//        for (int L = Lp + 1; L <= maxLevel; L++) {
//            if (!MathUtils.childCanHost(width, L, r)) break;
//
//            // parents that actually PASSED at previous level (conditional)
//            double parentsPassed = parentsVisited * qCondPrev;
//
//            // children visited at this level: we probe BOTH children of each passed parent
//            double nodesVisitedL = 2.0 * parentsPassed;
//
//            // unconditional stats for this level (what MathUtils already knows)
//            double betaL    = tree.getMembershipFpRate(L);
//            double qUncondL = MathUtils.fp_rate(probs, keySeq, width, L, betaL); // = βL + (1-βL) hL
//            double hL       = (qUncondL - betaL) / (1.0 - betaL);                // recover hL
//
//            // CONDITIONAL pass prob at level L given parent passed at L-1
//            double qCondL = betaL + (1.0 - betaL) * (hL / qCondPrev);
//            if (qCondL < 0.0) qCondL = 0.0;
//            if (qCondL > 1.0) qCondL = 1.0;
//
//            // Make expectedProbesPerNode(...) use qCondL as its positive fraction
//            // by passing an "effective" beta that solves: qCondL = betaEff + (1-betaEff)*hL
//            double denom   = 1.0 - hL;
//            double betaEff = denom > 0 ? (qCondL - hL) / denom : 1.0; // guard when hL ~ 1
//            if (betaEff < 0.0) betaEff = 0.0;
//            if (betaEff > 1.0) betaEff = 1.0;
//
//            double H_L = MathUtils.expectedProbesPerNode(probs, keySeq, betaEff, width, L);
//            C_vert += bloomProbeCost * H_L * nodesVisitedL;
//
//            // advance
//            parentsVisited = nodesVisitedL;
//            qCondPrev      = qCondL;
//        }
//
//        return C_hor + C_vert;
//    }


    public double costAtLevel(ImplicitTree<?> tree,
                              double[] probs,     // probs[sym] = global empirical p̂ for symbol 'sym'
                              int[] keySeq,       // keySeq[i] = symbol id at position i
                              int Lp,
                              double _unusedBloomFp,   // keep signature
                              int _unusedStopLp) {

        final int width    = tree.baseIntervalSize();
        final int maxLevel = tree.maxDepth() - 1;     // levels are 0..maxLevel
        final int r        = probs.length;

        // --- Horizontal (Lp) as you had it ---
        double betaLp = tree.getMembershipFpRate(Lp);
        double H_lp   = MathUtils.expectedProbesPerNode(probs, keySeq, betaLp, width, Lp);
        double parentsVisited = 1 << Lp;                    // nodes we probe at Lp
        double C_hor = bloomProbeCost * H_lp * parentsVisited;

        // --- Carry per-symbol conditional YES probability down the tree.
        // Base (no conditioning from above at Lp): q_cond(Lp) = q_uncond(Lp).
        double qCondPrev; {
            double qUncondLp = MathUtils.fp_rate(probs, keySeq, width, Lp, betaLp);
            // IMPORTANT: qUncondLp here is an ALL-PASS factor. We want the PER-SYMBOL rate.
            // For uniform or when you treat all symbols at level Lp equivalently, use:
            //   hLp and qSymLp, then set qCondPrev = qSymLp (per-symbol).
            // We recover per-symbol q at Lp by averaging over first occurrences implicitly via the same beta:
            //   qSymLp = βLp + (1-βLp) * hLp, where hLp is computed from a representative symbol prob.
            //
            // To stay consistent with your current structure, we’ll approximate by taking
            // the *geometric mean per first-occurrence*, i.e., qSymLp ≈ (qAllPass)^(1/m),
            // where m is the number of first occurrences U at Lp.
            int bLp  = width >> Lp;
            int ellp = Math.min(keySeq.length, bLp);
            int mLp  = countFirstOccurrences(keySeq, ellp);
            double qSymLp = (mLp > 0) ? Math.pow(Math.max(qUncondLp, 1e-300), 1.0 / mLp) : 0.0;
            qCondPrev = clamp01(qSymLp);
        }

        double C_vert = 0.0;

        for (int L = Lp + 1; L <= maxLevel; L++) {
            if (!MathUtils.childCanHost(width, L, r)) break;

            // Parents that actually PASSED at previous level (conditional, set-level):
            // Use the same “geometric mean per first-occurrence” idea to lift the per-symbol qCondPrev
            // to a set-level pass prob for L-1.
            int bPrev  = width >> (L - 1);
            int ellPrev = Math.min(keySeq.length, bPrev);
            int mPrev   = countFirstOccurrences(keySeq, ellPrev);
            double FCondPrev = (mPrev > 0) ? Math.pow(qCondPrev, mPrev) : 0.0;

            double parentsPassed = parentsVisited * FCondPrev;

            // children visited at this level: we probe BOTH children of each passed parent
            double nodesVisitedL = 2.0 * parentsPassed;

            // --- Per-symbol unconditional at this level (for the same global p̂) ---
            double betaL = tree.getMembershipFpRate(L);

            // We need the per-symbol conditional YES prob at level L given parent passed:
            // qCondL = betaL + (1 - betaL) * (hL / qCondPrev)
            // where qCondPrev is the PER-SYMBOL conditional pass at L-1.
            // To compute hL, we need a representative per-symbol p̂; we get it from probs[sym].
            // Since your generator is uniform, p̂ ≈ 1/|Σ| for every symbol. We can pick any symbol id.
            // We’ll take the first symbol in the sequence as representative for h computations.

            double pSym = probs[0];

            double hL      = MathUtils.h_b(width, L,     pSym);
            double hPrev   = MathUtils.h_b(width, L - 1, pSym);
            // Per-symbol uncond at parent (for reference):
            double qSymPrevUncond = betaL + (1.0 - betaL) * hPrev; // only used when fallback is needed

            // Guard tiny denominators
            double denom = Math.max(qCondPrev, 1e-12);
            double qCondL = betaL + (1.0 - betaL) * (hL / denom);
            // Some extra safety clamping
            qCondL = clamp01(qCondL);

            // Make expectedProbesPerNode(...) use qCondL as the per-symbol YES rate by
            // passing an effective beta that satisfies: qCondL = betaEff + (1 - betaEff) * hL
            double betaEff;
            double oneMinusH = 1.0 - hL;
            if (oneMinusH <= 1e-12) {
                betaEff = 1.0; // degenerate: hL ~ 1 ⇒ any betaEff gives q≈1; choose 1 to be safe.
            } else {
                betaEff = (qCondL - hL) / oneMinusH;
            }
            betaEff = clamp01(betaEff);

            double H_L = MathUtils.expectedProbesPerNode(probs, keySeq, betaEff, width, L);
            C_vert    += bloomProbeCost * H_L * nodesVisitedL;

            // advance
            parentsVisited = nodesVisitedL;
            qCondPrev      = qCondL;
        }

        return C_hor + C_vert;
    }

    private static int countFirstOccurrences(int[] keySeq, int limit) {
        java.util.HashSet<Integer> seen = new java.util.HashSet<>();
        int m = 0;
        for (int i = 0; i < limit; i++) {
            if (seen.add(keySeq[i])) m++;
        }
        return m;
    }

    private static double clamp01(double x) {
        if (x < 0.0) return 0.0;
        if (x > 1.0) return 1.0;
        return x;
    }

//    public double costAtLevel(ImplicitTree<?> tree,
//                              double[] probs,
//                              int[] keySeq,
//                              int Lp,
//                              double _unused,              // keep signature
//                              int _unusedStopLp) {
//
//        final int width    = tree.baseIntervalSize();
//        final int maxLevel = tree.maxDepth() - 1;       // valid levels: 0..maxLevel
//        final int r        = probs.length;
//
//        // --- Horizontal at Lp (use β specific to Lp) ---
//        double betaLp = tree.getMembershipFpRate(Lp);
//        double H_lp   = MathUtils.expectedProbesPerNode(probs, keySeq, betaLp, width, Lp);
//        double C_hor  = bloomProbeCost * H_lp * (1 << Lp);
//
//        // --- Vertical: unconditional branching ---
//        double parents = 1 << Lp;   // nodes processed at Lp (already paid by C_hor)
//        double C_vert  = 0.0;
//
//        for (int L = Lp + 1; L <= maxLevel; L++) {
//            // children to visit at level L: 2 * parents_that_all_passed at (L-1)
//            double betaPrev = tree.getMembershipFpRate(L - 1);
//            double F_prev   = MathUtils.fp_rate(probs, keySeq, width, L - 1, betaPrev); // uses ℓ = min(r, b_{L-1})
//            double nodesL   = 2.0 * parents * F_prev;
//
//            double betaL = tree.getMembershipFpRate(L);
//            double H_L   = MathUtils.expectedProbesPerNode(probs, keySeq, betaL, width, L);
//            C_vert      += bloomProbeCost * H_L * nodesL;
//
//            if (!MathUtils.childCanHost(width, L, r)) break; // cannot go deeper
//            parents = nodesL;                                 // these become parents for next L
//        }
//
//        // leaf term is off in your experiment; keep it off
//        return C_hor + C_vert;
//    }

    /**
     * Returns (verticalCost, followedFinal).
     * followedFinal is the expected "1 + F*(N-1)" at the last processed level
     */
    private Pair<Double, Double> verticalCostBranching(double[] probs,
                                                       int[] keySeq,
                                                       int width,
                                                       int maxDepth,
                                                       int r,
                                                       int Lp,
                                                       double bloomFp, int stopLp, ImplicitTree<?> tree) {

//        double minProb =         Arrays.stream(probs).min().getAsDouble();
        double betaLp = tree.getMembershipFpRate(Lp);

        double F_lp    = MathUtils.fp_rate(probs, keySeq, width, Lp, betaLp);
        if (!MathUtils.childCanHost(width, Lp, r)) {
            // No deeper levels to visit; followed = 1 + F(Lp)*(2^Lp - 1)
//            double F_lp    = MathUtils.fp_rate(probs, width, Lp, bloomFp);
            double followed = 1.0 + F_lp * ((1 << Lp) - 1.0);
            return Pair.create(0.0, followed);
        }

        // Start from children of Lp
        double followed = 1.0 + F_lp * ((1 << Lp) - 1.0);
        double N        = 2.0 * followed;               // nodes processed at Lp+1



        double cost = 0.0;
        for (int L = Lp + 1; L <= stopLp; L++) {
            // Pay per-node probes at this level
            double betaL = tree.getMembershipFpRate(L);

            double H_L = MathUtils.expectedProbesPerNode(probs, keySeq, betaL, width, L);
            cost += bloomProbeCost * H_L * N;

            // Update branching if we can go deeper
            if (MathUtils.childCanHost(width, L, r)) {
                double F_L = MathUtils.fp_rate(probs, keySeq, width, L, betaL);
                followed   = 1.0 + F_L * (N - 1.0);
                N          = 2.0 * followed;
            } else {
                break;
            }
        }

        return Pair.create(cost, followed);
    }





}
