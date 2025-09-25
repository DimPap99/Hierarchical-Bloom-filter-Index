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
    double bloomProbeCost = 1;
    double leafSearchCost = 1;

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
        final int r          = probs.length;
        double minProb = Arrays.stream(probs).min().getAsDouble();



        //Thats the Lp level that is closer to the root
        int minLp =  MathUtils.pruningLevel(tree, 0.99, minProb);

        //Get the maximum Lp level (Deepest Lp) that can actually FIT the pattern. After that level
        //Our pattern matching algorithm switches to leaf probing so there is no reason to search further
        //e.g. if pattern length is 17 then that level is 5, since 5^2 = 32 --> fits 17. (4^2 = 16 doesnt fit it)
        //11
        //minLp + 2;//
        int maxLpLevel = Math.min(MathUtils.pruningLevel(tree, 0.05, minProb), (int) (tree.maxDepth()-1 - Math.ceil(Math.log(p.originalSz)/Math.log(2))));


        double bestCost = Double.POSITIVE_INFINITY;
        int    bestLp   = 0;

        // Evaluate cost for each candidate Lp (no α sweep needed)
        for (int Lp = minLp; Lp <= maxLpLevel; Lp++) {

            double cost = costAtLevel(tree, probs, p.effectiveNgramArr, Lp, tree.getMembershipFpRate(Lp), maxDepth);
            if (cost < bestCost) {
                bestCost = cost;
                bestLp   = Lp;
                this.predictedBloomProbeCost = cost;
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
//        //  Horizontal (Lp)
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


    private static double[] qCondChildGivenParent(double[] probs,
                                                  int width,
                                                  int level,     // child level
                                                  double betaPrev,
                                                  double betaL) {
        int bL   = width >> level;
        int ell  = Math.min(probs.length, bL);
        double[] q = new double[ell];
        for (int i = 0; i < ell; i++) {
            double p     = probs[i];
            double hPrev = MathUtils.h_b(width, level - 1, p);
            double hL    = MathUtils.h_b(width, level,     p);

            double numer = hL + betaL * (hPrev - hL) + betaL * betaPrev * (1.0 - hPrev);
            double denom = betaPrev + (1.0 - betaPrev) * hPrev;

            double qc = (denom > 0.0) ? (numer / denom) : 1.0; // if denom==0, parent YES impossible → child won’t be visited; 1.0 is harmless here
            q[i] = clamp01(qc);
        }
        return q;
    }
    private static double[] qUncondAtLevel(double[] probs, int width, int level, double beta) {
        int bL   = width >> level;
        int ell  = Math.min(probs.length, bL);
        double[] q = new double[ell];
        for (int i = 0; i < ell; i++) {
            double p = probs[i];
            double h = MathUtils.h_b(width, level, p);
            q[i] = clamp01(beta + (1.0 - beta) * h);
        }
        return q;
    }

    private static double expectedProbesFromQ(double[] q, int[] keySeq, int width, int level) {
        int bL   = width >> level;
        int ell  = Math.min(keySeq.length, bL);
        if (ell <= 0) return 0.0;
        if (ell == 1) return 1.0;

        java.util.HashSet<Integer> seen = new java.util.HashSet<>();
        double total = 1.0;    // P(N>=1) = 1
        double prod  = 1.0;    // running product over q's of first occurrences so far

        for (int i = 0; i < ell - 1; i++) {
            int sym = keySeq[i];
            if (seen.add(sym)) {
                prod *= q[i];
                if (prod == 0.0) {
                    // all subsequent tail terms are zero
                    return total;
                }
            }
            total += prod;
        }
        return total;
    }
    private static int deepestVisitedLevel(int width, int r) {
        int q = width / Math.max(r,1);
        if (q <= 0) return 0;
        return 31 - Integer.numberOfLeadingZeros(q); // floor(log2(width/r))
    }

    private static double productFirstOccurrences(double[] q, int[] keySeq, int width, int level) {
        int bL   = width >> level;
        int ell  = Math.min(keySeq.length, bL);
        if (ell <= 0) return 1.0;  // empty prefix passes trivially

        java.util.HashSet<Integer> seen = new java.util.HashSet<>();
        double prod = 1.0;
        for (int i = 0; i < ell; i++) {
            int sym = keySeq[i];
            if (seen.add(sym)) {
                prod *= q[i];
                if (prod == 0.0) return 0.0; // short-circuit
            }
        }
        return prod;
    }





    // #distinct symbols among the first ell positions (first-occurrence count)
    private static int countDistinctAtLevel(int[] keySeq, int ell) {
        java.util.HashSet<Integer> seen = new java.util.HashSet<>();
        int m = 0;
        for (int i = 0; i < ell; i++) if (seen.add(keySeq[i])) m++;
        return m;
    }
//IE Principle cost
//
//    @Override
//    public double costAtLevel(ImplicitTree<?> tree,
//                              double[] probs,   // per-position per-slot masses p_i
//                              int[]    keySeq,  // symbols in probe order
//                              int      Lp,
//                              double   _unusedBloomFp,
//                              int      _unusedStopLp) {
//
//        final int width = tree.baseIntervalSize();
//        final int r     = keySeq.length;
//        final int Ldesc = deepestVisitedLevel(width, r);
//
//        double parentsVisited = 1 << Lp;   // S(Lp): nodes probed at base level
//        double total = 0.0;
//
//        //  Base level: unconditional masses
//        {
//            final int b   = width >> Lp;
//            final int ell = Math.min(r, b);
//            final double beta = tree.getMembershipFpRate(Lp);
//
//            // per-position masses at this node: unconditional == empirical
//            double[] pPerPos = new double[ell];
//            System.arraycopy(probs, 0, pPerPos, 0, ell);
//
//            HfPair hf = hfInsideNodeIEorB2(pPerPos, keySeq, ell, b, beta);
//            total += hf.H * parentsVisited;
//
//            // #children visited at next level = 2 * S(Lp) * F(Lp)
//            if (Lp < Ldesc) {
//                parentsVisited = 2.0 * parentsVisited * hf.F;  // S(Lp+1)
//                if (parentsVisited <= 0.0) return total;
//            } else {
//                return total;
//            }
//        }
//
//        //  Deeper levels: conditional masses
//        for (int L = Lp + 1; L <= Ldesc; L++) {
//            if (!MathUtils.childCanHost(width, L - 1, r)) break;
//
//            final int b   = width >> L;
//            final int ell = Math.min(r, b);
//            final double betaPrev = tree.getMembershipFpRate(L - 1);
//            final double beta     = tree.getMembershipFpRate(L);
//
//            // Build child-level conditional YES via your existing routine
//            double[] qCond = qCondChildGivenParent(probs, width, L, betaPrev, beta);
//
//            // Map YES → true child-block presence g → per-slot mass pCond for this node
//            double[] pPerPos = new double[ell];
//            for (int i = 0; i < ell; i++) {
//                double q = clamp01(qCond[i]);
//                double omb = 1.0 - beta;
//                double g = (omb > 0.0) ? (q - beta) / omb : 1.0;    // strip FP
//                if (g < 0.0) g = 0.0; else if (g > 1.0) g = 1.0;
//                pPerPos[i] = 1.0 - Math.pow(1.0 - g, 1.0 / Math.max(1, b));
//            }
//
//            // Work INSIDE this level’s nodes
//            HfPair hf = hfInsideNodeIEorB2(pPerPos, keySeq, ell, b, beta);
//
//            // IMPORTANT: parentsVisited already equals S(L) (children that reached this level).
//            double nodesVisitedL = parentsVisited;    // NOT 2 * parentsVisited
//            total += hf.H * nodesVisitedL;
//
//            // Prepare for next level: S(L+1) = 2 * S(L) * F(L)
//            if (L < Ldesc) {
//                parentsVisited = 2.0 * nodesVisitedL * hf.F;  // children count for next iteration
//                if (parentsVisited <= 0.0) break;
//            }
//        }
//
//        return total;
//    }

    // ====== return type ======



    @Override
    public double costAtLevel(ImplicitTree<?> tree,
                              double[] probs,   // per-position per-slot masses p_i, aligned 1-1 with keySeq
                              int[]    keySeq,  // symbols in probe order
                              int      Lp,
                              double   _unusedBloomFp,
                              int      _unusedStopLp) {

        final int width = tree.baseIntervalSize();
        final int r     = keySeq.length;
        final int Ldesc = deepestVisitedLevel(width, r);

        double total = 0.0;
        double betaL = 0;//tree.getMembershipFpRate(Lp);
        double betaPrev;
        // ---------- Base level Lp: UNCONDITIONAL IE from probs ----------
        HF nsLp = HF_uncond_pos_beta(width, Lp, keySeq, probs, betaL); // no Bloom => β=0 inside helper
        double nodes = (1 << Lp);
        total += nsLp.H * nodes;

        if (Lp >= Ldesc) return total;

        // Children that reach Lp+1 come only from parents that passed at Lp
        nodes = 2.0 * nodes * nsLp.F;
        if (nodes <= 0.0) return total;

        // ---------- Deeper levels: CONDITIONAL IE using qCond (child | parent passed) ----------
        for (int L = Lp + 1; L <= Ldesc; L++) {
            if (!MathUtils.childCanHost(width, L - 1, r)) break;

            betaL = 0f;//tree.getMembershipFpRate(Lp);
            betaPrev = 0f;//tree.getMembershipFpRate(Lp - 1);
            // Your existing conditional YES per position (β=0)
            double[] qCondL = qCondChildGivenParent(probs, width, L, betaPrev, betaL);

            // IE inside node using effective per-slot masses derived from qCond
            HF ns = HF_cond_from_q_pos_beta(width, L, keySeq, qCondL, betaL);
            total += ns.H * nodes;

            if (L < Ldesc) {
                nodes = 2.0 * nodes * ns.F;   // only passing parents spawn children
                if (nodes <= 0.0) break;
            }
        }

        return total;
    }


    public double costAtLevel_biased(ImplicitTree<?> tree,
                              double[] probs,   // probs[i] = p_i (per-position)
                              int[]    keySeq,  // keySeq[i] = symbol at position i
                              int      Lp,
                              double   _unusedBloomFp,
                              int      _unusedStopLp) {

        final int width = tree.baseIntervalSize();
        final int r     = keySeq.length;
        final int Ldesc = deepestVisitedLevel(width, r); //Lp+2; //


        //  Horizontal cost at Lp (UNCONDITIONAL, exact per-position q)
        double betaLp      = tree.getMembershipFpRate(Lp);
        double[] qUncondLp = qUncondAtLevel(probs, width, Lp, betaLp);   // q_i(Lp)
        double H_lp        = expectedProbesFromQ(qUncondLp, keySeq, width, Lp);
        double parentsVisited = 1 << Lp;
        double C_hor          = H_lp * parentsVisited;

        // We'll carry the per-position CONDITIONAL array down the tree.
        double[] qCondPrev = qUncondLp;
        double C_vert = 0.0;

        for (int L = Lp + 1; L <= Ldesc; L++) {
            // We reach level L iff children of (L-1) can still host the full pattern (i.e., b_L >= r).
            if (!MathUtils.childCanHost(width, L - 1, r)) break;

            // Branching: parents that pass at (L-1). Use the SAME q-array to form the pass prob.
            double FcondPrev     = productFirstOccurrences(qCondPrev, keySeq, width, L - 1);
            double parentsPassed = parentsVisited * FcondPrev;
            if (parentsPassed <= 0.0) break; // early exit if nothing passes

            double nodesVisitedL = 2.0 * parentsPassed;

            // Build child-level conditional q_i given the parent level
            double betaPrev = tree.getMembershipFpRate(L - 1);
            double betaL    = tree.getMembershipFpRate(L);
            double[] qCondL = qCondChildGivenParent(probs, width, L, betaPrev, betaL);

            // Expected probes at this level
            double H_L = expectedProbesFromQ(qCondL, keySeq, width, L);
            C_vert += H_L * nodesVisitedL;

            // Advance to next level
            parentsVisited = nodesVisitedL;
            qCondPrev      = qCondL;
        }

        return C_hor + C_vert;
    }


    // Result container
    // Result container
    public static final class HF {
        public final double H;  // expected probes in this node
        public final double F;  // all-YES prob for distincts that fit this node
        public HF(double H, double F) { this.H = H; this.F = F; }
    }

    // ---------- Base level: unconditional IE with Bloom (βL) ----------
    public static HF HF_uncond_pos_beta(int width, int level,
                                        int[] keySeq, double[] probs, double betaL) {
        final int bL  = width >> level;
        final int ell = Math.min(keySeq.length, bL);
        if (ell <= 0) return new HF(0.0, 1.0);

        java.util.HashSet<Integer> seen = new java.util.HashSet<>(ell * 2);
        java.util.ArrayList<Integer> first = new java.util.ArrayList<>();
        for (int pos = 0; pos < ell; pos++) if (seen.add(keySeq[pos])) first.add(pos);
        final int M = first.size();
        if (M == 0) return new HF(1.0, 1.0);

        double[] p = new double[M];                        // per-slot p_i by POSITION
        for (int m = 0; m < M; m++) p[m] = clamp01(probs[first.get(m)]);

        double[] Fm = new double[M];
        for (int m = 1; m <= M; m++) Fm[m - 1] = IE_prefix_collapsed_beta(p, m, bL, betaL);

        double H = 1.0;
        for (int m = 0; m < M; m++) {
            int next = (m + 1 < M) ? first.get(m + 1) : (ell - 1);  // 0-based
            int mult = next - first.get(m);                          // k_{m+1} - k_m
            H += mult * Fm[m];
        }
        return new HF(H, Fm[M - 1]);
    }

    // ---------- Deeper levels: conditional IE from qCond[] with Bloom (βL) ----------
    public static HF HF_cond_from_q_pos_beta(int width, int level,
                                             int[] keySeq, double[] qCond, double betaL) {
        final int bL  = width >> level;
        final int ell = Math.min(keySeq.length, bL);
        if (ell <= 0) return new HF(0.0, 1.0);

        java.util.HashSet<Integer> seen = new java.util.HashSet<>(ell * 2);
        java.util.ArrayList<Integer> first = new java.util.ArrayList<>();
        for (int pos = 0; pos < ell; pos++) if (seen.add(keySeq[pos])) first.add(pos);
        final int M = first.size();
        if (M == 0) return new HF(1.0, 1.0);

        final double omb = 1.0 - betaL;

        // qCond -> presence g (strip FPs) -> per-slot mass pEff
        double[] pEff = new double[M];
        for (int m = 0; m < M; m++) {
            double q = clamp01(qCond[first.get(m)]);
            double g = (omb > 0.0) ? clamp01((q - betaL) / omb) : 1.0;   // presence prob in child
            double p = 1.0 - Math.pow(1.0 - g, 1.0 / Math.max(1, bL));   // per-slot mass
            pEff[m] = clamp01(p);
        }

        double[] Fm = new double[M];
        for (int m = 1; m <= M; m++) Fm[m - 1] = IE_prefix_collapsed_beta(pEff, m, bL, betaL);

        double H = 1.0;
        for (int m = 0; m < M; m++) {
            int next = (m + 1 < M) ? first.get(m + 1) : (ell - 1);
            int mult = next - first.get(m);
            H += mult * Fm[m];
        }
        return new HF(H, Fm[M - 1]);
    }

    // ---------- Collapsed IE with Bloom (β) for first m distincts ----------
    private static double IE_prefix_collapsed_beta(double[] pFirstOccur, int m, int bL, double betaL) {
        final int subsets = 1 << m;
        final double omb = 1.0 - betaL;
        double F = 0.0;
        for (int mask = 0; mask < subsets; mask++) {
            int bits = Integer.bitCount(mask);      // |U|
            double sumP = 0.0;
            for (int i = 0; i < m; i++) if ((mask & (1 << i)) != 0) sumP += pFirstOccur[i];
            double base  = clamp01(1.0 - sumP);     // Pr(no symbol from U in a slot)
            double coeff = ((bits & 1) == 0 ? 1.0 : -1.0) * Math.pow(omb, bits);
            F += coeff * Math.pow(base, bL);
        }
        return clamp01(F);
    }

    private static double clamp01(double x) { return (x <= 0.0) ? 0.0 : (x >= 1.0) ? 1.0 : x; }


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
//        //  Horizontal at Lp (use β specific to Lp)
//        double betaLp = tree.getMembershipFpRate(Lp);
//        double H_lp   = MathUtils.expectedProbesPerNode(probs, keySeq, betaLp, width, Lp);
//        double C_hor  = bloomProbeCost * H_lp * (1 << Lp);
//
//        //  Vertical: unconditional branching
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
