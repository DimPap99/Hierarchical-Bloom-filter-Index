package estimators;

import tree.ImplicitTree;
import utilities.MathUtils;

public class CostFunctionIE extends AbstractCostFunction {

    @Override
    public double costAtLevel(ImplicitTree<?> tree,
                              double[] probs,
                              long[] keySeq,
                              int Lp,
                              double bloomFp,
                              int stopLp) {
        final int width = tree.baseIntervalSize();
        final int r = keySeq.length;
        final int Ldesc = MathUtils.deepestVisitedLevel(width, r);

        double total = 0.0;
        double nodes = 1 << Lp;

        MathUtils.HF nsLp = MathUtils.HF_uncond_pos_beta(width, Lp, keySeq, probs, 0.0);
        total += nsLp.H * nodes;

        if (Lp >= Ldesc) {
            return total;
        }

        nodes = 2.0 * nodes * nsLp.F;
        if (nodes <= 0.0) {
            return total;
        }

        for (int L = Lp + 1; L <= Ldesc; L++) {
            if (!MathUtils.childCanHost(width, L - 1, r)) {
                break;
            }

            double[] qCondL = MathUtils.qCondChildGivenParent(probs, width, L, 0.0, 0.0);
            MathUtils.HF ns = MathUtils.HF_cond_from_q_pos_beta(width, L, keySeq, qCondL, 0.0);
            total += ns.H * nodes;

            if (L < Ldesc) {
                nodes = 2.0 * nodes * ns.F;
                if (nodes <= 0.0) {
                    break;
                }
            }
        }

        return total;
    }
}
