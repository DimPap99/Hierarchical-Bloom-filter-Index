package estimators;

import tree.ImplicitTree;
import utilities.MathUtils;

public class CostFunctionIE extends AbstractCostFunction {

    private int ieMaxOrder = Integer.MAX_VALUE;

    public CostFunctionIE() {
    }

    public CostFunctionIE(int ieMaxOrder) {
        setIEMaxOrder(ieMaxOrder);
    }

    public CostFunctionIE setIEMaxOrder(int ieMaxOrder) {
        if (ieMaxOrder < 0) {
            throw new IllegalArgumentException("ieMaxOrder must be >= 0");
        }
        this.ieMaxOrder = ieMaxOrder;
        return this;
    }

    public int ieMaxOrder() {
        return ieMaxOrder;
    }

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
        double nodesAtLevel = 1 << Lp;

        MathUtils.HF parentStats = MathUtils.HF_uncond_pos_beta(width, Lp, keySeq, probs, 0.0, ieMaxOrder);
        total += parentStats.H * nodesAtLevel;

        int level = Lp;
        while (level < Ldesc && MathUtils.childCanHost(width, level, r)) {
            int nextLevel = level + 1;
            double[] qCond = MathUtils.qCondChildGivenParent(probs, width, nextLevel, 0.0, 0.0);
            MathUtils.HF childStats = MathUtils.HF_cond_from_q_pos_beta(width,
                                                                        nextLevel,
                                                                        keySeq,
                                                                        qCond,
                                                                        0.0,
                                                                        ieMaxOrder);

            double fanout = SelectiveFanout.multiplier(
                    level,
                    Lp,
                    Ldesc,
                    parentStats.F,
                    SelectiveFanout.costEfficiencyScore(parentStats.H, childStats.H));
            nodesAtLevel = fanout * nodesAtLevel * parentStats.F;
            if (nodesAtLevel <= 0.0) {
                break;
            }

            total += childStats.H * nodesAtLevel;
            parentStats = childStats;
            level = nextLevel;
        }

        return total;
    }
}
