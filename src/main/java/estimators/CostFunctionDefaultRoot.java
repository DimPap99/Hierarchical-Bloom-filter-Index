package estimators;

import search.Pattern;
import tree.ImplicitTree;

public class CostFunctionDefaultRoot extends AbstractCostFunction{
    @Override
    public double costAtLevel(ImplicitTree<?> tree, double[] probs, long[] keySeq, int Lp, double bloomFp, int stopLp) {
        return 0;
    }
    @Override
    public int minCostLp(ImplicitTree tree,
                         double confInit,
                         Pattern pattern,
                         double bfCost,
                         double leafCost,
                         boolean strides) {
        return tree.effectiveRoot(); //return the root
    }

}
