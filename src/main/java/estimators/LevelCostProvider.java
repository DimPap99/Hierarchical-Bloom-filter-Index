package estimators;

import tree.ImplicitTree;

/**
 * Supplies cost estimates for a given tree level.
 */
public interface LevelCostProvider {

    double costAtLevel(ImplicitTree<?> tree,
                       double[] probs,
                       long[] keySeq,
                       int Lp,
                       double bloomFp,
                       int stopLp);
}

