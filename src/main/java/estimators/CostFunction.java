package estimators;

import search.Pattern;
import tree.ImplicitTree;

public interface CostFunction {

    int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost);

    double getAlpha();

    double getEstimatedProbeCost();

    double costAtLevel(ImplicitTree<?> tree,  double[] probs, int[] keySeq, int Lp, double bloomFp, int stopLp);
//    double costAtLevel(int level);
}
