package estimators;

import search.Pattern;
import tree.ImplicitTree;

public interface CostFunction {

    int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost);

    double getAlpha();
}
