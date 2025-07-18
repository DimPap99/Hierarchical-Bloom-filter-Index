package estimators;

import search.Pattern;
import tree.ImplicitTree;

public interface CostFunction {

    int minCostLp(ImplicitTree tree,
                         double   bfFalsePosRate,   // Bloom collision per symbol (β)
                         double confInit,
                         Pattern p,
                        double bfCost, double leafCost);


}
