package estimators;

import org.apache.commons.math3.util.Pair;
import search.Pattern;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.Arrays;

import static utilities.MathUtils.pruningLevel;

//Dummy CF used to manually set confidence values for Lp
public class CostFunctionSetConf  implements CostFunction {
    @Override
    public int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost) {
        double pMax = Arrays.stream(tree.estimator.estimateALl(p)).min().getAsDouble();

        return pruningLevel(tree, confInit, pMax);
    }

    @Override
    public double getAlpha() {
        return 0;
    }
}







