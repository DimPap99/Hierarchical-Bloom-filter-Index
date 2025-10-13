package estimators;

import PMIndex.NgramModel;
import search.Pattern;
import tree.ImplicitTree;

import java.util.Arrays;

import static utilities.MathUtils.pruningLevel;

//Dummy CF used to manually set confidence values for Lp
public class CostFunctionSetConf  implements CostFunction {
    @Override
    public int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost, boolean strides) {
        double pMax = Arrays.stream(tree.estimator.estimateALl(p, strides)).min().getAsDouble();

        return pruningLevel(tree, confInit, pMax);
    }

    @Override
    public double getAlpha() {
        return 0;
    }

    @Override
    public double getEstimatedProbeCost() {
        return 0;
    }

    @Override
    public double costAtLevel(ImplicitTree<?> tree, double[] probs, int[] keySeq, int Lp, double bloomFp, int stopLp) {
        return 0;
    }

    @Override
    public void setModel(NgramModel.Model bigramModel) {

    }

    @Override
    public NgramModel.Model getBigramModel() {
        return null;
    }

//    @Override
//    public void setModel(NgramModel.Model bigramModel) {
//
//    }

//    @Override
//    public NgramModel.Model getBigramModel() {
//        return null;
//    }

}







