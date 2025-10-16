package estimators;

import PMIndex.NgramModel;
import search.Pattern;
import tree.ImplicitTree;

public interface CostFunction {

    int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost, boolean strides);

    double getAlpha();

    double getEstimatedProbeCost();

    double costAtLevel(ImplicitTree<?> tree,  double[] probs, long[] keySeq, int Lp, double bloomFp, int stopLp);

    void setModel(NgramModel.Model bigramModel);

    NgramModel.Model getBigramModel();
//    double costAtLevel(int level);
}
