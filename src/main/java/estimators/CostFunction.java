package estimators;

import PMIndex.NgramModel;
import search.Pattern;
import tree.ImplicitTree;

public interface CostFunction extends LevelCostProvider {

    int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost, boolean strides);

    double getAlpha();

    double getEstimatedProbeCost();

    void setModel(NgramModel.Model bigramModel);

    NgramModel.Model getBigramModel();
//    double costAtLevel(int level);
}
