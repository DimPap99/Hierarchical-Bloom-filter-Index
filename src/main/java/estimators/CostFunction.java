package estimators;

import PMIndex.NgramModel;
import search.Pattern;
import tree.ImplicitTree;

public interface CostFunction extends LevelCostProvider {

    int minCostLp(ImplicitTree tree, double confInit, Pattern p, double bfCost, double leafCost, boolean strides);

    // Predicted membership tests when Block Search starts at level lp. Cost functions that need the
    // whole pattern (its character length, the stride layout) override this; the default keeps the
    // key-sequence call the experiment harness has always made.
    default double costAtLevel(ImplicitTree<?> tree, Pattern p, int lp, boolean strides) {
        double[] probs = tree.estimator.estimateALl(p, strides);
        int stopLp = (int) (tree.maxDepth() - 1 - Math.ceil(Math.log(p.nGramToLong.length) / Math.log(2)));
        return costAtLevel(tree, probs, p.nGramToLong, lp, 0.0, stopLp);
    }

    double getAlpha();

    double getEstimatedProbeCost();

    void setModel(NgramModel.Model bigramModel);

    NgramModel.Model getBigramModel();
//    double costAtLevel(int level);
}
