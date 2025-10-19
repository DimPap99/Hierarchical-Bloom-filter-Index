package estimators;

import PMIndex.NgramModel;
import search.Pattern;
import tree.ImplicitTree;

/**
 * Shared scaffolding for cost functions that evaluate level costs through {@link LevelCostMinimizer}.
 */
public abstract class AbstractCostFunction implements CostFunction {

    protected double bloomProbeCost = 1.0;
    protected double leafSearchCost = 1.0;
    protected double predictedBloomProbeCost = 0.0;
    protected double alpha = 0.0;
    protected NgramModel.Model bigramModel;

    @Override
    public int minCostLp(ImplicitTree tree,
                         double confInit,
                         Pattern pattern,
                         double bfCost,
                         double leafCost,
                         boolean strides) {
        if (bfCost > 0.0) {
            this.bloomProbeCost = bfCost;
        }
        if (leafCost > 0.0) {
            this.leafSearchCost = leafCost;
        }

        LevelCostMinimizer.Result result = LevelCostMinimizer.minimize(this, tree, confInit, pattern, strides);
        this.predictedBloomProbeCost = result.bestCost();
        return result.bestLp();
    }

    @Override
    public double getAlpha() {
        return this.alpha;
    }

    @Override
    public double getEstimatedProbeCost() {
        return this.predictedBloomProbeCost;
    }

    @Override
    public void setModel(NgramModel.Model bigramModel) {
        this.bigramModel = bigramModel;
    }

    @Override
    public NgramModel.Model getBigramModel() {
        return this.bigramModel;
    }
}

