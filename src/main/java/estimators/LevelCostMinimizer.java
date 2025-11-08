package estimators;

import search.Pattern;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.Arrays;

/**
 * Iterates over candidate pruning levels and selects the minimum-cost option using a {@link LevelCostProvider}.
 */
public final class LevelCostMinimizer {
    public static boolean isDebugging = false;
    private LevelCostMinimizer() {
        // utility
    }

    public static Result minimize(LevelCostProvider provider,
                                  ImplicitTree tree,
                                  double confInit,
                                  Pattern pattern,
                                  boolean strides) {
        double[] probs = tree.estimator.estimateALl(pattern, strides);
        if (probs.length == 0) {
            return new Result(0, 0.0);
        }

        int width = tree.baseIntervalSize();
        int maxDepth = tree.maxDepth();
        int r = probs.length;
        double minProb = Arrays.stream(probs).min().orElse(0.0);
        int _maxlp = MathUtils.pruningLevel(tree, 0.99, minProb);
        int minLp = 0;//MathUtils.pruningLevel(tree, 0.99, minProb);
        int patternFitDepth = computePatternFitDepth(tree, pattern);
        int maxLpLevel = MathUtils.pruningLevel(tree, 0.95, minProb);

        if (maxLpLevel < minLp) {
            maxLpLevel = minLp;
        }

        long[] keySeq = strides ? pattern.effectiveNgramArr : pattern.nGramToLong;

        double bestCost = Double.POSITIVE_INFINITY;
        int bestLp = minLp;

        for (int Lp = minLp; Lp <= maxLpLevel; Lp++) {
            double cost = provider.costAtLevel(tree, probs, keySeq, Lp, tree.getMembershipFpRate(Lp), maxDepth);
            if (cost < bestCost) {
                bestCost = cost;
                bestLp = Lp;
            }
        }
        if(bestLp < _maxlp && !isDebugging) {
            bestCost = _maxlp;
        }
        return new Result(bestLp, bestCost);
    }

    private static int computePatternFitDepth(ImplicitTree tree, Pattern pattern) {
        int maxDepth = tree.maxDepth() - 1;
        double patternLength = Math.max(1, pattern.originalSz);
        int depthOffset = (int) Math.ceil(Math.log(patternLength) / Math.log(2.0));
        int fitDepth = maxDepth - depthOffset;
        return Math.max(0, fitDepth);
    }

    public static final class Result {
        private final int bestLp;
        private final double bestCost;

        public Result(int bestLp, double bestCost) {
            this.bestLp = bestLp;
            this.bestCost = bestCost;
        }

        public int bestLp() {
            return bestLp;
        }

        public double bestCost() {
            return bestCost;
        }
    }
}

