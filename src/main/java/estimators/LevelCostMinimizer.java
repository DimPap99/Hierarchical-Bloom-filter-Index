package estimators;

import search.Pattern;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.Arrays;

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
            return new Result(tree.effectiveRoot(), 0.0);
        }

        int width = tree.baseIntervalSize();
        int maxDepth = tree.maxDepth();
        int r = probs.length;
        double minProb = Arrays.stream(probs).min().orElse(0.0);
        int _maxlp = MathUtils.pruningLevel(tree, 0.99, minProb);
        // Levels above the effective root were dropped by a provisioning policy and cannot be searched.
        int minLp = tree.effectiveRoot();
        int maxLpLevel = (int) ((tree.maxDepth()-1) - Math.ceil(Math.log(pattern.nGramToLong.length)/Math.log(2)));

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
        // Also apply this clamp when SelectiveFanout is enabled, regardless of debugging flag
        if (bestLp > _maxlp && (SelectiveFanout.isSelectiveRegimeEnabled())) {
            bestLp = Math.max(minLp, _maxlp);
        }
        return new Result(bestLp, bestCost);
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
