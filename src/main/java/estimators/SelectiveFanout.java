package estimators;

/**
 * Centralised helper for adapting the branching multiplier when we operate in a selective regime.
 * The toggle is static on purpose so experiments can opt in without plumbing additional parameters.
 */
public final class SelectiveFanout {
    private static final double MIN_MULT = 1.0;
    private static final double MAX_MULT = 2.0;
    private static final double DEPTH_WEIGHT = 0.4;
    private static final double FEAS_WEIGHT = 0.4;
    private static final double COST_WEIGHT = 0.2;
    private static volatile boolean selectiveRegimeEnabled = Boolean.getBoolean("hbi.selectiveRegime");

    private SelectiveFanout() {}

    /** Enable or disable the adaptive selective regime heuristics. */
    public static void setSelectiveRegimeEnabled(boolean enabled) {
        selectiveRegimeEnabled = enabled;
    }

    /** Returns whether the selective regime is active. */
    public static boolean isSelectiveRegimeEnabled() {
        return selectiveRegimeEnabled;
    }

    /**
     * Computes the branching multiplier for the given parent level.
     * When selective mode is disabled we fall back to the historical "always two children" assumption.
     */
    static double multiplier(int parentLevel,
                             int startLevel,
                             int descLimit,
                             double feasibilityScore,
                             double costEfficiencyScore) {
        if (!selectiveRegimeEnabled) {
            return MAX_MULT;
        }
        double normalizedDepth = 0.0;
        if (descLimit > startLevel) {
            normalizedDepth = (double) Math.max(0, parentLevel - startLevel)
                    / (double) (descLimit - startLevel);
        }
        double probScore = clamp01(feasibilityScore);
        double costScore = clamp01(costEfficiencyScore);
        double score = (DEPTH_WEIGHT * normalizedDepth)
                + (FEAS_WEIGHT * probScore)
                + (COST_WEIGHT * costScore);
        score = clamp01(score);
        return MIN_MULT + (MAX_MULT - MIN_MULT) * score;
    }

    /**
     * Normalises the relative expected probe cost between parent and child so that cheaper children push the score up.
     */
    static double costEfficiencyScore(double parentCost, double childCost) {
        if (childCost <= 0.0) {
            return 0.0;
        }
        if (parentCost <= 0.0) {
            return 1.0;
        }
        double ratio = parentCost / childCost;
        return clamp01(ratio);
    }

    private static double clamp01(double value) {
        if (value <= 0.0) {
            return 0.0;
        }
        if (value >= 1.0) {
            return 1.0;
        }
        return value;
    }
}

