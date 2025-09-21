package PMIndex;

import utilities.PatternResult;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Collects optional runtime statistics for {@link HBI} without touching the hot path
 * when instrumentation is disabled.
 */
public final class HbiStats {

    private final List<Integer> lpLevels = new ArrayList<>();
    private final List<Double> alphas = new ArrayList<>();
    private PatternResult latestPatternResult;
    private boolean collectStats;
    private boolean experimentMode;

    public HbiStats(boolean collectStats, boolean experimentMode) {
        this.collectStats = collectStats;
        this.experimentMode = experimentMode;
    }

    public boolean isCollecting() {
        return collectStats;
    }

    public void setCollecting(boolean collectStats) {
        this.collectStats = collectStats;
        if (!collectStats) {
            lpLevels.clear();
            alphas.clear();
        }
    }

    public boolean isExperimentMode() {
        return experimentMode;
    }

    public void setExperimentMode(boolean experimentMode) {
        this.experimentMode = experimentMode;
    }

    public void recordLp(int lp) {
        if (collectStats) {
            lpLevels.add(lp);
        }
    }

    public void recordAlpha(double alpha) {
        if (collectStats) {
            alphas.add(alpha);
        }
    }

    public List<Integer> lpLevels() {
        return Collections.unmodifiableList(lpLevels);
    }

    public List<Double> alphas() {
        return Collections.unmodifiableList(alphas);
    }

    public PatternResult latestPatternResult() {
        return latestPatternResult;
    }

    public void setLatestPatternResult(PatternResult latestPatternResult) {
        this.latestPatternResult = latestPatternResult;
    }
}
