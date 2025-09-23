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
    private long totalQueryTimeNanos;
    private long totalLpTimeNanos;
    private long queryCount;

    public HbiStats(boolean collectStats, boolean experimentMode) {
        this.collectStats = collectStats;
        this.experimentMode = experimentMode;
        resetTiming();
    }

    public boolean isCollecting() {
        return collectStats;
    }

    public void setCollecting(boolean collectStats) {
        this.collectStats = collectStats;
        if (!collectStats) {
            lpLevels.clear();
            alphas.clear();
            resetTiming();
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

    public void recordQueryTiming(long queryTimeNanos, long lpTimeNanos) {
        if (!collectStats) {
            return;
        }
        totalQueryTimeNanos += queryTimeNanos;
        totalLpTimeNanos += lpTimeNanos;
        queryCount++;
    }

    public double averageQueryTimeMillis() {
        if (queryCount == 0) {
            return 0.0;
        }
        return (totalQueryTimeNanos / 1_000_000.0) / queryCount;
    }

    public double averageLpTimeMillis() {
        if (queryCount == 0) {
            return 0.0;
        }
        return (totalLpTimeNanos / 1_000_000.0) / queryCount;
    }

    public double lpShareOfQuery() {
        if (totalQueryTimeNanos == 0) {
            return 0.0;
        }
        return (double) totalLpTimeNanos / totalQueryTimeNanos;
    }

    public long totalQueryCount() {
        return queryCount;
    }

    public long totalQueryTimeNanos() {
        return totalQueryTimeNanos;
    }

    public long totalLpTimeNanos() {
        return totalLpTimeNanos;
    }

    private void resetTiming() {
        totalQueryTimeNanos = 0L;
        totalLpTimeNanos = 0L;
        queryCount = 0L;
    }
}
