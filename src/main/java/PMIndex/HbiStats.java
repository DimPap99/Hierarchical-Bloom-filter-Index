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
    private final List<Integer> cfLpLevels = new ArrayList<>();
    private final List<Double> alphas = new ArrayList<>();
    private final List<Long> minCostLpTimesNanos = new ArrayList<>();
    private PatternResult latestPatternResult;
    private boolean collectStats;
    private boolean experimentMode;
    private long totalQueryTimeNanos;
    private long totalLpTimeNanos;
    private long queryCount;
    private long totalMinCostLpTimeNanos;

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
            cfLpLevels.clear();
            resetExperimentStats();
            resetTiming();
        }
    }

    public boolean isExperimentMode() {
        return experimentMode;
    }

    public void setExperimentMode(boolean experimentMode) {
        this.experimentMode = experimentMode;
        if (!experimentMode) {
            resetExperimentStats();
        }
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

    public void recordCfLp(int cfLp) {
        if (collectStats) {
            cfLpLevels.add(cfLp);
        }
    }

    public List<Integer> cfLpLevels() {
        return Collections.unmodifiableList(cfLpLevels);
    }

    public List<Double> alphas() {
        return Collections.unmodifiableList(alphas);
    }

    public List<Long> minCostLpTimesNanos() {
        return Collections.unmodifiableList(minCostLpTimesNanos);
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

    public void recordMinCostLpTime(long durationNanos) {
        // Record CF minCost Lp timing when stats collection is enabled,
        // regardless of experiment mode. This keeps overhead opt‑in via collectStats.
        if (!collectStats) {
            return;
        }
        minCostLpTimesNanos.add(durationNanos);
        totalMinCostLpTimeNanos += durationNanos;
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

    /** Average chosen Lp level across queries (current algorithm). */
    public double averageChosenLp() {
        if (lpLevels.isEmpty()) return 0.0;
        long sum = 0;
        for (int v : lpLevels) sum += v;
        return sum * 1.0 / lpLevels.size();
    }

    /** Average chosen Lp level for CF min-cost (if recorded). */
    public double averageCfChosenLp() {
        if (cfLpLevels.isEmpty()) return 0.0;
        long sum = 0;
        for (int v : cfLpLevels) sum += v;
        return sum * 1.0 / cfLpLevels.size();
    }

    public double lpShareOfQuery() {
        if (totalQueryTimeNanos == 0) {
            return 0.0;
        }
        return (double) totalLpTimeNanos / totalQueryTimeNanos;
    }

    public double averageMinCostLpTimeMillis() {
        if (minCostLpTimesNanos.isEmpty()) {
            return 0.0;
        }
        return (totalMinCostLpTimeNanos / 1_000_000.0) / minCostLpTimesNanos.size();
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

    public long totalMinCostLpTimeNanos() {
        return totalMinCostLpTimeNanos;
    }

    private void resetTiming() {
        totalQueryTimeNanos = 0L;
        totalLpTimeNanos = 0L;
        queryCount = 0L;
    }

    private void resetExperimentStats() {
        minCostLpTimesNanos.clear();
        totalMinCostLpTimeNanos = 0L;
    }
}
