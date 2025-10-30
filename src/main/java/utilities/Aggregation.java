package utilities;

import utilities.BenchmarkEnums.IndexType;

import java.util.EnumMap;

/**
 * Aggregation helpers for benchmark stats to keep HBIDatasetBenchmarkMulti lean.
 */
public final class Aggregation {
    private Aggregation() {}

    public static final class AggregateStats {
        private final EnumMap<IndexType, StatsSum> sums = new EnumMap<>(IndexType.class);

        public void accumulate(IndexType indexType,
                               double avgInsertMsPerSymbol,
                               double totalInsertMs,
                               double totalQueryMs) {
            StatsSum sum = sums.computeIfAbsent(indexType, ignored -> new StatsSum());
            sum.avgInsertMsPerSymbolSum += avgInsertMsPerSymbol;
            sum.totalInsertMsSum += totalInsertMs;
            sum.totalQueryMsSum += totalQueryMs;
            sum.count++;
        }

        public StatsSnapshot snapshot(IndexType indexType) {
            StatsSum sum = sums.get(indexType);
            if (sum == null || sum.count == 0) {
                return null;
            }
            return new StatsSnapshot(
                    sum.avgInsertMsPerSymbolSum / sum.count,
                    sum.totalInsertMsSum / sum.count,
                    sum.totalQueryMsSum / sum.count,
                    sum.count);
        }

        public int count(IndexType indexType) {
            StatsSum sum = sums.get(indexType);
            return sum == null ? 0 : sum.count;
        }

        public boolean hasData() {
            return sums.values().stream().anyMatch(sum -> sum.count > 0);
        }

        private static final class StatsSum {
            double avgInsertMsPerSymbolSum;
            double totalInsertMsSum;
            double totalQueryMsSum;
            int count;
        }
    }

    public record StatsSnapshot(double avgInsertMsPerSymbol,
                                double avgInsertMs,
                                double avgQueryMs,
                                int count) {}
}
