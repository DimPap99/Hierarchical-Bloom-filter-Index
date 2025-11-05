package utilities;

import utilities.BenchmarkEnums.IndexType;
import utilities.BenchmarkEnums.QueryType;

import java.io.IOException;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Console and CSV reporting for aggregated benchmark results.
 */
public final class BenchmarkReporter {
    private BenchmarkReporter() {}

    public static void printConsoleSummary(Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated,
                                           List<IndexType> activeIndexes,
                                           List<QueryType> queryTypes,
                                           double currentFp) {
        System.out.println();
        aggregated.forEach((type, byNgram) -> {
            int datasetCount = activeIndexes.stream()
                    .mapToInt(index -> byNgram.values().stream()
                            .flatMap(map -> map.values().stream())
                            .mapToInt(stats -> stats.count(index))
                            .max().orElse(0))
                    .max().orElse(0);

            System.out.printf(Locale.ROOT,
                    "Aggregated averages for %s across %d dataset(s) at FPR=%.6f%n",
                    type.fileToken(), datasetCount, currentFp);

            byNgram.forEach((ng, statsByPattern) -> {
                System.out.printf(Locale.ROOT, "  ngram %d%n", ng);
                statsByPattern.forEach((pl, stats) -> {
                    System.out.printf(Locale.ROOT, "    Pattern %d%n", pl);
                    for (IndexType index : activeIndexes) {
                        Aggregation.StatsSnapshot snap = stats.snapshot(index);
                        if (snap == null || snap.count() == 0) {
                            System.out.printf(Locale.ROOT, "      %s -> no data%n", index.displayName());
                            continue;
                        }
                        System.out.printf(Locale.ROOT,
                                "      %s -> avgInsertMs/sym=%.6f, avgInsertMs=%.3f, avgQueryMs=%.3f, count=%d%n",
                                index.displayName(), snap.avgInsertMsPerSymbol(), snap.avgInsertMs(), snap.avgQueryMs(), snap.count());
                    }
                });
            });
        });
    }

    public static Path writeCsv(MultiBenchmarkOptions options,
                                Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated,
                                List<IndexType> activeIndexes,
                                List<QueryType> queryTypes,
                                String algo,
                                double currentFp) throws IOException {

        List<Object> header = new ArrayList<>();
        header.add("pattern_len");
        header.add("ngram");
        for (QueryType type : queryTypes) {
            for (IndexType index : activeIndexes) {
                header.add("avg_insert_ms_per_symbol_" + type.fileToken() + "_" + index.csvLabel());
                header.add("avg_insert_ms_" + type.fileToken() + "_" + index.csvLabel());
                header.add("avg_query_ms_" + type.fileToken() + "_" + index.csvLabel());
                header.add("count_" + type.fileToken() + "_" + index.csvLabel());
            }
        }
        header.add("algo");
        header.add("policy");
        header.add("alphabet");
        header.add("fpr");
        header.add("delta");
        header.add("window_power");
        header.add("tree_power");
        header.add("rank_eps");
        header.add("delta_q");
        header.add("delta_samp");
        header.add("quantile");
        header.add("policy_buckets");
        header.add("policy_n_req");
        header.add("policy_lb");
        header.add("policy_impossible");

        List<List<Object>> csvRows = new ArrayList<>();
        csvRows.add(header);

        for (int ng : options.ngrams()) {
            // Union of pattern lengths across types
            Set<Integer> allPatternLengths = new TreeSet<>();
            for (QueryType type : queryTypes) {
                Map<Integer, Map<Integer, Aggregation.AggregateStats>> byN = aggregated.get(type);
                Map<Integer, Aggregation.AggregateStats> statsByPl = (byN != null) ? byN.get(ng) : null;
                if (statsByPl != null) allPatternLengths.addAll(statsByPl.keySet());
            }

            for (Integer pl : allPatternLengths) {
                List<Object> row = new ArrayList<>();
                row.add(pl);
                row.add(ng);
                for (QueryType type : queryTypes) {
                    Map<Integer, Map<Integer, Aggregation.AggregateStats>> byN = aggregated.get(type);
                    Map<Integer, Aggregation.AggregateStats> byPl = (byN != null) ? byN.getOrDefault(ng, Map.of()) : Map.of();
                    Aggregation.AggregateStats stats = byPl.get(pl);
                    for (IndexType index : activeIndexes) {
                        Aggregation.StatsSnapshot snap = (stats != null) ? stats.snapshot(index) : null;
                        if (snap == null || snap.count() == 0) {
                            row.add(null); row.add(null); row.add(null); row.add(0);
                        } else {
                            row.add(snap.avgInsertMsPerSymbol());
                            row.add(snap.avgInsertMs());
                            row.add(snap.avgQueryMs());
                            row.add(snap.count());
                        }
                    }
                }
                row.add(algo);
                row.add(options.memPolicy());
                row.add(options.alphabetSizeFor(ng));
                row.add(currentFp);
                int deltaVal = options.reinsertPerWorkload() ? pl : options.suffixDelta(); // delayed suffix default delta
                row.add(deltaVal);
                row.add(options.windowPower());
                row.add(options.treePower());
                row.add(options.rankEpsTarget());
                row.add(options.deltaQ());
                row.add(options.deltaSamp());
                row.add(options.quantile());
                if (options.memPolicy() != Utils.MemPolicy.NONE) {
                    Utils.HopsDesignResult design = options.policyDesignFor(ng);
                    row.add(design.suggestedBuckets);
                    row.add(design.requiredSampleSize);
                    row.add(design.occupancyLowerBound);
                    row.add(design.impossible);
                } else {
                    row.add(null);
                    row.add(null);
                    row.add(null);
                    row.add(null);
                }
                csvRows.add(row);
            }
        }

        String fpToken = String.format(Locale.ROOT, "%.6f", currentFp).replace('.', 'p');
        String ngToken = (options.ngrams().size() == 1) ? String.valueOf(options.ngrams().get(0)) : "ngmulti";
        String treeToken = "t" + options.treePower();
        String fileName = "%s_%s_%s_%s_%s_fp%s.csv"
                .formatted(options.window(), treeToken, options.queryTypeLabel(), ngToken, algo, fpToken);
        Path csvPath = Path.of(fileName);
        CsvUtil.writeRows(csvPath, csvRows);
        return csvPath;
    }
}
