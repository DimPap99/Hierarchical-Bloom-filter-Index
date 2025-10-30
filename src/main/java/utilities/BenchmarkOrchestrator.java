package utilities;

import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import utilities.BenchmarkEnums.IndexType;
import utilities.BenchmarkEnums.QueryType;
import utilities.MultiQueryExperiment.InsertStats;
import utilities.MultiQueryExperiment.MultiRunResult;
import utilities.MultiQueryExperiment.QueryWorkload;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Encapsulates the outer benchmark loops (per FPR x n-gram x dataset) so
 * HBIDatasetBenchmarkMulti remains a thin CLI.
 */
public final class BenchmarkOrchestrator {
    private BenchmarkOrchestrator() {}
    public static void run(MultiBenchmarkOptions options) throws Exception {

        List<QueryType> queryTypes = (options.queryType() != null)
                ? List.of(options.queryType())
                : List.of(QueryType.UNIFORM, QueryType.MISSING, QueryType.RARE);
        List<IndexType> activeIndexes = new ArrayList<>();
        if (options.runHbi()) activeIndexes.add(IndexType.HBI);
        if (options.runSuffix()) activeIndexes.add(IndexType.SUFFIX);
        if (options.runRegexBaseline()) activeIndexes.add(IndexType.REGEX);

        System.out.printf(Locale.ROOT,
                "Running window %s (%s) with datasets under %s%n",
                options.window(), options.queryTypeLabel(), options.dataRoot());

        List<Path> datasetDirs = BenchmarkIO.listDatasetDirectories(options.dataRoot());
        if (datasetDirs.isEmpty()) {
            throw new IllegalStateException("No dataset folders found under " + options.dataRoot());
        }

        for (double currentFp : options.fpRates()) {
            Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated = new EnumMap<>(QueryType.class);
            queryTypes.forEach(type -> aggregated.put(type, new TreeMap<>()));

            System.out.printf(Locale.ROOT, "%n==== FPR = %.6f ====%n", currentFp);

            for (int ng : options.ngrams()) {
                System.out.printf(Locale.ROOT, "---- ngram = %d ----%n", ng);

                for (Path datasetDir : datasetDirs) {
                    Path datasetFile = BenchmarkIO.findSingleDatasetFile(datasetDir);
                    Path queryDir = options.queryRoot().resolve(datasetDir.getFileName());
                    if (!Files.isDirectory(queryDir)) {
                        System.out.printf(Locale.ROOT,
                                "Skipping dataset %s (missing query folder %s)%n",
                                datasetDir.getFileName(), queryDir);
                        continue;
                    }

                    System.out.printf(Locale.ROOT, "Dataset %s -> %s%n",
                            datasetDir.getFileName(), datasetFile.getFileName());

                    Map<QueryType, List<QueryWorkload>> workloadsByType = new EnumMap<>(QueryType.class);
                    boolean hasQueries = false;
                    for (QueryType type : queryTypes) {
                        List<Path> queryFiles = BenchmarkIO.findQueryFiles(queryDir, type);
                        if (queryFiles.isEmpty()) {
                            System.out.printf(Locale.ROOT,
                                    "No %s queries in %s, skipping.%n",
                                    type.fileToken(), queryDir);
                            workloadsByType.put(type, List.of());
                            continue;
                        }
                        hasQueries = true;
                        System.out.printf(Locale.ROOT, "  Query type %s%n", type.fileToken());
                        List<QueryWorkload> workloads = queryFiles.stream()
                                .map(path -> new QueryWorkload(
                                        BenchmarkIO.parsePatternLength(path.getFileName().toString()),
                                        path))
                                .sorted(Comparator.comparingInt(QueryWorkload::patternLength))
                                .collect(Collectors.toList());
                        workloadsByType.put(type, workloads);
                    }

                    if (!hasQueries) continue;

                    // Warmups
                    for (int warmIndex = 0; warmIndex < options.warmupRuns(); warmIndex++) {
                        if (options.runHbi()) {
                            HBI warm = IndexFactory.createHbi(
                                    options.windowLength(), options.treeLength(), options.alphabetSizeFor(ng), currentFp,
                                    options.runConfidence(), options.memPolicy(), ng);
                            warm.strides = true;
                            warm.stats().setCollecting(false);
                            warm.stats().setExperimentMode(false);
                            InsertStats warmIns = isSegments(options)
                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), warm, ng)
                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), warm, ng);
                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;
                                if (isSegments(options)) {
                                    SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, warm, ng, warmIns, false, false);
                                } else {
                                    MultiQueryExperiment.runQueries(workloads, warm, ng, warmIns, false, false);
                                }
                            }
                            warm = null;
                        }
                        if (options.runSuffix()) {
                            // Use delayed SSWSI; warmups don't target a specific pattern length, use configured default delta
                            //dont use ngrams for suffix cause it hurts performance
                            IPMIndexing suffixWarm = IndexFactory.createSuffixIndex(
                                    options.windowLength());
                            InsertStats suffixIns = isSegments(options)
                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffixWarm, 1)
                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffixWarm, 1);
                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;
                                if (isSegments(options)) {
                                    SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, suffixWarm, ng, suffixIns, false, false);
                                } else {
                                    MultiQueryExperiment.runQueries(workloads, suffixWarm, ng, suffixIns, false, false);
                                }
                            }
                            suffixWarm = null;
                        }
                        if (options.runRegexBaseline()) {
                            IPMIndexing regexWarm = new RegexIndex();
                            InsertStats regexIns = isSegments(options)
                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), regexWarm, 1)
                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), regexWarm, 1);
                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;
                                if (isSegments(options)) {
                                    SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, regexWarm, 1, regexIns, false, false);
                                } else {
                                    MultiQueryExperiment.runQueries(workloads, regexWarm, 1, regexIns, false, false);
                                }
                            }
                            regexWarm = null;
                        }
                    }

                    // Measured runs
                    for (int runIndex = 0; runIndex < options.runs(); runIndex++) {
                        if (!options.reinsertPerWorkload()) {
                            HBI hbi = null; InsertStats hbiIns = null;
                            if (options.runHbi()) {
                                hbi = IndexFactory.createHbi(
                                        options.windowLength(), options.treeLength(), options.alphabetSizeFor(ng), currentFp,
                                        options.runConfidence(), options.memPolicy(), ng);
                                hbi.strides = true;
                                hbi.stats().setCollecting(false);
                                hbi.stats().setExperimentMode(false);
                                hbiIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), hbi, ng)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, ng);
                            }
                            IPMIndexing suffix = null; InsertStats suffixIns = null;
                            if (options.runSuffix()) {
                                // Non-reinsert path: single index for all workloads. Use configured default delta.
                                suffix = IndexFactory.createSuffixIndex(
                                        options.windowLength());
                                suffixIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffix, ng)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, ng);
                            }
                            IPMIndexing regex = null; InsertStats regexIns = null;
                            if (options.runRegexBaseline()) {
                                regex = new RegexIndex();
                                regexIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), regex, 1)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), regex, 1);
                            }

                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;

                                MultiRunResult hbiResults = null;
                                if (hbi != null) {
                                    hbiResults = isSegments(options)
                                            ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, hbi, ng, hbiIns, false, false)
                                            : MultiQueryExperiment.runQueries(workloads, hbi, ng, hbiIns, false, false);
                                }
                                MultiRunResult suffixResults = null;
                                if (suffix != null) {
                                    suffixResults = isSegments(options)
                                            ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, suffix, ng, suffixIns, false, false)
                                            : MultiQueryExperiment.runQueries(workloads, suffix, ng, suffixIns, false, false);
                                }
                                MultiRunResult regexResults = null;
                                if (regex != null) {
                                    regexResults = isSegments(options)
                                            ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, regex, 1, regexIns, false, false)
                                            : MultiQueryExperiment.runQueries(workloads, regex, 1, regexIns, false, false);
                                }

                                accumulate(aggregated, activeIndexes, ng, type, hbiResults, suffixResults, regexResults);
                            }

                            hbi = null; suffix = null; regex = null;
                        } else {
                            // Reinsert only for delayed suffix index. HBI and Regex are built once per dataset.
                            HBI hbiSingle = null; InsertStats hbiSingleIns = null;
                            IPMIndexing regexSingle = null; InsertStats regexSingleIns = null;
                            if (options.runHbi()) {
                                hbiSingle = IndexFactory.createHbi(
                                        options.windowLength(), options.treeLength(), options.alphabetSizeFor(ng), currentFp,
                                        options.runConfidence(), options.memPolicy(), ng);
                                hbiSingle.strides = true;
                                hbiSingle.stats().setCollecting(false);
                                hbiSingle.stats().setExperimentMode(false);
                                hbiSingleIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), hbiSingle, ng)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), hbiSingle, ng);
                            }
                            if (options.runRegexBaseline()) {
                                regexSingle = new RegexIndex();
                                regexSingleIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), regexSingle, 1)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), regexSingle, 1);
                            }

                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;
                                for (QueryWorkload workload : workloads) {
                                    // Infer pattern length from query file name (e.g., "<len>.uniform.txt")
                                    int pl = BenchmarkIO.parsePatternLength(workload.queryFile().getFileName().toString());
                                    if (hbiSingle != null) {
                                        MultiRunResult res = isSegments(options)
                                                ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), hbiSingle, ng, hbiSingleIns, false, false)
                                                : MultiQueryExperiment.runQueries(List.of(workload), hbiSingle, ng, hbiSingleIns, false, false);
                                        addOne(aggregated, IndexType.HBI, ng, type, pl, res, workload);
                                    }
                                    if (options.runSuffix()) {
                                        // Reinsert per workload: set delta to inferred pattern length from query filename
                                        IPMIndexing suffix = IndexFactory.createSuffixIndex(
                                                options.windowLength());
                                        InsertStats ins = isSegments(options)
                                                ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffix, ng)
                                                : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, ng);
                                        MultiRunResult res = isSegments(options)
                                                ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), suffix, ng, ins, false, false)
                                                : MultiQueryExperiment.runQueries(List.of(workload), suffix, ng, ins, false, false);
                                        addOne(aggregated, IndexType.SUFFIX, ng, type, pl, res, workload);
                                    }
                                    if (regexSingle != null) {
                                        MultiRunResult res = isSegments(options)
                                                ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), regexSingle, 1, regexSingleIns, false, false)
                                                : MultiQueryExperiment.runQueries(List.of(workload), regexSingle, 1, regexSingleIns, false, false);
                                        addOne(aggregated, IndexType.REGEX, ng, type, pl, res, workload);
                                    }
                                }
                            }
                            hbiSingle = null; regexSingle = null;
                        }
                    }
                }
            }

            boolean anyResults = aggregated.values().stream()
                    .flatMap(m -> m.values().stream())
                    .flatMap(m2 -> m2.values().stream())
                    .anyMatch(Aggregation.AggregateStats::hasData);
            if (!anyResults) {
                System.out.println("No runs executed for FPR " + currentFp);
            } else {
                BenchmarkReporter.printConsoleSummary(aggregated, activeIndexes, queryTypes, currentFp);
                Path csvPath = BenchmarkReporter.writeCsv(options, aggregated, activeIndexes, queryTypes, "bs", currentFp);
                System.out.printf(Locale.ROOT, "Wrote aggregated results to %s%n", csvPath);
            }
        }
    }

    private static boolean isSegments(MultiBenchmarkOptions options) {
        return "segments".equalsIgnoreCase(options.mode());
    }

    private static void accumulate(Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated,
                                   List<IndexType> activeIndexes,
                                   int ng,
                                   QueryType type,
                                   MultiRunResult hbiRes,
                                   MultiRunResult suffixRes,
                                   MultiRunResult regexRes) {
        Map<Integer, Map<Integer, Aggregation.AggregateStats>> perType = aggregated.get(type);
        Map<Integer, Aggregation.AggregateStats> byPattern = perType.computeIfAbsent(ng, _k -> new TreeMap<>());
        if (hbiRes != null) hbiRes.results().forEach((wl, r) -> byPattern.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                .accumulate(IndexType.HBI, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs()));
        if (suffixRes != null) suffixRes.results().forEach((wl, r) -> byPattern.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                .accumulate(IndexType.SUFFIX, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs()));
        if (regexRes != null) regexRes.results().forEach((wl, r) -> byPattern.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                .accumulate(IndexType.REGEX, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs()));
    }

    private static void addOne(Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated,
                               IndexType indexType,
                               int ng,
                               QueryType type,
                               int pl,
                               MultiRunResult res,
                               QueryWorkload workload) {
        Map<Integer, Map<Integer, Aggregation.AggregateStats>> perType = aggregated.get(type);
        Map<Integer, Aggregation.AggregateStats> byPattern = perType.computeIfAbsent(ng, _k -> new TreeMap<>());
        var r = res.results().get(workload);
        byPattern.computeIfAbsent(pl, _k -> new Aggregation.AggregateStats())
                .accumulate(indexType, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs());
    }
}
