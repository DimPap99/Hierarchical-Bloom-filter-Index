package utilities;

import PMIndex.HBI;
import PMIndex.IPMIndexing;
import estimators.SelectiveFanout;
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

// Runs the outer benchmark loops across FPR, n-gram, and datasets.
public final class BenchmarkOrchestrator {
    private BenchmarkOrchestrator() {}
    public static void run(MultiBenchmarkOptions options) throws Exception {

        // Cache for Suffix baseline results so we don't rerun identical work
        // across FPR or n-gram loops when its configuration is stable.
        // We populate this during the first encountered (fp, ngram) combo and
        // reuse for subsequent iterations. Multiple 'runs' within the very first
        // (fp, ngram) still execute to preserve averaging semantics.
        Map<SuffixCacheKey, MultiRunResult> suffixResultsCache = new HashMap<>();
        Map<SuffixTreeCacheKey, MultiRunResult> suffixTreeResultsCache = new HashMap<>();

        List<QueryType> queryTypes = (options.queryType() != null)
                ? List.of(options.queryType())
                : List.of(QueryType.UNIFORM, QueryType.MISSING, QueryType.RARE);
        List<IndexType> activeIndexes = new ArrayList<>();
        if (options.runHbi()) activeIndexes.add(IndexType.HBI);
        if (options.runSuffix()) activeIndexes.add(IndexType.SUFFIX);
        if (options.runSuffixTreeBaseline()) activeIndexes.add(IndexType.SUFFIX_TREE);

        System.out.printf(Locale.ROOT,
                "Running window %s (%s) with datasets under %s%n",
                options.window(), options.queryTypeLabel(), options.dataRoot());

        List<Path> datasetDirs = BenchmarkIO.listDatasetDirectories(options.dataRoot());
        if (datasetDirs.isEmpty()) {
            throw new IllegalStateException("No dataset folders found under " + options.dataRoot());
        }

        int fpIndex = 0;
        for (double currentFp : options.fpRates()) {
            Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated = new EnumMap<>(QueryType.class);
            queryTypes.forEach(type -> aggregated.put(type, new TreeMap<>()));

            System.out.printf(Locale.ROOT, "%n==== FPR = %.6f ====%n", currentFp);

            int ngIndex = 0;
            for (int ng : options.ngrams()) {
                System.out.printf(Locale.ROOT, "---- ngram = %d ----%n", ng);

                int suffixNgForLoop = options.resolveSuffixNgram(ngIndex, ng);
                IndexFactory.setSuffixNgram(suffixNgForLoop);
                int suffixTreeNgForLoop = options.suffixTreeNgramFor(ngIndex, ng);
                IndexFactory.setSuffixTreeNgram(suffixTreeNgForLoop);

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
                    // Concise printout of effective settings for this loop (per FPR x ngram x dataset)
                    int effAlphabet = options.alphabetSizeFor(ng);
                    System.out.printf(Locale.ROOT,
                            "  Settings -> windowLen=%d, treeLen=%d, ngram=%d, alphabetBase=%d, alphabet=%d, suffixNgram=%d, suffixTreeNgram=%d, fpr=%.6f, runs=%d, algo=%s, strides=%b, mode=%s, policy=%s, reinsert=%b, eps=%.6f, deltaQ=%.3f, deltaSamp=%.3f, p=%.3f%n",
                            options.windowLength(), options.treeLength(), ng,
                            options.alphabetBase(), effAlphabet, suffixNgForLoop, suffixTreeNgForLoop, currentFp,
                            options.runs(), options.algorithm(), options.strides(), options.mode(), options.memPolicy(), options.reinsertPerWorkload(),
                            options.rankEpsTarget(), options.deltaQ(), options.deltaSamp(), options.quantile());

                    final double policyQuantile = options.quantile();
                    int policyBuckets = 0;
                    Utils.HopsDesignResult policyDesign = null;
                    if (options.memPolicy() != Utils.MemPolicy.NONE) {
                        int distinctEstimate = Math.max(1, options.alphabetSizeFor(ng));
                        policyDesign = Utils.designBucketsForRankTargetChebyshev(distinctEstimate,
                                options.rankEpsTarget(), options.deltaQ(), options.deltaSamp());
                        policyBuckets = policyDesign.suggestedBuckets;
                        System.out.printf(Locale.ROOT,
                                "  Policy auto-design -> Dhat=%d, buckets=%d, n_req=%d, LB=%d%s%n",
                                distinctEstimate,
                                policyDesign.suggestedBuckets,
                                policyDesign.requiredSampleSize,
                                policyDesign.occupancyLowerBound,
                                policyDesign.impossible ? " (target unattainable with Dhat)" : "");
                    }

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

                    workloadsByType = filterWorkloadsByMinLength(
                            workloadsByType,
                            options.minQueryLength(),
                            "configured minimum query length");
                    workloadsByType = filterWorkloadsByMaxLength(
                            workloadsByType,
                            options.maxQueryLength(),
                            "configured maximum query length");
                    workloadsByType = filterWorkloadsByMinLength(
                            workloadsByType,
                            ng,
                            "ngram");
                    boolean hasEligibleWorkloads = workloadsByType.values().stream().anyMatch(list -> list != null && !list.isEmpty());
                    if (!hasEligibleWorkloads) {
                        System.out.printf(Locale.ROOT,
                                "  No query workloads remain after filtering; skipping dataset %s for this configuration.%n",
                                datasetDir.getFileName());
                        continue;
                    }

                    // Warmups
                    for (int warmIndex = 0; warmIndex < options.warmupRuns(); warmIndex++) {
                        if (options.runHbi()) {
                            //testing to improve CF. Keep it as false to have the original behavior
                            boolean enableSelective = false;//(ng >= 6) || isCaidaWindow(options.window());
                            SelectiveFanout.setSelectiveRegimeEnabled(false);
                            HBI warm = IndexFactory.createHbi(
                                    options.windowLength(), options.treeLength(), options.alphabetSizeFor(ng), currentFp,
                                    options.runConfidence(), options.memPolicy(), ng, options.algorithm(), policyQuantile, policyBuckets);
                            warm.strides = options.strides();
                            warm.stats().setCollecting(true);
                            warm.stats().setExperimentMode(false);
                            InsertStats warmIns = isSegments(options)
                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), warm, ng)
                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), warm, ng);
                            forcePolicyIfActive(warm, options);
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
                            final int suffixNg = IndexFactory.getSuffixNgram();
                            // Use delayed SSWSI; warmups don't target a specific pattern length, use configured default delta
                            //dont use ngrams for suffix cause it hurts performance
                            IPMIndexing suffixWarm = IndexFactory.createSuffixIndex(
                                    options.windowLength(), options.alphabetSizeFor(suffixNg));
                            InsertStats suffixIns = isSegments(options)
                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffixWarm, suffixNg)
                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffixWarm, suffixNg);
                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;
                                if (isSegments(options)) {
                                    SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, suffixWarm, suffixNg, suffixIns, false, false);
                                } else {
                                    MultiQueryExperiment.runQueries(workloads, suffixWarm, suffixNg, suffixIns, false, false);
                                }
                            }
                            suffixWarm = null;
                        }
                        if (options.runSuffixTreeBaseline()) {
                            IPMIndexing suffixTreeWarm = IndexFactory.createSuffixTreeIndex(options.alphabetSizeFor(suffixTreeNgForLoop));
                            InsertStats suffixTreeIns = isSegments(options)
                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffixTreeWarm, suffixTreeNgForLoop)
                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffixTreeWarm, suffixTreeNgForLoop);
                            for (QueryType type : queryTypes) {
                                List<QueryWorkload> workloads = workloadsByType.get(type);
                                if (workloads == null || workloads.isEmpty()) continue;
                                if (isSegments(options)) {
                                    SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, suffixTreeWarm, suffixTreeNgForLoop, suffixTreeIns, false, false);
                                } else {
                                    MultiQueryExperiment.runQueries(workloads, suffixTreeWarm, suffixTreeNgForLoop, suffixTreeIns, false, false);
                                }
                            }
                            suffixTreeWarm = null;
                        }
                    }

                    // Measured runs
                    // For suffix reuse: during the first (fp, ng) combination only, we accumulate
                    // per-run suffix results so we can cache the averaged result for reuse.
                    Map<SuffixCacheKey, MRAccumulator> suffixAvgBuilders = new HashMap<>();
                    Map<SuffixTreeCacheKey, MRAccumulator> suffixTreeAvgBuilders = new HashMap<>();
                    for (int runIndex = 0; runIndex < options.runs(); runIndex++) {
                        if (!options.reinsertPerWorkload()) {
                            HBI hbi = null; InsertStats hbiIns = null;
                            if (options.runHbi()) {

                                boolean enableSelective = (ng >= 6) || isCaidaWindow(options.window());
                                SelectiveFanout.setSelectiveRegimeEnabled(enableSelective);
                                hbi = IndexFactory.createHbi(
                                        options.windowLength(), options.treeLength(), options.alphabetSizeFor(ng), currentFp,
                                        options.runConfidence(), options.memPolicy(), ng, options.algorithm(), policyQuantile, policyBuckets);
                                hbi.strides = options.strides();
                                hbi.stats().setCollecting(options.collectStats());
                                hbi.stats().setExperimentMode(false);
                                hbiIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), hbi, ng)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, ng);
                                forcePolicyIfActive(hbi, options);
                            }
                            // For Suffix baseline, decide whether to reuse cached results.
                            final int suffixNg = IndexFactory.getSuffixNgram();
                            final boolean canReuseSuffix = options.runSuffix() && options.reuseSuffixResults() && (fpIndex > 0 || ngIndex > 0);
                            IPMIndexing suffix = null; InsertStats suffixIns = null;
                            if (options.runSuffix() && !canReuseSuffix) {
                                // Build the suffix index only if we're not reusing results.
                                suffix = IndexFactory.createSuffixIndex(
                                        options.windowLength(), options.alphabetSizeFor(suffixNg));
                                suffixIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffix, suffixNg)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, suffixNg);
                            }
                            final boolean canReuseSuffixTree = options.runSuffixTreeBaseline() && options.reuseSuffixResults() && (fpIndex > 0 || ngIndex > 0);
                            IPMIndexing suffixTree = null; InsertStats suffixTreeIns = null;
                            if (options.runSuffixTreeBaseline() && !canReuseSuffixTree) {
                                suffixTree = IndexFactory.createSuffixTreeIndex(options.alphabetSizeFor(suffixTreeNgForLoop));
                                suffixTreeIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffixTree, suffixTreeNgForLoop)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffixTree, suffixTreeNgForLoop);
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
                                if (options.runSuffix()) {
                                    // Build the cache key per dataset + type workload set
                                    SuffixCacheKey cacheKey = SuffixCacheKey.forWorkloads(
                                            datasetFile.toString(), suffixNg, options.mode(), false, workloads);
                                    if (canReuseSuffix) {
                                        suffixResults = suffixResultsCache.get(cacheKey);
                                        if (suffixResults == null) {
                                            // Fallback: compute and cache now if missing
                                            IPMIndexing tmpSuffix = IndexFactory.createSuffixIndex(
                                                    options.windowLength(), options.alphabetSizeFor(suffixNg));
                                            InsertStats tmpIns = isSegments(options)
                                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), tmpSuffix, suffixNg)
                                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), tmpSuffix, suffixNg);
                                            suffixResults = isSegments(options)
                                                    ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, tmpSuffix, suffixNg, tmpIns, false, false)
                                                    : MultiQueryExperiment.runQueries(workloads, tmpSuffix, suffixNg, tmpIns, false, false);
                                            if (options.reuseSuffixResults()) suffixResultsCache.put(cacheKey, suffixResults);
                                        }
                                    } else if (suffix != null) {
                                        suffixResults = isSegments(options)
                                                ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, suffix, suffixNg, suffixIns, false, false)
                                                : MultiQueryExperiment.runQueries(workloads, suffix, suffixNg, suffixIns, false, false);
                                        // For averaging across runs, collect during first (fp,ng)
                                        if (options.reuseSuffixResults() && fpIndex == 0 && ngIndex == 0) {
                                            suffixAvgBuilders.computeIfAbsent(cacheKey, _k -> new MRAccumulator()).add(suffixResults);
                                        }
                                    }
                                }
                                // suffixResults averaged or reused per existing logic
                                MultiRunResult suffixTreeResults = null;
                                if (options.runSuffixTreeBaseline()) {
                                    SuffixTreeCacheKey treeKey = SuffixTreeCacheKey.forWorkloads(
                                            datasetFile.toString(), suffixTreeNgForLoop, options.mode(), false, workloads);
                                    if (canReuseSuffixTree) {
                                        suffixTreeResults = suffixTreeResultsCache.get(treeKey);
                                        if (suffixTreeResults == null) {
                                            IPMIndexing tmpTree = IndexFactory.createSuffixTreeIndex(options.alphabetSizeFor(suffixTreeNgForLoop));
                                            InsertStats tmpTreeIns = isSegments(options)
                                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), tmpTree, suffixTreeNgForLoop)
                                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), tmpTree, suffixTreeNgForLoop);
                                            suffixTreeResults = isSegments(options)
                                                    ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, tmpTree, suffixTreeNgForLoop, tmpTreeIns, false, false)
                                                    : MultiQueryExperiment.runQueries(workloads, tmpTree, suffixTreeNgForLoop, tmpTreeIns, false, false);
                                            if (options.reuseSuffixResults()) suffixTreeResultsCache.put(treeKey, suffixTreeResults);
                                        }
                                    } else if (suffixTree != null) {
                                        suffixTreeResults = isSegments(options)
                                                ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), workloads, suffixTree, suffixTreeNgForLoop, suffixTreeIns, false, false)
                                                : MultiQueryExperiment.runQueries(workloads, suffixTree, suffixTreeNgForLoop, suffixTreeIns, false, false);
                                        if (options.reuseSuffixResults() && fpIndex == 0 && ngIndex == 0) {
                                            suffixTreeAvgBuilders.computeIfAbsent(treeKey, _k -> new MRAccumulator()).add(suffixTreeResults);
                                        }
                                    }
                                    // suffixTreeResults averaged or reused per existing logic
                                }

                                // Aggregate all results (averages/reuse) for this run
                                accumulate(aggregated, activeIndexes, ng, type, hbiResults, suffixResults, suffixTreeResults);
                            }

                            hbi = null; suffix = null; suffixTree = null;
                        } else {
                            // Reinsert only for suffix-style baselines. HBI is built once per dataset.
                            HBI hbiSingle = null; InsertStats hbiSingleIns = null;
                            if (options.runHbi()) {
                                hbiSingle = IndexFactory.createHbi(
                                        options.windowLength(), options.treeLength(), options.alphabetSizeFor(ng), currentFp,
                                        options.runConfidence(), options.memPolicy(), ng, options.algorithm(), policyQuantile, policyBuckets);
                                hbiSingle.strides = options.strides();
                                // Honor --collect-stats/--stats even when not reinserting per workload
                                hbiSingle.stats().setCollecting(options.collectStats());
                                hbiSingle.stats().setExperimentMode(false);
                                hbiSingleIns = isSegments(options)
                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), hbiSingle, ng)
                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), hbiSingle, ng);
                                forcePolicyIfActive(hbiSingle, options);
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
                                        final int sfxNg = IndexFactory.getSuffixNgram();
                                        boolean canReuse = options.reuseSuffixResults() && (fpIndex > 0 || ngIndex > 0);
                                        // Build a cache key per single-workload when reinserting per workload
                                        SuffixCacheKey cacheKey = SuffixCacheKey.forWorkloads(
                                                datasetFile.toString(), sfxNg, options.mode(), true, List.of(workload));
                                        MultiRunResult res;
                                        if (canReuse) {
                                            res = suffixResultsCache.get(cacheKey);
                                            if (res == null) {
                                                IPMIndexing suffix = IndexFactory.createSuffixIndex(
                                                        options.windowLength(), options.alphabetSizeFor(sfxNg));
                                                InsertStats ins = isSegments(options)
                                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffix, sfxNg)
                                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, sfxNg);
                                                res = isSegments(options)
                                                        ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), suffix, sfxNg, ins, false, false)
                                                        : MultiQueryExperiment.runQueries(List.of(workload), suffix, sfxNg, ins, false, false);
                                                if (options.reuseSuffixResults()) suffixResultsCache.put(cacheKey, res);
                                            }
                                        } else {
                                            IPMIndexing suffix = IndexFactory.createSuffixIndex(
                                                    options.windowLength(), options.alphabetSizeFor(sfxNg));
                                            InsertStats ins = isSegments(options)
                                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffix, sfxNg)
                                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, sfxNg);
                                            res = isSegments(options)
                                                    ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), suffix, sfxNg, ins, false, false)
                                                    : MultiQueryExperiment.runQueries(List.of(workload), suffix, sfxNg, ins, false, false);
                                            if (options.reuseSuffixResults() && fpIndex == 0 && ngIndex == 0) {
                                                suffixAvgBuilders.computeIfAbsent(cacheKey, _k -> new MRAccumulator()).add(res);
                                            }
                                        }
                                        addOne(aggregated, IndexType.SUFFIX, sfxNg, type, pl, res, workload);
                                    }
                                    if (options.runSuffixTreeBaseline()) {
                                        boolean canReuseTree = options.reuseSuffixResults() && (fpIndex > 0 || ngIndex > 0);
                                        SuffixTreeCacheKey treeKey = SuffixTreeCacheKey.forWorkloads(
                                                datasetFile.toString(), suffixTreeNgForLoop, options.mode(), true, List.of(workload));
                                        MultiRunResult res;
                                        if (canReuseTree) {
                                            res = suffixTreeResultsCache.get(treeKey);
                                            if (res == null) {
                                                IPMIndexing tree = IndexFactory.createSuffixTreeIndex(options.alphabetSizeFor(suffixTreeNgForLoop));
                                                InsertStats ins = isSegments(options)
                                                        ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), tree, suffixTreeNgForLoop)
                                                        : MultiQueryExperiment.populateIndex(datasetFile.toString(), tree, suffixTreeNgForLoop);
                                                res = isSegments(options)
                                                        ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), tree, suffixTreeNgForLoop, ins, false, false)
                                                        : MultiQueryExperiment.runQueries(List.of(workload), tree, suffixTreeNgForLoop, ins, false, false);
                                                if (options.reuseSuffixResults()) suffixTreeResultsCache.put(treeKey, res);
                                            }
                                        } else {
                                            IPMIndexing tree = IndexFactory.createSuffixTreeIndex(options.alphabetSizeFor(suffixTreeNgForLoop));
                                            InsertStats ins = isSegments(options)
                                                    ? SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), tree, suffixTreeNgForLoop)
                                                    : MultiQueryExperiment.populateIndex(datasetFile.toString(), tree, suffixTreeNgForLoop);
                                            res = isSegments(options)
                                                    ? SegmentModeRunner.runQueriesSegments(datasetFile.toString(), List.of(workload), tree, suffixTreeNgForLoop, ins, false, false)
                                                    : MultiQueryExperiment.runQueries(List.of(workload), tree, suffixTreeNgForLoop, ins, false, false);
                                            if (options.reuseSuffixResults() && fpIndex == 0 && ngIndex == 0) {
                                                suffixTreeAvgBuilders.computeIfAbsent(treeKey, _k -> new MRAccumulator()).add(res);
                                            }
                                        }
                                        addOne(aggregated, IndexType.SUFFIX_TREE, suffixTreeNgForLoop, type, pl, res, workload);
                                    }
                                }
                            }
                            hbiSingle = null;
                        }
                        // If we are building cache averages, after finishing all runs for this dataset
                        // synthesize averaged results and store in the cache for reuse by later n-grams/FPRs.
                        if (options.runSuffix() && options.reuseSuffixResults() && fpIndex == 0 && ngIndex == 0
                                && runIndex == options.runs() - 1 && !suffixAvgBuilders.isEmpty()) {
                            for (Map.Entry<SuffixCacheKey, MRAccumulator> e : suffixAvgBuilders.entrySet()) {
                                suffixResultsCache.put(e.getKey(), e.getValue().buildAveraged());
                            }
                            suffixAvgBuilders.clear();
                        }
                        if (options.runSuffixTreeBaseline() && options.reuseSuffixResults() && fpIndex == 0 && ngIndex == 0
                                && runIndex == options.runs() - 1 && !suffixTreeAvgBuilders.isEmpty()) {
                            for (Map.Entry<SuffixTreeCacheKey, MRAccumulator> e : suffixTreeAvgBuilders.entrySet()) {
                                suffixTreeResultsCache.put(e.getKey(), e.getValue().buildAveraged());
                            }
                            suffixTreeAvgBuilders.clear();
                        }

                        // Encourage collection between iterations to reduce cross-run interference.
                        // This happens outside the measured sections and should not skew per-run timings.
                        System.gc();
                    }
                    // No special best-run aggregation; averages/reuse already reflected via MRAccumulator and caches.
                }
                ngIndex++;
            }
            fpIndex++;

            boolean anyResults = aggregated.values().stream()
                    .flatMap(m -> m.values().stream())
                    .flatMap(m2 -> m2.values().stream())
                    .anyMatch(Aggregation.AggregateStats::hasData);
            if (!anyResults) {
                System.out.println("No runs executed for FPR " + currentFp);
            } else {
                BenchmarkReporter.printConsoleSummary(aggregated, activeIndexes, queryTypes, currentFp);
                Path csvPath = BenchmarkReporter.writeCsv(options, aggregated, activeIndexes, queryTypes, options.algorithm(), currentFp);
                System.out.printf(Locale.ROOT, "Wrote aggregated results to %s%n", csvPath);
            }
        }
    }

    private static Map<QueryType, List<QueryWorkload>> filterWorkloadsByMinLength(
            Map<QueryType, List<QueryWorkload>> workloadsByType,
            int minThreshold,
            String thresholdLabel) {
        if (minThreshold <= 0) {
            return workloadsByType;
        }
        Map<QueryType, List<QueryWorkload>> filtered = new EnumMap<>(QueryType.class);
        for (Map.Entry<QueryType, List<QueryWorkload>> entry : workloadsByType.entrySet()) {
            QueryType type = entry.getKey();
            List<QueryWorkload> workloads = entry.getValue();
            if (workloads == null || workloads.isEmpty()) {
                filtered.put(type, List.of());
                continue;
            }
            List<QueryWorkload> eligible = workloads.stream()
                    .filter(workload -> workload.patternLength() >= minThreshold)
                    .collect(Collectors.toList());
            if (eligible.size() != workloads.size()) {
                System.out.printf(Locale.ROOT,
                        "  Query type %s: ignoring %d workloads (pattern length < %s %d)%n",
                        type.fileToken(), workloads.size() - eligible.size(), thresholdLabel, minThreshold);
            }
            filtered.put(type, eligible);
        }
        return filtered;
    }

    private static Map<QueryType, List<QueryWorkload>> filterWorkloadsByMaxLength(
            Map<QueryType, List<QueryWorkload>> workloadsByType,
            int maxThreshold,
            String thresholdLabel) {
        if (maxThreshold <= 0) {
            return workloadsByType;
        }
        Map<QueryType, List<QueryWorkload>> filtered = new EnumMap<>(QueryType.class);
        for (Map.Entry<QueryType, List<QueryWorkload>> entry : workloadsByType.entrySet()) {
            QueryType type = entry.getKey();
            List<QueryWorkload> workloads = entry.getValue();
            if (workloads == null || workloads.isEmpty()) {
                filtered.put(type, List.of());
                continue;
            }
            List<QueryWorkload> eligible = workloads.stream()
                    .filter(workload -> workload.patternLength() <= maxThreshold)
                    .collect(Collectors.toList());
            if (eligible.size() != workloads.size()) {
                System.out.printf(Locale.ROOT,
                        "  Query type %s: ignoring %d workloads (pattern length > %s %d)%n",
                        type.fileToken(), workloads.size() - eligible.size(), thresholdLabel, maxThreshold);
            }
            filtered.put(type, eligible);
        }
        return filtered;
    }

    private static boolean isCaidaWindow(String windowToken) {
        if (windowToken == null) return false;
        String w = windowToken.toLowerCase(java.util.Locale.ROOT);
        return w.contains("caida");
    }

    private static boolean isSegments(MultiBenchmarkOptions options) {
        return "segments".equalsIgnoreCase(options.mode());
    }

    private static void forcePolicyIfActive(HBI hbi, MultiBenchmarkOptions options) {
        if (hbi != null && options.memPolicy() != Utils.MemPolicy.NONE) {
//            hbi.forceApplyMemoryPolicy();
        }
    }

    private static void accumulate(Map<QueryType, Map<Integer, Map<Integer, Aggregation.AggregateStats>>> aggregated,
                                   List<IndexType> activeIndexes,
                                   int ng,
                                   QueryType type,
                                   MultiRunResult hbiRes,
                                   MultiRunResult suffixRes,
                                   MultiRunResult suffixTreeRes) {
        Map<Integer, Map<Integer, Aggregation.AggregateStats>> perType = aggregated.get(type);
        // HBI aggregated under its actual ng
        if (hbiRes != null) {
            Map<Integer, Aggregation.AggregateStats> byPatternHbi = perType.computeIfAbsent(ng, _k -> new TreeMap<>());
            hbiRes.results().forEach((wl, r) -> byPatternHbi.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                    .accumulate(IndexType.HBI, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen()));
        }
        // Suffix baseline: also aggregate under the active n-gram so it prints next to HBI
        // This ensures summaries like "Pattern X" show both HBI and SuffixTree together.
        if (suffixRes != null) {
            // Aggregate under current n-gram group for reporting
            Map<Integer, Aggregation.AggregateStats> byPatternSuffixAtNg = perType.computeIfAbsent(ng, _k -> new TreeMap<>());
            suffixRes.results().forEach((wl, r) -> byPatternSuffixAtNg.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                    .accumulate(IndexType.SUFFIX, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen()));

            // Retain aggregation under the suffix's configured n-gram for completeness
            int suffixNg = IndexFactory.getSuffixNgram();
            if (suffixNg != ng) {
                Map<Integer, Aggregation.AggregateStats> byPatternSuffixStrict = perType.computeIfAbsent(suffixNg, _k -> new TreeMap<>());
                suffixRes.results().forEach((wl, r) -> byPatternSuffixStrict.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                        .accumulate(IndexType.SUFFIX, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen()));
            }
        }
        if (suffixTreeRes != null) {
            Map<Integer, Aggregation.AggregateStats> byPatternSuffixTree = perType.computeIfAbsent(ng, _k -> new TreeMap<>());
            suffixTreeRes.results().forEach((wl, r) -> byPatternSuffixTree.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                    .accumulate(IndexType.SUFFIX_TREE, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen()));
            int suffixTreeNg = IndexFactory.getSuffixTreeNgram();
            if (suffixTreeNg != ng) {
                Map<Integer, Aggregation.AggregateStats> byPatternSuffixTreeStrict = perType.computeIfAbsent(suffixTreeNg, _k -> new TreeMap<>());
                suffixTreeRes.results().forEach((wl, r) -> byPatternSuffixTreeStrict.computeIfAbsent(wl.patternLength(), _k -> new Aggregation.AggregateStats())
                        .accumulate(IndexType.SUFFIX_TREE, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen()));
            }
        }
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
                .accumulate(indexType, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen());

        // Special-case: if we are adding Suffix results and ng differs from its configured n-gram,
        // also keep a copy under the suffix's native n-gram for completeness.
        if (indexType == IndexType.SUFFIX) {
            int suffixNg = IndexFactory.getSuffixNgram();
            if (suffixNg != ng) {
                Map<Integer, Aggregation.AggregateStats> byPatternSuffixStrict = perType.computeIfAbsent(suffixNg, _k -> new TreeMap<>());
                byPatternSuffixStrict.computeIfAbsent(pl, _k -> new Aggregation.AggregateStats())
                        .accumulate(indexType, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen());
            }
        } else if (indexType == IndexType.SUFFIX_TREE) {
            int suffixTreeNg = IndexFactory.getSuffixTreeNgram();
            if (suffixTreeNg != ng) {
                Map<Integer, Aggregation.AggregateStats> byPatternSuffixTreeStrict = perType.computeIfAbsent(suffixTreeNg, _k -> new TreeMap<>());
                byPatternSuffixTreeStrict.computeIfAbsent(pl, _k -> new Aggregation.AggregateStats())
                        .accumulate(indexType, r.avgInsertMsPerSymbol(), r.totalInsertTimeMs(), r.totalRunTimeMs(), r.avgLpMs(), r.avgCfLpMs(), r.avgLpChosen(), r.avgCfLpChosen());
            }
        }
    }
}

// Local cache key for suffix results reuse
record SuffixCacheKey(String dataset,
                      int suffixNgram,
                      String mode,
                      boolean reinsertPerWorkload,
                      String workloadKey) {

    static SuffixCacheKey forWorkloads(String dataset,
                                       int suffixNgram,
                                       String mode,
                                       boolean reinsertPerWorkload,
                                       List<QueryWorkload> workloads) {
        // Build a stable key from the query file paths in order
        String joined = workloads.stream()
                .map(w -> w.queryFile().toString())
                .collect(Collectors.joining("|"));
        return new SuffixCacheKey(dataset, suffixNgram, mode, reinsertPerWorkload, joined);
    }
}

record SuffixTreeCacheKey(String dataset,
                          int ngram,
                          String mode,
                          boolean reinsertPerWorkload,
                          String workloadKey) {

    static SuffixTreeCacheKey forWorkloads(String dataset,
                                           int ngram,
                                           String mode,
                                           boolean reinsertPerWorkload,
                                           List<QueryWorkload> workloads) {
        String joined = workloads.stream()
                .map(w -> w.queryFile().toString())
                .collect(Collectors.joining("|"));
        return new SuffixTreeCacheKey(dataset, ngram, mode, reinsertPerWorkload, joined);
    }
}

// Helper to average MultiRunResult across multiple runs for reuse.
final class MRAccumulator {
    private static final class Sum {
        double sumAvgInsertMsPerSymbol;
        double sumTotalInsertMs;
        double sumTotalQueryMs;
        double sumAvgQuerySize;
        double sumAvgLpMs;
        double sumAvgCfLpMs;
        double sumAvgLpChosen;
        double sumAvgCfLpChosen;
        int count;
    }

    private final Map<QueryWorkload, Sum> sums = new LinkedHashMap<>();

    void add(MultiRunResult res) {
        for (Map.Entry<QueryWorkload, ExperimentRunResult> e : res.results().entrySet()) {
            QueryWorkload wl = e.getKey();
            ExperimentRunResult r = e.getValue();
            Sum s = sums.computeIfAbsent(wl, _k -> new Sum());
            s.sumAvgInsertMsPerSymbol += r.avgInsertMsPerSymbol();
            s.sumTotalInsertMs += r.totalInsertTimeMs();
            s.sumTotalQueryMs += r.totalRunTimeMs();
            s.sumAvgQuerySize += r.avgQuerySize();
            s.sumAvgLpMs += r.avgLpMs();
            s.sumAvgCfLpMs += r.avgCfLpMs();
            s.sumAvgLpChosen += r.avgLpChosen();
            s.sumAvgCfLpChosen += r.avgCfLpChosen();
            s.count++;
        }
    }

    MultiRunResult buildAveraged() {
        Map<QueryWorkload, ExperimentRunResult> avg = new LinkedHashMap<>();
        for (Map.Entry<QueryWorkload, Sum> e : sums.entrySet()) {
            Sum s = e.getValue();
            double c = Math.max(1, s.count);
            double avgInsSym = s.sumAvgInsertMsPerSymbol / c;
            double avgInsMs = s.sumTotalInsertMs / c;
            double avgQueryMs = s.sumTotalQueryMs / c;
            double avgQSize = s.sumAvgQuerySize / c;
            avg.put(e.getKey(), new ExperimentRunResult(
                    avgQueryMs,
                    avgInsMs,
                    null,
                    avgQSize,
                    avgInsSym,
                    avgQueryMs,
                    s.sumAvgLpMs / c,
                    s.sumAvgCfLpMs / c,
                    s.sumAvgLpChosen / c,
                    s.sumAvgCfLpChosen / c,
                    null
            ));
        }
        return new MultiRunResult(avg);
    }
}
