import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import PMIndex.SuffixTreeIndex;
import estimators.CSEstimator;
import estimators.CostFunctionDefaultRoot;
import estimators.CostFunctionMaxProb;
import estimators.Estimator;
import membership.BloomFilter;
import membership.Membership;
import search.*;

import utilities.MultiQueryExperiment;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.EnumMap;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import utilities.CsvUtil;

public final class HBIDatasetBenchmarkPerDataset {

    private static final boolean USE_STRIDES = true;
    private static String algo = "root";

    private record BenchmarkOptions(Path dataRoot,
                                    Path queryRoot,
                                    String window,
                                    QueryType queryType,
                                    int ngram,
                                    int windowLength,
                                    int treeLength,
                                    int alphabetBase,
                                    double fpRate,
                                    double runConfidence,
                                    int warmupRuns,
                                    int runs,
                                    boolean runRegexBaseline,
                                    boolean runHbi,
                                    boolean runSuffix) {

        static BenchmarkOptions parse(String[] args) {
            Path dataRoot = Path.of("data");
            Path queryRoot = Path.of("queries");
            String window = "wzipf20_e1";
            QueryType queryType = QueryType.UNIFORM;
            Integer ngram = 2;
            Integer windowLength = null;
            Integer treeLength = null;
            Integer alphabetBase = 128;
            double fpRate = 0.001;
            double runConfidence = 0.99;
            int warmupRuns = 1;
            int runs = 3;
            boolean runRegexBaseline = false;
            boolean runHbi = true;
            boolean runSuffix = false;

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (!arg.startsWith("--")) {
                    continue;
                }
                String key;
                String value;
                int eq = arg.indexOf('=');
                if (eq >= 0) {
                    key = arg.substring(2, eq);
                    value = arg.substring(eq + 1);
                } else {
                    key = arg.substring(2);
                    if (i + 1 >= args.length) {
                        throw new IllegalArgumentException("Missing value for option --" + key);
                    }
                    value = args[++i];
                }
                switch (key) {
                    case "window" -> window = value;
                    case "type" -> queryType = "all".equalsIgnoreCase(value)
                            ? null
                            : QueryType.fromString(value);
                    case "data-root" -> dataRoot = Path.of(value);
                    case "query-root" -> queryRoot = Path.of(value);
                    case "ngrams" -> ngram = Integer.parseInt(value);
                    case "window-length" -> windowLength = Integer.parseInt(value);
                    case "tree-length" -> treeLength = Integer.parseInt(value);
                    case "alphabet-base" -> alphabetBase = Integer.parseInt(value);
                    case "fp" -> fpRate = Double.parseDouble(value);
                    case "confidence" -> runConfidence = Double.parseDouble(value);
                    case "warmup" -> warmupRuns = Integer.parseInt(value);
                    case "runs" -> runs = Integer.parseInt(value);
                    case "regex", "run-regex" -> runRegexBaseline = Boolean.parseBoolean(value);
                    case "run-hbi" -> runHbi = Boolean.parseBoolean(value);
                    case "run-suffix" -> runSuffix = Boolean.parseBoolean(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            if (window == null) {
                throw new IllegalArgumentException("--window is required");
            }

            if (ngram == null) {
                ngram = 4;
            }
            if (windowLength == null) {
                windowLength = deriveWindowLength(window, 1 << 20);
            }
            if (treeLength == null) {
                treeLength = windowLength;
            }

            Path resolvedDataRoot = dataRoot.resolve(window);
            Path resolvedQueryRoot = queryRoot.resolve(window);

            if (!Files.isDirectory(resolvedDataRoot)) {
                throw new IllegalArgumentException("Data directory not found: " + resolvedDataRoot);
            }
            if (!Files.isDirectory(resolvedQueryRoot)) {
                throw new IllegalArgumentException("Query directory not found: " + resolvedQueryRoot);
            }

            if (!runRegexBaseline && !runHbi && !runSuffix) {
                throw new IllegalArgumentException("At least one index must be enabled (HBI, suffix, or regex).");
            }

            return new BenchmarkOptions(
                    resolvedDataRoot,
                    resolvedQueryRoot,
                    window,
                    queryType,
                    ngram,
                    windowLength,
                    treeLength,
                    alphabetBase,
                    fpRate,
                    runConfidence,
                    warmupRuns,
                    runs,
                    runRegexBaseline,
                    runHbi,
                    runSuffix);
        }

        int alphabetSize() {
            double pow = Math.pow(alphabetBase, ngram);

            double sigma = Math.min(pow, windowLength);
            if (sigma >= Integer.MAX_VALUE) {
                return Integer.MAX_VALUE;
            }
            return (int) sigma;
        }

        String queryTypeLabel() {
            return queryType != null ? queryType.fileToken() : "all";
        }

        List<QueryType> orderedQueryTypes() {
            EnumSet<QueryType> selected = queryType != null
                    ? EnumSet.of(queryType)
                    : EnumSet.allOf(QueryType.class);

            List<QueryType> preferredOrder = List.of(QueryType.UNIFORM, QueryType.MISSING, QueryType.RARE);
            List<QueryType> entries = new ArrayList<>();
            for (QueryType candidate : preferredOrder) {
                if (selected.contains(candidate)) {
                    entries.add(candidate);
                }
            }
            return entries;
        }

        List<IndexType> activeIndexes() {
            List<IndexType> indexes = new ArrayList<>(3);
            if (runHbi) {
                indexes.add(IndexType.HBI);
            }
            if (runSuffix) {
                indexes.add(IndexType.SUFFIX);
            }
            if (runRegexBaseline) {
                indexes.add(IndexType.REGEX);
            }
            return indexes;
        }

        Path csvOutputFile() {
            String fileName = "%s_%s_%d_%s.csv".formatted(window, queryTypeLabel(), ngram, algo);
            return Path.of(fileName);
        }

        static int deriveWindowLength(String window, int defaultValue) {
            if (window == null || window.isEmpty()) {
                return defaultValue;
            }
            for (int i = 0; i < window.length(); i++) {
                if (Character.isDigit(window.charAt(i))) {
                    String numeric = window.substring(i);
                    try {
                        int power = Integer.parseInt(numeric);
                        return 1 << power;
                    } catch (NumberFormatException ignored) {
                        return defaultValue;
                    }
                }
            }
            return defaultValue;
        }
    }

    public static void main(String[] args) throws IOException {
        BenchmarkOptions options = BenchmarkOptions.parse(args);
        List<QueryType> queryTypes = options.orderedQueryTypes();
        List<IndexType> activeIndexes = options.activeIndexes();

        // New aggregation: per dataset number -> per query type -> per pattern length
        Map<Integer, Map<QueryType, Map<Integer, AggregateStats>>> aggregatedByDataset = new TreeMap<>();

        System.out.printf(Locale.ROOT,
                "Running window %s (%s) with datasets under %s%n",
                options.window(), options.queryTypeLabel(), options.dataRoot());

        List<Path> datasetDirs = listDatasetDirectories(options.dataRoot());
        if (datasetDirs.isEmpty()) {
            throw new IllegalStateException("No dataset folders found under " + options.dataRoot());
        }

        for (int dsIdx = 0; dsIdx < datasetDirs.size(); dsIdx++) {
            int datasetNumber = dsIdx + 1; // one-based for reporting
            Path datasetDir = datasetDirs.get(dsIdx);

            // Initialize buckets for this dataset
            Map<QueryType, Map<Integer, AggregateStats>> perTypeMap = new EnumMap<>(QueryType.class);
            for (QueryType qt : queryTypes) {
                perTypeMap.put(qt, new TreeMap<>());
            }
            aggregatedByDataset.put(datasetNumber, perTypeMap);

            Path datasetFile = findSingleDatasetFile(datasetDir);
            Path queryDir = options.queryRoot().resolve(datasetDir.getFileName());
            if (!Files.isDirectory(queryDir)) {
                System.out.printf(Locale.ROOT,
                        "Skipping dataset %s (missing query folder %s)%n",
                        datasetDir.getFileName(), queryDir);
                continue;
            }

            System.out.printf(Locale.ROOT, "Dataset #%d %s -> %s%n",
                    datasetNumber, datasetDir.getFileName(), datasetFile.getFileName());

            Map<QueryType, List<MultiQueryExperiment.QueryWorkload>> workloadsByType = new EnumMap<>(QueryType.class);
            boolean hasQueries = false;
            for (QueryType type : queryTypes) {
                List<Path> queryFiles = findQueryFiles(queryDir, type);
                if (queryFiles.isEmpty()) {
                    System.out.printf(Locale.ROOT,
                            "No %s queries in %s, skipping.%n",
                            type.fileToken(), queryDir);
                    workloadsByType.put(type, List.of());
                    continue;
                }
                hasQueries = true;
                System.out.printf(Locale.ROOT, "  Query type %s%n", type.fileToken());
                List<MultiQueryExperiment.QueryWorkload> workloads = queryFiles.stream()
                        .map(path -> new MultiQueryExperiment.QueryWorkload(
                                parsePatternLength(path.getFileName().toString()),
                                path))
                        .sorted(Comparator.comparingInt(MultiQueryExperiment.QueryWorkload::patternLength))
                        .collect(Collectors.toList());
                workloadsByType.put(type, workloads);
            }

            if (!hasQueries) {
                continue;
            }

            // Warmups
            for (int warmIndex = 0; warmIndex < options.warmupRuns(); warmIndex++) {
                if (options.runHbi()) {
                    HBI warmupIndex = newHbi(options);
                    warmupIndex.strides = USE_STRIDES;
                    warmupIndex.stats().setCollecting(false);
                    warmupIndex.stats().setExperimentMode(false);
                    MultiQueryExperiment.InsertStats warmupInsert =
                            MultiQueryExperiment.populateIndex(datasetFile.toString(), warmupIndex, options.ngram());
                    for (QueryType type : queryTypes) {
                        List<MultiQueryExperiment.QueryWorkload> workloads = workloadsByType.get(type);
                        if (workloads == null || workloads.isEmpty()) {
                            continue;
                        }
                        MultiQueryExperiment.runQueries(workloads, warmupIndex, options.ngram(), warmupInsert, false, false);
                    }
                }

                if (options.runSuffix()) {
                    SuffixTreeIndex suffixWarm = newSuffixTree(options);
                    MultiQueryExperiment.InsertStats suffixWarmInsert =
                            MultiQueryExperiment.populateIndex(datasetFile.toString(), suffixWarm, options.ngram());
                    for (QueryType type : queryTypes) {
                        List<MultiQueryExperiment.QueryWorkload> workloads = workloadsByType.get(type);
                        if (workloads == null || workloads.isEmpty()) {
                            continue;
                        }
                        MultiQueryExperiment.runQueries(workloads, suffixWarm, options.ngram(), suffixWarmInsert, false, false);
                    }
                }

                if (options.runRegexBaseline()) {
                    IPMIndexing regexWarm = new RegexIndex();
                    MultiQueryExperiment.InsertStats regexWarmInsert =
                            MultiQueryExperiment.populateIndex(datasetFile.toString(), regexWarm, 1);
                    for (QueryType type : queryTypes) {
                        List<MultiQueryExperiment.QueryWorkload> workloads = workloadsByType.get(type);
                        if (workloads == null || workloads.isEmpty()) {
                            continue;
                        }
                        MultiQueryExperiment.runQueries(workloads, regexWarm, 1, regexWarmInsert, false, false);
                    }
                }
            }

            // Measured runs for this dataset
            for (int runIndex = 0; runIndex < options.runs(); runIndex++) {
                HBI hbi = null;
                MultiQueryExperiment.InsertStats hbiInsertStats = null;
                if (options.runHbi()) {
                    hbi = newHbi(options);
                    hbi.strides = USE_STRIDES;
                    hbi.stats().setCollecting(false);
                    hbi.stats().setExperimentMode(false);
                    hbiInsertStats = MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, options.ngram());
                }

                IPMIndexing suffix = null;
                MultiQueryExperiment.InsertStats suffixInsertStats = null;
                if (options.runSuffix()) {
                    suffix = newSuffixTree(options);
                    suffixInsertStats = MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, options.ngram());
                }

                IPMIndexing regex = null;
                MultiQueryExperiment.InsertStats regexInsertStats = null;
                if (options.runRegexBaseline()) {
                    regex = new RegexIndex();
                    regexInsertStats = MultiQueryExperiment.populateIndex(datasetFile.toString(), regex, 1);
                }

                for (QueryType type : queryTypes) {
                    List<MultiQueryExperiment.QueryWorkload> workloads = workloadsByType.get(type);
                    if (workloads == null || workloads.isEmpty()) {
                        continue;
                    }

                    MultiQueryExperiment.MultiRunResult hbiResults = null;
                    if (hbi != null) {
                        hbiResults = MultiQueryExperiment.runQueries(workloads, hbi, options.ngram(), hbiInsertStats, false, false);
                    }

                    MultiQueryExperiment.MultiRunResult suffixResults = null;
                    if (suffix != null) {
                        suffixResults = MultiQueryExperiment.runQueries(workloads, suffix, options.ngram(), suffixInsertStats, false, false);
                    }

                    MultiQueryExperiment.MultiRunResult regexResults = null;
                    if (regex != null) {
                        regexResults = MultiQueryExperiment.runQueries(workloads, regex, 1, regexInsertStats, false, false);
                    }

                    Map<Integer, AggregateStats> perTypeAggregated = aggregatedByDataset.get(datasetNumber).get(type);

                    if (hbiResults != null) {
                        hbiResults.results().forEach((workload, result) -> {
                            int patternLength = workload.patternLength();
                            perTypeAggregated.computeIfAbsent(patternLength, ignored -> new AggregateStats())
                                    .accumulate(IndexType.HBI,
                                            result.avgInsertMsPerSymbol(),
                                            result.totalInsertTimeMs(),
                                            result.totalRunTimeMs());
                        });
                    }

                    if (suffixResults != null) {
                        suffixResults.results().forEach((workload, result) -> {
                            int patternLength = workload.patternLength();
                            perTypeAggregated.computeIfAbsent(patternLength, ignored -> new AggregateStats())
                                    .accumulate(IndexType.SUFFIX,
                                            result.avgInsertMsPerSymbol(),
                                            result.totalInsertTimeMs(),
                                            result.totalRunTimeMs());
                        });
                    }

                    if (regexResults != null) {
                        regexResults.results().forEach((workload, result) -> {
                            int patternLength = workload.patternLength();
                            perTypeAggregated.computeIfAbsent(patternLength, ignored -> new AggregateStats())
                                    .accumulate(IndexType.REGEX,
                                            result.avgInsertMsPerSymbol(),
                                            result.totalInsertTimeMs(),
                                            result.totalRunTimeMs());
                        });
                    }
                }
            }
        }

        boolean anyResults = aggregatedByDataset.values().stream()
                .flatMap(map -> map.values().stream())
                .flatMap(m -> m.values().stream())
                .anyMatch(AggregateStats::hasData);

        if (!anyResults) {
            System.out.println("No runs executed.");
            return;
        }

        System.out.println();

        // Human readable summaries per dataset
        for (Map.Entry<Integer, Map<QueryType, Map<Integer, AggregateStats>>> dsEntry : aggregatedByDataset.entrySet()) {
            int datasetNumber = dsEntry.getKey();
            Map<QueryType, Map<Integer, AggregateStats>> perType = dsEntry.getValue();

            System.out.printf(Locale.ROOT, "Aggregated averages for dataset #%d:%n", datasetNumber);

            for (QueryType type : perType.keySet()) {
                Map<Integer, AggregateStats> statsByPattern = perType.get(type);
                if (statsByPattern == null || statsByPattern.isEmpty()) {
                    System.out.printf(Locale.ROOT, "  %s -> No data available.%n", type.fileToken());
                    continue;
                }
                System.out.printf(Locale.ROOT, "  Query type %s%n", type.fileToken());

                statsByPattern.forEach((patternLength, stats) -> {
                    System.out.printf(Locale.ROOT, "    Pattern %d%n", patternLength);
                    for (IndexType index : activeIndexes) {
                        StatsSnapshot snapshot = stats.snapshot(index);
                        if (snapshot == null || snapshot.count() == 0) {
                            System.out.printf(Locale.ROOT, "      %s -> no data%n", index.displayName());
                            continue;
                        }
                        System.out.printf(Locale.ROOT,
                                "      %s -> avgInsert=%.6f ms, insert=%.3f ms, query=%.3f ms over %d run(s)%n",
                                index.displayName(),
                                snapshot.avgInsertMsPerSymbol(),
                                snapshot.avgInsertMs(),
                                snapshot.avgQueryMs(),
                                snapshot.count());
                    }
                });
            }
            System.out.println();
        }

        // CSV export per dataset
        List<Object> header = new ArrayList<>();
        header.add("dataset");
        header.add("patternLength");

        for (QueryType type : queryTypes) {
            for (IndexType index : activeIndexes) {
                String suffix = type.fileToken() + "_" + index.csvLabel();
                header.add("avgInsertMsPerSymbol_" + suffix);
                header.add("insertMs_" + suffix);
                header.add("avgMs_" + suffix);
                header.add("datasetCount_" + suffix); // now equals number of runs contributing for this dataset
            }
        }
        header.add("algorithm");

        List<List<?>> csvRows = new ArrayList<>();
        csvRows.add(header);

        for (Map.Entry<Integer, Map<QueryType, Map<Integer, AggregateStats>>> dsEntry : aggregatedByDataset.entrySet()) {
            int datasetNumber = dsEntry.getKey();
            Map<QueryType, Map<Integer, AggregateStats>> perType = dsEntry.getValue();

            // Collect all pattern lengths present in this dataset
            SortedSet<Integer> patternLengths = new TreeSet<>();
            perType.values().forEach(map -> patternLengths.addAll(map.keySet()));

            for (Integer patternLength : patternLengths) {
                List<Object> row = new ArrayList<>();
                row.add(datasetNumber);
                row.add(patternLength);

                for (QueryType type : queryTypes) {
                    Map<Integer, AggregateStats> statsByPattern = perType.getOrDefault(type, Map.of());
                    AggregateStats stats = statsByPattern.get(patternLength);

                    for (IndexType index : activeIndexes) {
                        StatsSnapshot snapshot = (stats != null) ? stats.snapshot(index) : null;
                        if (snapshot == null || snapshot.count() == 0) {
                            row.add(null);
                            row.add(null);
                            row.add(null);
                            row.add(0);
                        } else {
                            row.add(snapshot.avgInsertMsPerSymbol());
                            row.add(snapshot.avgInsertMs());
                            row.add(snapshot.avgQueryMs());
                            row.add(snapshot.count());
                        }
                    }
                }
                row.add(algo);
                csvRows.add(row);
            }
        }

        Path csvPath = options.csvOutputFile();
        CsvUtil.writeRows(csvPath, csvRows);
        System.out.printf(Locale.ROOT, "Wrote per-dataset results to %s%n", csvPath);
    }

    private static HBI newHbi(BenchmarkOptions options) {
        int alphabetSize = options.alphabetSize();
        Supplier<Estimator> estFactory = () -> new CSEstimator(options.treeLength(), 5, 16384);
        Supplier<Membership> memFactory = BloomFilter::new;
        Supplier<PruningPlan> prFactory = () -> new MostFreqPruning(options.runConfidence());
        Verifier verifier = new VerifierLinearLeafProbe();

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(options.windowLength())
                .fpRate(options.fpRate())
                .alphabetSize(alphabetSize)
                .treeLength(options.treeLength())
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(verifier)
                .costFunction(new CostFunctionDefaultRoot())
                .confidence(options.runConfidence())
                .experimentMode(false)
                .collectStats(false)
                .nGram(options.ngram())
                .build();

        return new HBI(configuration);
    }

    private static SuffixTreeIndex newSuffixTree(BenchmarkOptions options) {
        long expectedDistinct = Math.max(1L, (long) options.windowLength());
        double epsilon = 0.0001;
        int totalTokens = options.treeLength();
        return new SuffixTreeIndex(expectedDistinct, epsilon, totalTokens);
    }

    private static List<Path> listDatasetDirectories(Path root) throws IOException {
        try (Stream<Path> stream = Files.list(root)) {
            return stream
                    .filter(Files::isDirectory)
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), HBIDatasetBenchmarkPerDataset::compareNaturally))
                    .collect(Collectors.toList());
        }
    }

    private static Path findSingleDatasetFile(Path datasetDir) throws IOException {
        try (Stream<Path> stream = Files.list(datasetDir)) {
            return stream
                    .filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().endsWith(".txt"))
                    .findFirst()
                    .orElseThrow(() -> new IllegalStateException(
                            "No dataset .txt file found under " + datasetDir));
        }
    }

    private static List<Path> findQueryFiles(Path queryDir, QueryType type) throws IOException {
        String suffix = "." + type.fileToken() + ".txt";
        try (Stream<Path> stream = Files.list(queryDir)) {
            return stream
                    .filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().endsWith(suffix))
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), HBIDatasetBenchmarkPerDataset::compareNaturally))
                    .collect(Collectors.toList());
        }
    }

    private static int parsePatternLength(String fileName) {
        int dot = fileName.indexOf('.');
        if (dot < 0) {
            throw new IllegalArgumentException("Cannot determine pattern length from " + fileName);
        }
        String numeric = fileName.substring(0, dot);
        try {
            return Integer.parseInt(numeric);
        } catch (NumberFormatException ex) {
            throw new IllegalArgumentException("Invalid pattern length in " + fileName, ex);
        }
    }

    private static int compareNaturally(String left, String right) {
        try {
            int leftInt = Integer.parseInt(left);
            int rightInt = Integer.parseInt(right);
            return Integer.compare(leftInt, rightInt);
        } catch (NumberFormatException ignored) {
            return left.compareTo(right);
        }
    }

    private static final class AggregateStats {
        private final EnumMap<IndexType, StatsSum> sums = new EnumMap<>(IndexType.class);

        void accumulate(IndexType indexType,
                        double avgInsertMsPerSymbol,
                        double totalInsertMs,
                        double totalQueryMs) {
            StatsSum sum = sums.computeIfAbsent(indexType, ignored -> new StatsSum());
            sum.avgInsertMsPerSymbolSum += avgInsertMsPerSymbol;
            sum.totalInsertMsSum += totalInsertMs;
            sum.totalQueryMsSum += totalQueryMs;
            sum.count++;
        }

        StatsSnapshot snapshot(IndexType indexType) {
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

        boolean hasData() {
            return sums.values().stream().anyMatch(sum -> sum.count > 0);
        }

        private static final class StatsSum {
            private double avgInsertMsPerSymbolSum;
            private double totalInsertMsSum;
            private double totalQueryMsSum;
            private int count;
        }
    }

    private record StatsSnapshot(double avgInsertMsPerSymbol,
                                 double avgInsertMs,
                                 double avgQueryMs,
                                 int count) {
    }

    private enum IndexType {
        HBI("HBI", "hbi"),
        SUFFIX("SuffixTree", "suffix"),
        REGEX("Regex", "regex");

        private final String displayName;
        private final String csvLabel;

        IndexType(String displayName, String csvLabel) {
            this.displayName = displayName;
            this.csvLabel = csvLabel;
        }

        String displayName() {
            return displayName;
        }

        String csvLabel() {
            return csvLabel;
        }
    }

    private enum QueryType {
        MISSING("missing"),
        RARE("rare"),
        UNIFORM("uniform");

        private final String fileToken;

        QueryType(String fileToken) {
            this.fileToken = fileToken;
        }

        String fileToken() {
            return fileToken;
        }

        static QueryType fromString(String value) {
            return EnumSet.allOf(QueryType.class).stream()
                    .filter(type -> type.fileToken.equalsIgnoreCase(value))
                    .findFirst()
                    .orElseThrow(() -> new IllegalArgumentException(
                            "Unknown query type: " + value));
        }
    }
}
