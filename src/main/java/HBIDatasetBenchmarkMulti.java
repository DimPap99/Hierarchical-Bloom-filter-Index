import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import estimators.CostFunctionMaxProb;
import estimators.CSEstimator;
import estimators.Estimator;
import membership.BloomFilter;
import membership.Membership;
import search.*;

import utilities.ExperimentRunResult;
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

public final class HBIDatasetBenchmarkMulti {

    private static final boolean USE_STRIDES = true;



    public static void main(String[] args) throws IOException {
        BenchmarkOptions options = BenchmarkOptions.parse(args);
        List<QueryType> queryTypes = options.orderedQueryTypes();
        Map<QueryType, Map<Integer, AggregateStats>> aggregated = new EnumMap<>(QueryType.class);
        queryTypes.forEach(type -> aggregated.put(type, new TreeMap<>()));

        System.out.printf(Locale.ROOT,
                "Running window %s (%s) with datasets under %s%n",
                options.window(), options.queryTypeLabel(), options.dataRoot());

        List<Path> datasetDirs = listDatasetDirectories(options.dataRoot());
        if (datasetDirs.isEmpty()) {
            throw new IllegalStateException("No dataset folders found under " + options.dataRoot());
        }

        for (Path datasetDir : datasetDirs) {
            Path datasetFile = findSingleDatasetFile(datasetDir);
            Path queryDir = options.queryRoot().resolve(datasetDir.getFileName());
            if (!Files.isDirectory(queryDir)) {
                System.out.printf(Locale.ROOT,
                        "Skipping dataset %s (missing query folder %s)%n",
                        datasetDir.getFileName(), queryDir);
                continue;
            }

            System.out.printf(Locale.ROOT, "Dataset %s -> %s%n",
                    datasetDir.getFileName(), datasetFile.getFileName());

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

            for (int warmIndex = 0; warmIndex < options.warmupRuns(); warmIndex++) {
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

            for (int runIndex = 0; runIndex < options.runs(); runIndex++) {
                HBI hbi = newHbi(options);
                hbi.strides = USE_STRIDES;
                hbi.stats().setCollecting(false);
                hbi.stats().setExperimentMode(false);

                MultiQueryExperiment.InsertStats hbiInsertStats =
                        MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, options.ngram());

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

                    MultiQueryExperiment.MultiRunResult hbiResults =
                            MultiQueryExperiment.runQueries(workloads, hbi, options.ngram(), hbiInsertStats, false, false);

                    MultiQueryExperiment.MultiRunResult regexResults;
                    if (regex != null) {
                        regexResults = MultiQueryExperiment.runQueries(workloads, regex, 1, regexInsertStats, false, false);
                    } else {
                        regexResults = null;
                    }

                    Map<Integer, AggregateStats> perTypeAggregated = aggregated.get(type);

                    hbiResults.results().forEach((workload, result) -> {
                        int patternLength = workload.patternLength();
                        ExperimentRunResult regexResult =
                                (regexResults != null) ? regexResults.results().get(workload) : null;

                        RunAverages averages = new RunAverages(
                                result.avgInsertMsPerSymbol(),
                                result.totalInsertTimeMs(),
                                result.totalRunTimeMs(),
                                regexResult != null ? regexResult.totalInsertTimeMs() : 0.0,
                                regexResult != null ? regexResult.totalRunTimeMs() : 0.0);

                        perTypeAggregated.computeIfAbsent(patternLength, ignored -> new AggregateStats())
                                .accumulate(averages);
                    });
                }
            }
        }

        boolean anyResults = aggregated.values().stream().anyMatch(map -> !map.isEmpty());
        if (!anyResults) {
            System.out.println("No runs executed.");
            return;
        }

        System.out.println();

        boolean includeRegex = aggregated.values().stream()
                .flatMap(map -> map.values().stream())
                .anyMatch(stats -> stats.regexCount > 0);

        aggregated.forEach((type, statsByPattern) -> {
            int datasetCount = statsByPattern.values().stream()
                    .mapToInt(AggregateStats::count)
                    .max()
                    .orElse(0);

            System.out.printf(Locale.ROOT,
                    "Aggregated averages for %s across %d dataset(s):%n",
                    type.fileToken(), datasetCount);

            if (statsByPattern.isEmpty()) {
                System.out.println("  No data available.");
                System.out.println();
                return;
            }

            statsByPattern.forEach((patternLength, stats) -> {
                double avgInsertPerSymbol = stats.avgInsertMsPerSymbolSum / stats.count();
                double avgInsertMs = stats.totalInsertMsSum / stats.count();
                double avgQueryMs = stats.totalQueryMsSum / stats.count();
                System.out.printf(Locale.ROOT,
                        "  Pattern %d -> avgInsert=%.6f ms, insert=%.3f ms, query=%.3f ms over %d dataset(s)%n",
                        patternLength,
                        avgInsertPerSymbol,
                        avgInsertMs,
                        avgQueryMs,
                        stats.count());

                if (includeRegex && stats.regexCount > 0) {
                    double avgRegexInsert = stats.regexInsertMsSum / stats.regexCount;
                    double avgRegexQuery = stats.regexQueryMsSum / stats.regexCount;
                    System.out.printf(Locale.ROOT,
                            "         Regex -> insert=%.3f ms, query=%.3f ms over %d dataset(s)%n",
                            avgRegexInsert,
                            avgRegexQuery,
                            stats.regexCount);
                }
            });

            System.out.println();
        });

        List<Object> header = new ArrayList<>();
        header.add("patternLength");

        queryTypes.forEach(type -> {
            header.add("avgInsertMsPerSymbol_" + type.fileToken());
            header.add("insertMs_" + type.fileToken());
            header.add("avgMs_" + type.fileToken());
            header.add("datasetCount_" + type.fileToken());
            if (includeRegex) {
                header.add("regexAvgInsertMs_" + type.fileToken());
                header.add("regexAvgQueryMs_" + type.fileToken());
                header.add("regexDatasetCount_" + type.fileToken());
            }
        });

        List<List<?>> csvRows = new ArrayList<>();
        csvRows.add(header);

        SortedSet<Integer> allPatternLengths = new TreeSet<>();
        aggregated.values().forEach(map -> allPatternLengths.addAll(map.keySet()));

        for (Integer patternLength : allPatternLengths) {
            List<Object> row = new ArrayList<>();
            row.add(patternLength);

            for (QueryType type : queryTypes) {
                AggregateStats stats = aggregated.get(type).get(patternLength);
                if (stats == null || stats.count() == 0) {
                    row.add(null);
                    row.add(null);
                    row.add(null);
                    row.add(0);
                    if (includeRegex) {
                        row.add(null);
                        row.add(null);
                        row.add(0);
                    }
                    continue;
                }

                double avgInsertMsPerSymbol = stats.avgInsertMsPerSymbolSum / stats.count();
                double insertMs = stats.totalInsertMsSum / stats.count();
                double avgQueryMs = stats.totalQueryMsSum / stats.count();

                row.add(avgInsertMsPerSymbol);
                row.add(insertMs);
                row.add(avgQueryMs);
                row.add(stats.count());

                if (includeRegex) {
                    Double avgRegexInsert = stats.regexCount > 0
                            ? stats.regexInsertMsSum / stats.regexCount
                            : null;
                    Double avgRegexQuery = stats.regexCount > 0
                            ? stats.regexQueryMsSum / stats.regexCount
                            : null;
                    row.add(avgRegexInsert);
                    row.add(avgRegexQuery);
                    row.add(stats.regexCount);
                }
            }

            csvRows.add(row);
        }

        Path csvPath = options.csvOutputFile();
        CsvUtil.writeRows(csvPath, csvRows);
        System.out.printf(Locale.ROOT, "Wrote aggregated results to %s%n", csvPath);
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
                .costFunction(new CostFunctionMaxProb())
                .confidence(options.runConfidence())
                .experimentMode(false)
                .collectStats(false)
                .nGram(options.ngram())
                .build();

        return new HBI(configuration);
    }

    private static List<Path> listDatasetDirectories(Path root) throws IOException {
        try (Stream<Path> stream = Files.list(root)) {
            return stream
                    .filter(Files::isDirectory)
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), HBIDatasetBenchmarkMulti::compareNaturally))
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
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), HBIDatasetBenchmarkMulti::compareNaturally))
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
        private double avgInsertMsPerSymbolSum;
        private double totalInsertMsSum;
        private double totalQueryMsSum;
        private int count;
        private double regexInsertMsSum;
        private double regexQueryMsSum;
        private int regexCount;

        void accumulate(RunAverages averages) {
            avgInsertMsPerSymbolSum += averages.avgInsertMsPerSymbol();
            totalInsertMsSum += averages.totalInsertMs();
            totalQueryMsSum += averages.totalQueryMs();
            count++;

            if (averages.regexInsertMs() > 0 || averages.regexQueryMs() > 0) {
                regexInsertMsSum += averages.regexInsertMs();
                regexQueryMsSum += averages.regexQueryMs();
                regexCount++;
            }
        }

        int count() {
            return count;
        }
    }

    private record RunAverages(double avgInsertMsPerSymbol,
                               double totalInsertMs,
                               double totalQueryMs,
                               double regexInsertMs,
                               double regexQueryMs) {
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
                                    boolean runRegexBaseline) {

        static BenchmarkOptions parse(String[] args) {
            Path dataRoot = Path.of("data");
            Path queryRoot = Path.of("queries");
            String window = "w23";
            QueryType queryType = null;
            Integer ngram = 2;
            Integer windowLength = null;
            Integer treeLength = null;
            Integer alphabetBase = 130;
            double fpRate = 0.001;
            double runConfidence = 0.99;
            int warmupRuns = 1;
            int runs = 3;
            boolean runRegexBaseline = true;

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
                    case "regex" -> runRegexBaseline = Boolean.parseBoolean(value);
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
                windowLength = deriveWindowLength(window, 1 << 21);
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
                    runRegexBaseline);
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

        Path csvOutputFile() {
            String fileName = "%s_%s_%d.csv".formatted(window, queryTypeLabel(), ngram);
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
}
