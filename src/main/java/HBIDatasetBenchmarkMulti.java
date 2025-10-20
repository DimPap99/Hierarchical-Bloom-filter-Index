import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import estimators.CostFunctionMaxProb;
import estimators.CSEstimator;
import estimators.Estimator;
import membership.BloomFilter;
import membership.Membership;
import search.BlockSearch;
import search.MostFreqPruning;
import search.PruningPlan;
import search.Verifier;
import search.VerifierLinearLeafProbe;

import utilities.ExperimentRunResult;
import utilities.MultiQueryExperiment;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TreeMap;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import utilities.CsvUtil;

public final class HBIDatasetBenchmarkMulti {

    private static final boolean USE_STRIDES = true;

    public static void main(String[] args) throws IOException {
        BenchmarkOptions options = BenchmarkOptions.parse(args);
        Map<Integer, AggregateStats> aggregated = new TreeMap<>();

        System.out.printf(Locale.ROOT,
                "Running window %s (%s) with datasets under %s%n",
                options.window(), options.queryType().fileToken(), options.dataRoot());

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

            List<Path> queryFiles = findQueryFiles(queryDir, options.queryType());
            if (queryFiles.isEmpty()) {
                System.out.printf(Locale.ROOT,
                        "No %s queries in %s, skipping.%n",
                        options.queryType().fileToken(), queryDir);
                continue;
            }

            List<MultiQueryExperiment.QueryWorkload> workloads = queryFiles.stream()
                    .map(path -> new MultiQueryExperiment.QueryWorkload(
                            parsePatternLength(path.getFileName().toString()),
                            path))
                    .sorted(Comparator.comparingInt(MultiQueryExperiment.QueryWorkload::patternLength))
                    .collect(Collectors.toList());

            Map<Integer, RunAverages> resultsByPattern = runExperiments(options, datasetFile, workloads);

            resultsByPattern.forEach((patternLength, averages) -> {
                aggregated.computeIfAbsent(patternLength, ignored -> new AggregateStats())
                        .accumulate(averages);

                System.out.printf(Locale.ROOT,
                        "  Pattern %d -> avgInsert=%.6f ms, insert=%.3f ms, query=%.3f ms%n",
                        patternLength,
                        averages.avgInsertMsPerSymbol(),
                        averages.totalInsertMs(),
                        averages.totalQueryMs());

                if (options.runRegexBaseline() && (averages.regexInsertMs() > 0 || averages.regexQueryMs() > 0)) {
                    System.out.printf(Locale.ROOT,
                            "    Regex -> insert=%.3f ms, query=%.3f ms%n",
                            averages.regexInsertMs(),
                            averages.regexQueryMs());
                }
            });
        }

        if (aggregated.isEmpty()) {
            System.out.println("No runs executed.");
            return;
        }

        System.out.println();
        int datasetCount = aggregated.values().stream()
                .mapToInt(AggregateStats::count)
                .max()
                .orElse(0);

        System.out.printf(Locale.ROOT,
                "Aggregated averages across %d dataset(s):%n",
                datasetCount);

        boolean includeRegex = aggregated.values().stream()
                .anyMatch(stats -> stats.regexCount > 0);

        List<List<?>> csvRows = new ArrayList<>();
        if (includeRegex) {
            csvRows.add(List.of("patternLength",
                    "avgInsertMsPerSymbol",
                    "avgInsertMs",
                    "avgQueryMs",
                    "datasetCount",
                    "regexAvgInsertMs",
                    "regexAvgQueryMs",
                    "regexDatasetCount"));
        } else {
            csvRows.add(List.of("patternLength",
                    "avgInsertMsPerSymbol",
                    "avgInsertMs",
                    "avgQueryMs",
                    "datasetCount"));
        }

        aggregated.forEach((patternLength, stats) -> {
            double avgInsertPerSymbol = stats.avgInsertMsPerSymbolSum / stats.count;
            double avgInsertMs = stats.totalInsertMsSum / stats.count;
            double avgQueryMs = stats.totalQueryMsSum / stats.count;
            double avgRegexInsert = stats.regexCount > 0
                    ? stats.regexInsertMsSum / stats.regexCount
                    : 0.0;
            double avgRegexQuery = stats.regexCount > 0
                    ? stats.regexQueryMsSum / stats.regexCount
                    : 0.0;
            System.out.printf(Locale.ROOT,
                    "Pattern %d -> avgInsert=%.6f ms, insert=%.3f ms, query=%.3f ms over %d dataset(s)%n",
                    patternLength,
                    avgInsertPerSymbol,
                    avgInsertMs,
                    avgQueryMs,
                    stats.count());

            if (includeRegex && stats.regexCount > 0) {
                System.out.printf(Locale.ROOT,
                        "           Regex -> insert=%.3f ms, query=%.3f ms over %d dataset(s)%n",
                        avgRegexInsert,
                        avgRegexQuery,
                        stats.regexCount);
            }

            if (includeRegex) {
                csvRows.add(List.of(
                        patternLength,
                        avgInsertPerSymbol,
                        avgInsertMs,
                        avgQueryMs,
                        stats.count(),
                        avgRegexInsert,
                        avgRegexQuery,
                        stats.regexCount));
            } else {
                csvRows.add(List.of(
                        patternLength,
                        avgInsertPerSymbol,
                        avgInsertMs,
                        avgQueryMs,
                        stats.count()));
            }
        });

        Path csvPath = options.csvOutputFile();
        CsvUtil.writeRows(csvPath, csvRows);
        System.out.printf(Locale.ROOT, "Wrote aggregated results to %s%n", csvPath);
    }

    private static Map<Integer, RunAverages> runExperiments(BenchmarkOptions options,
                                                            Path datasetFile,
                                                            List<MultiQueryExperiment.QueryWorkload> workloads) throws IOException {
        if (workloads.isEmpty()) {
            return Map.of();
        }

        for (int i = 0; i < options.warmupRuns(); i++) {
            HBI warmupIndex = newHbi(options);
            warmupIndex.strides = USE_STRIDES;
            warmupIndex.stats().setCollecting(false);
            warmupIndex.stats().setExperimentMode(false);
            MultiQueryExperiment.run(datasetFile.toString(), workloads, warmupIndex, options.ngram(), false, false);
            if (options.runRegexBaseline()) {
                IPMIndexing regexWarmup = new RegexIndex();
                MultiQueryExperiment.run(datasetFile.toString(), workloads, regexWarmup, 1, false, false);
            }
        }

        Map<Integer, SampleAccumulator> hbiSamples = new TreeMap<>();
        Map<Integer, SampleAccumulator> regexSamples = new TreeMap<>();

        for (int i = 0; i < options.runs(); i++) {
            HBI hbi = newHbi(options);
            hbi.strides = USE_STRIDES;
            hbi.stats().setCollecting(false);
            hbi.stats().setExperimentMode(false);
            MultiQueryExperiment.MultiRunResult hbiResult = MultiQueryExperiment.run(
                    datasetFile.toString(),
                    workloads,
                    hbi,
                    options.ngram(),
                    false,
                    false);

            hbiResult.results().forEach((workload, result) ->
                    hbiSamples.computeIfAbsent(workload.patternLength(), ignored -> new SampleAccumulator())
                            .add(result));

            if (options.runRegexBaseline()) {
                IPMIndexing regex = new RegexIndex();
                MultiQueryExperiment.MultiRunResult regexResult = MultiQueryExperiment.run(
                        datasetFile.toString(),
                        workloads,
                        regex,
                        1,
                        false,
                        false);

                regexResult.results().forEach((workload, result) ->
                        regexSamples.computeIfAbsent(workload.patternLength(), ignored -> new SampleAccumulator())
                                .add(result));
            }
        }

        Map<Integer, RunAverages> averages = new TreeMap<>();
        hbiSamples.forEach((patternLength, sample) -> {
            if (sample.count == 0) {
                return;
            }
            double avgInsertMsPerSymbol = sample.avgInsertMsPerSymbolSum / sample.count;
            double avgInsertMs = sample.totalInsertMsSum / sample.count;
            double avgQueryMs = sample.totalQueryMsSum / sample.count;

            SampleAccumulator regexSample = regexSamples.get(patternLength);
            double regexInsert = regexSample == null || regexSample.count == 0
                    ? 0.0
                    : regexSample.totalInsertMsSum / regexSample.count;
            double regexQuery = regexSample == null || regexSample.count == 0
                    ? 0.0
                    : regexSample.totalQueryMsSum / regexSample.count;

            averages.put(patternLength, new RunAverages(
                    avgInsertMsPerSymbol,
                    avgInsertMs,
                    avgQueryMs,
                    regexInsert,
                    regexQuery));
        });

        return averages;
    }

    private static HBI newHbi(BenchmarkOptions options) {
        int alphabetSize = options.alphabetSize(options.windowLength);
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

    private static final class SampleAccumulator {
        private double avgInsertMsPerSymbolSum;
        private double totalInsertMsSum;
        private double totalQueryMsSum;
        private int count;

        void add(ExperimentRunResult result) {
            avgInsertMsPerSymbolSum += result.avgInsertMsPerSymbol();
            totalInsertMsSum += result.totalInsertTimeMs();
            totalQueryMsSum += result.totalRunTimeMs();
            count++;
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
            String window = "w21";
            QueryType queryType = QueryType.UNIFORM;
            Integer ngram = 2;
            Integer windowLength = null;
            Integer treeLength = null;
            Integer alphabetBase = 95;
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
                    case "type" -> queryType = QueryType.fromString(value);
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
            if (queryType == null) {
                throw new IllegalArgumentException("--type is required");
            }

            if (ngram == null) {
                ngram = 5;
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

        int alphabetSize(int windowLength) {
            double pow = Math.pow(alphabetBase, ngram);
            if (pow >= Integer.MAX_VALUE) {
                return Integer.MAX_VALUE;
            }
            double sigma = Math.min(pow, windowLength);
            return (int) sigma;
        }

        Path csvOutputFile() {
            String fileName = "%s_%s_%d.csv".formatted(window, queryType.fileToken(), ngram);
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
