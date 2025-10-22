import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.RegexIndex;
import PMIndex.SuffixTreeIndex;

import estimators.CSEstimator;
import estimators.CostFunctionMaxProb;
import estimators.Estimator;

import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;

import org.openjdk.jol.info.GraphLayout;

import search.BlockSearch;
import search.MostFreqPruning;
import search.PruningPlan;
import search.Verifier;
import search.VerifierLinearLeafProbe;

import utilities.CsvUtil;
import utilities.MemUtil;
import utilities.MemoryUsageReport;
import utilities.MultiQueryExperiment;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * MemoryBenchmark
 *
 * Build each index exactly once on a single dataset file and record total retained memory in MiB.
 * Columns: ngram, hbi_mem_mib, suffix_mem_mib, regexp_mem_mib, hbifp_rate.
 *
 * Dataset selection follows HBIDatasetBenchmarkMulti helpers:
 * it lists dataset subdirectories under --data-root/<window> and picks the first by natural order.
 */
public final class MemoryBenchmark {

    private static final boolean USE_STRIDES = true;

    public static void main(String[] args) throws IOException {
        MemoryOptions options = MemoryOptions.parse(args);
        MemUtil memUtil = new MemUtil();

        System.out.printf(Locale.ROOT,
                "Running memory benchmark on window %s under %s%n",
                options.window(), options.dataRoot());

        List<Path> datasetDirs = listDatasetDirectories(options.dataRoot());
        if (datasetDirs.isEmpty()) {
            throw new IllegalStateException("No dataset folders found under " + options.dataRoot());
        }

        // Pick the first dataset directory by natural order
        Path pickedDir = datasetDirs.get(0);
        Path datasetFile = findSingleDatasetFile(pickedDir);
        System.out.printf(Locale.ROOT, "Using dataset %s -> %s%n",
                pickedDir.getFileName(), datasetFile.getFileName());

        // CSV header in MiB
        List<List<?>> csvRows = new ArrayList<>();
        csvRows.add(List.of("ngram", "hbi_mem_mib", "suffix_mem_mib", "regexp_mem_mib", "hbifp_rate"));

        // ===== HBI =====
        HBI hbi = newHbi(options);
        hbi.strides = USE_STRIDES;
        MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, options.ngram());

        MemoryUsageReport hbiReport = memUtil.jolMemoryReportPartitionedWithTotal(hbi);
        System.out.println("\n=== HBI JOL partitioned report ===");
        System.out.println(hbiReport.report());
        double hbiMiB = hbiReport.totalMiB();
        System.out.printf(Locale.ROOT, "HBI total MiB: %.3f%n", hbiMiB);

        // ===== Suffix tree =====
        SuffixTreeIndex suffix = newSuffixTree(options);
        MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, options.ngram());

        MemoryUsageReport suffixReport = suffix.jolMemoryReportPartitionedWithTotal();
        System.out.println("\n=== SuffixTree JOL partitioned report ===");
        System.out.println(suffixReport.report());
        double suffixMiB = suffixReport.totalMiB();
        System.out.printf(Locale.ROOT, "SuffixTree total MiB: %.3f%n", suffixMiB);

        // ===== Regex baseline =====
        RegexIndex regex = new RegexIndex();
        // Regex baseline indexes the raw character stream, so use ngram = 1
        MultiQueryExperiment.populateIndex(datasetFile.toString(), regex, 1);

        long regexBytes = GraphLayout.parseInstance(regex).totalSize();
        double regexMiB = regexBytes / 1_048_576d;
        System.out.println("\n=== Regex baseline JOL footprint ===");
        System.out.println(GraphLayout.parseInstance(regex).toFootprint());
        System.out.printf(Locale.ROOT, "Regex total MiB: %.3f%n", regexMiB);

        // ===== CSV row =====
        csvRows.add(List.of(options.ngram(), hbiMiB, suffixMiB, regexMiB, options.fpRate()));
        Path out = options.csvOutputFile();
        CsvUtil.writeRows(out, csvRows);
        System.out.printf(Locale.ROOT, "Wrote memory results to %s%n", out);
    }

    // ===== Helpers mirroring HBIDatasetBenchmarkMulti =====

    private static HBI newHbi(MemoryOptions options) {
        int alphabetSize = options.alphabetSize();
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(options.treeLength());//new CSEstimator(options.treeLength(), 5, 16384);
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

    private static SuffixTreeIndex newSuffixTree(MemoryOptions options) {
        long expectedDistinct = Math.max(1L, (long) options.windowLength());
        double epsilon = 0.0001;
        int totalTokens = options.treeLength();
        return new SuffixTreeIndex(expectedDistinct, epsilon, totalTokens);
    }

    private static List<Path> listDatasetDirectories(Path root) throws IOException {
        try (Stream<Path> stream = Files.list(root)) {
            return stream
                    .filter(Files::isDirectory)
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), MemoryBenchmark::compareNaturally))
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

    private static int compareNaturally(String left, String right) {
        try {
            int leftInt = Integer.parseInt(left);
            int rightInt = Integer.parseInt(right);
            return Integer.compare(leftInt, rightInt);
        } catch (NumberFormatException ignored) {
            return left.compareTo(right);
        }
    }

    // ===== Options focused on memory CSV =====

    private record MemoryOptions(Path dataRoot,
                                 String window,
                                 int ngram,
                                 int windowLength,
                                 int treeLength,
                                 int alphabetBase,
                                 double fpRate,
                                 double runConfidence) {

        static MemoryOptions parse(String[] args) {
            Path dataRoot = Path.of("data");
            String window = "w21";
            Integer ngram = 5;
            Integer windowLength = null;
            Integer treeLength = null;
            Integer alphabetBase = 130;
            double fpRate = 0.001;
            double runConfidence = 0.99;

            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (!arg.startsWith("--")) continue;
                String key;
                String value;
                int eq = arg.indexOf('=');
                if (eq >= 0) {
                    key = arg.substring(2, eq);
                    value = arg.substring(eq + 1);
                } else {
                    key = arg.substring(2);
                    if (i + 1 >= args.length) throw new IllegalArgumentException("Missing value for option --" + key);
                    value = args[++i];
                }
                switch (key) {
                    case "window" -> window = value;
                    case "data-root" -> dataRoot = Path.of(value);
                    case "ngrams" -> ngram = Integer.parseInt(value);
                    case "window-length" -> windowLength = Integer.parseInt(value);
                    case "tree-length" -> treeLength = Integer.parseInt(value);
                    case "alphabet-base" -> alphabetBase = Integer.parseInt(value);
                    case "fp" -> fpRate = Double.parseDouble(value);
                    case "confidence" -> runConfidence = Double.parseDouble(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            if (ngram == null) ngram = 5;
            if (windowLength == null) windowLength = deriveWindowLength(window, 1 << 21);
            if (treeLength == null) treeLength = windowLength;

            Path resolvedDataRoot = dataRoot.resolve(window);
            if (!Files.isDirectory(resolvedDataRoot)) {
                throw new IllegalArgumentException("Data directory not found: " + resolvedDataRoot);
            }

            return new MemoryOptions(
                    resolvedDataRoot,
                    window,
                    ngram,
                    windowLength,
                    treeLength,
                    alphabetBase,
                    fpRate,
                    runConfidence);
        }

        int alphabetSize() {
            double pow = Math.pow(alphabetBase, ngram);
            double sigma = Math.min(pow, windowLength);
            if (sigma >= Integer.MAX_VALUE) return Integer.MAX_VALUE;
            return (int) sigma;
        }

        Path csvOutputFile() {
            String fileName = "memory_%s_%d.csv".formatted(window, ngram);
            return Path.of(fileName);
        }

        static int deriveWindowLength(String window, int defaultValue) {
            if (window == null || window.isEmpty()) return defaultValue;
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
