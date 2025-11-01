import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.RegexIndex;
import PMIndex.StreamingSlidingWindowIndex;

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
import utilities.MultiQueryExperiment;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.function.Supplier;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * MemoryBenchmark
 *
 * Build each index once on one dataset file and record retained memory in mebibytes (MiB),
 * measured using Java Object Layout, through GraphLayout.parseInstance(obj).totalSize().
 *
 * Output CSV columns:
 *
 *   window,
 *   window_power,
 *   window_length,
 *   tree_power,
 *   tree_length,
 *   ngram,
 *   fp_rate,
 *   hbi_mem_mib,
 *   suffix_core_mem_mib,
 *   suffix_total_mem_mib,
 *   regex_mem_mib
 *
 * hbi_mem_mib is the retained size of the Hierarchical Bloom filter Index (HBI).
 *
 * suffix_mem_mib is the retained size of the StreamingSlidingWindowIndex after subtracting
 * its AlphabetMapper dictionary map, so you measure just the rolling suffix segment hierarchy
 * and not the global token to identifier HashMap.
 *
 * regexp_mem_mib is the retained size of a trivial RegexIndex baseline.
 *
 * Important
 * We no longer pass lambda expressions or method references into HbiConfiguration.
 * Instead we pass small named classes that implement Supplier. This avoids hidden
 * lambda classes in the object graph. Java Object Layout in modern Java crashes
 * when it tries to inspect hidden lambda classes. With these named classes, it
 * will see only ordinary classes and it can compute the total size.
 */
public final class MemoryBenchmark {

    private static final boolean USE_STRIDES = true;

    public static void main(String[] args) throws IOException {
        MemoryOptions options = MemoryOptions.parse(args);

        System.out.printf(Locale.ROOT,
                "Running memory benchmark on window %s under %s%n",
                options.window(), options.dataRoot());

        List<Path> datasetDirs = listDatasetDirectories(options.dataRoot());
        if (datasetDirs.isEmpty()) {
            throw new IllegalStateException("No dataset folders found under " + options.dataRoot());
        }

        // pick first dataset directory by natural order
        Path pickedDir = datasetDirs.get(0);
        Path datasetFile = findSingleDatasetFile(pickedDir);
        System.out.printf(Locale.ROOT, "Using dataset %s -> %s%n",
                pickedDir.getFileName(), datasetFile.getFileName());

        // CSV header
        List<List<?>> csvRows = new ArrayList<>();
        csvRows.add(List.of(
                "window",
                "window_power",
                "window_length",
                "tree_power",
                "tree_length",
                "ngram",
                "fp_rate",
                "hbi_mem_mib",
                "suffix_core_mem_mib",
                "suffix_total_mem_mib",
                "regex_mem_mib"));

        Map<Integer, SuffixMeasurement> suffixCache = new HashMap<>();
        RegexMeasurement regexCache = null;

        for (double fpRate : options.fpRates()) {
            System.out.printf(Locale.ROOT, "%n==== FPR = %.6f ====%n", fpRate);
            for (int ngram : options.ngrams()) {
                System.out.printf(Locale.ROOT, "---- ngram = %d ----%n", ngram);

                HbiMeasurement hbiStats = buildAndMeasureHbi(datasetFile, options, ngram, fpRate);

                SuffixMeasurement suffixStats;
                if (options.reuseSuffixResults() && suffixCache.containsKey(ngram)) {
                    suffixStats = suffixCache.get(ngram);
                    System.out.printf(Locale.ROOT, "Reusing suffix measurement for ngram %d%n", ngram);
                    System.out.printf(Locale.ROOT,
                            "SuffixIndex total (reuse): %d B (%.3f MiB) [ngram=%d]%n",
                            suffixStats.totalBytes(),
                            suffixStats.totalMiB(),
                            ngram);
                    System.out.printf(Locale.ROOT,
                            "SuffixIndex core (reuse): %d B (%.3f MiB) [ngram=%d]%n",
                            suffixStats.coreBytes(),
                            suffixStats.coreMiB(),
                            ngram);
                } else {
                    suffixStats = buildAndMeasureSuffix(datasetFile, options, ngram);
                    if (options.reuseSuffixResults()) {
                        suffixCache.put(ngram, suffixStats);
                    }
                }

                if (regexCache == null) {
                    regexCache = buildAndMeasureRegex(datasetFile);
                }

                csvRows.add(List.of(
                        options.window(),
                        options.windowPower(),
                        options.windowLength(),
                        options.treePower(),
                        options.treeLength(),
                        ngram,
                        fpRate,
                        hbiStats.mib(),
                        suffixStats.coreMiB(),
                        suffixStats.totalMiB(),
                        regexCache.mib()));
            }
        }

        Path out = options.csvOutputFile();
        CsvUtil.writeRows(out, csvRows);
        System.out.printf(Locale.ROOT, "Wrote memory results to %s%n", out);
    }

    /**
     * Build an HBI for the requested (ngram, fpRate) combination using named Supplier classes,
     * so Java Object Layout never encounters hidden lambda implementations.
     */
    private static HBI newHbi(MemoryOptions options, int nGram, double fpRate) {
        final int windowLen  = options.windowLength();
        final int sigma      = options.alphabetSizeFor(nGram);
        final int treeLen    = options.treeLength();
        final double runConf = options.runConfidence();

        Supplier<Estimator> estFactory =
                new EstimatorFactory(treeLen);
        Supplier<Membership> memFactory =
                new MembershipFactory();
        Supplier<PruningPlan> prFactory =
                new PruningPlanFactory(runConf, fpRate);

        Verifier verifier = new VerifierLinearLeafProbe();

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(windowLen)
                .fpRate(fpRate)
                .alphabetSize(sigma)
                .treeLength(treeLen)
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(verifier)
                .costFunction(new CostFunctionMaxProb())
                .confidence(runConf)
                .experimentMode(false)
                .collectStats(false)
                .nGram(nGram)
                .build();

        return new HBI(configuration);
    }

    /**
     * Create the streaming suffix baseline using the same window size and expected alphabet size.
     */
    private static StreamingSlidingWindowIndex newSuffixIndex(MemoryOptions options, int nGram) {
        int windowSize = options.windowLength();
        int expectedAlphabetSize = options.alphabetSizeFor(nGram);
        return new StreamingSlidingWindowIndex(windowSize, expectedAlphabetSize);
    }

    private static HbiMeasurement buildAndMeasureHbi(Path datasetFile,
                                                     MemoryOptions options,
                                                     int ngram,
                                                     double fpRate) throws IOException {
        HBI hbi = newHbi(options, ngram, fpRate);
        hbi.strides = USE_STRIDES;
        MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, ngram);

        GraphLayout layout = GraphLayout.parseInstance(hbi);
        long bytes = layout.totalSize();
        double mib = bytes / 1_048_576d;

        System.out.println("\n=== HBI JOL footprint ===");
        System.out.println(layout.toFootprint());
        System.out.printf(Locale.ROOT,
                "HBI total: %d B (%.3f MiB) [ngram=%d, fp=%.6f]%n",
                bytes,
                mib,
                ngram,
                fpRate);

        return new HbiMeasurement(bytes, mib);
    }

    private static SuffixMeasurement buildAndMeasureSuffix(Path datasetFile,
                                                           MemoryOptions options,
                                                           int ngram) throws IOException {
        StreamingSlidingWindowIndex suffix = newSuffixIndex(options, ngram);
        MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, ngram);

        GraphLayout layout = GraphLayout.parseInstance(suffix);
        long totalBytes = layout.totalSize();
        long coreBytes = suffix.estimateRetainedBytesWithoutAlphabet();
        double coreMiB = coreBytes / 1_048_576d;
        double totalMiB = totalBytes / 1_048_576d;

        System.out.println("\n=== StreamingSlidingWindowIndex JOL footprint ===");
        System.out.println(layout.toFootprint());
        System.out.printf(Locale.ROOT,
                "SuffixIndex total (full, incl. AlphabetMapper): %d B (%.3f MiB) [ngram=%d]%n",
                totalBytes,
                totalMiB,
                ngram);
        System.out.printf(Locale.ROOT,
                "SuffixIndex core (excluding AlphabetMapper): %d B (%.3f MiB) [ngram=%d]%n",
                coreBytes,
                coreMiB,
                ngram);

        return new SuffixMeasurement(totalBytes, coreBytes, totalMiB, coreMiB);
    }

    private static RegexMeasurement buildAndMeasureRegex(Path datasetFile) throws IOException {
        RegexIndex regex = new RegexIndex();
        MultiQueryExperiment.populateIndex(datasetFile.toString(), regex, 1);

        GraphLayout layout = GraphLayout.parseInstance(regex);
        long bytes = layout.totalSize();
        double mib = bytes / 1_048_576d;

        System.out.println("\n=== Regex baseline JOL footprint ===");
        System.out.println(layout.toFootprint());
        System.out.printf(Locale.ROOT,
                "Regex total: %d B (%.3f MiB)%n",
                bytes,
                mib);

        return new RegexMeasurement(bytes, mib);
    }

    private record HbiMeasurement(long bytes, double mib) {}

    private record SuffixMeasurement(long totalBytes, long coreBytes, double totalMiB, double coreMiB) {}

    private record RegexMeasurement(long bytes, double mib) {}

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

    /**
     * Supplier<Estimator> without lambdas, so JOL will not encounter a hidden lambda class.
     */
    private static final class EstimatorFactory implements Supplier<Estimator> {
        private final int treeLen;

        private EstimatorFactory(int treeLen) {
            this.treeLen = treeLen;
        }

        @Override
        public Estimator get() {
            return new HashMapEstimator(treeLen);
        }
    }

    /**
     * Supplier<Membership> without lambdas.
     */
    private static final class MembershipFactory implements Supplier<Membership> {

        private MembershipFactory() {
        }

        @Override
        public Membership get() {
            return new BloomFilter();
        }
    }

    /**
     * Supplier<PruningPlan> without lambdas.
     * Stores the run time confidence and Bloom filter false positive rate as final fields.
     */
    private static final class PruningPlanFactory implements Supplier<PruningPlan> {
        private final double runConf;
        private final double fpRate;

        private PruningPlanFactory(double runConf, double fpRate) {
            this.runConf = runConf;
            this.fpRate = fpRate;
        }

        @Override
        public PruningPlan get() {
            return new MostFreqPruning(runConf, fpRate);
        }
    }

    // ===== Options and parsing logic =====

    private record MemoryOptions(Path dataRoot,
                                 String window,
                                 int windowPower,
                                 int treePower,
                                 int windowLength,
                                 int treeLength,
                                 int alphabetBase,
                                 List<Integer> ngrams,
                                 List<Double> fpRates,
                                 double runConfidence,
                                 boolean reuseSuffix) {

        static MemoryOptions parse(String[] args) {
            Path dataRoot = Path.of("data");
            String window = "w21";

            Integer windowLength = null;
            Integer treeLength = null;
            Integer windowPower = 21;
            Integer treePower = 21;
            Integer alphabetBase = 150;

            Integer defaultNgram = 8;
            List<Integer> ngramList = null;

            double fpRate = 0.15;
            List<Double> fpRates = null;
            String fpListArg = null;
            String fpGridArg = null;

            double runConfidence = 0.99;
            boolean reuseSuffix = true;

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
                    case "ngrams" -> {
                        if (value.contains(",")) {
                            ngramList = parseIntCsv(value);
                            if (ngramList.isEmpty()) throw new IllegalArgumentException("Empty --ngrams list");
                            defaultNgram = ngramList.get(0);
                        } else {
                            defaultNgram = Integer.parseInt(value);
                        }
                    }
                    case "window-length" -> windowLength = Integer.parseInt(value);
                    case "tree-length" -> treeLength = Integer.parseInt(value);
                    case "window-power", "wpow", "wpower" -> windowPower = Integer.parseInt(value);
                    case "tree-power", "tpow", "tpower" -> treePower = Integer.parseInt(value);
                    case "alphabet-base", "alphabet" -> alphabetBase = Integer.parseInt(value);
                    case "fp" -> {
                        if (value.contains(",")) {
                            fpListArg = value;
                        } else {
                            fpRate = Double.parseDouble(value);
                        }
                    }
                    case "fp-list" -> fpListArg = value;
                    case "fp-grid" -> fpGridArg = value;
                    case "confidence" -> runConfidence = Double.parseDouble(value);
                    case "reuse-suffix", "reuse-suffix-results" -> reuseSuffix = Boolean.parseBoolean(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            if (defaultNgram == null) defaultNgram = 8;
            if (ngramList == null) {
                ngramList = List.of(defaultNgram);
            } else {
                ngramList = ngramList.stream().distinct().sorted().collect(Collectors.toList());
            }

            if (fpListArg != null || fpGridArg != null) {
                List<Double> list = new ArrayList<>();
                if (fpListArg != null) list.addAll(parseFpList(fpListArg));
                if (fpGridArg != null) list.addAll(parseFpGrid(fpGridArg));
                fpRates = list.stream().distinct().sorted().collect(Collectors.toList());
                if (fpRates.isEmpty()) throw new IllegalArgumentException("Empty false positive rate list after parsing");
            } else {
                fpRates = List.of(fpRate);
            }

            if (windowLength == null) {
                int derived = deriveWindowLength(window, 1 << 21);
                windowLength = (windowPower != null) ? 1 << windowPower : derived;
            } else if (windowPower == null) {
                windowPower = 31 - Integer.numberOfLeadingZeros(windowLength);
            }

            if (treeLength == null) {
                treeLength = (treePower != null) ? 1 << treePower : windowLength;
            } else if (treePower == null) {
                treePower = 31 - Integer.numberOfLeadingZeros(treeLength);
            }

            if (windowLength < treeLength) {
                throw new IllegalArgumentException("Window length must be >= tree length");
            }
            if (windowPower != null && treePower != null && windowPower < treePower) {
                throw new IllegalArgumentException("Window power must be >= tree power");
            }

            Path resolvedDataRoot = dataRoot.resolve(window);
            if (!Files.isDirectory(resolvedDataRoot)) {
                throw new IllegalArgumentException("Data directory not found: " + resolvedDataRoot);
            }

            return new MemoryOptions(
                    resolvedDataRoot,
                    window,
                    windowPower != null ? windowPower : 21,
                    treePower != null ? treePower : 21,
                    windowLength,
                    treeLength,
                    alphabetBase,
                    ngramList,
                    fpRates,
                    runConfidence,
                    reuseSuffix);
        }

        int alphabetSizeFor(int ngram) {
            double pow = Math.pow(alphabetBase, ngram);
            double sigma = Math.min(pow, windowLength);
            if (sigma >= Integer.MAX_VALUE) return Integer.MAX_VALUE;
            return (int) sigma;
        }

        int alphabetSize() {
            return alphabetSizeFor(ngram());
        }

        int ngram() {
            return ngrams.get(0);
        }

        double fpRate() {
            return fpRates.get(0);
        }

        List<Integer> ngramList() {
            return ngrams;
        }

        List<Double> fpRateList() {
            return fpRates;
        }

        boolean reuseSuffixResults() {
            return reuseSuffix;
        }

        Path csvOutputFile() {
            boolean multiCombos = (ngrams.size() > 1) || (fpRates.size() > 1);
            if (multiCombos) {
                return Path.of("memory_%s_multi.csv".formatted(window));
            }
            return Path.of("memory_%s_%d.csv".formatted(window, ngram()));
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

        private static List<Integer> parseIntCsv(String csv) {
            String[] parts = csv.split(",");
            List<Integer> out = new ArrayList<>(parts.length);
            for (String p : parts) {
                if (!p.isBlank()) out.add(Integer.parseInt(p.trim()));
            }
            return out;
        }

        private static List<Double> parseFpList(String csv) {
            String[] parts = csv.split(",");
            List<Double> out = new ArrayList<>(parts.length);
            for (String p : parts) {
                if (!p.isBlank()) out.add(Double.parseDouble(p.trim()));
            }
            return out;
        }

        private static List<Double> parseFpGrid(String grid) {
            String[] parts = grid.split(":");
            if (parts.length != 3) throw new IllegalArgumentException("--fp-grid must be start:end:step");
            double start = Double.parseDouble(parts[0]);
            double end = Double.parseDouble(parts[1]);
            double step = Double.parseDouble(parts[2]);
            if (step <= 0.0) throw new IllegalArgumentException("Step must be positive");
            List<Double> out = new ArrayList<>();
            for (double x = start; x <= end + 1e-12; x += step) {
                double val = Math.max(start, Math.min(x, end));
                out.add(val);
            }
            return out;
        }
    }
}
