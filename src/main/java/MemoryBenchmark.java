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
import java.util.List;
import java.util.Locale;
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
 *   ngram,
 *   hbi_mem_mib,
 *   suffix_mem_mib,
 *   regexp_mem_mib,
 *   hbifp_rate
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
        csvRows.add(List.of("ngram", "hbi_mem_mib", "suffix_mem_mib", "regexp_mem_mib", "hbifp_rate"));

        // ===== HBI =====
        HBI hbi = newHbi(options);
        hbi.strides = USE_STRIDES;
        MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, options.ngram());

        long hbiBytes = GraphLayout.parseInstance(hbi).totalSize();
        double hbiMiB = hbiBytes / 1_048_576d;
        System.out.println("\n=== HBI JOL footprint ===");
        System.out.println(GraphLayout.parseInstance(hbi).toFootprint());
        System.out.printf(Locale.ROOT, "HBI total: %d B (%.3f MiB)%n", hbiBytes, hbiMiB);

        // ===== StreamingSlidingWindowIndex suffix-like baseline =====
        StreamingSlidingWindowIndex suffix = newSuffixIndex(options);
        MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, options.ngram());

        long suffixTotalBytes = GraphLayout.parseInstance(suffix).totalSize();
        long suffixCoreBytes = suffix.estimateRetainedBytesWithoutAlphabet();
        double suffixMiB = suffixCoreBytes / 1_048_576d;

        System.out.println("\n=== StreamingSlidingWindowIndex JOL footprint ===");
        System.out.println(GraphLayout.parseInstance(suffix).toFootprint());
        System.out.printf(Locale.ROOT,
                "SuffixIndex total (full, incl. AlphabetMapper): %d B (%.3f MiB)%n",
                suffixTotalBytes,
                suffixTotalBytes / 1_048_576d);
        System.out.printf(Locale.ROOT,
                "SuffixIndex core (excluding AlphabetMapper): %d B (%.3f MiB)%n",
                suffixCoreBytes,
                suffixMiB);

        // ===== Regex baseline =====
        RegexIndex regex = new RegexIndex();
        // Regex baseline indexes the raw token stream with ngram = 1
        MultiQueryExperiment.populateIndex(datasetFile.toString(), regex, 1);

        long regexBytes = GraphLayout.parseInstance(regex).totalSize();
        double regexMiB = regexBytes / 1_048_576d;
        System.out.println("\n=== Regex baseline JOL footprint ===");
        System.out.println(GraphLayout.parseInstance(regex).toFootprint());
        System.out.printf(Locale.ROOT, "Regex total: %d B (%.3f MiB)%n", regexBytes, regexMiB);

        // ===== CSV row =====
        csvRows.add(List.of(options.ngram(), hbiMiB, suffixMiB, regexMiB, options.fpRate()));
        Path out = options.csvOutputFile();
        CsvUtil.writeRows(out, csvRows);
        System.out.printf(Locale.ROOT, "Wrote memory results to %s%n", out);
    }

    /**
     * Build an HBI configured with small named Supplier classes, so there are no lambdas.
     * We snapshot primitive configuration fields out of MemoryOptions and pass those
     * into the Supplier implementations as final fields.
     */
    private static HBI newHbi(MemoryOptions options) {
        final int windowLen  = options.windowLength();
        final double fpRate  = options.fpRate();
        final int sigma      = options.alphabetSize();
        final int treeLen    = options.treeLength();
        final double runConf = options.runConfidence();
        final int nGram      = options.ngram();

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
    private static StreamingSlidingWindowIndex newSuffixIndex(MemoryOptions options) {
        int windowSize = options.windowLength();
        int expectedAlphabetSize = options.alphabetSize();
        return new StreamingSlidingWindowIndex(windowSize, expectedAlphabetSize);
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
                                 int ngram,
                                 int windowLength,
                                 int treeLength,
                                 int alphabetBase,
                                 double fpRate,
                                 double runConfidence) {

        static MemoryOptions parse(String[] args) {
            int ng = 8;
            Path dataRoot = Path.of("data");
            String window = "w21";
            Integer ngram = ng;
            Integer windowLength = null;
            Integer treeLength = null;
            Integer alphabetBase = 150;
            double fpRate = 0.15;
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

            if (ngram == null) ngram = ng;
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
