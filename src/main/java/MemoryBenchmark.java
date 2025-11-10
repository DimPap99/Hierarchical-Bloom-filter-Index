import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.RegexIndex;
import PMIndex.StreamingSlidingWindowIndex;
import PMIndex.IPMIndexing;

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
import utilities.MultiQueryExperiment;
import utilities.Utils;
import utilities.SegmentModeRunner;
import utilities.IndexFactory;

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
 *   suffix_tree_mem_mib,
 *   regex_mem_mib
 *
 * hbi_mem_mib is the retained size of the Hierarchical Bloom filter Index (HBI).
 *
 * suffix_total_mem_mib now keeps every byte reported by JOL, including the per-tree token remap
 * structures. This aligns the footprint with the actual implementation used in benchmarks.
 *
 * suffix_core_mem_mib subtracts the dictionary portion only (should be near-zero today) so legacy
 * comparisons that tracked “core” structures remain available.
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

    private record TreeSetting(int power, int length) {}

    private static final boolean USE_STRIDES = true;

    public static void main(String[] args) throws IOException {
        MemoryOptions options = MemoryOptions.parse(args);

        System.out.printf(Locale.ROOT,
                "Running memory benchmark on window %s under %s (mode=%s)%n",
                options.window(), options.dataRoot(), options.mode());

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
                "policy",
                "policy_buckets",
                "hbi_mem_mib",
                "suffix_core_mem_mib",
                "suffix_total_mem_mib",
                "suffix_tree_mem_mib",
                "regex_mem_mib"));

        Map<Integer, SuffixMeasurement> suffixCache = new HashMap<>();
        Map<Integer, SuffixTreeMeasurement> suffixTreeCache = new HashMap<>();
        RegexMeasurement regexCache = null;

        for (double fpRate : options.fpRates()) {
            System.out.printf(Locale.ROOT, "%n==== FPR = %.6f ====%n", fpRate);
            for (int ngram : options.ngrams()) {
                System.out.printf(Locale.ROOT, "---- ngram = %d ----%n", ngram);

                SuffixMeasurement suffixStats;
                if (options.reuseSuffixResults() && suffixCache.containsKey(ngram)) {
                    suffixStats = suffixCache.get(ngram);
                    System.out.printf(Locale.ROOT, "Reusing suffix measurement for suffix ngram %d%n", suffixStats.ngram());
                    System.out.printf(Locale.ROOT,
                            "SuffixIndex total (reuse, including token remaps): %d B (%.3f MiB) [ngram=%d]%n",
                            suffixStats.totalBytes(),
                            suffixStats.totalMiB(),
                            suffixStats.ngram());
                    System.out.printf(Locale.ROOT,
                            "SuffixIndex core (reuse, excluding token dictionary only): %d B (%.3f MiB) [ngram=%d]%n",
                            suffixStats.coreBytes(),
                            suffixStats.coreMiB(),
                            suffixStats.ngram());
                    System.out.printf(Locale.ROOT,
                            "  Token remap structures account for %d B (%.3f MiB) [ngram=%d]%n",
                            suffixStats.remapBytes(),
                            suffixStats.remapMiB(),
                            suffixStats.ngram());
                } else {
                    suffixStats = buildAndMeasureSuffix(datasetFile, options, ngram);
                    if (options.reuseSuffixResults()) {
                        suffixCache.put(ngram, suffixStats);
                    }
                }

                if (regexCache == null) {
                    regexCache = buildAndMeasureRegex(datasetFile, options);
                }

                SuffixTreeMeasurement suffixTreeStats;
                if (options.reuseSuffixResults() && suffixTreeCache.containsKey(ngram)) {
                    suffixTreeStats = suffixTreeCache.get(ngram);
                    System.out.printf(Locale.ROOT, "Reusing suffix-tree measurement for ngram %d%n", ngram);
                    System.out.printf(Locale.ROOT,
                            "SuffixTree total (reuse): %d B (%.3f MiB) [ngram=%d]%n",
                            suffixTreeStats.bytes(), suffixTreeStats.mib(), ngram);
                } else {
                    suffixTreeStats = buildAndMeasureSuffixTree(datasetFile, options, ngram);
                    if (options.reuseSuffixResults()) {
                        suffixTreeCache.put(ngram, suffixTreeStats);
                    }
                }

                for (TreeSetting treeSetting : options.treeSettings()) {
                    HbiMeasurement hbiStats = buildAndMeasureHbi(datasetFile, options, treeSetting, ngram, fpRate);

                    csvRows.add(List.of(
                            options.window(),
                            options.windowPower(),
                            options.windowLength(),
                            treeSetting.power(),
                            treeSetting.length(),
                            ngram,
                            fpRate,
                            options.memPolicy().name(),
                            options.policyBuckets(),
                            hbiStats.mib(),
                            suffixStats.coreMiB(),
                            suffixStats.totalMiB(),
                            suffixTreeStats.mib(),
                            regexCache.mib()));
                }
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
    private static HBI newHbi(MemoryOptions options, TreeSetting treeSetting, int nGram, double fpRate) {
        final int windowLen  = options.windowLength();
        final int sigma      = options.alphabetSizeFor(nGram);
        final int treeLen    = treeSetting.length();
        final double runConf = options.runConfidence();

        Supplier<Estimator> estFactory =
                new EstimatorFactory(treeLen);
        Supplier<Membership> memFactory =
                new MembershipFactory();
        Supplier<PruningPlan> prFactory =
                new PruningPlanFactory(runConf, fpRate);

        Verifier verifier = new VerifierLinearLeafProbe();

        int derivedBuckets = 0;
        if (options.memPolicy() != Utils.MemPolicy.NONE) {
            if (options.hasExplicitPolicyBuckets()) {
                derivedBuckets = options.policyBuckets();
            } else {
                Utils.HopsDesignResult design = options.policyDesignFor(nGram);
                if (design != null) {
                    derivedBuckets = Math.max(1, design.suggestedBuckets);
                    System.out.printf(Locale.ROOT,
                            "  Policy auto-design -> Dhat=%d, buckets=%d, n_req=%d, LB=%d%s%n",
                            Math.max(1, options.alphabetSizeFor(nGram)),
                            design.suggestedBuckets,
                            design.requiredSampleSize,
                            design.occupancyLowerBound,
                            design.impossible ? " (target unattainable with Dhat)" : "");
                }
            }
        }

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(windowLen)
                .fpRate(fpRate)
                .alphabetSize(sigma)
                .treeLength(treeLen)
                .memPolicy(options.memPolicy())
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(verifier)
                .costFunction(new CostFunctionMaxProb())
                .confidence(runConf)
                .experimentMode(false)
                .collectStats(false)
                .activateEstim(options.shouldActivateEstimators())
                .buckets(derivedBuckets)
                .nGram(nGram)
                .build();

        HBI hbi = new HBI(configuration);
        hbi.QUANTILE = options.quantile();
        return hbi;
    }

    /**
     * Create the streaming suffix baseline using the same window size and expected alphabet size.
     */
    private static StreamingSlidingWindowIndex newSuffixIndex(MemoryOptions options, int nGram) {
        int windowSize = options.windowLength();
        // Suffix index always operates on unigrams regardless of benchmarked n-gram setting.
        int expectedAlphabetSize = options.alphabetSizeFor(1);
        return new StreamingSlidingWindowIndex(windowSize, expectedAlphabetSize);
    }

    private static HbiMeasurement buildAndMeasureHbi(Path datasetFile,
                                                     MemoryOptions options,
                                                     TreeSetting treeSetting,
                                                     int ngram,
                                                     double fpRate) throws IOException {
        HBI hbi = newHbi(options, treeSetting, ngram, fpRate);
        hbi.strides = USE_STRIDES;
        if (options.isSegmentsMode()) {
            try {
                SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), hbi, ngram);
            } catch (Exception e) {
                throw new IOException("Failed to populate HBI in segments mode", e);
            }
        } else {
            MultiQueryExperiment.populateIndex(datasetFile.toString(), hbi, ngram);
        }
        if (options.memPolicy() != Utils.MemPolicy.NONE) {
            hbi.forceApplyMemoryPolicy();
        }

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
        int suffixNgram = 1;
        StreamingSlidingWindowIndex suffix = newSuffixIndex(options, ngram);
        if (options.isSegmentsMode()) {
            try {
                SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffix, suffixNgram);
            } catch (Exception e) {
                throw new IOException("Failed to populate suffix index in segments mode", e);
            }
        } else {
            MultiQueryExperiment.populateIndex(datasetFile.toString(), suffix, suffixNgram);
        }

        GraphLayout layout = GraphLayout.parseInstance(suffix);
        long layoutBytes = layout.totalSize();
        long dictBytes = suffix.estimateTokenDictionaryBytes();
        long remapBytes = suffix.estimateSuffixTreeRemapBytes();
        long totalBytes = layoutBytes;
        long coreBytes = Math.max(0L, layoutBytes - dictBytes);
        double coreMiB = coreBytes / 1_048_576d;
        double totalMiB = totalBytes / 1_048_576d;

        System.out.println("\n=== StreamingSlidingWindowIndex JOL footprint ===");
        System.out.println(layout.toFootprint());
        System.out.printf(Locale.ROOT,
                "SuffixIndex total (including token remaps): %d B (%.3f MiB) [ngram=%d]%n",
                totalBytes,
                totalMiB,
                suffixNgram);
        System.out.printf(Locale.ROOT,
                "SuffixIndex core (excluding token dictionary only): %d B (%.3f MiB) [ngram=%d]%n",
                coreBytes,
                coreMiB,
                suffixNgram);
        System.out.printf(Locale.ROOT,
                "  Token remap structures account for %d B (%.3f MiB) [ngram=%d]%n",
                remapBytes,
                remapBytes / 1_048_576d,
                suffixNgram);

        return new SuffixMeasurement(totalBytes, coreBytes, totalMiB, coreMiB, remapBytes, remapBytes / 1_048_576d, suffixNgram);
    }

    private static RegexMeasurement buildAndMeasureRegex(Path datasetFile,
                                                         MemoryOptions options) throws IOException {
        RegexIndex regex = new RegexIndex();
        if (options.isSegmentsMode()) {
            try {
                SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), regex, 1);
            } catch (Exception e) {
                throw new IOException("Failed to populate regex baseline in segments mode", e);
            }
        } else {
            MultiQueryExperiment.populateIndex(datasetFile.toString(), regex, 1);
        }

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

    private static SuffixTreeMeasurement buildAndMeasureSuffixTree(Path datasetFile,
                                                                   MemoryOptions options,
                                                                   int ngram) throws IOException {
        IPMIndexing suffixTree = IndexFactory.createSuffixTreeIndex(options.alphabetSizeFor(ngram));
        if (options.isSegmentsMode()) {
            try {
                SegmentModeRunner.insertDatasetSegments(datasetFile.toString(), suffixTree, ngram);
            } catch (Exception e) {
                throw new IOException("Failed to populate suffix-tree index in segments mode", e);
            }
        } else {
            MultiQueryExperiment.populateIndex(datasetFile.toString(), suffixTree, ngram);
        }

        GraphLayout layout = GraphLayout.parseInstance(suffixTree);
        long bytes = layout.totalSize();
        double mib = bytes / 1_048_576d;

        System.out.println("\n=== SuffixTree JOL footprint ===");
        System.out.println(layout.toFootprint());
        System.out.printf(Locale.ROOT,
                "SuffixTree total: %d B (%.3f MiB) [ngram=%d]%n",
                bytes,
                mib,
                ngram);

        return new SuffixTreeMeasurement(bytes, mib);
    }

    private record HbiMeasurement(long bytes, double mib) {}

    private record SuffixMeasurement(long totalBytes,
                                     long coreBytes,
                                     double totalMiB,
                                     double coreMiB,
                                     long remapBytes,
                                     double remapMiB,
                                     int ngram) {}

    private record RegexMeasurement(long bytes, double mib) {}
    private record SuffixTreeMeasurement(long bytes, double mib) {}

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
            return new CSEstimator(2048, 8, 1);//HashMapEstimator(treeLen);
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
                                 List<TreeSetting> treeSettings,
                                 int windowLength,
                                 int alphabetBase,
                                 String mode,
                                 List<Integer> ngrams,
                                 List<Double> fpRates,
                                 double runConfidence,
                                 double rankEpsTarget,
                                 double deltaQ,
                                 double deltaSamp,
                                 double quantile,
                                 boolean reuseSuffix,
                                 Utils.MemPolicy memPolicy,
                                 int policyBuckets) {

        static MemoryOptions parse(String[] args) {
            Path dataRoot = Path.of("data");
            String window = "w21";

            Integer windowLength = null;
            Integer treeLength = null;
            Integer windowPower = 21;
            Integer treePower = 21;
            List<Integer> treePowerList = null;
            Integer alphabetBase = 150;
            String mode = "chars";

            Integer defaultNgram = 8;
            List<Integer> ngramList = null;

            double fpRate = 0.15;
            List<Double> fpRates = null;
            String fpListArg = null;
            String fpGridArg = null;

            double runConfidence = 0.99;
            double rankEpsTarget = 0.025;
            double deltaQ = 0.05;
            double deltaSamp = 0.05;
            double quantile = 0.05;
            boolean reuseSuffix = true;
            Utils.MemPolicy memPolicy = Utils.MemPolicy.NONE;
            Integer policyBuckets = null;

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
                    case "tree-power", "tpow", "tpower" -> {
                        if (value.contains(",")) {
                            treePowerList = parseIntCsv(value);
                        } else {
                            treePower = Integer.parseInt(value);
                        }
                    }
                    case "tree-powers" -> treePowerList = parseIntCsv(value);
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
                    case "epsilon", "eps", "rank-eps", "eps-target" -> rankEpsTarget = Double.parseDouble(value);
                    case "delta-q" -> deltaQ = Double.parseDouble(value);
                    case "delta-samp", "delta-sample" -> deltaSamp = Double.parseDouble(value);
                    case "p", "quantile" -> quantile = Double.parseDouble(value);
                    case "mode" -> mode = value;
                    case "reuse-suffix", "reuse-suffix-results" -> reuseSuffix = Boolean.parseBoolean(value);
                    case "policy", "mem-policy", "memory-policy" -> memPolicy = parseMemPolicy(value);
                    case "policy-buckets", "mem-policy-buckets" -> policyBuckets = Integer.parseInt(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            if (!("chars".equalsIgnoreCase(mode) || "segments".equalsIgnoreCase(mode))) {
                throw new IllegalArgumentException("mode must be 'chars' or 'segments'");
            }

            if (!(rankEpsTarget > 0.0 && rankEpsTarget < 1.0)) {
                throw new IllegalArgumentException("epsilon must be in (0,1)");
            }
            if (!(deltaQ > 0.0 && deltaQ < 1.0)) {
                throw new IllegalArgumentException("delta-q must be in (0,1)");
            }
            if (!(deltaSamp > 0.0 && deltaSamp < 1.0)) {
                throw new IllegalArgumentException("delta-samp must be in (0,1)");
            }
            if (!(quantile > 0.0 && quantile < 1.0)) {
                throw new IllegalArgumentException("p must be in (0,1)");
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

            List<TreeSetting> treeSettings;
            if (treePowerList != null) {
                treePowerList = treePowerList.stream().distinct().sorted().collect(Collectors.toList());
                if (treePowerList.isEmpty()) {
                    throw new IllegalArgumentException("Empty --tree-power list after parsing");
                }
                if (treeLength != null && treePowerList.size() > 1) {
                    throw new IllegalArgumentException("Cannot specify --tree-length when multiple tree powers are provided");
                }
                treeSettings = new ArrayList<>(treePowerList.size());
                for (int tp : treePowerList) {
                    if (tp < 0) {
                        throw new IllegalArgumentException("Tree power must be non-negative: " + tp);
                    }
                    if (tp >= 31) {
                        throw new IllegalArgumentException("Tree power too large: " + tp);
                    }
                    int derivedLength = (treeLength != null) ? treeLength : 1 << tp;
                    if (windowLength < derivedLength) {
                        throw new IllegalArgumentException(
                                "Window length must be >= tree length for tree power " + tp);
                    }
                    if (windowPower != null && windowPower < tp) {
                        throw new IllegalArgumentException("Window power must be >= tree power (" + tp + ")");
                    }
                    treeSettings.add(new TreeSetting(tp, derivedLength));
                }
            } else {
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

                int resolvedPower = (treePower != null) ? treePower : 21;
                treeSettings = List.of(new TreeSetting(resolvedPower, treeLength));
            }

            Path resolvedDataRoot = dataRoot.resolve(window);
            if (!Files.isDirectory(resolvedDataRoot)) {
                throw new IllegalArgumentException("Data directory not found: " + resolvedDataRoot);
            }

            return new MemoryOptions(
                    resolvedDataRoot,
                    window,
                    windowPower != null ? windowPower : 21,
                    List.copyOf(treeSettings),
                    windowLength,
                    alphabetBase,
                    mode,
                    ngramList,
                    fpRates,
                    runConfidence,
                    rankEpsTarget,
                    deltaQ,
                    deltaSamp,
                    quantile,
                    reuseSuffix,
                    memPolicy,
                    policyBuckets != null ? Math.max(0, policyBuckets) : 0);
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

        boolean isSegmentsMode() {
            return "segments".equalsIgnoreCase(mode);
        }

        private TreeSetting primaryTree() {
            return treeSettings.get(0);
        }

        int treePower() {
            return primaryTree().power();
        }

        int treeLength() {
            return primaryTree().length();
        }

        boolean shouldActivateEstimators() {
            return memPolicy != Utils.MemPolicy.NONE;
        }

        public double rankEpsTarget() {
            return rankEpsTarget;
        }

        public double deltaQ() {
            return deltaQ;
        }

        public double deltaSamp() {
            return deltaSamp;
        }

        public double quantile() {
            return quantile;
        }

        public int policyBuckets() {
            return policyBuckets;
        }

        boolean hasExplicitPolicyBuckets() {
            return policyBuckets > 0;
        }

        Utils.HopsDesignResult policyDesignFor(int ngram) {
            if (memPolicy == Utils.MemPolicy.NONE) return null;
            int distinctEstimate = Math.max(1, alphabetSizeFor(ngram));
            return Utils.designBucketsForRankTargetChebyshev(distinctEstimate, rankEpsTarget, deltaQ, deltaSamp);
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
            boolean multiCombos = (ngrams.size() > 1) || (fpRates.size() > 1) || (treeSettings.size() > 1);
            String policyToken = memPolicy.name().toLowerCase(Locale.ROOT);
            if (multiCombos) {
                return Path.of("memory_%s_%s_multi.csv".formatted(window, policyToken));
            }
            return Path.of("memory_%s_%d_%s.csv".formatted(window, ngram(), policyToken));
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

        private static Utils.MemPolicy parseMemPolicy(String token) {
            if (token == null || token.isBlank()) return Utils.MemPolicy.NONE;
            return switch (token.toUpperCase(Locale.ROOT)) {
                case "NONE" -> Utils.MemPolicy.NONE;
                case "REACTIVE" -> Utils.MemPolicy.REACTIVE;
                case "PREDICTIVE" -> Utils.MemPolicy.PREDICTIVE;
                default -> throw new IllegalArgumentException("Unknown memory policy: " + token);
            };
        }
    }
}
//can
