import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;

import estimators.CostFunctionMaxProb;
import estimators.Estimator;
import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;
import search.BlockSearch;
import search.MostFreqPruning;
import search.PruningPlan;
import search.Verifier;
import search.VerifierLinearLeafProbe;

import java.io.IOException;
import java.nio.file.Path;
import java.util.List;
import java.util.Locale;
import java.util.function.Supplier;

/**
 * Simple driver that benchmarks both the refactored HBI and the legacy RegexIndex
 * on the same data / query workload. Command-line switches allow overriding the
 * default dataset, query file, and tuning knobs without editing the source.
 */
public final class Main {

    private static final Path DEFAULT_DATA_FILE = Path.of("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/pg2701.txt");
    private static final Path DEFAULT_QUERIES_FILE = Path.of("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/pg2701/unique_substrings_pg2701_10.txt");
    private static final int DEFAULT_WINDOW_LEN = 1 << 21;
    private static final int DEFAULT_TREE_LEN = 1 << 20;
    private static final int DEFAULT_ALPHABET_BASE = 75;
    private static final double DEFAULT_FP_RATE = 0.001;
    private static final int DEFAULT_RUNS = 1;
    private static final int DEFAULT_MIN_NGRAM = 1;
    private static final int DEFAULT_MAX_NGRAM = 3;
    private static final double WARMUP_CONFIDENCE = 0.999;
    private static final double RUN_CONFIDENCE = 0.99;

    public static void main(String[] args) throws IOException {
        CliOptions options = CliOptions.parse(args);

        System.out.printf(Locale.ROOT,
                "Data file: %s%nQueries file: %s%nWindow: %d  Tree: %d  Alphabet base: %d  FP: %.3g%n",
                options.dataFile, options.queriesFile, options.windowLength, options.treeLength, options.alphabetBase, options.fpRate);

        for (int nGram = options.minNgram; nGram <= options.maxNgram; nGram++) {
            double hbiTotalMs = 0;
            double ipmTotalMs = 0;
            double hbiTotalMsInsert = 0;
            double ipmTotalMsInsert = 0;
            double avgLp = 0;
            double avgAlpha = 0;

            int alphabetSize = options.alphabetSizeFor(nGram);

            System.out.println("N-gram: " + nGram);
            System.out.println("Window Size: " + options.windowLength);
            System.out.println("Tree Length: " + options.treeLength);
            System.out.println("Alphabet: " + alphabetSize);

            for (int i = 0; i < options.runs; i++) {
                HBI hbi = newHbi(options, alphabetSize, WARMUP_CONFIDENCE);
                hbi.stats().setCollecting(true);

                Experiment.run(options.dataFile.toString(), options.queriesFile.toString(), hbi, nGram, false, false);

                IPMIndexing ipm = new RegexIndex();
                Experiment.run(options.dataFile.toString(), options.queriesFile.toString(), ipm, 1, false, false);

                avgLp = average(hbi.stats().lpLevels());
                avgAlpha = average(hbi.stats().alphas());
            }

            if (options.runs > 0) {
                System.out.printf(Locale.ROOT, "Avg LP: %.3f%n", avgLp);
                System.out.printf(Locale.ROOT, "Avg Alpha: %.3f%n", avgAlpha);
            }
        }
    }

    private static HBI newHbi(CliOptions options, int alphabetSize, double conf) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(options.treeLength);
        Supplier<Membership> memFactory = BloomFilter::new;
        Supplier<PruningPlan> prFactory = () -> new MostFreqPruning(conf);
        Verifier verifier = new VerifierLinearLeafProbe();

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(options.windowLength)
                .fpRate(options.fpRate)
                .alphabetSize(alphabetSize)
                .treeLength(options.treeLength)
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(verifier)
                .costFunction(new CostFunctionMaxProb())
                .confidence(conf)
                .build();

        return new HBI(configuration);
    }

    private static double average(List<? extends Number> values) {
        if (values.isEmpty()) {
            return 0.0;
        }
        double sum = 0.0;
        for (Number value : values) {
            sum += value.doubleValue();
        }
        return sum / values.size();
    }

    private static final class CliOptions {
        final Path dataFile;
        final Path queriesFile;
        final int windowLength;
        final int treeLength;
        final int alphabetBase;
        final double fpRate;
        final int runs;
        final int minNgram;
        final int maxNgram;

        private CliOptions(Path dataFile,
                           Path queriesFile,
                           int windowLength,
                           int treeLength,
                           int alphabetBase,
                           double fpRate,
                           int runs,
                           int minNgram,
                           int maxNgram) {
            this.dataFile = dataFile;
            this.queriesFile = queriesFile;
            this.windowLength = windowLength;
            this.treeLength = treeLength;
            this.alphabetBase = alphabetBase;
            this.fpRate = fpRate;
            this.runs = runs;
            this.minNgram = minNgram;
            this.maxNgram = maxNgram;
        }

        static CliOptions parse(String[] args) {
            Path data = DEFAULT_DATA_FILE;
            Path queries = DEFAULT_QUERIES_FILE;
            int window = DEFAULT_WINDOW_LEN;
            int tree = DEFAULT_TREE_LEN;
            int alphabet = DEFAULT_ALPHABET_BASE;
            double fp = DEFAULT_FP_RATE;
            int runs = DEFAULT_RUNS;
            int minNgram = DEFAULT_MIN_NGRAM;
            int maxNgram = DEFAULT_MAX_NGRAM;

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
                    case "data" -> data = Path.of(value);
                    case "queries" -> queries = Path.of(value);
                    case "window" -> window = Integer.parseInt(value);
                    case "tree" -> tree = Integer.parseInt(value);
                    case "alphabet-base" -> alphabet = Integer.parseInt(value);
                    case "fp" -> fp = Double.parseDouble(value);
                    case "runs" -> runs = Integer.parseInt(value);
                    case "ngrams" -> {
                        String[] bounds = value.split("[-:]");
                        if (bounds.length == 2) {
                            minNgram = Integer.parseInt(bounds[0]);
                            maxNgram = Integer.parseInt(bounds[1]);
                        } else {
                            int parsed = Integer.parseInt(value);
                            minNgram = parsed;
                            maxNgram = parsed;
                        }
                    }
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            return new CliOptions(data, queries, window, tree, alphabet, fp, runs, minNgram, maxNgram);
        }

        int alphabetSizeFor(int nGram) {
            return (int) Math.pow(alphabetBase, nGram);
        }
    }
}
