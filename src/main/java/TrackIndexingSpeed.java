import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import estimators.CostFunctionMarkov;
import estimators.Estimator;
import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;
import search.BlockSearch;
import search.Verifier;
import search.VerifierLinearLeafProbe;
import utilities.CharRingBuffer;
import utilities.DatasetReader;
import utilities.RingBuffer;

import java.nio.file.Path;
import java.util.Locale;
import java.util.function.Supplier;

/**
 * Indexing-only benchmark: streams a dataset through the HBI index and measures
 * the time to “slide” the window without running any queries.
 */
public final class TrackIndexingSpeed {

    private static final Path DEFAULT_DATA_FILE = Path.of("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/pg2701.txt");
    private static final int DEFAULT_WINDOW_LEN = 1 << 21;
    private static final int DEFAULT_TREE_LEN = 1 << 20;
    private static final int DEFAULT_ALPHABET_BASE = 75;
    private static final double DEFAULT_FP_RATE = 0.001;
    private static final int DEFAULT_NGRAM = 3;
    private static final int DEFAULT_RUNS = 1;

    public static void main(String[] args) {
        CliOptions options = CliOptions.parse(args);

        int alphabetSize = options.alphabetSizeFor(options.nGram);
        System.out.printf(Locale.ROOT,
                "Data: %s%nN-gram: %d  Window: %d  Tree: %d  Alphabet: %d  FP: %.3g  Runs: %d%n",
                options.dataFile, options.nGram, options.windowLength, options.treeLength,
                alphabetSize, options.fpRate, options.runs);

        double totalMs = 0.0;
        long totalChars = 0L;
        long totalGrams = 0L;

        for (int i = 0; i < options.runs; i++) {
            HBI hbi = newHbi(options, alphabetSize);
            RunStats s = indexDataset(hbi, options.dataFile.toString(), options.nGram);
            totalMs += s.millis;
            totalChars += s.charsProcessed;
            totalGrams += s.gramsInserted;
        }

        if (options.runs > 0) {
            double avgMs = totalMs / options.runs;
            double avgChars = totalChars / (double) options.runs;
            double avgGrams = totalGrams / (double) options.runs;
            double gramsPerSec = (avgGrams / (avgMs / 1000.0));
            double charsPerSec = (avgChars / (avgMs / 1000.0));
            System.out.printf(Locale.ROOT,
                    "Indexing avg: %.3f ms | grams: %.0f | chars: %.0f | grams/s: %.2f | chars/s: %.2f%n",
                    avgMs, avgGrams, avgChars, gramsPerSec, charsPerSec);
        }
    }

    private static RunStats indexDataset(HBI hbi, String datasetPath, int nGram) {
        long startMs = System.currentTimeMillis();
        long chars = 0L;
        long grams = 0L;

        RingBuffer<Character> window = new CharRingBuffer(nGram);
        try (DatasetReader reader = new DatasetReader(datasetPath, datasetPath, true)) { // include line breaks, mirror Experiment
            for (char c : reader) {
                window.append(c);
                chars++;
                if (window.isFilled()) {
                    hbi.insert(window.snapshot().toString());
                    grams++;
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        long durMs = System.currentTimeMillis() - startMs;
        return new RunStats(durMs, chars, grams);
    }

    private static HBI newHbi(CliOptions options, int alphabetSize) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(options.treeLength);
        Supplier<Membership> memFactory = BloomFilter::new;
        Verifier verifier = new VerifierLinearLeafProbe();

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(options.windowLength)
                .fpRate(options.fpRate)
                .alphabetSize(alphabetSize)
                .treeLength(options.treeLength)
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(() -> null)
                .verifier(verifier)
                .costFunction(new CostFunctionMarkov())
                .confidence(0.99)
                .collectStats(false)
                .experimentMode(false)
                .build();

        return new HBI(configuration);
    }

    private static final class RunStats {
        final long millis;
        final long charsProcessed;
        final long gramsInserted;

        RunStats(long millis, long charsProcessed, long gramsInserted) {
            this.millis = millis;
            this.charsProcessed = charsProcessed;
            this.gramsInserted = gramsInserted;
        }
    }

    private static final class CliOptions {
        final Path dataFile;
        final int windowLength;
        final int treeLength;
        final int alphabetBase;
        final double fpRate;
        final int nGram;
        final int runs;

        private CliOptions(Path dataFile,
                           int windowLength,
                           int treeLength,
                           int alphabetBase,
                           double fpRate,
                           int nGram,
                           int runs) {
            this.dataFile = dataFile;
            this.windowLength = windowLength;
            this.treeLength = treeLength;
            this.alphabetBase = alphabetBase;
            this.fpRate = fpRate;
            this.nGram = nGram;
            this.runs = runs;
        }

        static CliOptions parse(String[] args) {
            Path data = DEFAULT_DATA_FILE;
            int window = DEFAULT_WINDOW_LEN;
            int tree = DEFAULT_TREE_LEN;
            int alphabet = DEFAULT_ALPHABET_BASE;
            double fp = DEFAULT_FP_RATE;
            int ngram = DEFAULT_NGRAM;
            int runs = DEFAULT_RUNS;

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
                    case "window" -> window = Integer.parseInt(value);
                    case "tree" -> tree = Integer.parseInt(value);
                    case "alphabet-base" -> alphabet = Integer.parseInt(value);
                    case "fp" -> fp = Double.parseDouble(value);
                    case "ngrams", "ngram" -> ngram = Integer.parseInt(value);
                    case "runs" -> runs = Integer.parseInt(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            return new CliOptions(data, window, tree, alphabet, fp, ngram, runs);
        }

        int alphabetSizeFor(int nGram) {
            return (int) Math.pow(alphabetBase, nGram);
        }
    }
}
