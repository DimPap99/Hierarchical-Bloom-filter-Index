package datagenerators;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Locale;

// CLI entry point for generating Markov-chain datasets.
// Uses MarkovChainGenerator and writes the stream to disk.
public final class GenerateMarkovDataset {

    private GenerateMarkovDataset() {
        // no instances
    }

    public static void main(String[] args) throws IOException {
        CliOptions options = CliOptions.parse(args);

        int minCodePoint = options.startCodePoint;
        int maxCodePoint = minCodePoint + options.alphabetSize;

        MarkovChainGenerator.MarkovDataset dataset = MarkovChainGenerator.generate(
                options.length,
                minCodePoint,
                maxCodePoint,
                options.order,
                options.seed);

        Path absolute = options.outputFile.toAbsolutePath();
        Path parent = absolute.getParent();
        if (parent != null) {
            Files.createDirectories(parent);
        }
        Files.writeString(absolute, dataset.data(), StandardCharsets.UTF_8);

        System.out.printf(Locale.ROOT,
                "Generated %d characters (order %d, alphabet %d) -> %s%n",
                options.length,
                options.order,
                options.alphabetSize,
                absolute);
    }

    private static final class CliOptions {
        final int length;
        final int alphabetSize;
        final int startCodePoint;
        final int order;
        final long seed;
        final Path outputFile;

        private CliOptions(int length,
                           int alphabetSize,
                           int startCodePoint,
                           int order,
                           long seed,
                           Path outputFile) {
            this.length = length;
            this.alphabetSize = alphabetSize;
            this.startCodePoint = startCodePoint;
            this.order = order;
            this.seed = seed;
            this.outputFile = outputFile;
        }

        static CliOptions parse(String[] args) {
            int length = 1 << 16;
            int alphabetSize = 40;
            int startCodePoint = 70;
            int order = 1;
            long seed = System.nanoTime();
            Path output = Path.of("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/wmarkov/1/"+order+"_Wmarkov.txt");

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
                    case "length" -> length = parseInt(value, "length");
                    case "alphabet" -> alphabetSize = parseInt(value, "alphabet");
                    case "start" -> startCodePoint = parseCodePoint(value);
                    case "order" -> order = parseInt(value, "order");
                    case "seed" -> seed = Long.parseLong(value);
                    case "output" -> output = Path.of(value);
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }

            if (alphabetSize <= 0) {
                throw new IllegalArgumentException("alphabet size must be positive");
            }
            if ((long) startCodePoint + alphabetSize > Character.MAX_VALUE + 1L) {
                throw new IllegalArgumentException("alphabet exceeds character range");
            }
            return new CliOptions(length, alphabetSize, startCodePoint, order, seed, output);
        }

        private static int parseInt(String value, String label) {
            try {
                return Integer.parseInt(value);
            } catch (NumberFormatException ex) {
                throw new IllegalArgumentException("Invalid value for " + label + ": " + value, ex);
            }
        }

        private static int parseCodePoint(String value) {
            if (value.length() == 1) {
                return value.charAt(0);
            }
            return parseInt(value, "start");
        }
    }
}
