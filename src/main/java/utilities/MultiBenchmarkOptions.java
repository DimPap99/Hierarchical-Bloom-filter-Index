package utilities;

import utilities.BenchmarkEnums.QueryType;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Parsed options for HBIDatasetBenchmarkMulti, extracted for reuse.
 */
public record MultiBenchmarkOptions(
        Path dataRoot,
        Path queryRoot,
        String window,
        QueryType queryType,
        String mode,                 // "chars" or "segments"
        Integer defaultNgram,        // primary n-gram
        List<Integer> ngrams,        // list to iterate
        int windowLength,
        int treeLength,
        int alphabetBase,
        double fpRate,
        List<Double> fpRates,
        double runConfidence,
        int warmupRuns,
        int runs,
        String algorithm,
        boolean runRegexBaseline,
        boolean runHbi,
        boolean runSuffix,
        Utils.MemPolicy memPolicy,
        int suffixDelta,
        boolean reinsertPerWorkload) {

    public static MultiBenchmarkOptions parse(String[] args) {
        Path dataRoot = Path.of("data");
        Path queryRoot = Path.of("queries");
        String window = "w21";
        QueryType queryType = QueryType.UNIFORM;
        String mode = "chars";

        Integer defaultNgram = 4;
        List<Integer> ngramList = null;
        Integer windowLength = 1 << 21;
        Integer treeLength = 1 << 21;
        Integer alphabetBase = 150;

        double fpRate = 0.15;
        List<Double> fpRates = null;
        String fpListArg = null;
        String fpGridArg = null;

        double runConfidence = 0.99;
        int warmupRuns = 0;
        int runs = 1;
        boolean runRegexBaseline = false;
        boolean runHbi = true;
        boolean runSuffix = true;
        Utils.MemPolicy policy = Utils.MemPolicy.NONE;
        int suffixDelta = 160; // default delta for delayed suffix when not reinserting per-workload
        boolean reinsertPerWorkload = false;
        String algorithm = "bs"; // default algorithm

        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if (!arg.startsWith("--")) continue;
            String key; String value;
            int eq = arg.indexOf('=');
            if (eq >= 0) { key = arg.substring(2, eq); value = arg.substring(eq + 1);} else {
                key = arg.substring(2);
                if (i + 1 >= args.length) throw new IllegalArgumentException("Missing value for option --" + key);
                value = args[++i];
            }

            switch (key) {
                case "window" -> window = value;
                case "type" -> queryType = "all".equalsIgnoreCase(value) ? null : QueryType.fromString(value);
                case "data-root" -> dataRoot = Path.of(value);
                case "query-root" -> queryRoot = Path.of(value);
                case "mode" -> mode = value;
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
                case "alphabet-base" -> alphabetBase = Integer.parseInt(value);
                case "fp" -> { if (value.contains(",")) fpListArg = value; else fpRate = Double.parseDouble(value);}                
                case "fp-list" -> fpListArg = value;
                case "fp-grid" -> fpGridArg = value; // start:end:step
                case "confidence" -> runConfidence = Double.parseDouble(value);
                case "warmup" -> warmupRuns = Integer.parseInt(value);
                case "runs" -> runs = Integer.parseInt(value);
                case "algorithm", "algo" -> algorithm = value;
                case "regex", "run-regex" -> runRegexBaseline = Boolean.parseBoolean(value);
                case "run-hbi" -> runHbi = Boolean.parseBoolean(value);
                case "run-suffix" -> runSuffix = Boolean.parseBoolean(value);
                case "reinsert-per-workload" -> reinsertPerWorkload = Boolean.parseBoolean(value);
                case "suffix-delta" -> suffixDelta = Integer.parseInt(value);
                default -> throw new IllegalArgumentException("Unknown option --" + key);
            }
        }

        if (window == null) throw new IllegalArgumentException("--window is required");
        if (defaultNgram == null) defaultNgram = 4;
        if (windowLength == null) windowLength = deriveWindowLength(window, 1 << 21);
        if (treeLength == null) treeLength = windowLength;

        Path resolvedDataRoot = dataRoot.resolve(window);
        Path resolvedQueryRoot = queryRoot.resolve(window);
        if (!Files.isDirectory(resolvedDataRoot)) throw new IllegalArgumentException("Data directory not found: " + resolvedDataRoot);
        if (!Files.isDirectory(resolvedQueryRoot)) throw new IllegalArgumentException("Query directory not found: " + resolvedQueryRoot);
        if (!runRegexBaseline && !runHbi && !runSuffix) throw new IllegalArgumentException("At least one index must be enabled");

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

        fpRate = fpRates.get(0);

        // normalize algorithm token
        algorithm = (algorithm == null) ? "bs" : algorithm.toLowerCase(Locale.ROOT).trim();

        return new MultiBenchmarkOptions(
                resolvedDataRoot,
                resolvedQueryRoot,
                window,
                queryType,
                mode,
                defaultNgram,
                ngramList,
                windowLength,
                treeLength,
                alphabetBase,
                fpRate,
                fpRates,
                runConfidence,
                warmupRuns,
                runs,
                algorithm,
                runRegexBaseline,
                runHbi,
                runSuffix,
                policy,
                suffixDelta,
                reinsertPerWorkload
        );
    }

    public int alphabetSizeFor(int ngram) {
        double pow = Math.pow(alphabetBase, ngram);
        double sigma = Math.min(pow, windowLength);
        return (sigma >= Integer.MAX_VALUE) ? Integer.MAX_VALUE : (int) sigma;
    }

    public String queryTypeLabel() { return queryType != null ? queryType.fileToken() : "all"; }

    public static int deriveWindowLength(String window, int def) {
        if (window == null || window.isEmpty()) return def;
        for (int i = 0; i < window.length(); i++) {
            if (Character.isDigit(window.charAt(i))) {
                String numeric = window.substring(i);
                try {
                    int power = Integer.parseInt(numeric);
                    return 1 << power;
                } catch (NumberFormatException ignored) { return def; }
            }
        }
        return def;
    }

    // helpers
    private static List<Integer> parseIntCsv(String csv) {
        String[] parts = csv.split(",");
        List<Integer> out = new ArrayList<>(parts.length);
        for (String p : parts) if (!p.isBlank()) out.add(Integer.parseInt(p.trim()));
        return out;
    }
    private static List<Double> parseFpList(String csv) {
        String[] parts = csv.split(",");
        List<Double> out = new ArrayList<>(parts.length);
        for (String p : parts) if (!p.isBlank()) out.add(Double.parseDouble(p.trim()));
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
