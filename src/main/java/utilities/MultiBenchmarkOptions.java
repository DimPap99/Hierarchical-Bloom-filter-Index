package utilities;

import utilities.BenchmarkEnums.QueryType;

import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.stream.Collectors;

// Parsed options for HBIDatasetBenchmarkMulti.
public record MultiBenchmarkOptions(
        Path dataRoot,
        Path queryRoot,
        String window,
        int windowPower,
        int treePower,
        QueryType queryType,
        String mode,                 // "chars" or "segments"
        Integer defaultNgram,        // primary n-gram
        List<Integer> ngrams,        // list to iterate
        int windowLength,
        int treeLength,
        int alphabetBase,
        int minQueryLength,
        int maxQueryLength,
        double fpRate,
        List<Double> fpRates,
        double runConfidence,
        int warmupRuns,
        int runs,
        String algorithm,
        boolean strides,
        boolean runSuffixTreeBaseline,
        boolean runHbi,
        boolean runSuffix,
        int suffixNgram,
        List<Integer> suffixNgramList,
        boolean suffixNgramFollowsPrimary,
        boolean baselinesMatchPrimaryNgram,
        Utils.MemPolicy memPolicy,
        int suffixDelta,
        boolean reuseSuffixResults,
        boolean reinsertPerWorkload,
        double rankEpsTarget,
        double deltaQ,
        double deltaSamp,
        double quantile,
        boolean collectStats) {

    public static MultiBenchmarkOptions parse(String[] args) {
        Path dataRoot = Path.of("data");
        Path queryRoot = Path.of("queries");
        String window = "w21";
        QueryType queryType = QueryType.UNIFORM;
        String mode = "chars";

        Integer defaultNgram = 4;
        List<Integer> ngramList = null;
        Integer windowLength = null; // derive after parsing from power or token
        Integer treeLength = null;   // derive after parsing from power or explicit
        Integer windowPower = 21;
        Integer treePower = 21;
        Integer alphabetBase = 150;
        int minQueryLength = 0;
        int maxQueryLength = 0;

        double fpRate = 0.15;
        List<Double> fpRates = null;
        String fpListArg = null;
        String fpGridArg = null;

        double runConfidence = 0.99;
        double rankEpsTarget = 0.025;
        double deltaQ = 0.05;
        double deltaSamp = 0.05;
        double quantile = 0.05;
        int warmupRuns = 0;
        int runs = 1;
        boolean runSuffixTreeBaseline = false;
        boolean runHbi = true;
        boolean runSuffix = true;
        int suffixNgram = 1;
        List<Integer> suffixNgramList = null;
        boolean suffixNgramFollowsPrimary = false;
        Utils.MemPolicy policy = Utils.MemPolicy.NONE;
        int suffixDelta = 160; // default delta for delayed suffix when not reinserting per-workload
        boolean reuseSuffixResults = true; // default: reuse suffix results across FPR/ng loops
        boolean reinsertPerWorkload = false;
        boolean baselinesMatchPrimaryNgram = false;
        String algorithm = "bs"; // default algorithm
        boolean strides = true; // default to strides enabled
        boolean collectStats = false; // default off; enable via --collect-stats=true

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
                case "policy", "mem-policy", "memory-policy" -> policy = parseMemPolicy(value);
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
                case "alphabet", "alphabet-base" -> alphabetBase = Integer.parseInt(value);
                case "min-query-length", "min-query-len", "query-length-start" -> minQueryLength = Integer.parseInt(value);
                case "max-query-length", "max-query-len", "query-length-end" -> maxQueryLength = Integer.parseInt(value);
                case "fp" -> { if (value.contains(",")) fpListArg = value; else fpRate = Double.parseDouble(value);}                
                case "fp-list" -> fpListArg = value;
                case "fp-grid" -> fpGridArg = value; // start:end:step
                case "confidence" -> runConfidence = Double.parseDouble(value);
                case "warmup" -> warmupRuns = Integer.parseInt(value);
                case "runs" -> runs = Integer.parseInt(value);
                case "algorithm", "algo" -> algorithm = value;
                case "strides", "use-strides" -> strides = Boolean.parseBoolean(value);
                case "suffix-tree", "run-suffix-tree", "regex", "run-regex" -> runSuffixTreeBaseline = Boolean.parseBoolean(value);
                case "run-hbi" -> runHbi = Boolean.parseBoolean(value);
                case "run-suffix" -> runSuffix = Boolean.parseBoolean(value);
                case "suffix-ngram" -> {
                    if ("match".equalsIgnoreCase(value) || "reuse".equalsIgnoreCase(value)) {
                        suffixNgramFollowsPrimary = true;
                        suffixNgramList = null;
                    } else if (value.contains(",")) {
                        suffixNgramList = parseIntCsv(value);
                        if (suffixNgramList.isEmpty()) {
                            throw new IllegalArgumentException("Empty --suffix-ngram list");
                        }
                        suffixNgramFollowsPrimary = false;
                    } else {
                        suffixNgram = Integer.parseInt(value);
                        suffixNgramList = null;
                        suffixNgramFollowsPrimary = false;
                    }
                }
                case "reinsert-per-workload" -> reinsertPerWorkload = Boolean.parseBoolean(value);
                case "suffix-delta" -> suffixDelta = Integer.parseInt(value);
                case "reuse-suffix", "reuse-suffix-results" -> reuseSuffixResults = Boolean.parseBoolean(value);
                case "baseline-match-ngrams", "match-baseline-ngrams", "baselines-match" -> baselinesMatchPrimaryNgram = Boolean.parseBoolean(value);
                case "epsilon", "eps", "rank-eps", "eps-target" -> rankEpsTarget = Double.parseDouble(value);
                case "delta-q" -> deltaQ = Double.parseDouble(value);
                case "delta-samp", "delta-sample" -> deltaSamp = Double.parseDouble(value);
                case "p", "quantile" -> quantile = Double.parseDouble(value);
                case "collect-stats", "stats" -> collectStats = Boolean.parseBoolean(value);
                default -> throw new IllegalArgumentException("Unknown option --" + key);
            }
        }

        if (!(rankEpsTarget > 0.0 && rankEpsTarget < 1.0)) {
            throw new IllegalArgumentException("--epsilon must be in (0,1)");
        }
        if (!(deltaQ > 0.0 && deltaQ < 1.0)) {
            throw new IllegalArgumentException("--delta-q must be in (0,1)");
        }
        if (!(deltaSamp > 0.0 && deltaSamp < 1.0)) {
            throw new IllegalArgumentException("--delta-samp must be in (0,1)");
        }
        if (!(quantile > 0.0 && quantile <= 1.0)) {
            throw new IllegalArgumentException("--p must be in (0,1]");
        }

        if (window == null) throw new IllegalArgumentException("--window is required");
        if (defaultNgram == null) defaultNgram = 4;
        // derive lengths if not provided explicitly
        if (windowLength == null) {
            // Prefer explicit power if provided; else derive from window token
            int derivedFromToken = deriveWindowLength(window, 1 << 21);
            windowLength = (windowPower != null) ? (1 << windowPower) : derivedFromToken;
        } else {
            // keep a best-effort power for reporting if not overridden
            if (windowPower == null) {
                windowPower = 31 - Integer.numberOfLeadingZeros(windowLength);
            }
        }
        if (treeLength == null) {
            treeLength = (treePower != null) ? (1 << treePower) : windowLength;
        } else if (treePower == null) {
            treePower = 31 - Integer.numberOfLeadingZeros(treeLength);
        }

        if (minQueryLength < 0) {
            throw new IllegalArgumentException("--min-query-length must be >= 0");
        }
        if (maxQueryLength < 0) {
            throw new IllegalArgumentException("--max-query-length must be >= 0");
        }
        if (maxQueryLength > 0 && minQueryLength > 0 && maxQueryLength < minQueryLength) {
            throw new IllegalArgumentException("--max-query-length must be >= --min-query-length when both are set");
        }

        // Validate relationship: window must be >= tree
        if (windowLength < treeLength) {
            throw new IllegalArgumentException(
                    "Invalid sizes: window length (" + windowLength + ") must be >= tree length (" + treeLength + ")");
        }
        if (windowPower != null && treePower != null && windowPower < treePower) {
            throw new IllegalArgumentException(
                    "Invalid powers: window power (" + windowPower + ") must be >= tree power (" + treePower + ")");
        }

        Path resolvedDataRoot = dataRoot.resolve(window);
        Path resolvedQueryRoot = queryRoot.resolve(window);
        if (!Files.isDirectory(resolvedDataRoot)) throw new IllegalArgumentException("Data directory not found: " + resolvedDataRoot);
        if (!Files.isDirectory(resolvedQueryRoot)) throw new IllegalArgumentException("Query directory not found: " + resolvedQueryRoot);
        if (!runSuffixTreeBaseline && !runHbi && !runSuffix) throw new IllegalArgumentException("At least one index must be enabled");

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
        suffixNgram = Math.max(1, suffixNgram);
        if (suffixNgramList != null) {
            suffixNgramList = suffixNgramList.stream()
                    .map(n -> Math.max(1, n))
                    .collect(Collectors.toUnmodifiableList());
        } else {
            suffixNgramList = List.of();
        }

        // normalize algorithm token
        algorithm = (algorithm == null) ? "bs" : algorithm.toLowerCase(Locale.ROOT).trim();

        if (baselinesMatchPrimaryNgram) {
            suffixNgramFollowsPrimary = true;
        }

        return new MultiBenchmarkOptions(
                resolvedDataRoot,
                resolvedQueryRoot,
                window,
                windowPower != null ? windowPower : 21,
                treePower != null ? treePower : 21,
                queryType,
                mode,
                defaultNgram,
                ngramList,
                windowLength,
                treeLength,
                alphabetBase,
                minQueryLength,
                maxQueryLength,
                fpRate,
                fpRates,
                runConfidence,
                warmupRuns,
                runs,
                algorithm,
                strides,
                runSuffixTreeBaseline,
                runHbi,
                runSuffix,
                suffixNgram,
                suffixNgramList,
                suffixNgramFollowsPrimary,
                baselinesMatchPrimaryNgram,
                policy,
                suffixDelta,
                reuseSuffixResults,
                reinsertPerWorkload,
                rankEpsTarget,
                deltaQ,
                deltaSamp,
                quantile,
                collectStats
        );
    }

    public int resolveSuffixNgram(int primaryNgramIndex, int primaryNgram) {
        if (suffixNgramFollowsPrimary) {
            return Math.max(1, primaryNgram);
        }
        if (!suffixNgramList.isEmpty()) {
            int idx = Math.min(Math.max(primaryNgramIndex, 0), suffixNgramList.size() - 1);
            return suffixNgramList.get(idx);
        }
        return suffixNgram;
    }

    public int suffixTreeNgramFor(int primaryNgramIndex, int primaryNgram) {
        if (baselinesMatchPrimaryNgram) {
            return Math.max(1, primaryNgram);
        }
        return resolveSuffixNgram(primaryNgramIndex, primaryNgram);
    }

    public int alphabetSizeFor(int ngram) {
        double pow = Math.pow(alphabetBase, ngram);
        double sigma = Math.min(pow, windowLength);
        return (sigma >= Integer.MAX_VALUE) ? Integer.MAX_VALUE : (int) sigma;
    }

    public Utils.HopsDesignResult policyDesignFor(int ngram) {
        int distinctEstimate = Math.max(1, alphabetSizeFor(ngram));
        return Utils.designBucketsForRankTargetChebyshev(distinctEstimate, rankEpsTarget, deltaQ, deltaSamp);
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

    private static Utils.MemPolicy parseMemPolicy(String token) {
        if (token == null || token.isBlank()) return Utils.MemPolicy.NONE;
        String t = token.trim().toUpperCase();
        if ("NONE".equals(t)) return Utils.MemPolicy.NONE;
        if ("REACTIVE".equals(t)) return Utils.MemPolicy.REACTIVE;
        if ("PREDICTIVE".equals(t)) return Utils.MemPolicy.PREDICTIVE;
        throw new IllegalArgumentException("Unknown mem policy '" + token + "' (expected NONE, REACTIVE, PREDICTIVE)");
    }
}
