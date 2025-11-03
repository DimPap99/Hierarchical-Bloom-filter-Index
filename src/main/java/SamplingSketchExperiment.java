import estimators.BottomKSampler;
import estimators.HOPS;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;
import utilities.CsvUtil;
import utilities.TokenHasher;
import utilities.Utils;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.Random;

import static estimators.HOPS.mix64;

    public class SamplingSketchExperiment {

        // =========================
        //  Experiment configuration
        // =========================

        private enum Dist {UNIFORM, ZIPF}
        private enum DatasetMode { SEGMENTS, WORDS, CHARS }

        private static final List<Double> DEFAULT_DELTA_QS = List.of(0.05, 0.075, 0.01);
        private static final List<Double> DEFAULT_DELTA_SAMPS = List.of(0.05, 0.025, 0.09);
        private static final String DEFAULT_CSV_FILENAME = "sampling_results.csv";

        private static double EXPONENT;
        private static double SUMMARY_DELTA = Double.NaN;

        private static final class Options {
            long streamLen = 10_000_000L;
            boolean streamLenProvided = false;
            int alphabet = 10_000;
            boolean alphabetProvided = false;
            int buckets = 2_500;
            Dist dist = Dist.ZIPF;
            double zipfExponent = 1.5;
            List<Double> quantiles = List.of(0.05);
            boolean quantilesProvided = false;
            long seed = 123_456_789L;
            int trials = 50;
            boolean verbose = false;
            boolean autoDesign = true;
            Integer distinctGuess = null;
            boolean distinctProvided = false;
            List<Double> rankEpsTargets = List.of(0.05);
            boolean rankEpsProvided = false;
            List<Double> deltaQs = DEFAULT_DELTA_QS;
            boolean deltaQsProvided = false;
            List<Double> deltaSamps = DEFAULT_DELTA_SAMPS;
            boolean deltaSampsProvided = false;
            Path csvPath = Paths.get(DEFAULT_CSV_FILENAME);
            boolean csvAppend = true;
            Path datasetPath = null;
            DatasetMode datasetMode = DatasetMode.SEGMENTS;
            int datasetNgram = 1;
            double summaryDelta = Double.NaN;
            boolean summaryDeltaProvided = false;
            String datasetLabel = "synthetic";
        }

        private record DatasetStream(ArrayList<Long> tokens, int alphabetSize) {}

        private static final class TokenEncoder {
            long encode(String token) {
                return TokenHasher.hashToPositiveLong(token);
            }
        }

        // =========================
        //  Core utilities
        // =========================

        /**
         * Percent error on values using max(1, truth) to avoid division by zero.
         */
        private static double percentValueError(int est, int truth) {
            int denom = Math.max(1, truth);
            return Math.abs(est - truth) / (double) denom;
        }
        /** Suite-A/B CSV header in the exact order you’ve been using. */
        // --- replace suiteACsvHeader() with this ---
        private static List<?> suiteACsvHeader() {
            return List.of(
                    "dataset","dataset_mode","stream_len","alphabet","buckets_configured",
                    "source_dist","zipf_exponent","quantile_target",
                    "eps_target","delta_q","delta_samp","delta_tot","Dhat",
                    "B_suggested","n_req","E_nb","Var_nb","LB_nb",
                    "Sactive","Bused","nb_MPQ","eps_DKW_MPQ",
                    "p_lo_MPQ","p_hi_MPQ","p_band_width_MPQ",
                    "x_lo_MPQ","x_hi_MPQ","x_band_width_MPQ",
                    "MPQ_xhat","MPQ_value_error_pct","MPQ_achieved_percentile","MPQ_rank_error",
                    "inside_band_MPQ","occupancy_ok","union_success",
                    "ns_per_insert_MPQ","ns_close_MPQ",
                    "nb_BK","eps_DKW_BK","x_lo_BK","x_hi_BK",
                    "BK_xhat","BK_value_error_pct","BK_achieved_percentile","BK_rank_error",
                    "ns_per_insert_BK","ns_close_BK",
                    "truth_x_p","mpq_time_ms","bk_time_ms","exponent"
            );
        }

        // --- replace suiteACsvRow(...) with this ---
        private static List<?> suiteACsvRow(TrialResult r, double exponent) {
            // design fields (when auto-design is ON)
            int    B_suggested = (r.design != null) ? r.design.suggestedBuckets    : r.Bused;
            int    n_req       = (r.design != null) ? r.design.requiredSampleSize : -1;
            double E_nb        = (r.design != null) ? r.design.expectedNonEmpty    : Double.NaN;
            double Var_nb      = (r.design != null) ? r.design.variance           : Double.NaN;
            int    LB_nb       = (r.design != null) ? r.design.occupancyLowerBound : -1;

            boolean unionSuccess = r.mpqInDKWValueBand && r.occLBMet;
            double  deltaTot     = r.delta + r.deltaSampUsed;
            double  pBandWidth   = r.pHiMPQ - r.pLoMPQ;
            int     xBandWidth   = r.xHiMPQ - r.xLoMPQ;

            return List.of(
                    r.datasetLabel, r.datasetMode, r.streamLength, r.alphabetConfigured, r.bucketsConfigured,
                    r.distributionLabel, r.zipfExponent, r.pTarget,
                    // knobs
                    r.epsTargetUsed, r.delta, r.deltaSampUsed, deltaTot, r.DhatUsed,
                    // design
                    B_suggested, n_req, E_nb, Var_nb, LB_nb,
                    // sizes
                    r.Sactive, r.Bused, r.nbMPQ, r.epsMPQ,
                    // bands (MPQ)
                    r.pLoMPQ, r.pHiMPQ, pBandWidth,
                    r.xLoMPQ, r.xHiMPQ, xBandWidth,
                    // MPQ estimate & checks
                    r.xMPQ, r.pctErrMPQ, r.pAchMPQ, r.rankErrMPQ,
                    r.mpqInDKWValueBand ? "TRUE" : "FALSE",
                    r.occLBMet ? "TRUE" : "FALSE",
                    unionSuccess ? "TRUE" : "FALSE",
                    // costs (MPQ)
                    r.nsPerInsertMPQ, r.closeNsPerItemMPQ,
                    // BK side
                    r.nbBK, r.epsBK, r.xLoBK, r.xHiBK,
                    r.xBK, r.pctErrBK, r.pAchBK, r.rankErrBK,
                    r.nsPerInsertBK, r.closeNsPerItemBK,
                    // truth + times
                    r.xTrue, (int) r.tMPQMs, (int) r.tBKMs, exponent
            );
        }

        private static Options parseArgs(String[] args) {
            if (args.length > 0 && !args[0].startsWith("--")) {
                return parseLegacyArgs(args);
            }

            Options opts = new Options();

            for (int i = 0; i < args.length; i++) {
                String flag = args[i];
                String value;
                switch (flag) {
                    case "--stream-len" -> {
                        value = requireValue(flag, args, ++i);
                        opts.streamLen = Long.parseLong(value);
                        opts.streamLenProvided = true;
                    }
                    case "--alphabet" -> {
                        value = requireValue(flag, args, ++i);
                        opts.alphabet = Integer.parseInt(value);
                        opts.alphabetProvided = true;
                    }
                    case "--buckets" -> opts.buckets = Integer.parseInt(requireValue(flag, args, ++i));
                    case "--dist" -> opts.dist = Dist.valueOf(requireValue(flag, args, ++i).toUpperCase(Locale.ROOT));
                    case "--zipf" -> opts.zipfExponent = Double.parseDouble(requireValue(flag, args, ++i));
                    case "--quantile" -> {
                        value = requireValue(flag, args, ++i);
                        opts.quantiles = parseDoubleList(value);
                        if (opts.quantiles.isEmpty()) {
                            throw new IllegalArgumentException("--quantile list may not be empty");
                        }
                        opts.quantilesProvided = true;
                    }
                    case "--seed" -> opts.seed = Long.parseLong(requireValue(flag, args, ++i));
                    case "--trials" -> opts.trials = Integer.parseInt(requireValue(flag, args, ++i));
                    case "--verbose" -> opts.verbose = Boolean.parseBoolean(requireValue(flag, args, ++i));
                    case "--auto-design" -> opts.autoDesign = Boolean.parseBoolean(requireValue(flag, args, ++i));
                    case "--distinct" -> {
                        value = requireValue(flag, args, ++i);
                        opts.distinctGuess = Integer.parseInt(value);
                        opts.distinctProvided = true;
                    }
                    case "--rank-eps" -> {
                        value = requireValue(flag, args, ++i);
                        opts.rankEpsTargets = parseDoubleList(value);
                        if (opts.rankEpsTargets.isEmpty()) {
                            throw new IllegalArgumentException("--rank-eps list may not be empty");
                        }
                        opts.rankEpsProvided = true;
                    }
                    case "--delta" -> {
                        value = requireValue(flag, args, ++i);
                        opts.summaryDelta = Double.parseDouble(value);
                        opts.summaryDeltaProvided = true;
                    }
                    case "--delta-q" -> {
                        value = requireValue(flag, args, ++i);
                        opts.deltaQs = parseDoubleList(value);
                        opts.deltaQsProvided = true;
                    }
                    case "--delta-samp" -> {
                        value = requireValue(flag, args, ++i);
                        opts.deltaSamps = parseDoubleList(value);
                        opts.deltaSampsProvided = true;
                    }
                    case "--csv" -> {
                        value = requireValue(flag, args, ++i);
                        opts.csvPath = value.equalsIgnoreCase("none") ? null : Paths.get(value);
                    }
                    case "--csv-append" -> opts.csvAppend = Boolean.parseBoolean(requireValue(flag, args, ++i));
                    case "--dataset" -> opts.datasetPath = Paths.get(requireValue(flag, args, ++i));
                    case "--dataset-mode" -> opts.datasetMode = parseDatasetMode(requireValue(flag, args, ++i));
                    case "--dataset-ngram" -> opts.datasetNgram = parsePositiveInt(flag, requireValue(flag, args, ++i));
                    case "--help" -> {
                        printUsage();
                        System.exit(0);
                    }
                    default -> throw new IllegalArgumentException("Unknown option '" + flag + "'. Use --help for usage.");
                }
            }

            if (opts.summaryDeltaProvided && !opts.deltaQsProvided) {
                opts.deltaQs = List.of(opts.summaryDelta);
                opts.deltaQsProvided = true;
            }
            if (!opts.deltaSampsProvided && opts.deltaQsProvided && opts.deltaQs.size() == 1) {
                opts.deltaSamps = List.of(DEFAULT_DELTA_SAMPS.get(0));
            }

            if (opts.deltaQs.isEmpty()) {
                opts.deltaQs = DEFAULT_DELTA_QS;
            }
            if (opts.deltaSamps.isEmpty()) {
                opts.deltaSamps = DEFAULT_DELTA_SAMPS;
            }

            if (opts.deltaQs.size() != opts.deltaSamps.size()) {
                throw new IllegalArgumentException("delta_q and delta_samp must have the same number of entries");
            }

            if (opts.summaryDeltaProvided) {
                SUMMARY_DELTA = opts.summaryDelta;
            } else if (!opts.deltaQs.isEmpty()) {
                SUMMARY_DELTA = opts.deltaQs.size() == 1 ? opts.deltaQs.get(0) : Double.NaN;
            }

            return opts;
        }

        private static Options parseLegacyArgs(String[] args) {
            Options opts = new Options();

            long N = opts.streamLen;
            int A = opts.alphabet;
            int B = opts.buckets;
            Dist dist = opts.dist;
            double s = opts.zipfExponent;
            double p = opts.quantiles.get(0);
            long seed = opts.seed;
            double delta = 0.01;
            int trials = opts.trials;
            boolean verbose = opts.verbose;
            boolean autoDesignB = opts.autoDesign;
            Integer distinctGuess = null;
            Double rankEpsTarget = opts.rankEpsTargets.get(0);
            Double deltaSamp = DEFAULT_DELTA_SAMPS.get(0);
            String csvPath = DEFAULT_CSV_FILENAME;
            boolean csvAppend = true;

            if (args.length >= 1) { N = Long.parseLong(args[0]); opts.streamLenProvided = true; }
            if (args.length >= 2) { A = Integer.parseInt(args[1]); opts.alphabetProvided = true; }
            if (args.length >= 3) { B = Integer.parseInt(args[2]); }
            if (args.length >= 4) { dist = Dist.valueOf(args[3].toUpperCase(Locale.ROOT)); }
            if (args.length >= 5) { s = Double.parseDouble(args[4]); }
            if (args.length >= 6) { p = Double.parseDouble(args[5]); opts.quantilesProvided = true; }
            if (args.length >= 7) { seed = Long.parseLong(args[6]); }
            if (args.length >= 8) { delta = Double.parseDouble(args[7]); opts.summaryDelta = delta; opts.summaryDeltaProvided = true; }
            if (args.length >= 9) { trials = Integer.parseInt(args[8]); }
            if (args.length >= 10) { verbose = Boolean.parseBoolean(args[9]); }
            if (args.length >= 11) { autoDesignB = Boolean.parseBoolean(args[10]); }
            if (args.length >= 12) { distinctGuess = Integer.parseInt(args[11]); opts.distinctProvided = true; }
            if (args.length >= 13) { rankEpsTarget = Double.parseDouble(args[12]); opts.rankEpsProvided = true; }
            if (args.length >= 14) { deltaSamp = Double.parseDouble(args[13]); opts.deltaSamps = List.of(deltaSamp); opts.deltaSampsProvided = true; }
            if (args.length >= 15) { csvPath = args[14]; }
            if (args.length >= 16) { csvAppend = Boolean.parseBoolean(args[15]); }
            if (args.length >= 17) { opts.datasetNgram = parsePositiveInt("--dataset-ngram", args[16]); }

            opts.streamLen = N;
            opts.alphabet = A;
            opts.buckets = B;
            opts.dist = dist;
            opts.zipfExponent = s;
            opts.quantiles = List.of(p);
            opts.quantilesProvided = true;
            opts.seed = seed;
            opts.trials = trials;
            opts.verbose = verbose;
            opts.autoDesign = autoDesignB;
            opts.distinctGuess = distinctGuess;
            opts.rankEpsTargets = List.of(rankEpsTarget);
            opts.deltaQs = List.of(opts.summaryDeltaProvided ? opts.summaryDelta : delta);
            opts.deltaQsProvided = true;
            if (!opts.deltaSampsProvided) {
                opts.deltaSamps = List.of(deltaSamp);
            }
            if (opts.deltaQs.isEmpty()) {
                opts.deltaQs = DEFAULT_DELTA_QS;
            }
            if (opts.deltaSamps.isEmpty()) {
                opts.deltaSamps = DEFAULT_DELTA_SAMPS;
            }
            opts.csvPath = csvPath != null && csvPath.equalsIgnoreCase("none") ? null : Paths.get(csvPath);
            opts.csvAppend = csvAppend;

            SUMMARY_DELTA = opts.summaryDeltaProvided ? opts.summaryDelta : delta;

            return opts;
        }

        private static DatasetMode parseDatasetMode(String value) {
            String norm = value.trim().toUpperCase(Locale.ROOT);
            return switch (norm) {
                case "SEGMENTS" -> DatasetMode.SEGMENTS;
                case "WORDS" -> DatasetMode.WORDS;
                case "CHARS", "CHARACTERS" -> DatasetMode.CHARS;
                default -> throw new IllegalArgumentException("Unknown dataset mode '" + value + "' (expected segments, words, chars)");
            };
        }

        private static List<Double> parseDoubleList(String csv) {
            if (csv == null || csv.isBlank()) {
                return List.of();
            }
            String[] parts = csv.split(",");
            List<Double> out = new ArrayList<>(parts.length);
            for (String part : parts) {
                if (part.isBlank()) continue;
                out.add(Double.parseDouble(part.trim()));
            }
            return out.isEmpty() ? List.of() : List.copyOf(out);
        }

        private static String requireValue(String flag, String[] args, int index) {
            if (index >= args.length) {
                throw new IllegalArgumentException("Missing value for " + flag);
            }
            return args[index];
        }

        private static int parsePositiveInt(String flag, String value) {
            try {
                int parsed = Integer.parseInt(value);
                if (parsed <= 0) {
                    throw new IllegalArgumentException(flag + " must be > 0 (was " + value + ")");
                }
                return parsed;
            } catch (NumberFormatException ex) {
                throw new IllegalArgumentException(flag + " must be an integer (was " + value + ")", ex);
            }
        }

        private static void printUsage() {
            System.out.println("SamplingSketchExperiment usage:\n" +
                    "  --stream-len <long>      Number of stream items (default 10M)\n" +
                    "  --alphabet <int>         Alphabet size for synthetic streams (default 10k)\n" +
                    "  --buckets <int>          HOPS/BK buckets (default 2500)\n" +
                    "  --dist <uniform|zipf>    Synthetic distribution (default zipf)\n" +
                    "  --zipf <double>          Zipf exponent if dist=zipf (default 1.5)\n" +
                    "  --quantile <list>        Target quantile(s) p (comma-separated, default 0.05)\n" +
                    "  --seed <long>            Base RNG seed (default 123456789)\n" +
                    "  --trials <int>           Number of trials (default 50)\n" +
                    "  --verbose <bool>         Per-trial logging (default false)\n" +
                    "  --auto-design <bool>     Enable Chebyshev auto-design (default true)\n" +
                    "  --distinct <int>         Override distinct guess D̂ (default depends on stream)\n" +
                    "  --rank-eps <list>        Target rank epsilon(s) (default 0.05)\n" +
                    "  --delta <double>         Convenience shorthand when using a single delta_q\n" +
                    "  --delta-q <list>         Comma-separated delta_q values (default 0.05,0.075,0.01)\n" +
                    "  --delta-samp <list>      Comma-separated delta_samp values (default 0.05,0.025,0.09)\n" +
                    "  --csv <path>             Output CSV path (default sampling_results.csv)\n" +
                    "  --csv-append <bool>      Append to CSV (default true)\n" +
                    "  --dataset <path>         Use dataset file instead of synthetic stream\n" +
                    "  --dataset-mode <mode>    segments: each line, words: whitespace tokens, chars: per character\n" +
                    "  --dataset-ngram <int>    N-gram width for dataset tokens (default 1)\n" +
                    "  --help                   Show this message\n");
        }


        /** Write rows to CSV with append semantics; auto-writes header if file is new/empty. */
        private static void writeSuiteACsv(Path out, List<TrialResult> rs) throws IOException {
            if (rs == null || rs.isEmpty()) return;

            boolean needHeader = true;
            if (Files.exists(out) && Files.size(out) > 0) {
                needHeader = false;
            }

            List<List<?>> rows = new ArrayList<>(rs.size() + 1);
            if (needHeader) rows.add(suiteACsvHeader());
            for (TrialResult r : rs) rows.add(suiteACsvRow(r, EXPONENT));

            try (OutputStream os = Files.newOutputStream(out, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
                CsvUtil.writeRows(os, rows, CsvUtil.Config.defaults());
            }
            System.out.printf("CSV: wrote %d trial row(s) to %s (append)\n", rs.size(), out.toAbsolutePath());
        }

        private static DatasetStream loadDataset(Options opts) throws IOException {
            Objects.requireNonNull(opts.datasetPath, "datasetPath");

            ArrayList<Long> tokens = new ArrayList<>();
            LongOpenHashSet distinct = new LongOpenHashSet();
            TokenEncoder encoder = new TokenEncoder();

            int ngram = Math.max(1, opts.datasetNgram);
            long limit = opts.streamLenProvided ? opts.streamLen : Long.MAX_VALUE;

            ArrayDeque<String> window = new ArrayDeque<>(ngram);

            try (BufferedReader reader = Files.newBufferedReader(opts.datasetPath, StandardCharsets.UTF_8)) {
                String line;
                outer: while ((line = reader.readLine()) != null) {
                    if (opts.datasetMode == DatasetMode.SEGMENTS) {
                        String token = line.trim();
                        if (token.isEmpty()) continue;
                        if (emitNgram(token, window, ngram, encoder, tokens, distinct, limit)) break;
                    } else if (opts.datasetMode == DatasetMode.WORDS) {
                        String trimmed = line.trim();
                        if (trimmed.isEmpty()) continue;
                        String[] parts = trimmed.split("\\s+");
                        for (String part : parts) {
                            if (emitNgram(part, window, ngram, encoder, tokens, distinct, limit)) {
                                break outer;
                            }
                        }
                    } else { // CHARS
                        for (int idx = 0, len = line.length(); idx < len; idx++) {
                            String ch = String.valueOf(line.charAt(idx));
                            if (emitNgram(ch, window, ngram, encoder, tokens, distinct, limit)) {
                                break outer;
                            }
                        }
                    }
                }
            }

            if (tokens.isEmpty()) {
                throw new IllegalStateException("Dataset produced no tokens (after applying n-gram size " + ngram + ")");
            }

            return new DatasetStream(tokens, distinct.size());
        }

        private static boolean emitNgram(String rawToken,
                                         ArrayDeque<String> window,
                                         int ngram,
                                         TokenEncoder encoder,
                                         ArrayList<Long> outTokens,
                                         LongOpenHashSet distinct,
                                         long limit) {
            if (ngram <= 1) {
                long value = encoder.encode(rawToken);
                outTokens.add(value);
                distinct.add(value);
                return outTokens.size() >= limit;
            }

            window.addLast(rawToken);
            if (window.size() < ngram) {
                return false;
            }

            String ngramToken = buildNgram(window);
            long value = encoder.encode(ngramToken);
            outTokens.add(value);
            distinct.add(value);
            window.removeFirst();
            return outTokens.size() >= limit;
        }

        private static String buildNgram(ArrayDeque<String> window) {
            StringBuilder sb = new StringBuilder();
            boolean first = true;
            for (String token : window) {
                if (!first) {
                    sb.append('\u0001');
                }
                sb.append(token);
                first = false;
            }
            return sb.toString();
        }

        /**
         * Empirical CDF F_m(x) from an ascending array 'a': proportion of entries <= x.
         */
        private static double empiricalCdfLeft(int[] a, int x) {
            int idx = upperBound(a, x) - 1; // last index with a[i] <= x
            if (idx < 0) return 0.0;
            return (idx + 1) / (double) a.length;
        }

        /**
         * Upper bound: first index with a[i] > x in ascending 'a'.
         */
        private static int upperBound(int[] a, int x) {
            int lo = 0, hi = a.length;
            while (lo < hi) {
                int mid = (lo + hi) >>> 1;
                if (a[mid] <= x) lo = mid + 1;
                else hi = mid;
            }
            return lo;
        }

        /**
         * Holds per-item timing results in nanoseconds.
         */
        private static final class PerItemCost {
            final double nsPerInsertMPQ;
            final double nsPerInsertBK;
            final double closeNsPerItemMPQ;
            final double closeNsPerItemBK;

            PerItemCost(double a, double b, double c, double d) {
                this.nsPerInsertMPQ = a;
                this.nsPerInsertBK = b;
                this.closeNsPerItemMPQ = c;
                this.closeNsPerItemBK = d;
            }
        }

        /**
         * Measure per-item costs for MPQ (HOPS) and Bottom-k on the given stream.
         * We time only the insert loops for update cost, and we time the close-time
         * steps separately. We include a short warm-up to let the JIT compile.
         */
        private static PerItemCost measurePerItemCosts(
                ArrayList<Long> stream,
                int buckets,
                long seed
        ) {
            final int n = stream.size();

            // Seeds consistent with your runExperiment
            long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
            long prioSeed = mix64(seed ^ 0xf00dbabe42aa1357L);

            // -------- MPQ / HOPS update cost --------
            // Warm-up
            {
                HOPS warm = new HOPS(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);
                for (long k : stream) warm.insert(k);
            }
            // Timed
            HOPS mpq = new HOPS(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);
            long t0 = System.nanoTime();
            for (long k : stream) mpq.insert(k);
            long t1 = System.nanoTime();
            double nsPerInsertMPQ = (t1 - t0) / (double) n;

            // -------- Bottom-k update cost --------
            // Warm-up
            {
                BottomKSampler warm = new BottomKSampler(/*k=*/buckets, /*seed=*/seed ^ 0x5eed5eedL);
                for (long k : stream) warm.offer(k);
            }
            // Timed
            BottomKSampler bk = new BottomKSampler(/*k=*/buckets, /*seed=*/seed ^ 0x5eed5eedL);
            long t2 = System.nanoTime();
            for (long k : stream) bk.offer(k);
            long t3 = System.nanoTime();
            double nsPerInsertBK = (t3 - t2) / (double) n;

            // -------- Close-time cost (sketch-only and quantile computation) --------
            // MPQ close: get representatives and sort their frequency proxies if available later.
            long t4 = System.nanoTime();
            long[] repsMPQ = mpq.getRepresentatives();
            // Sorting representatives by their counts belongs to the *estimation* step.
            // If you later plug in your frequency sketch, include that lookup here.
            Arrays.sort(repsMPQ); // placeholder to avoid dead-code elimination
            long t5 = System.nanoTime();
            double closeNsPerItemMPQ = (t5 - t4) / (double) n;

            // Bottom-k close: extract and sort keys
            long t6 = System.nanoTime();
            long[] repsBK = bk.sampleKeys();
            Arrays.sort(repsBK);
            long t7 = System.nanoTime();
            double closeNsPerItemBK = (t7 - t6) / (double) n;

            return new PerItemCost(nsPerInsertMPQ, nsPerInsertBK, closeNsPerItemMPQ, closeNsPerItemBK);
        }

        /**
         * DKW rank epsilon for confidence 1-delta over n samples.
         */
        private static double dkwRankEpsilon(int n, double delta) {
            if (n <= 0) return 1.0;
            if (delta <= 0.0 || delta >= 1.0) throw new IllegalArgumentException("delta must be in (0,1)");
            return Math.sqrt(Math.log(2.0 / delta) / (2.0 * n));
        }

        /**
         * Sorted array (ascending) of exact per-key counts.
         */
        private static int[] sortedCounts(HashMap<Long, Integer> freq) {
            int[] a = new int[freq.size()];
            int i = 0;
            for (int v : freq.values()) a[i++] = v;
            Arrays.sort(a);
            return a;
        }

        /**
         * Value at quantile p from an ascending array 'a' using the left-continuous definition:
         * return the smallest x such that empirical CDF ≥ p.
         * p=0 returns a[0]; p in (0,1] maps to index ceil(p*m)-1 with m=a.length.
         */
        private static int valueAtQuantile(int[] a, double p) {
            if (a.length == 0) return 0;
            if (p <= 0.0) return a[0];
            if (p >= 1.0) return a[a.length - 1];
            int rank = (int) Math.ceil(p * a.length) - 1;
            if (rank < 0) rank = 0;
            if (rank >= a.length) rank = a.length - 1;
            return a[rank];
        }

        // =========================
        //  Occupancy math (exact)
        // =========================

        // =========================
        //  One trial (truth + samplers)
        // =========================

        /**
         * Run one MPQ-vs-Bottom-k vs Truth quantile experiment and return detailed results.
         * Reports both value errors and rank errors, plus guarantee checks.
         */
        /**
         * Run one MPQ-vs-Bottom-k vs Truth quantile experiment and return detailed results.
         * Reports both value errors and rank errors, plus guarantee checks.
         */
        private static TrialResult runExperiment(
                List<Long> prebuiltStream,
                long streamLen,
                int alphabetSize,
                int buckets,
                Dist dist,
                double zipfS,
                double p,
                long seed,
                double delta,
                boolean verbose,
                boolean autoDesignB,
                Integer distinctGuess,
                Double rankEpsTarget,
                Double deltaSamp,
                String datasetLabel,
                DatasetMode datasetMode,
                String distributionLabel,
                int datasetNgram
        ) {
            if (!(p > 0.0 && p <= 1.0)) {
                throw new IllegalArgumentException("p must be in (0,1]");
            }

            final ArrayList<Long> stream;
            if (prebuiltStream != null) {
                if (prebuiltStream instanceof ArrayList<?> arr && prebuiltStream.getClass() == ArrayList.class) {
                    @SuppressWarnings("unchecked")
                    ArrayList<Long> cast = (ArrayList<Long>) arr;
                    stream = cast;
                } else {
                    stream = new ArrayList<>(prebuiltStream);
                }
            } else {
                Random rng = new Random(seed);
                stream = new ArrayList<>((int) Math.min(streamLen, Integer.MAX_VALUE));
                ZipfDistribution z = (dist == Dist.ZIPF) ? makeZipfIntDist(alphabetSize, zipfS, seed) : null;
                for (long t = 0; t < streamLen; t++) {
                    long key = (dist == Dist.UNIFORM) ? rng.nextInt(alphabetSize) : sampleZipfInt(Objects.requireNonNull(z));
                    stream.add(key);
                }
            }

            streamLen = stream.size();

            // Exact truth
            long tTruth0 = System.currentTimeMillis();
            HashMap<Long, Integer> freq = new HashMap<>(Math.max(16, alphabetSize * 2));
            for (long key : stream) freq.merge(key, 1, Integer::sum);
            int[] truthSorted = sortedCounts(freq);
            long tTruthMs = System.currentTimeMillis() - tTruth0;

            int Sactive = freq.size();
            int xTrue = valueAtQuantile(truthSorted, p);

            // Auto-design B if requested (Chebyshev)
            int Bused = buckets;
            Utils.HopsDesignResult design = null;
            double epsT = (rankEpsTarget != null) ? rankEpsTarget : 0.01;
            double deltaS = (deltaSamp != null) ? deltaSamp : 0.05;
            int DhatUsed = Math.max(1, (distinctGuess != null && distinctGuess > 0) ? distinctGuess : Sactive);

            if (autoDesignB) {
                design = Utils.designBucketsForRankTargetChebyshev(DhatUsed, epsT, /*delta_q=*/delta, /*delta_samp=*/deltaS);
                Bused = design.suggestedBuckets;

                // Best-achievable rank error with realized sample size at design LB:
                int nAch = Math.max(1, Math.min(DhatUsed, design.occupancyLowerBound)); // cannot exceed Dhat
                double epsAch = Math.sqrt(Math.log(2.0 / delta) / (2.0 * nAch));
                if (verbose) {
                    System.out.printf("Auto-design Chebyshev: Dhat=%d, epsTarget=%.6f, delta_q=%.3f, delta_samp=%.3f%n",
                            DhatUsed, epsT, delta, deltaS);
                    System.out.printf("Suggested B=%d, n_req=%d, E[n_b]=%.2f, Var=%.2f, LB(n_b)=%d%s%n",
                            design.suggestedBuckets, design.requiredSampleSize, design.expectedNonEmpty, design.variance,
                            design.occupancyLowerBound,
                            design.impossible ? "  (WARNING: target impossible with Dhat)" : "");
                }

                if (design.impossible) {
                    System.out.printf("NOTE: With Dhat=%d, the tightest DKW rank band at delta_q=%.3f is eps≈%.5f%n",
                            DhatUsed, delta, epsAch);
                }
            }

            // Measure per-item costs for the chosen B (uses same stream and seeds)
            PerItemCost cost = measurePerItemCosts(stream, Bused, seed);

            // --- MPQ sampler (HOPS) ---
            long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
            long prioSeed = mix64(seed ^ 0xf00dbabe42aa1357L);
            HOPS mpq = new HOPS(Bused, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

            long tMPQ0 = System.currentTimeMillis();
            for (long key : stream) mpq.insert(key);
            int xMPQ = mpq.estimateQuantile(k -> freq.getOrDefault(k, 0), p);
            long tMPQMs = System.currentTimeMillis() - tMPQ0;
            int nbMPQ = mpq.sampleSize();
            double epsMPQ = dkwRankEpsilon(nbMPQ, delta);
            double pLoMPQ = Math.max(0.0, p - epsMPQ);
            double pHiMPQ = Math.min(1.0, p + epsMPQ);
            int xLoMPQ = valueAtQuantile(truthSorted, pLoMPQ);
            int xHiMPQ = valueAtQuantile(truthSorted, pHiMPQ);
            boolean mpqInDKWValueBand = (xMPQ >= xLoMPQ && xMPQ <= xHiMPQ);

            double pAchMPQ = empiricalCdfLeft(truthSorted, xMPQ);
            double rankErrMPQ = Math.abs(pAchMPQ - p);

            // --- Bottom-k sampler (k=Bused) ---
            BottomKSampler bk = new BottomKSampler(/*k=*/Bused, /*seed=*/seed ^ 0x5eed5eedL);

            long tBK0 = System.currentTimeMillis();
            for (long key : stream) bk.offer(key);
            long[] repsBK = bk.sampleKeys();
            int nbBK = repsBK.length;
            int[] sampleCountsBK = new int[nbBK];
            for (int i = 0; i < nbBK; i++) sampleCountsBK[i] = freq.getOrDefault(repsBK[i], 0);
            Arrays.sort(sampleCountsBK);
            long tBKMs = System.currentTimeMillis() - tBK0;

            int xBK = valueAtQuantile(sampleCountsBK, p);
            double epsBK = dkwRankEpsilon(nbBK, delta);
            double pLoBK = Math.max(0.0, p - epsBK);
            double pHiBK = Math.min(1.0, p + epsBK);
            int xLoBK = valueAtQuantile(truthSorted, pLoBK);
            int xHiBK = valueAtQuantile(truthSorted, pHiBK);
            boolean bkInDKWValueBand = (xBK >= xLoBK && xBK <= xHiBK);

            double pAchBK = empiricalCdfLeft(truthSorted, xBK);
            double rankErrBK = Math.abs(pAchBK - p);

            // Design checks
            int nLB = -1;
            boolean occLBMet = true;
            if (design != null) {
                nLB = design.occupancyLowerBound;                // Chebyshev LB at suggested B
                occLBMet = nbMPQ >= nLB;
            }

            if (verbose) {
                System.out.printf("Active distinct S=%d, Bused=%d%n", Sactive, Bused);
                System.out.printf("Truth time: %d ms   MPQ time: %d ms   Bottom-k time: %d ms%n",
                        tTruthMs, tMPQMs, tBKMs);
                System.out.printf("Per-item cost (ns): MPQ insert=%.1f, MPQ close=%.3f, BK insert=%.1f, BK close=%.3f%n",
                        cost.nsPerInsertMPQ, cost.closeNsPerItemMPQ, cost.nsPerInsertBK, cost.closeNsPerItemBK);

                System.out.printf("True  x_p = %d%n", xTrue);

                System.out.printf("MPQ   nb=%d, x̂_p=%d, value error=%.3f%%, achieved percentile=%.5f, rank error=%.5f%n",
                        nbMPQ, xMPQ, 100.0 * percentValueError(xMPQ, xTrue), pAchMPQ, rankErrMPQ);
                System.out.printf("MPQ   DKW eps=%.5f ⇒ rank band [%.5f, %.5f], value band [%d, %d] ⇒ %s%n",
                        epsMPQ, pLoMPQ, pHiMPQ, xLoMPQ, xHiMPQ, mpqInDKWValueBand ? "INSIDE" : "OUTSIDE");
                if (design != null) {
                    System.out.printf("MPQ   occupancy LB check: observed nb=%d, LB(n_b)=%d ⇒ %s%n",
                            nbMPQ, nLB, occLBMet ? "OK" : "LOW");
                }

                System.out.printf("BK    nb=%d, x̂_p=%d, value error=%.3f%%, achieved percentile=%.5f, rank error=%.5f%n",
                        nbBK, xBK, 100.0 * percentValueError(xBK, xTrue), pAchBK, rankErrBK);
                System.out.printf("BK    DKW eps=%.5f ⇒ rank band [%.5f, %.5f], value band [%d, %d] ⇒ %s%n",
                        epsBK, pLoBK, pHiBK, xLoBK, xHiBK, bkInDKWValueBand ? "INSIDE" : "OUTSIDE");
            }

            return new TrialResult(
                    datasetLabel,
                    datasetMode != null ? datasetMode.name().toLowerCase(Locale.ROOT) : "n/a",
                    streamLen,
                    alphabetSize,
                    buckets,
                    distributionLabel,
                    zipfS,
                    datasetNgram,
                    // sizes and occupancy
                    Sactive, Bused, nbMPQ, nbBK,
                    // truth and estimates
                    xTrue, xMPQ, xBK,
                    // value errors (percentage)
                    100.0 * percentValueError(xMPQ, xTrue),
                    100.0 * percentValueError(xBK, xTrue),
                    // rank metrics
                    p, pAchMPQ, pAchBK, rankErrMPQ, rankErrBK,
                    // DKW parameters and value bands and checks
                    delta, dkwRankEpsilon(nbMPQ, delta), dkwRankEpsilon(nbBK, delta),
                    pLoMPQ, pHiMPQ, xLoMPQ, xHiMPQ, mpqInDKWValueBand,
                    pLoBK, pHiBK, xLoBK, xHiBK, bkInDKWValueBand,
                    // occupancy LB check
                    nLB, occLBMet,
                    // timings
                    tTruthMs, tMPQMs, tBKMs,
                    // per-item costs (ns)
                    cost.nsPerInsertMPQ, cost.closeNsPerItemMPQ, cost.nsPerInsertBK, cost.closeNsPerItemBK,
                    // deltas and design knob actually used
                    epsT, deltaS, DhatUsed,
                    // design info
                    design
            );
        }

        // =========================
        //  Trial result and summary
        // =========================

        /**
         * Immutable result for one trial with both value/rank diagnostics, guarantees, costs, and design knobs.
         */
        private static final class TrialResult {
            // Metadata
            final String datasetLabel;
            final String datasetMode;
            final long streamLength;
            final int alphabetConfigured;
            final int bucketsConfigured;
            final String distributionLabel;
            final double zipfExponent;
            final int datasetNgram;

            // Sizes
            final int Sactive;
            final int Bused;
            final int nbMPQ;
            final int nbBK;

            // Truth and estimates
            final int xTrue;
            final int xMPQ;
            final int xBK;

            // Value errors in percent
            final double pctErrMPQ;
            final double pctErrBK;

            // Rank targets and achieved percentiles
            final double pTarget;
            final double pAchMPQ;
            final double pAchBK;
            final double rankErrMPQ;
            final double rankErrBK;

            // DKW configuration and bands
            final double delta;   // delta_q
            final double epsMPQ;
            final double epsBK;

            final double pLoMPQ, pHiMPQ;
            final int xLoMPQ, xHiMPQ;
            final boolean mpqInDKWValueBand;

            final double pLoBK, pHiBK;
            final int xLoBK, xHiBK;
            final boolean bkInDKWValueBand;

            // Occupancy lower bound check
            final int nLB;            // Chebyshev LB at Bused (or -1 if design disabled)
            final boolean occLBMet;

            // Timings in milliseconds
            final long tTruthMs, tMPQMs, tBKMs;

            // Per-item costs in nanoseconds
            final double nsPerInsertMPQ;
            final double closeNsPerItemMPQ;
            final double nsPerInsertBK;
            final double closeNsPerItemBK;

            // Design knobs actually used (for CSV)
            final double epsTargetUsed;
            final double deltaSampUsed;
            final int DhatUsed;

            // Optional design info
            final Utils.HopsDesignResult design;

            TrialResult(
                    String datasetLabel,
                    String datasetMode,
                    long streamLength,
                    int alphabetConfigured,
                    int bucketsConfigured,
                    String distributionLabel,
                    double zipfExponent,
                    int datasetNgram,
                    int Sactive, int Bused, int nbMPQ, int nbBK,
                    int xTrue, int xMPQ, int xBK,
                    double pctErrMPQ, double pctErrBK,
                    double pTarget, double pAchMPQ, double pAchBK, double rankErrMPQ, double rankErrBK,
                    double delta, double epsMPQ, double epsBK,
                    double pLoMPQ, double pHiMPQ, int xLoMPQ, int xHiMPQ, boolean mpqInDKWValueBand,
                    double pLoBK, double pHiBK, int xLoBK, int xHiBK, boolean bkInDKWValueBand,
                    int nLB, boolean occLBMet,
                    long tTruthMs, long tMPQMs, long tBKMs,
                    double nsPerInsertMPQ, double closeNsPerItemMPQ, double nsPerInsertBK, double closeNsPerItemBK,
                    double epsTargetUsed, double deltaSampUsed, int DhatUsed,
                    Utils.HopsDesignResult design
            ) {
                this.datasetLabel = datasetLabel;
                this.datasetMode = datasetMode;
                this.streamLength = streamLength;
                this.alphabetConfigured = alphabetConfigured;
                this.bucketsConfigured = bucketsConfigured;
                this.distributionLabel = distributionLabel;
                this.zipfExponent = zipfExponent;
                this.datasetNgram = datasetNgram;

                this.Sactive = Sactive;
                this.Bused = Bused;
                this.nbMPQ = nbMPQ;
                this.nbBK = nbBK;

                this.xTrue = xTrue;
                this.xMPQ = xMPQ;
                this.xBK = xBK;

                this.pctErrMPQ = pctErrMPQ;
                this.pctErrBK = pctErrBK;

                this.pTarget = pTarget;
                this.pAchMPQ = pAchMPQ;
                this.pAchBK = pAchBK;
                this.rankErrMPQ = rankErrMPQ;
                this.rankErrBK = rankErrBK;

                this.delta = delta;
                this.epsMPQ = epsMPQ;
                this.epsBK = epsBK;

                this.pLoMPQ = pLoMPQ;
                this.pHiMPQ = pHiMPQ;
                this.xLoMPQ = xLoMPQ;
                this.xHiMPQ = xHiMPQ;
                this.mpqInDKWValueBand = mpqInDKWValueBand;

                this.pLoBK = pLoBK;
                this.pHiBK = pHiBK;
                this.xLoBK = xLoBK;
                this.xHiBK = xHiBK;
                this.bkInDKWValueBand = bkInDKWValueBand;

                this.nLB = nLB;
                this.occLBMet = occLBMet;

                this.tTruthMs = tTruthMs;
                this.tMPQMs = tMPQMs;
                this.tBKMs = tBKMs;

                this.nsPerInsertMPQ = nsPerInsertMPQ;
                this.closeNsPerItemMPQ = closeNsPerItemMPQ;
                this.nsPerInsertBK = nsPerInsertBK;
                this.closeNsPerItemBK = closeNsPerItemBK;

                this.epsTargetUsed = epsTargetUsed;
                this.deltaSampUsed = deltaSampUsed;
                this.DhatUsed = DhatUsed;

                this.design = design;
            }
        }


        /** Print averages and guarantee satisfaction rates over many trials. */
        /**
         * Print averages and guarantee satisfaction rates over many trials, including per-item costs.
         */
        /**
         * Print averages and guarantee satisfaction rates over many trials, including per-item costs.
         * NOW ALSO prints the average achieved epsilon (DKW band half-width) for MPQ and BK.
         */
        /**
         * Print averages and guarantee satisfaction rates over many trials, including per-item costs.
         * Now also prints:
         *  - avg achieved percentile p̂ (MPQ and BK)
         *  - avg rank error |p̂ - p| (MPQ and BK) with clearer labels
         *  - avg achieved epsilon from DKW (epsMPQ, epsBK)
         */
        private static void summarize(List<TrialResult> rs, double p, double delta) {
            int T = rs.size();
            if (T == 0) return;

            double avgS = 0, avgB = 0, avgNbMPQ = 0, avgNbBK = 0;

            // value error (in count space, %)
            double avgValErrMPQ = 0, avgValErrBK = 0;

            // rank stuff
            double avgRankErrMPQ = 0, avgRankErrBK = 0;
            double avgPAchMPQ    = 0, avgPAchBK    = 0; // achieved percentiles p̂
            double avgEpsMPQ     = 0, avgEpsBK     = 0; // DKW half-widths ε

            double insideMPQ = 0, insideBK = 0;
            double occOK = 0, occTotal = 0;
            double unionOK = 0, unionTotal = 0;

            double avgTruthMs = 0, avgMPQMs = 0, avgBKMs = 0;

            // per-item costs
            double avgNsInsMPQ = 0, avgNsCloseMPQ = 0, avgNsInsBK = 0, avgNsCloseBK = 0;

            for (TrialResult r : rs) {
                avgS += r.Sactive;
                avgB += r.Bused;
                avgNbMPQ += r.nbMPQ;
                avgNbBK += r.nbBK;

                avgValErrMPQ += r.pctErrMPQ;
                avgValErrBK  += r.pctErrBK;

                avgRankErrMPQ += r.rankErrMPQ;
                avgRankErrBK  += r.rankErrBK;

                avgPAchMPQ += r.pAchMPQ;
                avgPAchBK  += r.pAchBK;

                avgEpsMPQ  += r.epsMPQ;
                avgEpsBK   += r.epsBK;

                insideMPQ += r.mpqInDKWValueBand ? 1 : 0;
                insideBK  += r.bkInDKWValueBand ? 1 : 0;

                if (r.nLB >= 0) {
                    occTotal += 1;
                    occOK += r.occLBMet ? 1 : 0;

                    unionTotal += 1;
                    boolean unionSuccess = r.occLBMet && r.mpqInDKWValueBand;
                    unionOK += unionSuccess ? 1 : 0;
                }

                avgTruthMs += r.tTruthMs;
                avgMPQMs   += r.tMPQMs;
                avgBKMs    += r.tBKMs;

                avgNsInsMPQ    += r.nsPerInsertMPQ;
                avgNsCloseMPQ  += r.closeNsPerItemMPQ;
                avgNsInsBK     += r.nsPerInsertBK;
                avgNsCloseBK   += r.closeNsPerItemBK;
            }

            // divide by T to get averages
            avgS /= T;
            avgB /= T;
            avgNbMPQ /= T;
            avgNbBK /= T;

            avgValErrMPQ /= T;
            avgValErrBK  /= T;

            avgRankErrMPQ /= T;
            avgRankErrBK  /= T;

            avgPAchMPQ /= T;
            avgPAchBK  /= T;

            avgEpsMPQ /= T;
            avgEpsBK  /= T;

            double rateMPQ = 100.0 * insideMPQ / T;
            double rateBK  = 100.0 * insideBK  / T;

            double rateOcc   = (occTotal   > 0) ? 100.0 * occOK   / occTotal   : Double.NaN;
            double rateUnion = (unionTotal > 0) ? 100.0 * unionOK / unionTotal : Double.NaN;

            avgTruthMs /= T;
            avgMPQMs   /= T;
            avgBKMs    /= T;

            avgNsInsMPQ    /= T;
            avgNsCloseMPQ  /= T;
            avgNsInsBK     /= T;
            avgNsCloseBK   /= T;

            System.out.println("\n=== Aggregate over " + T + " trials ===");
            String deltaText = Double.isNaN(delta) ? "varies" : String.format(Locale.ROOT, "%.3f", delta);
            System.out.printf("Target percentile p = %.5f (delta_q = %s)\n", p, deltaText);

            System.out.printf("Avg Sactive = %.1f, Avg Bused = %.1f\n", avgS, avgB);
            System.out.printf("Avg sample sizes: MPQ nb=%.1f, BK nb=%.1f\n", avgNbMPQ, avgNbBK);

            // Value error (count space)
            System.out.printf("Avg VALUE error vs truth x_p (%%):  MPQ=%.3f%%, BK=%.3f%%\n",
                    avgValErrMPQ, avgValErrBK);

            // Achieved percentile and rank error
            System.out.printf("Avg achieved percentile p̂:        MPQ=%.5f, BK=%.5f\n",
                    avgPAchMPQ, avgPAchBK);
            System.out.printf("Avg RANK error |p̂−p|:             MPQ=%.5f, BK=%.5f\n",
                    avgRankErrMPQ, avgRankErrBK);

            // Certified epsilon from DKW given nb per trial
            System.out.printf("Avg certified DKW ε (half-width):  MPQ=%.5f, BK=%.5f\n",
                    avgEpsMPQ, avgEpsBK);

            // Band success
            System.out.printf("Inside DKW→value band:             MPQ=%.1f%% of trials, BK=%.1f%% of trials\n",
                    rateMPQ, rateBK);

            if (!Double.isNaN(rateOcc)) {
                System.out.printf("Occupancy LB satisfied:            %.1f%% of trials\n", rateOcc);
            }
            if (!Double.isNaN(rateUnion)) {
                System.out.printf("Union success (occ LB && band):    %.1f%% of trials\n", rateUnion);
            }

            System.out.printf("Avg times (ms): Truth=%.1f, MPQ=%.1f, Bottom-k=%.1f\n",
                    avgTruthMs, avgMPQMs, avgBKMs);

            System.out.printf("Avg per-item cost (ns): MPQ insert=%.1f, MPQ close=%.3f, BK insert=%.1f, BK close=%.3f\n",
                    avgNsInsMPQ, avgNsCloseMPQ, avgNsInsBK, avgNsCloseBK);
        }


        /**
         * Print Suite-A CSV header (exactly the columns you asked for).
         */
        private static void printSuiteACsvHeader() {
            System.out.println(
                    "dataset,dataset_mode,stream_len,alphabet,buckets_configured,source_dist,zipf_exponent,quantile_target,dataset_ngram," +
                            "eps_target,delta_q,delta_samp,Dhat,B_suggested,n_req,E_nb,Var_nb,LB_nb," +
                            "Sactive,Bused,nb_MPQ,eps_DKW_MPQ,p_lo_MPQ,p_hi_MPQ,x_lo_MPQ,x_hi_MPQ," +
                            "MPQ_xhat,MPQ_value_error_pct,MPQ_achieved_percentile,MPQ_rank_error," +
                            "inside_band_MPQ,occupancy_ok,ns_per_insert_MPQ,ns_close_MPQ," +
                            "nb_BK,eps_DKW_BK,x_lo_BK,x_hi_BK,BK_xhat,BK_value_error_pct," +
                            "BK_achieved_percentile,BK_rank_error,ns_per_insert_BK,ns_close_BK," +
                            "truth_x_p,mpq_time_ms,bk_time_ms"
            );
        }

        /**
         * Print one Suite-A CSV row for a single trial result.
         */
        private static void printSuiteACsvRow(TrialResult r) {
            String inside = r.mpqInDKWValueBand ? "TRUE" : "FALSE";
            String occOK = r.occLBMet ? "TRUE" : "FALSE";

            // Design fields (present when auto-design is enabled)
            int B_suggested = (r.design != null) ? r.design.suggestedBuckets : r.Bused;
            int n_req = (r.design != null) ? r.design.requiredSampleSize : -1;
            double E_nb = (r.design != null) ? r.design.expectedNonEmpty : Double.NaN;
            double Var_nb = (r.design != null) ? r.design.variance : Double.NaN;
            int LB_nb = (r.design != null) ? r.design.occupancyLowerBound : -1;

            System.out.printf(
                    Locale.US,
                    "%s,%s,%d,%d,%d,%s,%.4f,%.5f,%d," +         // dataset metadata
                            "%.3f,%.2f,%.3f,%d,%d,%d,%.2f,%.2f,%d," +   // eps_target, delta_q, delta_samp, Dhat, B_suggested, n_req, E_nb, Var_nb, LB_nb
                            "%d,%d,%d,%.5f,%.5f,%.5f,%d,%d," +          // Sactive, Bused, nb_MPQ, eps_DKW_MPQ, p_lo_MPQ, p_hi_MPQ, x_lo_MPQ, x_hi_MPQ
                            "%d,%.3f,%.5f,%.5f,%s,%s,%.3f,%.3f," +      // MPQ_xhat, MPQ_value_error_pct, MPQ_achieved_percentile, MPQ_rank_error, inside_band_MPQ, occupancy_ok, ns_per_insert_MPQ, ns_close_MPQ
                            "%d,%.5f,%d,%d,%d,%.3f,%.5f,%.5f,%.3f,%.3f," + // nb_BK, eps_DKW_BK, x_lo_BK, x_hi_BK, BK_xhat, BK_value_error_pct, BK_achieved_percentile, BK_rank_error, ns_per_insert_BK, ns_close_BK
                            "%d,%d,%d%n",                                // truth_x_p, mpq_time_ms, bk_time_ms
                    r.datasetLabel, r.datasetMode, r.streamLength, r.alphabetConfigured, r.bucketsConfigured, r.distributionLabel,
                    r.zipfExponent, r.pTarget, r.datasetNgram,
                    r.epsTargetUsed, r.delta, r.deltaSampUsed, r.DhatUsed, B_suggested, n_req, E_nb, Var_nb, LB_nb,
                    r.Sactive, r.Bused, r.nbMPQ, r.epsMPQ, r.pLoMPQ, r.pHiMPQ, r.xLoMPQ, r.xHiMPQ,
                    r.xMPQ, r.pctErrMPQ, r.pAchMPQ, r.rankErrMPQ, inside, occOK, r.nsPerInsertMPQ, r.closeNsPerItemMPQ,
                    r.nbBK, r.epsBK, r.xLoBK, r.xHiBK, r.xBK, r.pctErrBK, r.pAchBK, r.rankErrBK, r.nsPerInsertBK, r.closeNsPerItemBK,
                    r.xTrue, r.tMPQMs, r.tBKMs
            );
        }


        // =========================
        //  Zipf helpers
        // =========================

        static ZipfDistribution makeZipfIntDist(int alphabetSize, double exponent, long seed) {
            if (alphabetSize <= 0) throw new IllegalArgumentException("alphabetSize must be > 0");
            if (exponent <= 0.0) throw new IllegalArgumentException("exponent must be > 0");
            RandomGenerator rg = new Well19937c(seed); // deterministic
            return new ZipfDistribution(rg, alphabetSize, exponent);
        }

        static int sampleZipfInt(ZipfDistribution dist) {
            return dist.sample() - 1; // map {1..N} -> {0..N-1}
        }

        // =========================
        //  Main: CLI compatible + optional auto-design
        // =========================

        public static void main(String[] args) throws IOException {
            Options opts = parseArgs(args);

            ArrayList<Long> datasetStream = null;
            int effectiveAlphabet = opts.alphabet;

            opts.datasetLabel = "synthetic-" + opts.dist.name().toLowerCase(Locale.ROOT);
            if (opts.datasetNgram > 1) {
                opts.datasetLabel = opts.datasetLabel + "-ng" + opts.datasetNgram;
            }

            if (opts.datasetPath != null) {
                DatasetStream ds = loadDataset(opts);
                datasetStream = ds.tokens();
                if (!opts.alphabetProvided) {
                    effectiveAlphabet = ds.alphabetSize();
                }
                if (!opts.streamLenProvided) {
                    opts.streamLen = datasetStream.size();
                } else if (opts.streamLen < datasetStream.size()) {
                    datasetStream = new ArrayList<>(datasetStream.subList(0, (int) opts.streamLen));
                }
                if (!opts.distinctProvided) {
                    opts.distinctGuess = null;
                }
                opts.streamLen = Math.min(opts.streamLen, datasetStream.size());
                EXPONENT = Double.NaN;
                opts.datasetLabel = opts.datasetPath.getFileName().toString();
                if (opts.datasetNgram > 1) {
                    opts.datasetLabel = opts.datasetLabel + "-ng" + opts.datasetNgram;
                }
            } else {
                effectiveAlphabet = opts.alphabet;
                EXPONENT = opts.zipfExponent;
                if (!opts.distinctProvided) {
                    opts.distinctGuess = effectiveAlphabet;
                }
            }
            opts.alphabet = effectiveAlphabet;

            if (!opts.deltaQsProvided) {
                opts.deltaQs = DEFAULT_DELTA_QS;
            }
            if (!opts.deltaSampsProvided) {
                opts.deltaSamps = DEFAULT_DELTA_SAMPS;
            }
            if (opts.deltaQs.size() != opts.deltaSamps.size()) {
                throw new IllegalArgumentException("delta_q and delta_samp must have the same number of entries");
            }

            if (Double.isNaN(SUMMARY_DELTA)) {
                SUMMARY_DELTA = opts.deltaQs.size() == 1 ? opts.deltaQs.get(0) : Double.NaN;
            }

            String distributionLabel = (opts.datasetPath != null) ? "dataset" : opts.dist.name().toLowerCase(Locale.ROOT);

            System.out.printf("Stream source: %s%n", opts.datasetPath != null ? opts.datasetPath + " (mode=" + opts.datasetMode + ")" : "synthetic (" + opts.dist + ")");
            if (opts.datasetPath != null) {
                String modeDesc = switch (opts.datasetMode) {
                    case SEGMENTS -> "each non-empty line treated as a token";
                    case WORDS -> "tokens split on whitespace";
                    case CHARS -> "individual UTF-16 characters";
                };
                System.out.printf("Tokenisation: %s%n", modeDesc);
            }
            System.out.printf("Stream length: %d, alphabet=%d, buckets=%d%n", opts.streamLen, effectiveAlphabet, opts.buckets);
            System.out.printf("Quantiles: %s  Rank eps targets: %s%n", opts.quantiles, opts.rankEpsTargets);
            System.out.printf("Dataset n-gram: %d%n", opts.datasetNgram);
            System.out.printf("delta_q values: %s  delta_samp values: %s%n", opts.deltaQs, opts.deltaSamps);

            List<TrialResult> allResults = new ArrayList<>(opts.trials * opts.deltaQs.size() * opts.quantiles.size());

            for (int qIdx = 0; qIdx < opts.quantiles.size(); qIdx++) {
                double quantile = opts.quantiles.get(qIdx);
                double epsTarget = opts.rankEpsTargets.size() == 1
                        ? opts.rankEpsTargets.get(0)
                        : opts.rankEpsTargets.get(Math.min(qIdx, opts.rankEpsTargets.size() - 1));
                double zipfForCall = opts.datasetPath != null ? Double.NaN : opts.zipfExponent;

                for (int z = 0; z < opts.deltaQs.size(); z++) {
                    double deltaQ = opts.deltaQs.get(z);
                    double deltaSamp = opts.deltaSamps.get(z);

                    List<TrialResult> comboResults = new ArrayList<>(opts.trials);

                    for (int i = 0; i < opts.trials; i++) {
                        long trialSeed = mix64(
                                opts.seed
                                        + i * 0x9E3779B97F4A7C15L
                                        + qIdx * 0x6A09E667F3BCC909L
                                        + z * 0xD1B54A32D192ED03L
                        );
                        TrialResult r = runExperiment(
                                datasetStream,
                                opts.streamLen,
                                effectiveAlphabet,
                                opts.buckets,
                                opts.dist,
                                zipfForCall,
                                quantile,
                                trialSeed,
                                deltaQ,
                                opts.verbose,
                                opts.autoDesign,
                                opts.distinctGuess,
                                epsTarget,
                                deltaSamp,
                                opts.datasetLabel,
                                opts.datasetPath != null ? opts.datasetMode : null,
                                distributionLabel,
                                opts.datasetNgram
                        );
                        comboResults.add(r);
                        allResults.add(r);
                    }

                    if (opts.autoDesign && !comboResults.isEmpty() && comboResults.get(0).design != null && opts.verbose) {
                        Utils.HopsDesignResult d = comboResults.get(0).design;
                        System.out.printf("Design summary (quantile %.5f, delta_q %.3f): suggested B=%d, n_req=%d, LB(n_b)=%d, E[n_b]=%.2f%s%n",
                                quantile,
                                deltaQ,
                                d.suggestedBuckets,
                                d.requiredSampleSize,
                                d.occupancyLowerBound,
                                d.expectedNonEmpty,
                                d.impossible ? " (WARNING: target impossible)" : "");
                    }

                    summarize(comboResults, quantile, deltaQ);
                }
            }

            if (opts.csvPath != null) {
                if (!opts.csvAppend) {
                    Files.deleteIfExists(opts.csvPath);
                }
                writeSuiteACsv(opts.csvPath, allResults);
            } else if (allResults.size() == 1) {
                printSuiteACsvHeader();
                printSuiteACsvRow(allResults.get(0));
            }
        }
    }
