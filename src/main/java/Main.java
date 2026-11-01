import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import PMIndex.SuffixTreeIndex;
import utilities.DataCatalog;
import utilities.ExperimentRunResult;
import utilities.IndexFactory;
import utilities.ResultVerifier;
import utilities.Workload;
import utilities.WorkloadRunner;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.function.IntFunction;

// Quick driver: runs HBI, and optionally a baseline, on one workload over a range of n-gram sizes.
//
//   java Main                                             pg2701, patterns of length 10, n-grams 1-3
//   java Main --dataset w19 --part 2 --pattern-length 40 --query-type rare --ngrams 4
//   java Main --dataset mawi19 --alphabet-base 800000     (read as segments: one record per line)
//   java Main --dataset adversarial/blocks/blocks_w19   (names can be nested folders)
//   java Main --data w19/1/1_W19.txt --queries w19/1/10.rare.txt   (explicit files, relative to the roots)
//   java Main --dataset w21 --schedule window             re-run the queries after every window
//   java Main --dataset w21 --baseline suffixtree --ngrams 2      compare answers with the online suffix tree
//   java Main --dataset w21 --skip-queries true           insert speed only
//
// Data selection flags are DataCatalog.Selection's (--dataset, --part, --pattern-length, --query-type,
// --reader, --data, --queries, --data-root, --query-root). The others:
//   --window <len> / --window-power <p>   default: the number in the dataset name (w19 -> 2^19), else 2^21
//   --tree <len>                          default: the window length
//   --alphabet-base <n>                   the alphabet for n-gram size n is base^n (default 75)
//   --fp <rate>  --conf <c>  --ngrams <n | lo-hi>
//   --runs <n>  --warmup <n>              timed runs and untimed warm-up runs per n-gram size
//   --schedule end|window                 query once after inserting everything (default), or every window
//   --baseline regex|suffixtree|none      default regex; RegexIndex reads characters only, so it is skipped
//                                         for segments. With the end schedule every index's answers are
//                                         checked against each other (all are exact).
//   --skip-queries true                   only insert, and report the insert speed
//   --verbose true                        print every query's matches
public final class Main {

    public static void main(String[] args) {
        Options o = Options.parse(args);
        Workload workload = o.selection.resolve();
        if (o.skipQueries) {
            workload = workload.withQueries(null);
        }

        System.out.println("Workload: " + workload);
        System.out.printf(Locale.ROOT,
                "Window: %d  Tree: %d  Alphabet base: %d  FP: %.3g  Confidence: %.3g  Schedule: %s  Baseline: %s%n",
                o.windowLength, o.treeLength, o.alphabetBase, o.fpRate, o.confidence, o.schedule, o.baseline);

        boolean verify = workload.hasQueries() && o.schedule.equals("end") && !o.baseline.equals("none");
        ResultVerifier.reset(verify);
        WorkloadRunner.Options runOptions = WorkloadRunner.Options.QUIET.verbose(o.verbose);

        for (int nGram = o.minNgram; nGram <= o.maxNgram; nGram++) {
            int alphabetSize = o.alphabetSizeFor(nGram);
            int n = nGram;
            System.out.printf(Locale.ROOT, "%nN-gram: %d  Alphabet: %d%n", nGram, alphabetSize);

            List<HBI> hbis = new ArrayList<>();
            List<ExperimentRunResult> hbiRuns = time(o, workload, nGram, runOptions, ignored -> {
                HBI hbi = IndexFactory.createMarkovHbi(o.windowLength, o.treeLength, alphabetSize, o.fpRate, o.confidence, n);
                hbi.stats().setCollecting(true);
                hbis.add(hbi);
                return hbi;
            });
            report("HBI", hbiRuns, workload.hasQueries());
            reportLp(hbis.subList(Math.min(o.warmup, hbis.size()), hbis.size()));

            IntFunction<IPMIndexing> baseline = o.baselineFor(workload);
            if (baseline != null) {
                int baselineNgram = o.baseline.equals("regex") ? 1 : nGram;
                report(o.baseline, time(o, workload, baselineNgram, runOptions, baseline), workload.hasQueries());
            }
        }

        if (verify) {
            System.out.println();
            System.out.println(ResultVerifier.summary());
        }
    }

    // Warm-up runs first (results dropped), then the timed runs.
    private static List<ExperimentRunResult> time(Options o, Workload workload, int nGram,
                                                  WorkloadRunner.Options runOptions, IntFunction<IPMIndexing> newIndex) {
        List<ExperimentRunResult> timed = new ArrayList<>();
        for (int run = 0; run < o.warmup + o.runs; run++) {
            ExperimentRunResult r = o.run(workload, newIndex.apply(run), nGram, runOptions);
            if (run >= o.warmup) {
                timed.add(r);
            }
        }
        return timed;
    }

    private static void report(String label, List<ExperimentRunResult> runs, boolean withQueries) {
        if (runs.isEmpty()) {
            return;
        }
        double insertMs = 0, perSymbolMs = 0, queryMs = 0, loadMs = 0;
        for (ExperimentRunResult r : runs) {
            insertMs += r.totalInsertTimeMs();
            perSymbolMs += r.avgInsertMsPerSymbol();
            queryMs += r.totalRunTimeMs();
            loadMs += r.avgQueryLoadMs();
        }
        int n = runs.size();
        double perSymbol = perSymbolMs / n;
        System.out.printf(Locale.ROOT, "%-10s insert %.3f ms (%.6f ms per n-gram, %.0f n-grams/s)",
                label, insertMs / n, perSymbol, perSymbol == 0 ? 0.0 : 1000.0 / perSymbol);
        if (withQueries) {
            System.out.printf(Locale.ROOT, ", queries %.3f ms total, %.3f ms per query load, avg query %.2f n-grams",
                    queryMs / n, loadMs / n, runs.get(0).avgQuerySize());
        }
        System.out.println();
    }

    private static void reportLp(List<HBI> timedRuns) {
        double lp = 0.0, alpha = 0.0;
        int samples = 0;
        for (HBI hbi : timedRuns) {
            if (hbi.stats().lpLevels().isEmpty()) continue;
            lp += average(hbi.stats().lpLevels());
            alpha += average(hbi.stats().alphas());
            samples++;
        }
        if (samples > 0) {
            System.out.printf(Locale.ROOT, "           avg LP level %.3f, avg alpha %.3f%n", lp / samples, alpha / samples);
        }
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

    private static final class Options {
        final DataCatalog.Selection selection = new DataCatalog.Selection("pg2701");
        int windowLength;
        int treeLength;
        int alphabetBase = 75;
        double fpRate = 0.001;
        double confidence = 0.999;
        int runs = 1;
        int warmup = 0;
        int minNgram = 1;
        int maxNgram = 3;
        String schedule = "end";
        String baseline = "regex";
        boolean skipQueries = false;
        boolean verbose = false;

        static Options parse(String[] args) {
            Options o = new Options();
            Integer window = null, windowPower = null, tree = null;
            for (int i = 0; i < args.length; i++) {
                String arg = args[i];
                if (!arg.startsWith("--")) {
                    throw new IllegalArgumentException("Unexpected argument " + arg);
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
                if (o.selection.accept(key, value)) {
                    continue;
                }
                switch (key) {
                    case "window" -> window = Integer.parseInt(value);
                    case "window-power" -> windowPower = Integer.parseInt(value);
                    case "tree" -> tree = Integer.parseInt(value);
                    case "alphabet-base" -> o.alphabetBase = Integer.parseInt(value);
                    case "fp" -> o.fpRate = Double.parseDouble(value);
                    case "conf" -> o.confidence = Double.parseDouble(value);
                    case "runs" -> o.runs = Integer.parseInt(value);
                    case "warmup" -> o.warmup = Integer.parseInt(value);
                    case "schedule" -> o.schedule = value.toLowerCase(Locale.ROOT);
                    case "baseline" -> o.baseline = value.toLowerCase(Locale.ROOT);
                    case "skip-queries" -> o.skipQueries = Boolean.parseBoolean(value);
                    case "verbose" -> o.verbose = Boolean.parseBoolean(value);
                    case "ngrams" -> {
                        String[] bounds = value.split("[-:]");
                        o.minNgram = Integer.parseInt(bounds[0]);
                        o.maxNgram = Integer.parseInt(bounds[bounds.length - 1]);
                    }
                    default -> throw new IllegalArgumentException("Unknown option --" + key);
                }
            }
            if (!o.schedule.equals("end") && !o.schedule.equals("window")) {
                throw new IllegalArgumentException("--schedule must be end or window");
            }
            if (!List.of("regex", "suffixtree", "none").contains(o.baseline)) {
                throw new IllegalArgumentException("--baseline must be regex, suffixtree or none");
            }
            if (o.runs < 1 || o.warmup < 0) {
                throw new IllegalArgumentException("--runs must be at least 1 and --warmup at least 0");
            }
            if (window != null) {
                o.windowLength = window;
            } else if (windowPower != null) {
                o.windowLength = 1 << windowPower;
            } else {
                o.windowLength = o.selection.windowLength(1 << 21);
            }
            o.treeLength = (tree != null) ? tree : o.windowLength;
            return o;
        }

        // A fresh baseline index per run, or null when there is none for this workload.
        IntFunction<IPMIndexing> baselineFor(Workload workload) {
            return switch (baseline) {
                case "regex" -> workload.mode() == Workload.ReaderMode.CHARS ? run -> new RegexIndex() : null;
                case "suffixtree" -> run -> new SuffixTreeIndex();
                default -> null;
            };
        }

        ExperimentRunResult run(Workload workload, IPMIndexing index, int nGram, WorkloadRunner.Options options) {
            return schedule.equals("window")
                    ? WorkloadRunner.stream(workload, index, nGram, windowLength, options)
                    : WorkloadRunner.run(workload, index, nGram, options);
        }

        int alphabetSizeFor(int nGram) {
            return (int) Math.min(Integer.MAX_VALUE, Math.pow(alphabetBase, nGram));
        }
    }
}
