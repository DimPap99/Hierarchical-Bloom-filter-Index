import PMIndex.IPMIndexing;
import PMIndex.StreamingSlidingWindowIndex;
import search.Pattern;
import tree.ssws.SuffixTree;

import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Random;

// Nanosecond-resolution benchmarks for the SSWSI baseline (docs/sswsi-results.md). The main harness times
// each query workload in whole milliseconds, which cannot resolve SSWSI's 1-30 ms workloads.
//
//   index <dataset file> <query dir> <n> <reps> [window] [legacy|amortized] [fastutil|fks] [chars|segments]
//       Streams the dataset like the harness into a fresh index per rep, after one warm-up rep, and prints the
//       median insert time. Then runs every <m>.<type>.txt workload once untimed, then each one best of 7, and
//       prints its hit total and a checksum over the sorted absolute positions of every query, so two builds or
//       two modes can be checked for identical answers. The options after <reps> can come in any order.
//       Defaults: window 2^20, legacy (StreamingSlidingWindowIndex), fastutil, chars.
//       chars: newlines dropped, one n-gram of characters per symbol (MultiQueryExperiment).
//       segments: one token per non-empty line, n consecutive tokens per symbol, queries split on spaces
//       (SegmentModeRunner).
//       amortized: the paper's structure (SlidingSuffixIndex), with the given child dictionary.
//   build [legacy|static] [fastutil|fks]
//       Times one suffix tree build on random 63-bit tokens, 2^4..2^20 of them, 100 or all distinct: legacy is
//       SuffixTree.build, static is Dc3SuffixTreeBuilder (the array tree of SlidingSuffixIndex).
//
// It uses reflection where the API changed after Phase 0, so it also compiles and runs against older checkouts
// (legacy only).
//
// Run it with a fixed, pre-touched heap: -Xms8g -Xmx8g -XX:+AlwaysPreTouch. Without pre-touching, queries that
// happen to allocate into heap pages the JVM has not written yet pay a page fault per fresh page; at n=1, m=160
// that made one mode 1.6x slower than another, depending only on how much the last build had allocated.
public final class SswsiNanoBenchmark {

    private static final int[] QUERY_LENGTHS = {10, 20, 40, 80, 160};
    private static final String[] QUERY_TYPES = {"uniform", "missing", "rare"};

    public static void main(String[] args) throws Exception {
        if (args.length >= 5 && args[0].equals("index")) {
            int window = 1 << 20;
            String mode = "legacy";
            String dictionary = "fastutil";
            boolean segments = false;
            for (int i = 5; i < args.length; i++) {
                String a = args[i].toLowerCase(Locale.ROOT);
                switch (a) {
                    case "legacy", "amortized" -> mode = a;
                    case "fastutil", "fks" -> dictionary = a;
                    case "chars" -> segments = false;
                    case "segments" -> segments = true;
                    default -> window = Integer.parseInt(a);
                }
            }
            index(Path.of(args[1]), Path.of(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), window,
                    mode, dictionary, segments);
        } else if (args.length >= 1 && args[0].equals("build")) {
            boolean legacy = args.length < 2 || args[1].equalsIgnoreCase("legacy");
            build(legacy, (args.length >= 3) ? args[2] : "fastutil");
        } else {
            System.err.println("usage: SswsiNanoBenchmark index <dataset file> <query dir> <n> <reps> [window]"
                    + " [legacy|amortized] [fastutil|fks] [chars|segments]");
            System.err.println("       SswsiNanoBenchmark build [legacy|static] [fastutil|fks]");
            System.exit(2);
        }
    }

    private static IPMIndexing newIndex(String mode, String dictionary, int window) throws Exception {
        if (mode.equals("legacy")) {
            return new StreamingSlidingWindowIndex(window);
        }
        Class<?> kind = Class.forName("tree.ssws.ChildDictionary$Kind");
        Object k = kind.getMethod("parse", String.class).invoke(null, dictionary);
        return (IPMIndexing) Class.forName("PMIndex.SlidingSuffixIndex")
                .getConstructor(int.class, kind).newInstance(window, k);
    }

    private static void index(Path dataset, Path queryDir, int n, int reps, int window, String mode,
                              String dictionary, boolean segments) throws Exception {
        String text = segments ? null : Files.readString(dataset).replace("\n", "").replace("\r", "");
        String[] tokens = segments ? nonEmptyLines(dataset).toArray(new String[0]) : null;
        int symbols = segments ? tokens.length - n + 1 : text.length() - n + 1;
        double[] insertMs = new double[reps];
        IPMIndexing index = null;
        for (int rep = -1; rep < reps; rep++) {
            index = null;
            System.gc();
            index = newIndex(mode, dictionary, window);
            long t0 = System.nanoTime();
            if (segments) {
                for (int i = 0; i + n <= tokens.length; i++) {
                    index.insert(Pattern.joinTokens(tokens, i, n));
                }
            } else {
                for (int i = 0; i + n <= text.length(); i++) {
                    index.insert(text.substring(i, i + n));
                }
            }
            if (rep >= 0) {
                insertMs[rep] = (System.nanoTime() - t0) / 1e6;
            }
        }
        Arrays.sort(insertMs);
        System.out.printf("%s%s %s n=%d window=%d symbols=%d insert ms: median %.0f min %.0f max %.0f (%.3f us/symbol median)%n",
                mode, mode.equals("legacy") ? "" : "/" + dictionary, segments ? "segments" : "chars", n, window,
                symbols, insertMs[reps / 2], insertMs[0], insertMs[reps - 1], insertMs[reps / 2] * 1e3 / symbols);

        long offset = 0;
        try {
            Method windowStart = index.getClass().getMethod("windowStartPosition");
            offset = (long) windowStart.invoke(index);       // report() is window-relative since Phase 0
        } catch (NoSuchMethodException olderCheckout) {
            // report() returned absolute positions
        }
        // Like the harness: segments-mode lines are split before the timer, the Pattern is built inside it.
        List<String> labels = new ArrayList<>();
        List<List<String>> workloads = new ArrayList<>();
        List<List<String[]>> splitWorkloads = new ArrayList<>();
        for (String type : QUERY_TYPES) {
            for (int m : QUERY_LENGTHS) {
                Path file = queryDir.resolve(m + "." + type + ".txt");
                if (!Files.exists(file)) {
                    continue;
                }
                List<String> queries = nonEmptyLines(file);
                List<String[]> split = new ArrayList<>();
                if (segments) {
                    for (String q : queries) {
                        split.add(q.split(" +"));
                    }
                }
                labels.add(String.format("%-7s m=%-3d", type, m));
                workloads.add(queries);
                splitWorkloads.add(split);
            }
        }
        // One untimed pass over every workload first, so all workloads are timed after the same warm-up. Cells
        // still differ between JVM runs (legacy most, by up to 2.4x on Gutenberg n=1): take medians over several.
        for (int w = 0; w < workloads.size(); w++) {
            runAll(index, workloads.get(w), splitWorkloads.get(w), n, segments);
        }
        for (int w = 0; w < workloads.size(); w++) {
            List<String> queries = workloads.get(w);
            List<String[]> split = splitWorkloads.get(w);
            long best = Long.MAX_VALUE;
            for (int rep = 0; rep < 7; rep++) {
                long t0 = System.nanoTime();
                runAll(index, queries, split, n, segments);
                best = Math.min(best, System.nanoTime() - t0);
            }
            long hits = 0, checksum = 0;
            for (int q = 0; q < queries.size(); q++) {
                final long base = offset;
                Pattern p = segments ? new Pattern(split.get(q), n) : new Pattern(queries.get(q), n);
                long[] positions = index.report(p).stream()
                        .mapToLong(x -> x + base).sorted().toArray();
                hits += positions.length;
                checksum = checksum * 31 + Arrays.hashCode(positions);
            }
            System.out.printf("  %s queries=%d best-of-7 %.3f ms (%.2f us/query) hits=%d checksum=%016x%n",
                    labels.get(w), queries.size(), best / 1e6, best / 1e3 / queries.size(), hits, checksum);
        }
    }

    private static void runAll(IPMIndexing index, List<String> queries, List<String[]> split, int n,
                               boolean segments) {
        for (int q = 0; q < queries.size(); q++) {
            index.report(segments ? new Pattern(split.get(q), n) : new Pattern(queries.get(q), n));
        }
    }

    // Non-empty lines; a trailing '\r' (CRLF files) is not part of the line.
    private static List<String> nonEmptyLines(Path file) throws Exception {
        List<String> out = new ArrayList<>();
        for (String line : Files.readAllLines(file, StandardCharsets.UTF_8)) {
            if (line.endsWith("\r")) {
                line = line.substring(0, line.length() - 1);
            }
            if (!line.isEmpty()) {
                out.add(line);
            }
        }
        return out;
    }

    private static void build(boolean legacy, String dictionary) throws Exception {
        Method oneArg = null, twoArg = null, staticBuild = null;
        Object builder = null, kind = null;
        if (legacy) {
            try {
                oneArg = SuffixTree.class.getMethod("build", long[].class);
            } catch (NoSuchMethodException olderCheckout) {
                twoArg = SuffixTree.class.getMethod("build", long[].class, int.class);
            }
        } else {
            Class<?> kindClass = Class.forName("tree.ssws.ChildDictionary$Kind");
            kind = kindClass.getMethod("parse", String.class).invoke(null, dictionary);
            Class<?> builderClass = Class.forName("tree.ssws.Dc3SuffixTreeBuilder");
            builder = builderClass.getConstructor().newInstance();
            staticBuild = builderClass.getMethod("build", long[].class, int.class, int.class, kindClass);
        }
        System.out.println(legacy ? "SuffixTree.build (legacy)" : "Dc3SuffixTreeBuilder, " + dictionary);
        Random rnd = new Random(2);
        for (int log : new int[]{4, 8, 12, 16, 18, 20}) {
            int n = 1 << log;
            for (int distinct : new int[]{100, n}) {
                if (distinct > n) {
                    continue;
                }
                long[] alphabet = new long[distinct];
                for (int i = 0; i < distinct; i++) {
                    alphabet[i] = rnd.nextLong() & Long.MAX_VALUE;
                }
                long[] tokens = new long[n];
                for (int i = 0; i < n; i++) {
                    tokens[i] = alphabet[(distinct == n) ? i : rnd.nextInt(distinct)];
                }
                int reps = Math.max(5, (1 << 22) / n);
                long best = Long.MAX_VALUE;
                for (int round = 0; round < 5; round++) {
                    long t0 = System.nanoTime();
                    for (int r = 0; r < reps; r++) {
                        if (staticBuild != null) {
                            staticBuild.invoke(builder, tokens, 0, n, kind);
                        } else if (oneArg != null) {
                            oneArg.invoke(null, (Object) tokens);
                        } else {
                            twoArg.invoke(null, tokens, 1 << 21);
                        }
                    }
                    best = Math.min(best, (System.nanoTime() - t0) / reps);
                }
                System.out.printf("n=2^%-2d distinct=%-8s build %10.1f us%n",
                        log, (distinct == n) ? "all" : distinct, best / 1e3);
            }
        }
    }
}
