import PMIndex.IPMIndexing;
import PMIndex.StreamingSlidingWindowIndex;
import org.openjdk.jol.info.GraphLayout;
import org.openjdk.jol.util.Multiset;
import search.Pattern;
import tree.ssws.SuffixTree;
import utilities.MemUtil;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;

// Nanosecond-resolution benchmarks for the SSWSI baseline (docs/sswsi-results.md). The main harness times
// each query workload in whole milliseconds, which cannot resolve SSWSI's 1-30 ms workloads.
//
//   index <dataset file> <query dir> <n> <reps> [window] [legacy|amortized] [fc|dc3] [fastutil|fks] [chars|segments]
//       Streams the dataset like the harness into a fresh index per rep, after one warm-up rep, and prints the
//       median insert time. Then runs every <m>.<type>.txt workload once untimed, then each one best of 7, and
//       prints its hit total and a checksum over the sorted absolute positions of every query, so two builds or
//       two modes can be checked for identical answers. The options after <reps> can come in any order.
//       Defaults: window 2^20, legacy (StreamingSlidingWindowIndex), dc3, fks, chars.
//       chars: newlines dropped, one n-gram of characters per symbol (MultiQueryExperiment).
//       segments: one token per non-empty line, n consecutive tokens per symbol, queries split on spaces
//       (SegmentModeRunner).
//       amortized: the paper's structure (SlidingSuffixIndex), with the given builder (Farach-Colton or DC3)
//       and child dictionary.
//   build [legacy|static] [fc|dc3] [fastutil|fks]
//       Times one suffix tree build on random 63-bit tokens, 2^4..2^21 of them, 100 or all distinct: legacy is
//       SuffixTree.build, static is the array tree of SlidingSuffixIndex with the given builder and dictionary.
//   memory <dataset file> <n> [window] [legacy|amortized] [fc|dc3] [fastutil|fks] [chars|segments]
//       Streams the dataset like `index` into one index, then prints what MemoryBenchmark measures: the JOL size of
//       everything reachable from the index, per stored symbol, and its largest classes; for amortized, also the
//       stored text, segment trees, boundary trees and range-maximum structures. Then the builders' per-thread
//       work arrays, which the index does not reference, so JOL's total leaves them out.
//       JOL on JDK 25 needs -Djdk.attach.allowAttachSelf=true -XX:+EnableDynamicAgentLoading.
//
// It uses reflection where the API changed after Phase 0, so it also compiles and runs against older checkouts
// (legacy only), compiled together with this checkout's utilities/MemUtil.java.
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
            String builder = "dc3";
            String dictionary = "fks";
            boolean segments = false;
            for (int i = 5; i < args.length; i++) {
                String a = args[i].toLowerCase(Locale.ROOT);
                switch (a) {
                    case "legacy", "amortized" -> mode = a;
                    case "fc", "dc3" -> builder = a;
                    case "fastutil", "fks" -> dictionary = a;
                    case "chars" -> segments = false;
                    case "segments" -> segments = true;
                    default -> window = Integer.parseInt(a);
                }
            }
            index(Path.of(args[1]), Path.of(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), window,
                    mode, builder, dictionary, segments);
        } else if (args.length >= 1 && args[0].equals("build")) {
            boolean legacy = true;
            String builder = "dc3";
            String dictionary = "fks";
            for (int i = 1; i < args.length; i++) {
                String a = args[i].toLowerCase(Locale.ROOT);
                switch (a) {
                    case "legacy" -> legacy = true;
                    case "static" -> legacy = false;
                    case "fc", "dc3" -> builder = a;
                    case "fastutil", "fks" -> dictionary = a;
                    default -> throw new IllegalArgumentException("unknown build option " + a);
                }
            }
            build(legacy, builder, dictionary);
        } else if (args.length >= 3 && args[0].equals("memory")) {
            int window = 1 << 20;
            String mode = "legacy";
            String builder = "dc3";
            String dictionary = "fks";
            boolean segments = false;
            for (int i = 3; i < args.length; i++) {
                String a = args[i].toLowerCase(Locale.ROOT);
                switch (a) {
                    case "legacy", "amortized" -> mode = a;
                    case "fc", "dc3" -> builder = a;
                    case "fastutil", "fks" -> dictionary = a;
                    case "chars" -> segments = false;
                    case "segments" -> segments = true;
                    default -> window = Integer.parseInt(a);
                }
            }
            memory(Path.of(args[1]), Integer.parseInt(args[2]), window, mode, builder, dictionary, segments);
        } else {
            System.err.println("usage: SswsiNanoBenchmark index <dataset file> <query dir> <n> <reps> [window]"
                    + " [legacy|amortized] [fc|dc3] [fastutil|fks] [chars|segments]");
            System.err.println("       SswsiNanoBenchmark build [legacy|static] [fc|dc3] [fastutil|fks]");
            System.err.println("       SswsiNanoBenchmark memory <dataset file> <n> [window] [legacy|amortized]"
                    + " [fc|dc3] [fastutil|fks] [chars|segments]");
            System.exit(2);
        }
    }

    private static final double MIB = 1024.0 * 1024.0;

    private static void memory(Path dataset, int n, int window, String mode, String builder, String dictionary,
                               boolean segments) throws Exception {
        IPMIndexing index = newIndex(mode, builder, dictionary, window);
        int symbols = 0;
        if (segments) {
            String[] tokens = nonEmptyLines(dataset).toArray(new String[0]);
            for (int i = 0; i + n <= tokens.length; i++, symbols++) {
                index.insert(Pattern.joinTokens(tokens, i, n));
            }
        } else {
            String text = Files.readString(dataset).replace("\n", "").replace("\r", "");
            for (int i = 0; i + n <= text.length(); i++, symbols++) {
                index.insert(text.substring(i, i + n));
            }
        }
        // Symbols stored: S for amortized (w to 2w - 1 of them once w have arrived), the segments for legacy.
        long stored;
        if (mode.equals("legacy")) {
            Field f = StreamingSlidingWindowIndex.class.getDeclaredField("totalStoredLength");
            f.setAccessible(true);
            stored = ((Number) f.get(index)).longValue();
        } else {
            stored = (long) index.getClass().getMethod("storedSymbols").invoke(index);
        }

        GraphLayout layout = GraphLayout.parseInstance(index);
        long total = layout.totalSize();
        System.out.printf("%s%s %s n=%d window=%d symbols=%d stored=%d: %.1f MiB, %.1f B per stored symbol%n",
                mode, mode.equals("legacy") ? "" : "/" + builder + "/" + dictionary, segments ? "segments" : "chars",
                n, window, symbols, stored, total / MIB, total / (double) stored);
        if (!mode.equals("legacy")) {
            @SuppressWarnings("unchecked")
            Map<String, Long> parts = (Map<String, Long>) index.getClass()
                    .getMethod("estimateRetainedBytesByPart").invoke(index);
            for (Map.Entry<String, Long> e : parts.entrySet()) {
                System.out.printf("  %-15s %8.1f MiB %6.1f B per stored symbol%n", e.getKey(), e.getValue() / MIB,
                        e.getValue() / (double) stored);
            }
        }
        Multiset<Class<?>> sizes = layout.getClassSizes();
        Multiset<Class<?>> counts = layout.getClassCounts();
        List<Class<?>> classes = new ArrayList<>(sizes.keys());
        classes.sort(Comparator.comparingLong((Class<?> c) -> sizes.count(c)).reversed());
        System.out.println("  largest classes:");
        for (Class<?> c : classes.subList(0, Math.min(8, classes.size()))) {
            System.out.printf("    %-60s %10d objects %8.1f MiB%n", c.getName(), counts.count(c), sizes.count(c) / MIB);
        }

        // Work arrays kept per thread between builds, sized to the largest tree built so far.
        long scratchTotal = 0;
        System.out.println("  per-thread work arrays (not in the total):");
        for (Map.Entry<String, Long> e : MemUtil.sswsiWorkArrayBytes().entrySet()) {
            scratchTotal += e.getValue();
            System.out.printf("    %-45s %8.1f MiB %6.1f B per window symbol%n", e.getKey(), e.getValue() / MIB,
                    e.getValue() / (double) window);
        }
        System.out.printf("    %-45s %8.1f MiB %6.1f B per window symbol%n", "total", scratchTotal / MIB,
                scratchTotal / (double) window);
    }

    private static IPMIndexing newIndex(String mode, String builder, String dictionary, int window)
            throws Exception {
        if (mode.equals("legacy")) {
            return new StreamingSlidingWindowIndex(window);
        }
        Class<?> dictKind = Class.forName("tree.ssws.ChildDictionary$Kind");
        Class<?> builderKind = Class.forName("tree.ssws.SuffixTreeBuilder$Kind");
        Class<?> builderType = Class.forName("tree.ssws.SuffixTreeBuilder");
        Object d = dictKind.getMethod("parse", String.class).invoke(null, dictionary);
        Object b = builderKind.getMethod("create").invoke(builderKind.getMethod("parse", String.class).invoke(null, builder));
        return (IPMIndexing) Class.forName("PMIndex.SlidingSuffixIndex")
                .getConstructor(int.class, builderType, dictKind).newInstance(window, b, d);
    }

    private static void index(Path dataset, Path queryDir, int n, int reps, int window, String mode,
                              String builder, String dictionary, boolean segments) throws Exception {
        String text = segments ? null : Files.readString(dataset).replace("\n", "").replace("\r", "");
        String[] tokens = segments ? nonEmptyLines(dataset).toArray(new String[0]) : null;
        int symbols = segments ? tokens.length - n + 1 : text.length() - n + 1;
        double[] insertMs = new double[reps];
        IPMIndexing index = null;
        for (int rep = -1; rep < reps; rep++) {
            index = null;
            System.gc();
            index = newIndex(mode, builder, dictionary, window);
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
                mode, mode.equals("legacy") ? "" : "/" + builder + "/" + dictionary, segments ? "segments" : "chars", n, window,
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

    private static void build(boolean legacy, String builderName, String dictionary) throws Exception {
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
            Class<?> builderKind = Class.forName("tree.ssws.SuffixTreeBuilder$Kind");
            builder = builderKind.getMethod("create")
                    .invoke(builderKind.getMethod("parse", String.class).invoke(null, builderName));
            staticBuild = Class.forName("tree.ssws.SuffixTreeBuilder")
                    .getMethod("build", long[].class, int.class, int.class, kindClass);
        }
        System.out.println(legacy ? "SuffixTree.build (legacy)" : "static tree, " + builderName + ", " + dictionary);
        Random rnd = new Random(2);
        for (int log : new int[]{4, 8, 10, 12, 14, 16, 18, 20, 21}) {
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
