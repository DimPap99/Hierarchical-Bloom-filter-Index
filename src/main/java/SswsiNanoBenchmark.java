import PMIndex.StreamingSlidingWindowIndex;
import search.Pattern;
import tree.ssws.SuffixTree;

import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

// Nanosecond-resolution benchmarks for the SSWSI baseline (docs/sswsi-results.md). The main harness times
// each query workload in whole milliseconds, which cannot resolve SSWSI's 1-30 ms workloads.
//
//   index <dataset file> <query dir> <n> <reps> [window]
//       Streams the dataset like the harness (newlines dropped, one n-gram per symbol) into a fresh index
//       per rep, after one warm-up rep, and prints the median insert time. Then runs every
//       <m>.<type>.txt workload best of 7 and prints its hit total and a checksum over the sorted
//       absolute positions of every query, so two builds can be checked for identical answers.
//   build
//       Times SuffixTree.build on random 63-bit tokens, 2^4..2^20 of them, 100 or all distinct.
//
// It uses reflection where the API changed in Phase 0, so it also runs against older checkouts.
public final class SswsiNanoBenchmark {

    private static final int[] QUERY_LENGTHS = {10, 20, 40, 80, 160};
    private static final String[] QUERY_TYPES = {"uniform", "missing", "rare"};

    public static void main(String[] args) throws Exception {
        if (args.length >= 5 && args[0].equals("index")) {
            int window = (args.length >= 6) ? Integer.parseInt(args[5]) : 1 << 20;
            index(Path.of(args[1]), Path.of(args[2]), Integer.parseInt(args[3]), Integer.parseInt(args[4]), window);
        } else if (args.length == 1 && args[0].equals("build")) {
            build();
        } else {
            System.err.println("usage: SswsiNanoBenchmark index <dataset file> <query dir> <n> <reps> [window]");
            System.err.println("       SswsiNanoBenchmark build");
            System.exit(2);
        }
    }

    private static void index(Path dataset, Path queryDir, int n, int reps, int window) throws Exception {
        String text = Files.readString(dataset).replace("\n", "").replace("\r", "");
        double[] insertMs = new double[reps];
        StreamingSlidingWindowIndex index = null;
        for (int rep = -1; rep < reps; rep++) {
            index = null;
            System.gc();
            index = new StreamingSlidingWindowIndex(window);
            long t0 = System.nanoTime();
            for (int i = 0; i + n <= text.length(); i++) {
                index.insert(text.substring(i, i + n));
            }
            if (rep >= 0) {
                insertMs[rep] = (System.nanoTime() - t0) / 1e6;
            }
        }
        Arrays.sort(insertMs);
        int symbols = text.length() - n + 1;
        System.out.printf("n=%d symbols=%d insert ms: median %.0f min %.0f max %.0f (%.3f us/symbol median)%n",
                n, symbols, insertMs[reps / 2], insertMs[0], insertMs[reps - 1], insertMs[reps / 2] * 1e3 / symbols);

        long offset = 0;
        try {
            Method windowStart = index.getClass().getMethod("windowStartPosition");
            offset = (long) windowStart.invoke(index);       // report() is window-relative since Phase 0
        } catch (NoSuchMethodException olderCheckout) {
            // report() returned absolute positions
        }
        for (String type : QUERY_TYPES) {
            for (int m : QUERY_LENGTHS) {
                Path file = queryDir.resolve(m + "." + type + ".txt");
                if (!Files.exists(file)) {
                    continue;
                }
                List<String> queries = new ArrayList<>();
                for (String q : Files.readAllLines(file)) {
                    if (!q.isEmpty()) {
                        queries.add(q);
                    }
                }
                long best = Long.MAX_VALUE;
                for (int rep = 0; rep < 7; rep++) {
                    long t0 = System.nanoTime();
                    for (String q : queries) {
                        index.report(new Pattern(q, n));
                    }
                    best = Math.min(best, System.nanoTime() - t0);
                }
                long hits = 0, checksum = 0;
                for (String q : queries) {
                    final long base = offset;
                    long[] positions = index.report(new Pattern(q, n)).stream()
                            .mapToLong(p -> p + base).sorted().toArray();
                    hits += positions.length;
                    checksum = checksum * 31 + Arrays.hashCode(positions);
                }
                System.out.printf("  %-7s m=%-3d queries=%d best-of-7 %.3f ms (%.2f us/query) hits=%d checksum=%016x%n",
                        type, m, queries.size(), best / 1e6, best / 1e3 / queries.size(), hits, checksum);
            }
        }
    }

    private static void build() throws Exception {
        Method oneArg = null, twoArg = null;
        try {
            oneArg = SuffixTree.class.getMethod("build", long[].class);
        } catch (NoSuchMethodException olderCheckout) {
            twoArg = SuffixTree.class.getMethod("build", long[].class, int.class);
        }
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
                        if (oneArg != null) {
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
