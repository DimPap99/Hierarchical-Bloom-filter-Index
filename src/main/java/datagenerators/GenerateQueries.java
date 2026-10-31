package datagenerators;

import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Generates the query workloads of one window, in the layout the benchmarks read:
 * {@code <query dir>/<r>.uniform.txt}, {@code <r>.missing.txt} and {@code <r>.rare.txt} for
 * r = 10, 20, 40, 80, 160 (or the lengths given), with 500 queries each.
 *
 * <p>Symbols are characters ({@code chars}: newlines dropped, as the harness reads text) or records
 * ({@code segments}: one per line; a query line joins its records with single spaces).
 * <ul>
 *   <li><b>uniform</b>: substrings at uniformly random positions (thesis Sec. 5.1). All occur.</li>
 *   <li><b>missing</b>: the uniform samples, of which exactly half are turned into non-occurrences:
 *       30% of their symbols (rounded up), at distinct random positions, are replaced by symbols of the
 *       window, redrawn until the pattern does not occur in the window (thesis Sec. 5.2.1). A replacement
 *       is the symbol at a uniformly random window position, so symbols are drawn as often as they occur:
 *       the near misses look like the window's own text, and an index cannot reject them just because
 *       they contain a new or rare symbol.</li>
 *   <li><b>rare</b>: substrings containing a rare symbol. Symbol frequencies are taken from the lowest
 *       up until the positions holding them cover at least 1% of the window; a position is drawn
 *       uniformly among those, and the substring covers it at a random offset. All occur.</li>
 * </ul>
 * Deterministic for a given seed. Absence is checked exactly: a pattern is accepted as missing only
 * if its rolling hash matches no substring of the window, and equal strings have equal hashes.
 *
 * <p>Usage: {@code GenerateQueries <window file> <query dir> <chars|segments> <seed> [lengths]},
 * where lengths is a comma-separated list such as {@code 10,20,40,80,160,320}.
 */
public final class GenerateQueries {

    static final int[] LENGTHS = {10, 20, 40, 80, 160};
    static final int QUERIES = 500;
    static final double MUTATED_SHARE = 0.5;
    static final double MUTATED_SYMBOLS = 0.3;
    static final double RARE_COVERAGE = 0.01;

    private GenerateQueries() {
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 4 || args.length > 5 || !(args[2].equals("chars") || args[2].equals("segments"))) {
            System.err.println("usage: GenerateQueries <window file> <query dir> <chars|segments> <seed> [lengths]");
            System.exit(2);
        }
        Window window = Window.read(Path.of(args[0]), args[2].equals("segments"));
        Path out = Path.of(args[1]);
        long seed = Long.parseLong(args[3]);
        int[] lengths = (args.length == 5)
                ? Arrays.stream(args[4].split(",")).mapToInt(Integer::parseInt).toArray()
                : LENGTHS;
        Files.createDirectories(out);
        System.out.println(generate(window, out, seed, lengths));
    }

    static String generate(Window win, Path out, long seed, int[] lengths) throws IOException {
        int n = win.seq.length;
        int[] freq = new int[win.sigma];
        for (int s : win.seq) {
            freq[s]++;
        }
        int[] rarePool = rarePositions(win.seq, freq);
        long singletons = Arrays.stream(freq).filter(f -> f == 1).count();
        StringBuilder report = new StringBuilder(String.format(
                "%s: %d symbols, %d distinct (%d occur once); rare pool %d positions (%.2f%%)%n",
                win.source, n, win.sigma, singletons, rarePool.length, 100.0 * rarePool.length / n));

        for (int r : lengths) {
            if (r > n) {
                continue;
            }
            Random rnd = new Random(seed * 1_000_003L + r);
            long base = rnd.nextLong() | 1L;
            LongOpenHashSet present = substringHashes(win.seq, r, base);

            int[][] uniform = new int[QUERIES][];
            for (int q = 0; q < QUERIES; q++) {
                int s = rnd.nextInt(n - r + 1);
                uniform[q] = Arrays.copyOfRange(win.seq, s, s + r);
            }

            int[][] missing = new int[QUERIES][];
            Integer[] order = new Integer[QUERIES];
            for (int q = 0; q < QUERIES; q++) {
                order[q] = q;
                missing[q] = uniform[q];
            }
            shuffle(order, rnd);
            int mutate = (int) Math.round(QUERIES * MUTATED_SHARE);
            int replaced = (int) Math.ceil(r * MUTATED_SYMBOLS);
            long attempts = 0;
            for (int k = 0; k < mutate; k++) {
                int q = order[k];
                int[] mutant;
                do {
                    if (++attempts > 1_000_000L) {
                        throw new IllegalStateException("cannot build a missing pattern of length " + r
                                + " for " + win.source + ": the window has too few distinct symbols");
                    }
                    mutant = uniform[q].clone();
                    int[] positions = new int[r];
                    for (int i = 0; i < r; i++) {
                        positions[i] = i;
                    }
                    for (int i = 0; i < replaced; i++) {
                        int j = i + rnd.nextInt(r - i);
                        int t = positions[i]; positions[i] = positions[j]; positions[j] = t;
                        mutant[positions[i]] = win.seq[rnd.nextInt(n)];
                    }
                } while (present.contains(hash(mutant, base)));
                missing[q] = mutant;
            }

            int[][] rare = new int[QUERIES][];
            for (int q = 0; q < QUERIES; q++) {
                int p = rarePool[rnd.nextInt(rarePool.length)];
                int s = Math.max(0, Math.min(p - rnd.nextInt(r), n - r));
                rare[q] = Arrays.copyOfRange(win.seq, s, s + r);
            }

            write(out.resolve(r + ".uniform.txt"), uniform, win);
            write(out.resolve(r + ".missing.txt"), missing, win);
            write(out.resolve(r + ".rare.txt"), rare, win);
            report.append(String.format("  r=%-3d uniform %d, missing %d (%d mutated, %.1f draws each), rare %d%n",
                    r, QUERIES, QUERIES, mutate, attempts / (double) mutate, QUERIES));
        }
        return report.toString();
    }

    /** Positions holding the least frequent symbols, lowest frequencies first, until they cover RARE_COVERAGE. */
    static int[] rarePositions(int[] seq, int[] freq) {
        int maxFreq = 0;
        for (int f : freq) {
            maxFreq = Math.max(maxFreq, f);
        }
        long[] positionsAtFreq = new long[maxFreq + 1];
        for (int f : freq) {
            positionsAtFreq[f] += f;
        }
        long need = (long) Math.ceil(seq.length * RARE_COVERAGE);
        int threshold = 0;
        long covered = 0;
        while (covered < need && threshold < maxFreq) {
            threshold++;
            covered += positionsAtFreq[threshold];
        }
        List<Integer> pool = new ArrayList<>();
        for (int p = 0; p < seq.length; p++) {
            if (freq[seq[p]] <= threshold) {
                pool.add(p);
            }
        }
        return pool.stream().mapToInt(Integer::intValue).toArray();
    }

    static LongOpenHashSet substringHashes(int[] seq, int r, long base) {
        LongOpenHashSet set = new LongOpenHashSet(seq.length);
        long power = 1;
        for (int i = 0; i < r; i++) {
            power *= base;
        }
        long h = 0;
        for (int i = 0; i < seq.length; i++) {
            h = h * base + (seq[i] + 1L);
            if (i >= r) {
                h -= (seq[i - r] + 1L) * power;
            }
            if (i >= r - 1) {
                set.add(h);
            }
        }
        return set;
    }

    static long hash(int[] pattern, long base) {
        long h = 0;
        for (int v : pattern) {
            h = h * base + (v + 1L);
        }
        return h;
    }

    private static void shuffle(Integer[] a, Random rnd) {
        for (int i = a.length - 1; i > 0; i--) {
            int j = rnd.nextInt(i + 1);
            Integer t = a[i]; a[i] = a[j]; a[j] = t;
        }
    }

    private static void write(Path file, int[][] patterns, Window win) throws IOException {
        try (BufferedWriter w = Files.newBufferedWriter(file, StandardCharsets.UTF_8)) {
            for (int[] p : patterns) {
                for (int i = 0; i < p.length; i++) {
                    if (win.segments && i > 0) {
                        w.write(' ');
                    }
                    w.write(win.symbols[p[i]]);
                }
                w.write('\n');
            }
        }
    }

    /** A window as a sequence of symbol ids 0..sigma-1, with the spelling of each id. */
    static final class Window {
        final Path source;
        final boolean segments;
        final int[] seq;
        final String[] symbols;
        final int sigma;

        private Window(Path source, boolean segments, int[] seq, String[] symbols) {
            this.source = source;
            this.segments = segments;
            this.seq = seq;
            this.symbols = symbols;
            this.sigma = symbols.length;
        }

        static Window read(Path file, boolean segments) throws IOException {
            Object2IntOpenHashMap<String> ids = new Object2IntOpenHashMap<>();
            ids.defaultReturnValue(-1);
            List<String> spelled = new ArrayList<>();
            int[] seq = new int[1 << 16];
            int n = 0;
            try (BufferedReader in = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
                if (segments) {
                    String line;
                    while ((line = in.readLine()) != null) {
                        if (line.isEmpty()) {
                            continue;
                        }
                        if (line.indexOf(' ') >= 0) {
                            throw new IOException("record contains a space, which query lines use as separator: " + line);
                        }
                        if (n == seq.length) {
                            seq = Arrays.copyOf(seq, n * 2);
                        }
                        seq[n++] = id(ids, spelled, line);
                    }
                } else {
                    int c;
                    while ((c = in.read()) != -1) {
                        if (c == '\n' || c == '\r') {
                            continue;
                        }
                        if (n == seq.length) {
                            seq = Arrays.copyOf(seq, n * 2);
                        }
                        seq[n++] = id(ids, spelled, String.valueOf((char) c));
                    }
                }
            }
            return new Window(file, segments, Arrays.copyOf(seq, n), spelled.toArray(new String[0]));
        }

        private static int id(Object2IntOpenHashMap<String> ids, List<String> spelled, String symbol) {
            int id = ids.getInt(symbol);
            if (id < 0) {
                id = spelled.size();
                ids.put(symbol, id);
                spelled.add(symbol);
            }
            return id;
        }
    }
}
