package datagenerators;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Samples experiment windows from a long record stream (one record per line), as in the thesis
 * (Sec. 5.1): for each window size w, pick a start row uniformly at random and take the next w rows.
 * A start is redrawn if its window would overlap another window of the same size, so the windows of
 * one size are independent samples rather than mostly the same rows.
 *
 * <p>Writes {@code <data root>/<name><log2 w>/<i>/<i>_W<log2 w>.txt} for i = 1..count, plus a
 * {@code SOURCE.txt} next to the window directories that records the input, its row count, the seed
 * and every start row, so the windows can be regenerated.
 *
 * <p>Usage: {@code SampleWindows <records.txt> <data root> <name> <seed> <count> <log2 w>...}
 * e.g. {@code SampleWindows data/traces/202501011400.tokens.txt data mawi 1 5 19 20 21}
 */
public final class SampleWindows {

    private SampleWindows() {
    }

    public static void main(String[] args) throws IOException {
        if (args.length < 6) {
            System.err.println("usage: SampleWindows <records.txt> <data root> <name> <seed> <count> <log2 w>...");
            System.exit(2);
        }
        Path records = Path.of(args[0]);
        Path dataRoot = Path.of(args[1]);
        String name = args[2];
        long seed = Long.parseLong(args[3]);
        int count = Integer.parseInt(args[4]);
        List<Integer> powers = new ArrayList<>();
        for (int i = 5; i < args.length; i++) {
            powers.add(Integer.parseInt(args[i]));
        }

        long rows = countLines(records);
        System.out.printf("%s: %d rows%n", records, rows);

        Random rnd = new Random(seed);
        List<Window> windows = new ArrayList<>();
        for (int power : powers) {
            long w = 1L << power;
            if (w > rows) {
                throw new IllegalArgumentException("window 2^" + power + " is longer than the input (" + rows + " rows)");
            }
            for (int i = 1; i <= count; i++) {
                long start;
                int draws = 0;
                do {
                    if (++draws > 10_000) {
                        throw new IllegalArgumentException("cannot place " + count + " disjoint windows of 2^"
                                + power + " rows in " + rows + " rows");
                    }
                    start = (long) (rnd.nextDouble() * (rows - w + 1));
                } while (overlaps(windows, power, start, w));
                Path dir = dataRoot.resolve(name + power).resolve(Integer.toString(i));
                Files.createDirectories(dir);
                windows.add(new Window(power, i, start, w, dir.resolve(i + "_W" + power + ".txt")));
            }
        }

        // One pass over the input feeds every window whose row range covers the current row.
        for (Window win : windows) {
            win.out = Files.newBufferedWriter(win.file, StandardCharsets.UTF_8);
        }
        try (BufferedReader in = Files.newBufferedReader(records, StandardCharsets.UTF_8)) {
            String line;
            long row = 0;
            while ((line = in.readLine()) != null) {
                if (line.isEmpty()) {
                    continue;
                }
                for (Window win : windows) {
                    if (row >= win.start && row < win.start + win.length) {
                        win.out.write(line);
                        win.out.write('\n');
                    }
                }
                row++;
            }
        }
        for (Window win : windows) {
            win.out.close();
        }

        for (int power : powers) {
            Path source = dataRoot.resolve(name + power).resolve("SOURCE.txt");
            StringBuilder sb = new StringBuilder();
            sb.append("Windows of 2^").append(power).append(" rows sampled by datagenerators.SampleWindows\n");
            sb.append("input: ").append(records).append(" (").append(rows).append(" rows)\n");
            sb.append("seed: ").append(seed).append(", windows: ").append(count).append(", sizes: ").append(powers).append('\n');
            for (Window win : windows) {
                if (win.power == power) {
                    sb.append(win.index).append(": rows ").append(win.start).append(" to ")
                            .append(win.start + win.length - 1).append('\n');
                }
            }
            Files.writeString(source, sb.toString(), StandardCharsets.UTF_8);
            System.out.print(sb);
        }
    }

    private static boolean overlaps(List<Window> windows, int power, long start, long w) {
        for (Window other : windows) {
            if (other.power == power && start < other.start + other.length && other.start < start + w) {
                return true;
            }
        }
        return false;
    }

    private static long countLines(Path file) throws IOException {
        long n = 0;
        try (BufferedReader in = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String line;
            while ((line = in.readLine()) != null) {
                if (!line.isEmpty()) {
                    n++;
                }
            }
        }
        return n;
    }

    private static final class Window {
        final int power;
        final int index;
        final long start;
        final long length;
        final Path file;
        BufferedWriter out;

        Window(int power, int index, long start, long length, Path file) {
            this.power = power;
            this.index = index;
            this.start = start;
            this.length = length;
            this.file = file;
        }
    }
}
