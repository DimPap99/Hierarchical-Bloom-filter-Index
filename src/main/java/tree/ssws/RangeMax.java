package tree.ssws;

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.Arrays;

/**
 * Range-maximum structure over a tree's suffix array, for reporting only the occurrences that start at or after a
 * threshold (Bille et al., Sec. 3.2): find the largest start in the range; if it is below the threshold, stop;
 * otherwise report it and continue on both sides of it. That costs O(1) queries per reported occurrence, plus one.
 *
 * <p>Entries at or above {@code limit} count as minus infinity. For a boundary tree, {@code limit} is the length of
 * its left half, which makes this the paper's BL' (the suffix's start if it starts in the left half, minus infinity
 * otherwise) without a second array.
 *
 * <p>Linear space: blocks of {@value #BLOCK} entries, a sparse table over the block maxima, and a scan inside a
 * block. A full sparse table over the entries would cost about 21 ints per entry at w = 2^21.
 */
public final class RangeMax {

    private static final int SHIFT = 5;
    private static final int BLOCK = 1 << SHIFT;
    /** Ranges of at most this many entries are scanned instead of split. */
    private static final int SCAN = 2 * BLOCK;

    private final int[] a;
    private final int limit;
    private final int[][] table;    // table[k][b]: position of the maximum of blocks b .. b + 2^k - 1

    /** Over all suffix starts of the tree. */
    public RangeMax(StaticSuffixTree tree) {
        this(tree.sa, Integer.MAX_VALUE);
    }

    /** Over the suffix starts below {@code limit}; the others never count as the maximum and are never reported. */
    public RangeMax(StaticSuffixTree tree, int limit) {
        this(tree.sa, limit);
    }

    RangeMax(int[] a, int limit) {
        this.a = a;
        this.limit = limit;
        int n = a.length;
        if (n <= SCAN) {
            table = null;   // every query scans
            return;
        }
        int blocks = (n + BLOCK - 1) >>> SHIFT;
        int levels = 32 - Integer.numberOfLeadingZeros(blocks);
        table = new int[levels][];
        int[] base = new int[blocks];
        for (int b = 0; b < blocks; b++) {
            base[b] = scanMax(b << SHIFT, Math.min(n, (b + 1) << SHIFT) - 1);
        }
        table[0] = base;
        for (int k = 1; k < levels; k++) {
            int half = 1 << (k - 1);
            int[] prev = table[k - 1];
            int[] cur = new int[blocks - (1 << k) + 1];
            for (int b = 0; b < cur.length; b++) {
                cur[b] = better(prev[b], prev[b + half]);
            }
            table[k] = cur;
        }
    }

    private int value(int i) {
        int v = a[i];
        return (v < limit) ? v : -1;
    }

    private int better(int i, int j) {
        return (value(i) >= value(j)) ? i : j;
    }

    private int scanMax(int l, int r) {
        int best = l;
        int bestValue = value(l);
        for (int i = l + 1; i <= r; i++) {
            int v = value(i);
            if (v > bestValue) {
                bestValue = v;
                best = i;
            }
        }
        return best;
    }

    /** Position of a maximum of {@code [l, r]}; entries at or above the limit count as -1. */
    int argMax(int l, int r) {
        int bl = l >>> SHIFT;
        int br = r >>> SHIFT;
        if (br - bl <= 1) {
            return scanMax(l, r);
        }
        int best = better(scanMax(l, ((bl + 1) << SHIFT) - 1), scanMax(br << SHIFT, r));
        int k = 31 - Integer.numberOfLeadingZeros(br - bl - 1);
        best = better(best, table[k][bl + 1]);
        return better(best, table[k][br - (1 << k)]);
    }

    /**
     * Appends {@code base + a[i]} for every {@code i} in {@code [l, r]} with {@code threshold <= a[i] < limit}, in
     * no particular order.
     */
    public void reportAtLeast(int l, int r, int threshold, IntArrayList out, int base) {
        if (threshold < 0) {
            threshold = 0;   // -1 stands for minus infinity
        }
        if (r - l < SCAN) {
            scanReport(l, r, threshold, out, base);
            return;
        }
        int[] stack = STACK.get();
        int sp = 0;
        stack[sp++] = l;
        stack[sp++] = r;
        while (sp > 0) {
            int hi = stack[--sp];
            int lo = stack[--sp];
            if (hi - lo < SCAN) {
                scanReport(lo, hi, threshold, out, base);
                continue;
            }
            int j = argMax(lo, hi);
            int v = value(j);
            if (v < threshold) {
                continue;
            }
            out.add(base + v);
            if (sp + 4 > stack.length) {
                stack = Arrays.copyOf(stack, 2 * stack.length);
                STACK.set(stack);
            }
            stack[sp++] = lo;            // j - lo >= 0 and hi - j >= 0 entries remain; empty ranges stop at once
            stack[sp++] = j - 1;
            stack[sp++] = j + 1;
            stack[sp++] = hi;
        }
    }

    private void scanReport(int l, int r, int threshold, IntArrayList out, int base) {
        for (int i = l; i <= r; i++) {
            int v = a[i];
            if (v >= threshold && v < limit) {
                out.add(base + v);
            }
        }
    }

    private static final ThreadLocal<int[]> STACK = ThreadLocal.withInitial(() -> new int[64]);
}
