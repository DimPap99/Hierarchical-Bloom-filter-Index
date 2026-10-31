package tree.ssws;

import java.util.Arrays;

/**
 * Rank-space reduction for suffix tree construction.
 *
 * <p>Renames the tokens to 1..sigma and appends the terminator 0, which is the input format of
 * {@link SuffixArrayBuilder}. Names are given in order of first occurrence. Any injective naming
 * gives the same suffix tree, but this one keeps the suffix array close to sequential when most
 * tokens are distinct (n-gram mode), which makes DC3, the LCP pass and the tree build several times
 * more cache-friendly than names in sorted order.
 *
 * <p>This replaces the hashing of Bille et al., Lemma 2, with a deterministic LSD radix sort that
 * groups equal tokens: it is linear, needs no randomness and never has to retry. Up to {@link #TINY}
 * tokens are named by direct comparison. The sort uses 8-bit digits below {@link #MEDIUM} tokens,
 * 11-bit digits below {@link #LARGE} and 16-bit digits above, which keeps the per-build cost of
 * clearing the digit counts small next to the tokens themselves; it skips a pass when every token
 * has the same digit there.
 */
final class RankReducer {

    static final int TINY = 32;
    static final int MEDIUM = 1 << 10;
    static final int LARGE = 1 << 16;

    private RankReducer() {
    }

    /** Returns {@code r} of length {@code n + 1}: {@code r[i]} names {@code text[i]}, and {@code r[n] = 0}. */
    static int[] rankWithTerminator(long[] text) {
        return rankWithTerminator(text, 0, text.length);
    }

    /** The same for the {@code n} tokens {@code text[offset, offset + n)}: {@code r[i]} names {@code text[offset + i]}. */
    static int[] rankWithTerminator(long[] text, int offset, int n) {
        return rankWithTerminator(text, offset, n, 0);
    }

    /** The same, with {@code extra} more zeros after the terminator (Farach-Colton pads to an even length). */
    static int[] rankWithTerminator(long[] text, int offset, int n, int extra) {
        int[] names = new int[n + 1 + extra];
        if (n <= TINY) {
            nameByComparison(text, offset, n, names);
        } else {
            nameBySorting(text, offset, n, names);
        }
        return names;
    }

    private static void nameByComparison(long[] text, int offset, int n, int[] names) {
        int next = 0;
        for (int i = 0; i < n; i++) {
            long v = text[offset + i];
            int j = 0;
            while (text[offset + j] != v) {
                j++;
            }
            names[i] = (j == i) ? ++next : names[j];
        }
    }

    private static void nameBySorting(long[] text, int offset, int n, int[] names) {
        Workspace ws = WORKSPACE.get();
        ws.ensure(n);
        System.arraycopy(text, offset, ws.keysA, 0, n);
        for (int i = 0; i < n; i++) {
            ws.idxA[i] = i;
        }
        boolean inB = radixSort(ws, n, (n < MEDIUM) ? 8 : (n < LARGE) ? 11 : 16);
        long[] keys = inB ? ws.keysB : ws.keysA;
        int[] idx = inB ? ws.idxB : ws.idxA;

        // The sort is stable, so each run of equal keys starts with the key's first position f.
        // Store f + 1 for every member of the run.
        int k = 0;
        while (k < n) {
            long v = keys[k];
            int first = idx[k] + 1;
            do {
                names[idx[k]] = first;
                k++;
            } while (k < n && keys[k] == v);
        }
        // Left to right, a first occurrence gets the next name; any later one copies the name already
        // given to its first occurrence.
        int next = 0;
        for (int i = 0; i < n; i++) {
            int f = names[i] - 1;
            names[i] = (f == i) ? ++next : names[f];
        }
    }

    /**
     * Stable LSD radix sort of {@code ws.keysA[0..n)} as unsigned values, carrying the positions in
     * {@code ws.idxA}. Returns true if the result ended up in the B buffers.
     */
    private static boolean radixSort(Workspace ws, int n, int bits) {
        int radix = 1 << bits;
        int mask = radix - 1;
        int digits = (Long.SIZE + bits - 1) / bits;
        int[] count = ws.count(digits * radix);

        long[] keysA = ws.keysA;
        for (int i = 0; i < n; i++) {
            long v = keysA[i];
            for (int d = 0; d < digits; d++) {
                count[d * radix + (int) ((v >>> (d * bits)) & mask)]++;
            }
        }

        long[] srcKeys = ws.keysA, dstKeys = ws.keysB;
        int[] srcIdx = ws.idxA, dstIdx = ws.idxB;
        boolean inB = false;
        for (int d = 0; d < digits; d++) {
            int base = d * radix;
            int shift = d * bits;
            // Every token has the same digit here: the pass would not move anything.
            if (count[base + (int) ((srcKeys[0] >>> shift) & mask)] == n) {
                continue;
            }
            int sum = 0;
            for (int b = 0; b < radix; b++) {
                int c = count[base + b];
                count[base + b] = sum;
                sum += c;
            }
            for (int i = 0; i < n; i++) {
                long v = srcKeys[i];
                int pos = count[base + (int) ((v >>> shift) & mask)]++;
                dstKeys[pos] = v;
                dstIdx[pos] = srcIdx[i];
            }
            long[] tk = srcKeys; srcKeys = dstKeys; dstKeys = tk;
            int[] ti = srcIdx; srcIdx = dstIdx; dstIdx = ti;
            inB = !inB;
        }
        return inB;
    }

    private static final ThreadLocal<Workspace> WORKSPACE = ThreadLocal.withInitial(Workspace::new);

    /** Buffers reused across builds on the same thread, so a build does not zero-fill fresh arrays. */
    private static final class Workspace {
        long[] keysA = new long[0];
        long[] keysB = new long[0];
        int[] idxA = new int[0];
        int[] idxB = new int[0];
        private int[] count = new int[0];

        void ensure(int n) {
            if (keysA.length < n) {
                int cap = Math.max(n, keysA.length + (keysA.length >> 1));
                keysA = new long[cap];
                keysB = new long[cap];
                idxA = new int[cap];
                idxB = new int[cap];
            }
        }

        int[] count(int size) {
            if (count.length < size) {
                count = new int[size];
            } else {
                Arrays.fill(count, 0, size, 0);
            }
            return count;
        }
    }
}
