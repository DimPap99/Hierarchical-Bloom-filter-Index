package tree.ssws;

/**
 * Range-minimum queries in O(1) over {@code a[0, n)}, in linear space. Used by {@link FarachColton} for LCP
 * queries over a suffix array and for lowest common ancestors.
 *
 * <p>Blocks of 32 entries. Inside a block, {@code mask[i]} marks the positions on the stack of suffix minima of
 * the block's prefix up to {@code i}: position {@code j <= i} is on it when no entry in {@code (j, i]} is smaller.
 * The minimum of {@code [l, i]} is then the lowest marked position at or after {@code l}, one
 * {@code numberOfTrailingZeros} away. Across blocks, a sparse table over the block minima.
 */
final class RangeMin {

    private static final int SHIFT = 5;
    private static final int BLOCK = 1 << SHIFT;

    private final int[] a;
    private final int[] mask;
    private final int[][] table;    // table[k][b]: position of the minimum of blocks b .. b + 2^k - 1

    RangeMin(int[] a, int n) {
        this.a = a;
        this.mask = new int[Math.max(n, 1)];
        int blocks = (n + BLOCK - 1) >>> SHIFT;
        int[] stack = new int[BLOCK];
        int[] base = new int[Math.max(blocks, 1)];
        for (int b = 0; b < blocks; b++) {
            int start = b << SHIFT;
            int end = Math.min(n, start + BLOCK);
            int top = 0;
            int bits = 0;
            for (int i = start; i < end; i++) {
                while (top > 0 && a[stack[top - 1]] > a[i]) {
                    bits &= ~(1 << (stack[--top] - start));
                }
                stack[top++] = i;
                bits |= 1 << (i - start);
                mask[i] = bits;
            }
            base[b] = start + Integer.numberOfTrailingZeros(mask[end - 1]);
        }
        int levels = (blocks <= 1) ? 1 : 32 - Integer.numberOfLeadingZeros(blocks);
        table = new int[levels][];
        table[0] = base;
        for (int k = 1; k < levels; k++) {
            int half = 1 << (k - 1);
            int[] prev = table[k - 1];
            int[] cur = new int[blocks - (1 << k) + 1];
            for (int b = 0; b < cur.length; b++) {
                cur[b] = smaller(prev[b], prev[b + half]);
            }
            table[k] = cur;
        }
    }

    private int smaller(int i, int j) {
        return (a[j] < a[i]) ? j : i;
    }

    private int inBlock(int l, int r) {
        return (r & -BLOCK) + Integer.numberOfTrailingZeros(mask[r] & (-1 << (l & (BLOCK - 1))));
    }

    /** Position of a minimum of {@code a[l..r]}, for {@code 0 <= l <= r < n}. */
    int argMin(int l, int r) {
        int bl = l >>> SHIFT;
        int br = r >>> SHIFT;
        if (bl == br) {
            return inBlock(l, r);
        }
        int best = smaller(inBlock(l, (bl << SHIFT) + BLOCK - 1), inBlock(br << SHIFT, r));
        if (br - bl > 1) {
            int k = 31 - Integer.numberOfLeadingZeros(br - bl - 1);
            best = smaller(best, smaller(table[k][bl + 1], table[k][br - (1 << k)]));
        }
        return best;
    }

    /** The minimum of {@code a[l..r]}; short ranges are scanned, which is cheaper than three table lookups. */
    int min(int l, int r) {
        return a[argMin(l, r)];
    }
}
