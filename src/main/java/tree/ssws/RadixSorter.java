package tree.ssws;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * RadixSorter
 *
 * Helper for:
 *  - Integer Least Significant Digit radix sort (stable)
 *  - Rank transform
 *  - Parallel radix sort of (key,index) where key is 64-bit unsigned
 *
 * The sortParallel version below has been optimized so that it only
 * processes as many 8-bit passes as are actually needed for the
 * observed key range, instead of unconditionally doing 8 passes.
 */
public final class RadixSorter {

    private static final int RADIX = 256;      // 8 bits per pass
    private static final int BITS_PER_BYTE = 8;

    private RadixSorter() {
        // utility
    }

    /**
     * Stable Least Significant Digit radix sort in-place on an int[].
     * Handles negative values by offsetting them so bucket indices stay >= 0.
     */
    public static void sortInPlace(int[] values) {
        if (values == null || values.length <= 1) {
            return;
        }
        int min = values[0];
        int max = values[0];
        for (int v : values) {
            if (v < min) {
                min = v;
            }
            if (v > max) {
                max = v;
            }
        }
        int offset = (min < 0) ? -min : 0;
        int maxShifted = max + offset;

        int passes = 0;
        while ((maxShifted >> (passes * 8)) > 0) {
            passes++;
        }
        if (passes == 0) {
            passes = 1;
        }

        int[] buffer = new int[values.length];
        int[] count = new int[RADIX];

        for (int pass = 0; pass < passes; pass++) {
            Arrays.fill(count, 0);
            int shift = pass * 8;
            for (int value : values) {
                int shifted = value + offset;
                int bucket = (shifted >> shift) & 0xFF;
                count[bucket]++;
            }
            for (int i = 1; i < RADIX; i++) {
                count[i] += count[i - 1];
            }
            for (int i = values.length - 1; i >= 0; i--) {
                int shifted = values[i] + offset;
                int bucket = (shifted >> shift) & 0xFF;
                buffer[--count[bucket]] = values[i];
            }
            System.arraycopy(buffer, 0, values, 0, values.length);
        }
    }

    /**
     * Returns a sorted copy of values using sortInPlace logic.
     */
    public static int[] sortedCopy(int[] values) {
        if (values == null) {
            return null;
        }
        int[] copy = Arrays.copyOf(values, values.length);
        sortInPlace(copy);
        return copy;
    }

    /**
     * Stable Least Significant Digit radix sort for (key,index) pairs.
     * Reorders both arrays in parallel so that keys[] is sorted as
     * unsigned 64-bit integers and indices[] follows the same order.
     *
     * Optimization:
     *   We detect how many bytes of significance there actually are
     *   across all keys and only do that many byte passes, instead of
     *   always doing 8 passes unconditionally.
     */
    public static void sortParallel(long[] keys, int[] indices) {
        if (keys == null || indices == null) {
            throw new IllegalArgumentException("keys and indices must be non-null");
        }
        if (keys.length != indices.length) {
            throw new IllegalArgumentException("keys and indices must have the same length");
        }

        final int n = keys.length;
        if (n <= 1) {
            return;
        }

        long[] keyBuffer = new long[n];
        int[] indexBuffer = new int[n];
        int[] count = new int[RADIX];

        // Figure out how many 8-bit passes we actually need.
        // We OR all keys together so maxBits has a 1 in any bit
        // position that was ever set in any key.
        long maxBits = 0L;
        for (long k : keys) {
            maxBits |= k;
        }

        // Find the highest non-zero byte.
        int highestByte = 7;
        while (highestByte > 0) {
            long mask = 0xFFL << (highestByte * BITS_PER_BYTE);
            if ((maxBits & mask) != 0L) {
                break;
            }
            highestByte--;
        }

        // Perform LSD passes up to highestByte inclusive.
        for (int bytePos = 0; bytePos <= highestByte; bytePos++) {
            Arrays.fill(count, 0);

            int shift = bytePos * BITS_PER_BYTE;

            // count
            for (long key : keys) {
                int bucket = (int) ((key >>> shift) & 0xFFL);
                count[bucket]++;
            }

            // exclusive prefix sum
            int total = 0;
            for (int i = 0; i < RADIX; i++) {
                int c = count[i];
                count[i] = total;
                total += c;
            }

            // scatter
            for (int i = 0; i < n; i++) {
                long key = keys[i];
                int bucket = (int) ((key >>> shift) & 0xFFL);
                int pos = count[bucket]++;
                keyBuffer[pos] = key;
                indexBuffer[pos] = indices[i];
            }

            // copy back
            System.arraycopy(keyBuffer, 0, keys, 0, n);
            System.arraycopy(indexBuffer, 0, indices, 0, n);
        }
    }

    /**
     * Map input values into ranks {0,...,distinct-1} in stable order.
     * Equal values get same rank.
     */
    public static int[] rankTransform(int[] values) {
        if (values == null) {
            return null;
        }
        if (values.length == 0) {
            return new int[0];
        }

        int[] sorted = sortedCopy(values);

        Map<Integer, Integer> ranks = new HashMap<>();
        int nextRank = 0;
        for (int v : sorted) {
            if (!ranks.containsKey(v)) {
                ranks.put(v, nextRank++);
            }
        }

        int[] mapped = new int[values.length];
        for (int i = 0; i < values.length; i++) {
            mapped[i] = ranks.get(values[i]);
        }
        return mapped;
    }
}
