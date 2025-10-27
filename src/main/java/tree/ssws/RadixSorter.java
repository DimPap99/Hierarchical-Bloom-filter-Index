package tree.ssws;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Utility class providing an LSD radix sort for integer alphabets and a helper to map
 * arbitrary integer sequences into a compact rank space. The implementation is stable and
 * allocates only two temporary arrays regardless of the number of passes.
 */
public final class RadixSorter {

    private static final int RADIX = 256; // process eight bits per pass
    private static final int BITS_PER_BYTE = 8;

    private RadixSorter() {
        // utility class
    }

    /**
     * Sorts {@code values} in-place using a stable LSD radix sort. The input may contain
     * negative values; internally we offset the values so counting buckets remain non-negative.
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
        // determine number of passes (at least one to handle zero correctly)
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
     * Returns a sorted copy of {@code values} using the same radix sort.
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
     * Stable LSD radix sort for parallel arrays. Reorders {@code indices} so they are sorted by the
     * associated {@code keys} interpreted as unsigned 64-bit integers. Both arrays are modified in
     * place.
     */
    public static void sortParallel(long[] keys, int[] indices) {
        if (keys == null || indices == null) {
            throw new IllegalArgumentException("keys and indices must be non-null");
        }
        if (keys.length != indices.length) {
            throw new IllegalArgumentException("keys and indices must have the same length");
        }
        int n = keys.length;
        if (n <= 1) {
            return;
        }
        long[] keyBuffer = new long[n];
        int[] indexBuffer = new int[n];
        int[] count = new int[RADIX];

        for (int shift = 0; shift < Long.SIZE; shift += BITS_PER_BYTE) {
            Arrays.fill(count, 0);
            for (long key : keys) {
                int bucket = (int) ((key >>> shift) & 0xFFL);
                count[bucket]++;
            }
            int total = 0;
            for (int i = 0; i < RADIX; i++) {
                int c = count[i];
                count[i] = total;
                total += c;
            }
            for (int i = 0; i < n; i++) {
                long key = keys[i];
                int bucket = (int) ((key >>> shift) & 0xFFL);
                int pos = count[bucket]++;
                keyBuffer[pos] = key;
                indexBuffer[pos] = indices[i];
            }
            System.arraycopy(keyBuffer, 0, keys, 0, n);
            System.arraycopy(indexBuffer, 0, indices, 0, n);
        }
    }

    /**
     * Maps {@code values} into rank space {0, 1, ..., distinct-1}. Relative ordering is preserved
     * and equal values receive identical ranks. Radix sort is used internally to avoid the need for
     * comparison-based sorting on potentially large alphabets.
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
