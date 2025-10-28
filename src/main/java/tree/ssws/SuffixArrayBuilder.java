package tree.ssws;

import java.util.Arrays;

/**
 * SuffixArrayBuilder
 *
 * This class builds:
 *
 *   - A suffix array in linear time using the DC3 algorithm
 *     which is also known as Skew. DC3 partitions suffixes
 *     by index modulo three, sorts two of the three classes
 *     with radix sort on triples (s[i], s[i+1], s[i+2]),
 *     recursively ranks them, induces the remaining class,
 *     and merges in linear time.
 *
 *   - A Longest Common Prefix (LCP) array in linear time
 *     using Kasai's algorithm.
 *

 */
final class SuffixArrayBuilder {

    private SuffixArrayBuilder() {
    }

    /**
     * Build the suffix array of "text" using the DC3 / Skew algorithm.
     * "text" must be a non-negative integer alphabet and should already
     * include the sentinel as the lexicographically smallest symbol.
     */
    static int[] buildSuffixArray(int[] text) {
        if (text == null) {
            throw new IllegalArgumentException("text must be non-null");
        }
        int n = text.length;
        if (n == 0) {
            return new int[0];
        }

        // DC3 requires padding with up to 3 zeros so that we can safely
        // read s[i+1] and s[i+2] without range checks. We copy once.
        int[] s = Arrays.copyOf(text, n + 3);

        // The Skew algorithm assumes non-negative integers. We already
        // enforced that during normalization in FarachColtonSuffixTree.build.
        int max = 0;
        for (int i = 0; i < n; i++) {
            int value = s[i];
            if (value < 0) {
                throw new IllegalArgumentException("text must be non-negative");
            }
            if (value > max) {
                max = value;
            }
        }

        int[] sa = new int[n];
        int alphabetSize = Math.max(max, n); // safe upper bound for counting arrays
        skew(s, sa, n, alphabetSize);
        return sa;
    }

    /**
     * The DC3 / Skew algorithm. This code has been tuned to:
     *   - Allocate working arrays once per call (sample, sa12, rank12, mod0, buffer, count)
     *   - Reuse a shared counting array "count" in radixPass and countingSort
     *     instead of allocating a new array each time
     *
     * This significantly reduces heap churn for large inputs, which
     * directly helps StreamingSlidingWindowIndex where we rebuild
     * suffix arrays for entire segments and boundaries.
     */
    private static void skew(int[] s, int[] sa, int n, int alphabetSize) {
        if (n == 1) {
            sa[0] = 0;
            return;
        }
        if (n == 2) {
            if (s[0] <= s[1]) {
                sa[0] = 0;
                sa[1] = 1;
            } else {
                sa[0] = 1;
                sa[1] = 0;
            }
            return;
        }

        // Split suffixes by index modulo 3.
        final int n0 = (n + 2) / 3;
        final int n1 = (n + 1) / 3;
        final int n2 = n / 3;
        final int n12 = n1 + n2;

        // Working arrays for the modulo 1 and modulo 2 suffixes.
        int[] sample = new int[n12 + 3];
        int[] sa12 = new int[n12 + 3];

        // Build the list of indices i where i % 3 != 0.
        int sampleIndex = 0;
        for (int i = 0; i < n + (n0 - n1); i++) {
            if (i % 3 != 0) {
                sample[sampleIndex++] = i;
            }
        }

        // Shared counting array. We size it once to handle any value we will see.
        // We know all keys are <= max(alphabetSize, n12, n) which is O(n),
        // and we pass maxKey bounds to zero-fill only the prefix we need.
        int[] count = new int[Math.max(Math.max(alphabetSize + 2, n12 + 2), n + 5)];

        // Radix sort the triples (s[i], s[i+1], s[i+2]) for i % 3 != 0.
        radixPass(sample, sa12, s, 2, n12, count);
        radixPass(sa12, sample, s, 1, n12, count);
        radixPass(sample, sa12, s, 0, n12, count);

        // Assign names (ranks) to those triples. Equal triples get the same name.
        int name = 0;
        int prev0 = -1;
        int prev1 = -1;
        int prev2 = -1;
        for (int i = 0; i < n12; i++) {
            int pos = sa12[i];
            int c0 = s[pos];
            int c1 = s[pos + 1];
            int c2 = s[pos + 2];

            if (c0 != prev0 || c1 != prev1 || c2 != prev2) {
                name++;
                prev0 = c0;
                prev1 = c1;
                prev2 = c2;
            }

            if (pos % 3 == 1) {
                sample[pos / 3] = name;
            } else {
                sample[pos / 3 + n1] = name;
            }
        }

        // If names are not yet unique, recurse. Otherwise we already have sorted order.
        if (name < n12) {
            int[] recursiveSa = new int[n12];
            // We recurse on the "sample" array of names, padded to length n12+3
            skew(Arrays.copyOf(sample, n12 + 3), recursiveSa, n12, name);

            // Map ranks back to actual indices in the original text.
            for (int i = 0; i < n12; i++) {
                if (recursiveSa[i] < n1) {
                    sa12[i] = recursiveSa[i] * 3 + 1;
                } else {
                    sa12[i] = (recursiveSa[i] - n1) * 3 + 2;
                }
            }
        }

        // rank12[pos] gives the rank (1-based) among the non-mod-0 suffixes.
        int[] rank12 = new int[n + 3];
        for (int i = 0; i < n12; i++) {
            rank12[sa12[i]] = i + 1;
        }

        // Sort the mod-0 suffixes using counting sort driven first by rank12 then by s[i].
        int[] mod0 = new int[n0];
        for (int i = 0, idx = 0; i < n; i += 3) {
            mod0[idx++] = i;
        }
        int[] buffer = new int[n0];
        countingSort(mod0, buffer, n0, rank12, 1, n12, count);
        countingSort(buffer, mod0, n0, s, 0, alphabetSize, count);

        // Merge the two sorted lists (sa12 and mod0) in linear time.
        int p = 0;
        int t = 0;
        int k = 0;
        while (p < n0 && t < n12) {
            int i0 = sa12[t];
            int j0 = mod0[p];
            if (suffixLessOrEqual(i0, j0, s, rank12)) {
                sa[k++] = i0;
                t++;
            } else {
                sa[k++] = j0;
                p++;
            }
        }
        while (p < n0) {
            sa[k++] = mod0[p++];
        }
        while (t < n12) {
            sa[k++] = sa12[t++];
        }
    }

    /**
     * Compare two suffixes in constant time using the precomputed ranks.
     *
     * If left % 3 == 1:
     *   Compare (s[left], rank12[left+1]) versus (s[right], rank12[right+1])
     *
     * Otherwise (left % 3 == 2):
     *   Compare (s[left], s[left+1], rank12[left+2])
     *   versus   (s[right], s[right+1], rank12[right+2])
     *
     * This is the standard DC3 merge comparator and runs in O(1) time.
     */
    private static boolean suffixLessOrEqual(int left, int right,
                                             int[] s,
                                             int[] rank12) {
        if (left % 3 == 1) {
            if (s[left] != s[right]) {
                return s[left] < s[right];
            }
            return rank12[left + 1] <= rank12[right + 1];
        } else {
            if (s[left] != s[right]) {
                return s[left] < s[right];
            }
            if (s[left + 1] != s[right + 1]) {
                return s[left + 1] < s[right + 1];
            }
            return rank12[left + 2] <= rank12[right + 2];
        }
    }

    /**
     * Stable counting sort by the key s[index + offset], reusing the shared "count" array.
     * This replaces the old version which allocated a fresh counting array for each pass.
     */
    private static void radixPass(int[] source,
                                  int[] dest,
                                  int[] s,
                                  int offset,
                                  int length,
                                  int[] count) {
        if (length == 0) {
            return;
        }

        // Find the maximum key value so we know how many buckets are actually used.
        int maxValue = 0;
        for (int i = 0; i < length; i++) {
            int value = s[source[i] + offset];
            if (value > maxValue) {
                maxValue = value;
            }
        }

        // Zero only the prefix we will write.
        Arrays.fill(count, 0, maxValue + 2, 0);

        for (int i = 0; i < length; i++) {
            int value = s[source[i] + offset];
            count[value + 1]++;
        }
        for (int i = 1; i < maxValue + 2; i++) {
            count[i] += count[i - 1];
        }
        for (int i = 0; i < length; i++) {
            int idx = source[i];
            int value = s[idx + offset];
            dest[count[value]++] = idx;
        }
    }

    /**
     * Counting sort for mod-0 suffixes. This is also now reusing the shared
     * "count" array rather than allocating a new one each call.
     *
     * We sort by "key[index + offset]". The "maxKey" argument is an upper bound
     * on that key so that we can bound how much of "count" we must clear.
     */
    private static void countingSort(int[] source,
                                     int[] dest,
                                     int length,
                                     int[] key,
                                     int offset,
                                     int maxKey,
                                     int[] count) {
        if (length == 0) {
            return;
        }

        Arrays.fill(count, 0, maxKey + 2, 0);

        for (int i = 0; i < length; i++) {
            int value = key[source[i] + offset];
            count[value + 1]++;
        }
        for (int i = 1; i < maxKey + 2; i++) {
            count[i] += count[i - 1];
        }
        for (int i = 0; i < length; i++) {
            int index = source[i];
            int value = key[index + offset];
            dest[count[value]++] = index;
        }
    }

    /**
     * Build Kasai's Longest Common Prefix array in O(n) time.
     *
     * Longest Common Prefix[i] is the length of the longest common prefix
     * between suffixes at sa[i] and sa[i+1]. We compute it with Kasai's algorithm,
     * which runs in linear time by remembering the match length h from the
     * previous step and decreasing h by one at each move.
     */
    static int[] buildLcpArray(int[] text, int[] sa) {
        if (text == null || sa == null) {
            throw new IllegalArgumentException("text and sa must be non-null");
        }
        int n = text.length;
        if (n == 0) {
            return new int[0];
        }

        int[] rank = new int[n];
        for (int i = 0; i < n; i++) {
            rank[sa[i]] = i;
        }

        int[] lcp = new int[n - 1];
        int h = 0;
        for (int i = 0; i < n; i++) {
            int r = rank[i];
            if (r == n - 1) {
                h = 0;
                continue;
            }
            int j = sa[r + 1];
            while (i + h < n && j + h < n && text[i + h] == text[j + h]) {
                h++;
            }
            lcp[r] = h;
            if (h > 0) {
                h--;
            }
        }
        return lcp;
    }
}
