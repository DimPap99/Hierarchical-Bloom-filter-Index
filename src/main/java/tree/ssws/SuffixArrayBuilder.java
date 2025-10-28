package tree.ssws;

import java.util.Arrays;

/**
 * SuffixArrayBuilder
 *
 * Provides:
 *   - buildSuffixArray(int[] text): DC3 / Skew linear time suffix array construction
 *   - buildLcpArray(int[] text, int[] sa): Kasai's Longest Common Prefix in linear time
 *
 * Performance improvements:
 *   We reuse a per-thread Dc3Workspace to avoid allocating fresh large
 *   scratch arrays (sample, sa12, rank12, mod0, buffer, count) on every call.
 *   This significantly reduces garbage collector pressure if you build
 *   many suffix arrays repeatedly (for example in a sliding window index).
 */
final class SuffixArrayBuilder {

    private SuffixArrayBuilder() {
    }

    /**
     * Small reusable workspace for DC3 / skew.
     * We store and reuse scratch arrays here instead of reallocating them every build.
     */
    private static final class Dc3Workspace {
        int[] sample;
        int[] sa12;
        int[] rank12;
        int[] mod0;
        int[] buffer;
        int[] count;

        Dc3Workspace(int cap, int countCap) {
            sample = new int[cap];
            sa12 = new int[cap];
            rank12 = new int[cap];
            mod0 = new int[cap];
            buffer = new int[cap];
            count = new int[countCap];
        }

        void ensureCapacity(int n, int alphabetSize) {
            // n is the full text length
            // n12 ~ 2/3 n in DC3, n0 ~ 1/3 n. We just allocate ~n for all.
            int cap = n + 5;
            if (sample.length < cap) {
                sample = Arrays.copyOf(sample, cap);
                sa12   = Arrays.copyOf(sa12,   cap);
                rank12 = Arrays.copyOf(rank12, cap);
                mod0   = Arrays.copyOf(mod0,   cap);
                buffer = Arrays.copyOf(buffer, cap);
            }
            int neededCount = Math.max(Math.max(alphabetSize + 2, n + 2), 16);
            if (count.length < neededCount) {
                count = Arrays.copyOf(count, neededCount);
            }
        }
    }

    /**
     * One workspace per thread. This saves us from frequent allocation
     * in streaming / repeated builds, but does not require changing the
     * public method signatures.
     */
    private static final ThreadLocal<Dc3Workspace> TLS_WORKSPACE =
            ThreadLocal.withInitial(() -> new Dc3Workspace(16, 32));

    /**
     * Build the suffix array of "text" using DC3 / Skew.
     * "text" must be a non-negative integer alphabet and must already
     * include a globally smallest sentinel at the end.
     */
    static int[] buildSuffixArray(int[] text) {
        if (text == null) {
            throw new IllegalArgumentException("text must be non-null");
        }
        int n = text.length;
        if (n == 0) {
            return new int[0];
        }

        // DC3 expects padding up to 3 zeros so we can read s[i+1], s[i+2] without bounds checks.
        int[] s = Arrays.copyOf(text, n + 3);

        // Confirm non-negative and track max symbol for counting sort bounds.
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
        int alphabetSize = Math.max(max, n);

        Dc3Workspace ws = TLS_WORKSPACE.get();
        ws.ensureCapacity(n, alphabetSize);

        skewWithWorkspace(s, sa, n, alphabetSize, ws);

        return sa;
    }

    /**
     * A lightly modified DC3 / Skew that uses (and reuses) a Dc3Workspace.
     * The core algorithm is the same:
     *
     *   1. Sort mod-1 and mod-2 suffixes by 3-character tuple
     *   2. Assign ranks (names)
     *   3. Recurse if ranks are not unique
     *   4. Sort mod-0 suffixes by (first char, rank of following suffix)
     *   5. Merge
     */
    private static void skewWithWorkspace(int[] s,
                                          int[] sa,
                                          int n,
                                          int alphabetSize,
                                          Dc3Workspace ws) {

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

        // Partition suffixes by index mod 3.
        final int n0 = (n + 2) / 3;
        final int n1 = (n + 1) / 3;
        final int n2 = (n) / 3;
        final int n12 = n1 + n2;

        int[] sample = ws.sample;  // length >= n + 5
        int[] sa12   = ws.sa12;
        int[] rank12 = ws.rank12;
        int[] mod0   = ws.mod0;
        int[] buffer = ws.buffer;
        int[] count  = ws.count;

        // Build list of indices i where i % 3 != 0 into 'sample'.
        int sampleIndex = 0;
        for (int i = 0; i < n + (n0 - n1); i++) {
            if (i % 3 != 0) {
                sample[sampleIndex++] = i;
            }
        }

        // Radix sort the triples (s[i], s[i+1], s[i+2]) for i % 3 != 0.
        radixPass(sample, sa12, s, 2, n12, count);
        radixPass(sa12, sample, s, 1, n12, count);
        radixPass(sample, sa12, s, 0, n12, count);

        // Assign names (ranks) to those triples
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

        // If names are not unique, recurse on the "sample" array of names.
        if (name < n12) {
            int[] recursiveSa = new int[n12];
            int[] sampleCopy = Arrays.copyOf(sample, n12 + 3);
            skewWithWorkspace(sampleCopy, recursiveSa, n12, name, ws);

            for (int i = 0; i < n12; i++) {
                if (recursiveSa[i] < n1) {
                    sa12[i] = recursiveSa[i] * 3 + 1;
                } else {
                    sa12[i] = (recursiveSa[i] - n1) * 3 + 2;
                }
            }
        }

        // rank12[pos] = rank (1-based) among mod-1/mod-2 suffixes
        Arrays.fill(rank12, 0, n + 3, 0);
        for (int i = 0; i < n12; i++) {
            rank12[sa12[i]] = i + 1;
        }

        // Sort the mod-0 suffixes (positions 0,3,6,...) into mod0.
        int p0 = 0;
        for (int i = 0; i < n; i += 3) {
            mod0[p0++] = i;
        }

        // We need a stable counting sort of mod0 by:
        //   key (rank12[i+1]) then key (s[i])
        // We do it in two passes like the original code but reuse 'count' and 'buffer'
        countingSort(mod0, buffer, n0, rank12, 1, n12, count);
        countingSort(buffer, mod0, n0, s, 0, alphabetSize, count);

        // Merge the two sorted lists (sa12 for mod-1/mod-2, mod0 for mod-0)
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
     * suffixLessOrEqual:
     * Compare two suffixes in O(1) using precomputed ranks.
     *
     * If left % 3 == 1:
     *   compare (s[left], rank12[left+1])
     * else (left % 3 == 2):
     *   compare (s[left], s[left+1], rank12[left+2])
     */
    private static boolean suffixLessOrEqual(int left,
                                             int right,
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
     * radixPass:
     * Stable counting sort of 'source' into 'dest' using key s[index + offset].
     * We reuse the provided "count" array to avoid new allocations.
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

        // Find maximum key to know how many buckets we actually need.
        int maxValue = 0;
        for (int i = 0; i < length; i++) {
            int value = s[source[i] + offset];
            if (value > maxValue) {
                maxValue = value;
            }
        }

        // Zero only the needed prefix of count.
        Arrays.fill(count, 0, maxValue + 2, 0);

        // Frequency
        for (int i = 0; i < length; i++) {
            int value = s[source[i] + offset];
            count[value + 1]++;
        }
        // Prefix sum
        for (int i = 1; i < maxValue + 2; i++) {
            count[i] += count[i - 1];
        }
        // Scatter
        for (int i = 0; i < length; i++) {
            int idx = source[i];
            int value = s[idx + offset];
            dest[count[value]++] = idx;
        }
    }

    /**
     * countingSort:
     * Stable counting sort for the mod-0 suffixes array.
     * Sorts "source" by key[index + offset], writes result to "dest".
     * Uses/reuses the shared "count" array.
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
     * Kasai's algorithm in linear time.
     * Longest Common Prefix[i] = LCP between suffix sa[i] and suffix sa[i+1].
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
