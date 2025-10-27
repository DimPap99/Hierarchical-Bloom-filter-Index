package tree.ssws;

import java.util.Arrays;

/**
 * Utility class providing linear-time construction of suffix and LCP arrays over integer
 * alphabets. The suffix array is built with the DC3 / Skew algorithm and therefore performs a
 * constant number of stable radix passes over O(n) sized arrays. The LCP array is produced with
 * Kasai's algorithm.
 */
final class SuffixArrayBuilder {

    private SuffixArrayBuilder() {
    }

    static int[] buildSuffixArray(int[] text) {
        if (text == null) {
            throw new IllegalArgumentException("text must be non-null");
        }
        int n = text.length;
        if (n == 0) {
            return new int[0];
        }
        int[] s = Arrays.copyOf(text, n + 3);
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
        skew(s, sa, n, max);
        return sa;
    }

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

        int n0 = (n + 2) / 3;
        int n1 = (n + 1) / 3;
        int n2 = n / 3;
        int n12 = n1 + n2;

        int[] sample = new int[n12 + 3];
        int[] sa12 = new int[n12 + 3];

        int sampleIndex = 0;
        for (int i = 0; i < n + (n0 - n1); i++) {
            if (i % 3 != 0) {
                sample[sampleIndex++] = i;
            }
        }

        radixPass(sample, sa12, s, 2, n12, alphabetSize);
        radixPass(sa12, sample, s, 1, n12, alphabetSize);
        radixPass(sample, sa12, s, 0, n12, alphabetSize);

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

        if (name < n12) {
            int[] recursiveSa = new int[n12];
            skew(Arrays.copyOf(sample, n12 + 3), recursiveSa, n12, name);
            for (int i = 0; i < n12; i++) {
                if (recursiveSa[i] < n1) {
                    sa12[i] = recursiveSa[i] * 3 + 1;
                } else {
                    sa12[i] = (recursiveSa[i] - n1) * 3 + 2;
                }
            }
        }

        int[] rank12 = new int[n + 3];
        for (int i = 0; i < n12; i++) {
            rank12[sa12[i]] = i + 1;
        }

        int[] mod0 = new int[n0];
        for (int i = 0, idx = 0; i < n; i += 3) {
            mod0[idx++] = i;
        }
        int[] buffer = new int[n0];
        countingSort(mod0, buffer, n0, rank12, 1, n12);
        countingSort(buffer, mod0, n0, s, 0, alphabetSize);

        int p = 0;
        int t = 0;
        int k = 0;
        while (p < n0 && t < n12) {
            int i = sa12[t];
            int j = mod0[p];
            if (suffixLessOrEqual(i, j, s, rank12)) {
                sa[k++] = i;
                t++;
            } else {
                sa[k++] = j;
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

    private static boolean suffixLessOrEqual(int left, int right, int[] s, int[] rank12) {
        if (left % 3 == 1) {
            if (s[left] != s[right]) {
                return s[left] < s[right];
            }
            return rank12[left + 1] <= rank12[right + 1];
        }
        if (s[left] != s[right]) {
            return s[left] < s[right];
        }
        if (s[left + 1] != s[right + 1]) {
            return s[left + 1] < s[right + 1];
        }
        return rank12[left + 2] <= rank12[right + 2];
    }

    private static void radixPass(int[] source, int[] dest, int[] s, int offset, int length, int alphabetSize) {
        if (length == 0) {
            return;
        }
        int[] count = new int[alphabetSize + 2];
        for (int i = 0; i < length; i++) {
            int value = s[source[i] + offset];
            count[value + 1]++;
        }
        for (int i = 1; i < count.length; i++) {
            count[i] += count[i - 1];
        }
        for (int i = 0; i < length; i++) {
            int index = source[i];
            int value = s[index + offset];
            dest[count[value]++] = index;
        }
    }

    private static void countingSort(int[] source, int[] dest, int length, int[] key, int offset, int maxKey) {
        if (length == 0) {
            return;
        }
        int[] count = new int[maxKey + 2];
        for (int i = 0; i < length; i++) {
            int value = key[source[i] + offset];
            count[value + 1]++;
        }
        for (int i = 1; i < count.length; i++) {
            count[i] += count[i - 1];
        }
        for (int i = 0; i < length; i++) {
            int index = source[i];
            int value = key[index + offset];
            dest[count[value]++] = index;
        }
    }

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
