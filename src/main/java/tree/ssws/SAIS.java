package tree.ssws;

import java.util.Arrays;

/**
 * SA-IS (Nong–Zhang–Chan) suffix array construction in O(n) time.
 *
 * Input: int[] s of length n with s[n-1] == 0 (unique smallest sentinel),
 *        and 0 <= s[i] <= K for some K (K can be >= n; still linear on a word RAM).
 * Output: suffix array sa of length n.
 */
final class SAIS {

    private SAIS() {}

    static int[] buildSuffixArray(int[] s) {
        int n = s.length;
        if (n == 0) return new int[0];
        if (n == 1) return new int[]{0};

        int K = 0;
        for (int v : s) if (v > K) K = v;

        int[] sa = new int[n];
        induceSA(s, sa, K);
        return sa;
    }

    // --- SA-IS core ---

    private static void induceSA(int[] s, int[] sa, int K) {
        final int n = s.length;
        final boolean[] isS = new boolean[n];
        final boolean[] isLMS = new boolean[n];

        // classify S/L
        isS[n - 1] = true; // sentinel is S-type
        for (int i = n - 2; i >= 0; --i) {
            if (s[i] < s[i + 1]) isS[i] = true;
            else if (s[i] > s[i + 1]) isS[i] = false;
            else isS[i] = isS[i + 1];
        }
        for (int i = 1; i < n; ++i) {
            isLMS[i] = isS[i] && !isS[i - 1];
        }

        // bucket sizes
        int[] bucketSize = new int[K + 1];
        for (int v : s) bucketSize[v]++;

        // 1) place LMS suffixes using bucket ends
        Arrays.fill(sa, -1);
        int[] bucketEnd = new int[K + 1];
        bucketEnd[0] = bucketSize[0] - 1;
        for (int i = 1; i <= K; i++) bucketEnd[i] = bucketEnd[i - 1] + bucketSize[i];

        for (int i = 0; i < n; i++) {
            if (!isLMS[i]) continue;
            int c = s[i];
            sa[bucketEnd[c]--] = i;
        }

        // 2) induce L-type using bucket starts
        induceL(s, sa, isS, bucketSize);

        // 3) induce S-type using bucket ends
        induceS(s, sa, isS, bucketSize);

        // 4) name LMS substrings
        int m = 0;
        for (boolean b : isLMS) if (b) m++;

        int[] lmsOrder = new int[m];
        {
            int p = 0;
            for (int i = 0; i < n; i++) if (isLMS[i]) lmsOrder[p++] = i;
        }
        // LMS order as they appear in SA
        int[] lmsInSa = new int[m];
        {
            int p = 0;
            for (int pos : sa) if (pos >= 0 && isLMS[pos]) lmsInSa[p++] = pos;
        }

        int[] lmsName = new int[n];
        Arrays.fill(lmsName, -1);

        int name = 0;
        lmsName[lmsInSa[0]] = name;
        for (int i = 1; i < m; i++) {
            int a = lmsInSa[i - 1];
            int b = lmsInSa[i];
            if (!lmsSubstrEqual(s, isLMS, a, b)) name++;
            lmsName[b] = name;
        }

        // 5) build reduced string s1 of LMS names in LMS order (not SA order!)
        int[] s1 = new int[m];
        int idx = 0;
        for (int pos : lmsOrder) s1[idx++] = lmsName[pos];

        int[] sa1;
        if (name + 1 == m) {
            // already unique → directly place ranks
            sa1 = new int[m];
            for (int i = 0; i < m; i++) sa1[s1[i]] = i;
        } else {
            sa1 = new int[m];
            induceSA(s1, sa1, name);
        }

        // 6) reorder LMS positions according to sa1
        int[] orderedLMS = new int[m];
        for (int i = 0; i < m; i++) orderedLMS[i] = lmsOrder[sa1[i]];

        // 7) place LMS suffixes again using bucket ends in this precise order, then induce
        Arrays.fill(sa, -1);
        Arrays.fill(bucketEnd, 0);
        bucketEnd[0] = bucketSize[0] - 1;
        for (int i = 1; i <= K; i++) bucketEnd[i] = bucketEnd[i - 1] + bucketSize[i];

        for (int i = m - 1; i >= 0; i--) {
            int p = orderedLMS[i];
            int c = s[p];
            sa[bucketEnd[c]--] = p;
        }

        induceL(s, sa, isS, bucketSize);
        induceS(s, sa, isS, bucketSize);
    }

    private static void induceL(int[] s, int[] sa, boolean[] isS, int[] bucketSize) {
        final int n = s.length;
        final int K = bucketSize.length - 1;

        int[] bucketStart = new int[K + 1];
        bucketStart[0] = 0;
        for (int i = 1; i <= K; i++) bucketStart[i] = bucketStart[i - 1] + bucketSize[i - 1];

        for (int i = 0; i < n; i++) {
            int j = sa[i] - 1;
            if (j >= 0 && !isS[j]) {
                int c = s[j];
                sa[bucketStart[c]++] = j;
            }
        }
    }

    private static void induceS(int[] s, int[] sa, boolean[] isS, int[] bucketSize) {
        final int n = s.length;
        final int K = bucketSize.length - 1;

        int[] bucketEnd = new int[K + 1];
        bucketEnd[0] = bucketSize[0] - 1;
        for (int i = 1; i <= K; i++) bucketEnd[i] = bucketEnd[i - 1] + bucketSize[i];

        for (int i = n - 1; i >= 0; i--) {
            int j = sa[i] - 1;
            if (j >= 0 && isS[j]) {
                int c = s[j];
                sa[bucketEnd[c]--] = j;
            }
        }
    }

    private static boolean lmsSubstrEqual(int[] s, boolean[] isLMS, int a, int b) {
        if (a == b) return true;
        int n = s.length;

        int i = a, j = b;
        while (true) {
            boolean aLMS = isLMS[i];
            boolean bLMS = isLMS[j];
            if (s[i] != s[j] || aLMS != bLMS) return false;
            if (aLMS && bLMS && i != a && j != b) return true; // both reached next LMS
            i++; j++;
            // safety: sentinel must be reached together if ever
            if (i >= n || j >= n) return i == j;
        }
    }
}
