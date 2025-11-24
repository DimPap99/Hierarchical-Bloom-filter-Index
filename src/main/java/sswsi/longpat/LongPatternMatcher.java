package sswsi.longpat;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

// Matcher for long streaming patterns using the three-phase logic from the paper.
public final class LongPatternMatcher {

    private LongPatternMatcher() {}

    public static List<Integer> matchLongPattern(int[] P, int[] S, int windowStartGlobal) {
        if (P == null || S == null || P.length == 0 || S.length == 0) {
            return Collections.emptyList();
        }
        int m = P.length;

        int q = Math.max(1, m / 4); // floor(m/4)

        // Phase 1: prefix KMP and period detection on P[0..q-1]
        int[] prefix = new int[q];
        System.arraycopy(P, 0, prefix, 0, q);
        int[] pi = failure(prefix);
        int period = q - pi[q - 1];
        boolean periodic = (q % period == 0) && (period <= q / 2);

        // Phase 2: KMP over S for prefix occurrences; build explicit list or chains
        List<Integer> starts = kmpAll(prefix, S);
        // In our windowed setting |S| = O(m), so starts is O(1) unless periodic

        ArrayList<Integer> result = new ArrayList<>();
        if (!periodic) {
            // Non-periodic: explicitly filter to full matches
            for (int s : starts) {
                if (s + m <= S.length && equalTail(P, S, s, q)) {
                    result.add(windowStartGlobal + s);
                }
            }
            return result;
        }

        // Periodic: compress into chains where differences equal period
        List<int[]> chains = new ArrayList<>(); // each chain as [a1, ak]
        int i = 0;
        while (i < starts.size()) {
            int a1 = starts.get(i);
            int ak = a1;
            int j = i + 1;
            while (j < starts.size() && starts.get(j) - ak == period) {
                ak = starts.get(j);
                j++;
            }
            chains.add(new int[]{a1, ak});
            i = j;
        }

        // Phase 3: extend occurrences or chains by checking remaining P[q..m-1]
        for (int[] c : chains) {
            int a1 = c[0], ak = c[1];
            // Check from ak backwards along the chain if needed
            // For each candidate a in [a1..ak step period], verify the tail
            for (int a = a1; a <= ak; a += period) {
                if (a + m <= S.length && equalTail(P, S, a, q)) {
                    result.add(windowStartGlobal + a);
                }
            }
        }
        return result;
    }

    private static boolean equalTail(int[] P, int[] S, int start, int from) {
        for (int i = from; i < P.length; i++) {
            if (S[start + i] != P[i]) return false;
        }
        return true;
    }

    private static List<Integer> kmpAll(int[] pattern, int[] text) {
        ArrayList<Integer> res = new ArrayList<>();
        if (pattern.length == 0) return res;
        int[] pi = failure(pattern);
        int j = 0;
        for (int i = 0; i < text.length; i++) {
            while (j > 0 && text[i] != pattern[j]) j = pi[j - 1];
            if (text[i] == pattern[j]) {
                j++;
                if (j == pattern.length) {
                    res.add(i - pattern.length + 1);
                    j = pi[j - 1];
                }
            }
        }
        return res;
    }

    private static int[] failure(int[] pattern) {
        int[] pi = new int[pattern.length];
        int j = 0;
        for (int i = 1; i < pattern.length; i++) {
            while (j > 0 && pattern[i] != pattern[j]) j = pi[j - 1];
            if (pattern[i] == pattern[j]) j++;
            pi[i] = j;
        }
        return pi;
    }
}
