package tree.ssws;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.junit.jupiter.api.Test;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class StaticSuffixTreeTest {

    private static final SuffixTreeBuilder DC3 = new Dc3SuffixTreeBuilder();

    // Small texts: every substring as a pattern. Random texts over 1, 2, 4, 26 and n symbols, periodic texts, a^n.
    @Test
    void everySubstringOfSmallTexts() {
        Random rnd = new Random(41);
        for (long[] text : texts(rnd, new int[]{0, 1, 2, 3, 4, 5, 7, 8, 13, 31, 32, 33, 40})) {
            Trees trees = new Trees(text, rnd);
            for (int i = 0; i < text.length; i++) {
                for (int j = i + 1; j <= text.length; j++) {
                    trees.check(Arrays.copyOfRange(text, i, j), rnd);
                }
            }
            trees.check(new long[]{-7L}, rnd);   // a symbol that does not occur
        }
    }

    // Larger texts: random present, mutated and absent patterns.
    @Test
    void randomPatternsOnLargerTexts() {
        Random rnd = new Random(43);
        for (long[] text : texts(rnd, new int[]{64, 100, 257, 1000, 4096, 5000})) {
            Trees trees = new Trees(text, rnd);
            for (int q = 0; q < 300; q++) {
                int m = 1 + rnd.nextInt(Math.min(text.length, (q % 3 == 0) ? text.length : 24));
                int s = rnd.nextInt(text.length - m + 1);
                long[] present = Arrays.copyOfRange(text, s, s + m);
                trees.check(present, rnd);
                long[] mutated = present.clone();
                mutated[rnd.nextInt(m)] = text[rnd.nextInt(text.length)];
                trees.check(mutated, rnd);
            }
        }
    }

    // A run of one symbol: the tree is as deep as it is long. Nothing may recurse on depth.
    @Test
    void longRunOfOneSymbol() {
        int n = 200_000;
        long[] text = new long[n];
        Arrays.fill(text, 5L);
        for (ChildDictionary.Kind kind : ChildDictionary.Kind.values()) {
            StaticSuffixTree tree = DC3.build(text, 0, n, kind);
            assertEquals(n, tree.internalNodeCount());     // one node per depth 0 .. n-1
            for (int m : new int[]{1, 2, 8, 1000, n - 1, n}) {
                long[] pattern = new long[m];
                Arrays.fill(pattern, 5L);
                int node = tree.locus(pattern);
                assertEquals(n - m + 1, tree.hi(node) - tree.lo(node) + 1, "a^" + m);
            }
            assertEquals(StaticSuffixTree.NONE, tree.locus(new long[n + 1]));
        }
    }

    // The tree reads text[offset, offset + length) and nothing outside it.
    @Test
    void offsetViewMatchesCopy() {
        Random rnd = new Random(47);
        long[] buffer = new long[700];
        for (int i = 0; i < buffer.length; i++) {
            buffer[i] = 1 + rnd.nextInt(3);
        }
        int offset = 123;
        int length = 400;
        StaticSuffixTree view = DC3.build(buffer, offset, length, ChildDictionary.Kind.FASTUTIL);
        long[] copy = Arrays.copyOfRange(buffer, offset, offset + length);
        StaticSuffixTree own = DC3.build(copy, 0, length, ChildDictionary.Kind.FASTUTIL);
        assertEquals(own.internalNodeCount(), view.internalNodeCount());
        for (int i = 0; i <= length; i++) {
            assertEquals(own.suffix(i), view.suffix(i));
        }
        for (int q = 0; q < 200; q++) {
            int m = 1 + rnd.nextInt(30);
            int s = rnd.nextInt(length - m + 1);
            long[] p = Arrays.copyOfRange(copy, s, s + m);
            assertEquals(own.locus(p), view.locus(p));
            // A pattern that would continue past the view must not match across its end.
            long[] tail = Arrays.copyOfRange(buffer, offset + length - m + 1, offset + length + 1);
            assertEquals(occurrences(copy, tail), occurrencesOf(view, view.locus(tail)));
        }
    }

    // Range maximum against a scan, including entries at or above the limit, which count as minus infinity.
    @Test
    void rangeMaxMatchesScan() {
        Random rnd = new Random(53);
        for (int n : new int[]{1, 2, 31, 32, 33, 63, 64, 65, 100, 1000, 4097, 10_000}) {
            int[] a = new int[n];
            for (int i = 0; i < n; i++) {
                a[i] = rnd.nextInt(2 * n);
            }
            for (int limit : new int[]{Integer.MAX_VALUE, n, n / 3}) {
                RangeMax rmq = new RangeMax(a, limit);
                for (int q = 0; q < 2000; q++) {
                    int l = rnd.nextInt(n);
                    int r = l + rnd.nextInt(n - l);
                    int best = -1;
                    for (int i = l; i <= r; i++) {
                        best = Math.max(best, (a[i] < limit) ? a[i] : -1);
                    }
                    int j = rmq.argMax(l, r);
                    assertTrue(j >= l && j <= r);
                    assertEquals(best, (a[j] < limit) ? a[j] : -1, "n=" + n + " [" + l + "," + r + "]");

                    int threshold = rnd.nextInt(2 * n + 2) - 1;
                    IntArrayList expected = new IntArrayList();
                    for (int i = l; i <= r; i++) {
                        if (a[i] >= Math.max(0, threshold) && a[i] < limit) {
                            expected.add(1000 + a[i]);
                        }
                    }
                    IntArrayList out = new IntArrayList();
                    rmq.reportAtLeast(l, r, threshold, out, 1000);
                    expected.sort(null);
                    out.sort(null);
                    assertEquals(expected, out, "every qualifying entry exactly once");
                }
            }
        }
    }

    /** The same text built with both dictionaries. */
    private static final class Trees {
        final long[] text;
        final StaticSuffixTree fastutil;
        final StaticSuffixTree fks;

        Trees(long[] text, Random rnd) {
            this.text = text;
            this.fastutil = DC3.build(text, 0, text.length, ChildDictionary.Kind.FASTUTIL);
            this.fks = DC3.build(text, 0, text.length, ChildDictionary.Kind.FKS);
            checkStructure(text, fastutil);
            checkStructure(text, fks);
            assertEquals(fastutil.internalNodeCount(), fks.internalNodeCount());
        }

        void check(long[] pattern, Random rnd) {
            TreeSet<Integer> expected = occurrences(text, pattern);
            int node = fastutil.locus(pattern);
            assertEquals(node, fks.locus(pattern), "FASTUTIL and FKS disagree");
            assertEquals(expected, occurrencesOf(fastutil, node), "locus range");
            if (node == StaticSuffixTree.NONE) {
                return;
            }
            int lo = fastutil.lo(node);
            int hi = fastutil.hi(node);
            assertTrue(lo >= 1, "the terminator's suffix never matches");

            // Threshold reporting over the suffix array (leftmost segment).
            RangeMax starts = new RangeMax(fastutil);
            for (int threshold : new int[]{0, 1, rnd.nextInt(text.length + 1), text.length / 2, text.length}) {
                IntArrayList out = new IntArrayList();
                starts.reportAtLeast(lo, hi, threshold, out, 0);
                assertEquals(expected.tailSet(threshold), sortedNoDuplicates(out));
            }

            // BL' reporting (boundary trees): starts in [threshold, half).
            int half = rnd.nextInt(text.length + 1);
            RangeMax leftStarts = new RangeMax(fastutil, half);
            for (int threshold : new int[]{0, half - pattern.length + 1, rnd.nextInt(half + 1)}) {
                IntArrayList out = new IntArrayList();
                leftStarts.reportAtLeast(lo, hi, threshold, out, -3);
                TreeSet<Integer> want = new TreeSet<>();
                for (int s : expected.subSet(Math.max(0, threshold), half)) {
                    want.add(s - 3);
                }
                assertEquals(want, sortedNoDuplicates(out));
            }
        }
    }

    private static TreeSet<Integer> sortedNoDuplicates(IntArrayList out) {
        TreeSet<Integer> set = new TreeSet<>(out);
        assertEquals(out.size(), set.size(), "reported twice");
        return set;
    }

    private static TreeSet<Integer> occurrencesOf(StaticSuffixTree tree, int node) {
        TreeSet<Integer> set = new TreeSet<>();
        if (node != StaticSuffixTree.NONE) {
            IntArrayList out = new IntArrayList();
            tree.reportAll(node, out, 0);
            set.addAll(out);
            assertEquals(out.size(), set.size());
        }
        return set;
    }

    private static TreeSet<Integer> occurrences(long[] text, long[] pattern) {
        TreeSet<Integer> hits = new TreeSet<>();
        outer:
        for (int i = 0; i + pattern.length <= text.length; i++) {
            for (int k = 0; k < pattern.length; k++) {
                if (text[i + k] != pattern[k]) {
                    continue outer;
                }
            }
            hits.add(i);
        }
        return hits;
    }

    /**
     * n+1 leaves visited in suffix-array order by a depth-first traversal; every internal node has at least two
     * children (the root of the empty text has one) with strictly increasing first-symbol ranks and contiguous
     * leaf ranges that partition the parent's; depths grow downwards; an internal node's depth is the LCP of the
     * first and last suffix below it; a leaf's depth is its suffix length with the terminator.
     */
    static void checkStructure(long[] text, StaticSuffixTree tree) {
        int n = text.length;
        int[] ranks = RankReducer.rankWithTerminator(text);
        int[] sa = bruteForceSuffixArray(ranks);
        assertEquals(n + 1, tree.leafCount());
        int[] leafOrder = new int[n + 1];
        int leaves = 0;
        int internal = 0;
        ArrayDeque<Integer> stack = new ArrayDeque<>();
        stack.push(0);
        while (!stack.isEmpty()) {
            int v = stack.pop();
            if (StaticSuffixTree.isLeaf(v)) {
                assertEquals(leaves, tree.lo(v));
                assertEquals(n + 1 - tree.suffix(~v), tree.depth(v));
                leafOrder[leaves++] = tree.suffix(~v);
                continue;
            }
            internal++;
            int k = tree.childCount(v);
            assertTrue(k >= 2 || (v == 0 && n == 0), "unary internal node");
            if (v != 0) {
                assertEquals(lcp(ranks, tree.suffix(tree.lo(v)), tree.suffix(tree.hi(v))), tree.depth(v));
            } else {
                assertEquals(0, tree.depth(v));
            }
            int nextLo = tree.lo(v);
            int prevRank = -1;
            for (int i = 0; i < k; i++) {
                int c = tree.child(v, i);
                assertEquals(nextLo, tree.lo(c), "children's leaf ranges must be contiguous");
                assertTrue(tree.hi(c) >= tree.lo(c));
                nextLo = tree.hi(c) + 1;
                assertTrue(tree.depth(c) > tree.depth(v));
                int rank = ranks[tree.suffix(tree.lo(c)) + tree.depth(v)];
                assertTrue(rank > prevRank, "first symbols must be distinct and sorted");
                prevRank = rank;
                if (rank > 0) {
                    assertEquals(c, tree.childBySymbol(v, tree.symbolAt(tree.suffix(tree.lo(c)) + tree.depth(v))));
                }
            }
            assertEquals(tree.hi(v) + 1, nextLo, "children's leaf ranges must cover the parent's");
            for (int i = k - 1; i >= 0; i--) {
                stack.push(tree.child(v, i));
            }
        }
        assertEquals(n + 1, leaves);
        assertEquals(tree.internalNodeCount(), internal);
        assertArrayEquals(sa, leafOrder, "depth-first leaf order must be the suffix array");
    }

    private static int lcp(int[] s, int a, int b) {
        int h = 0;
        while (a + h < s.length && b + h < s.length && s[a + h] == s[b + h]) {
            h++;
        }
        return h;
    }

    private static int[] bruteForceSuffixArray(int[] s) {
        Integer[] order = new Integer[s.length];
        for (int i = 0; i < s.length; i++) {
            order[i] = i;
        }
        Arrays.sort(order, (a, b) -> {
            while (a < s.length && b < s.length) {
                if (s[a] != s[b]) {
                    return Integer.compare(s[a], s[b]);
                }
                a++;
                b++;
            }
            return Integer.compare(s.length - a, s.length - b);
        });
        int[] sa = new int[s.length];
        for (int i = 0; i < s.length; i++) {
            sa[i] = order[i];
        }
        return sa;
    }

    /** For each length: random texts over 1, 2, 4, 26 and n symbols, then periodic texts and a^n. */
    private static List<long[]> texts(Random rnd, int[] lengths) {
        List<long[]> out = new ArrayList<>();
        for (int n : lengths) {
            for (int sigma : new int[]{1, 2, 4, 26, Math.max(1, n)}) {
                long[] alphabet = new long[sigma];
                for (int i = 0; i < sigma; i++) {
                    alphabet[i] = rnd.nextLong() & Long.MAX_VALUE;   // like TokenHasher output
                }
                long[] t = new long[n];
                for (int i = 0; i < n; i++) {
                    t[i] = alphabet[(sigma == n) ? i : rnd.nextInt(sigma)];
                }
                out.add(t);
            }
            out.add(periodic(n, 1, 2));
            out.add(periodic(n, 1, 2, 3));
            out.add(periodic(n, 1, 1, 2));
            out.add(fibonacci(n));
            out.add(thueMorse(n));
        }
        return out;
    }

    private static long[] periodic(int n, long... unit) {
        long[] t = new long[n];
        for (int i = 0; i < n; i++) {
            t[i] = unit[i % unit.length];
        }
        return t;
    }

    private static long[] fibonacci(int n) {
        StringBuilder a = new StringBuilder("a");
        StringBuilder b = new StringBuilder("ab");
        while (b.length() < n) {
            StringBuilder next = new StringBuilder(b).append(a);
            a = b;
            b = next;
        }
        long[] t = new long[n];
        for (int i = 0; i < n; i++) {
            t[i] = b.charAt(i);
        }
        return t;
    }

    private static long[] thueMorse(int n) {
        long[] t = new long[n];
        for (int i = 0; i < n; i++) {
            t[i] = Integer.bitCount(i) & 1;   // includes the symbol 0
        }
        return t;
    }
}
