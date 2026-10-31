package tree.ssws;

import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

// Farach-Colton must give exactly the suffix array and LCP array of DC3 + Kasai: the depth-first leaf order and the
// internal depths of the tree. Base 2 makes every step run down to tiny inputs; the default base is checked too.
class FarachColtonTest {

    @Test
    void allBinaryStringsUpToLength14() {
        for (int len = 1; len <= 14; len++) {
            int[] t = new int[len];
            for (int mask = 0; mask < (1 << len); mask++) {
                for (int i = 0; i < len; i++) {
                    t[i] = 1 + ((mask >> i) & 1);
                }
                assertSameAsDc3(t, 2, 2);
            }
        }
    }

    @Test
    void allTernaryStringsUpToLength9() {
        for (int len = 1; len <= 9; len++) {
            int total = 1;
            for (int i = 0; i < len; i++) {
                total *= 3;
            }
            int[] t = new int[len];
            for (int code = 0; code < total; code++) {
                for (int i = 0, c = code; i < len; i++, c /= 3) {
                    t[i] = 1 + c % 3;
                }
                assertSameAsDc3(t, 3, 2);
            }
        }
    }

    @Test
    void randomStrings() {
        Random rnd = new Random(61);
        int[] lengths = {1, 2, 3, 5, 31, 32, 33, 34, 63, 64, 65, 100, 255, 1000, 4097, 30_000, 200_000};
        for (int n : lengths) {
            for (int sigma : new int[]{1, 2, 4, 26, n}) {
                int[] t = new int[n];
                for (int i = 0; i < n; i++) {
                    t[i] = 1 + rnd.nextInt(sigma);
                }
                for (int base : new int[]{2, FarachColton.BASE}) {
                    if (n <= 30_000 || base == FarachColton.BASE) {
                        assertSameAsDc3(t, sigma, base);
                    }
                }
            }
        }
    }

    @Test
    void repetitiveAndCombinatorialWords() {
        for (int n : new int[]{1, 2, 7, 64, 1000, 65_536}) {
            for (int base : new int[]{2, FarachColton.BASE}) {
                assertSameAsDc3(periodic(n, 1), 1, base);              // a^n
                assertSameAsDc3(periodic(n, 1, 2), 2, base);           // (ab)^n
                assertSameAsDc3(periodic(n, 1, 1, 2), 2, base);
                assertSameAsDc3(fibonacci(n), 2, base);
                assertSameAsDc3(thueMorse(n), 2, base);
            }
        }
        for (int[] k : new int[][]{{2, 4}, {2, 10}, {2, 15}, {3, 6}, {4, 8}}) {
            int[] word = deBruijn(k[0], k[1]);
            assertEquals((int) Math.pow(k[0], k[1]), word.length);
            if (word.length <= 1024) {   // every k-mer appears exactly once, cyclically
                java.util.HashSet<String> kmers = new java.util.HashSet<>();
                for (int i = 0; i < word.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    for (int j = 0; j < k[1]; j++) {
                        sb.append(word[(i + j) % word.length]).append(',');
                    }
                    kmers.add(sb.toString());
                }
                assertEquals(word.length, kmers.size());
            }
            for (int base : new int[]{2, FarachColton.BASE}) {
                assertSameAsDc3(word, k[0], base);
            }
        }
    }

    // The tree built from Farach-Colton's arrays has the structure StaticSuffixTreeTest checks, and the same nodes
    // and answers as the DC3 tree.
    @Test
    void builderGivesTheSameTreesAsDc3() {
        Random rnd = new Random(67);
        SuffixTreeBuilder fc = new FarachColtonSuffixTreeBuilder();
        SuffixTreeBuilder dc3 = new Dc3SuffixTreeBuilder();
        for (int n : new int[]{0, 1, 2, 3, 33, 100, 1000, 3000}) {
            for (int sigma : new int[]{1, 2, 26, Math.max(1, n)}) {
                long[] alphabet = new long[sigma];
                for (int i = 0; i < sigma; i++) {
                    alphabet[i] = rnd.nextLong() & Long.MAX_VALUE;
                }
                long[] text = new long[n];
                for (int i = 0; i < n; i++) {
                    text[i] = alphabet[(sigma == n) ? i : rnd.nextInt(sigma)];
                }
                for (ChildDictionary.Kind kind : ChildDictionary.Kind.values()) {
                    StaticSuffixTree a = fc.build(text, 0, n, kind);
                    StaticSuffixTree b = dc3.build(text, 0, n, kind);
                    StaticSuffixTreeTest.checkStructure(text, a);
                    assertEquals(b.internalNodeCount(), a.internalNodeCount());
                    for (int i = 0; i <= n; i++) {
                        assertEquals(b.suffix(i), a.suffix(i));
                    }
                    for (int q = 0; q < 50 && n > 0; q++) {
                        int m = 1 + rnd.nextInt(Math.min(n, 30));
                        int s = rnd.nextInt(n - m + 1);
                        long[] p = Arrays.copyOfRange(text, s, s + m);
                        int na = a.locus(p);
                        int nb = b.locus(p);
                        assertEquals(b.lo(nb), a.lo(na));
                        assertEquals(b.hi(nb), a.hi(na));
                    }
                }
            }
        }
    }

    @Test
    void rangeMinMatchesScan() {
        Random rnd = new Random(71);
        for (int n : new int[]{1, 2, 31, 32, 33, 64, 65, 1000, 5000}) {
            for (int spread : new int[]{3, 1000}) {
                int[] a = new int[n];
                for (int i = 0; i < n; i++) {
                    a[i] = rnd.nextInt(spread);
                }
                RangeMin rmq = new RangeMin(a, n);
                for (int q = 0; q < 3000; q++) {
                    int l = rnd.nextInt(n);
                    int r = l + rnd.nextInt(n - l);
                    int min = Integer.MAX_VALUE;
                    for (int i = l; i <= r; i++) {
                        min = Math.min(min, a[i]);
                    }
                    int j = rmq.argMin(l, r);
                    assertTrue(j >= l && j <= r);
                    assertEquals(min, a[j], "n=" + n + " [" + l + ", " + r + "]");
                }
            }
        }
    }

    /** Farach-Colton and DC3 + Kasai on {@code text} (symbols in [1, sigma]) must agree exactly. */
    private static void assertSameAsDc3(int[] text, int sigma, int base) {
        int n = text.length;
        int[] s = Arrays.copyOf(text, n + 2);
        int[] sa = new int[n + 1];
        int[] lcp = new int[n];
        FarachColton.build(s, n, sigma, sa, lcp, base);
        int[] ranked = Arrays.copyOf(text, n + 1);
        int[] expectedSa = SuffixArrayBuilder.buildSuffixArray(ranked);
        int[] expectedLcp = SuffixArrayBuilder.buildLcpArray(ranked, expectedSa);
        String label = "n=" + n + " sigma=" + sigma + " base=" + base
                + (n <= 40 ? " text=" + Arrays.toString(text) : "");
        assertArrayEquals(expectedSa, sa, "suffix array, " + label);
        assertArrayEquals(expectedLcp, lcp, "LCP array, " + label);
    }

    private static int[] periodic(int n, int... unit) {
        int[] t = new int[n];
        for (int i = 0; i < n; i++) {
            t[i] = unit[i % unit.length];
        }
        return t;
    }

    private static int[] fibonacci(int n) {
        StringBuilder a = new StringBuilder("a");
        StringBuilder b = new StringBuilder("ab");
        while (b.length() < n) {
            StringBuilder next = new StringBuilder(b).append(a);
            a = b;
            b = next;
        }
        int[] t = new int[n];
        for (int i = 0; i < n; i++) {
            t[i] = b.charAt(i) - 'a' + 1;
        }
        return t;
    }

    private static int[] thueMorse(int n) {
        int[] t = new int[n];
        for (int i = 0; i < n; i++) {
            t[i] = 1 + (Integer.bitCount(i) & 1);
        }
        return t;
    }

    /** A de Bruijn sequence of order {@code k} over {@code [1, sigma]} (the standard Lyndon-word construction). */
    private static int[] deBruijn(int sigma, int k) {
        List<Integer> out = new ArrayList<>();
        int[] a = new int[k + 1];
        int[] stackT = new int[k + 2];
        int[] stackP = new int[k + 2];
        int top = 0;
        stackT[top] = 1;
        stackP[top] = 1;
        top++;
        // Iterative form of db(t, p): if t > k, emit a[1..p] when k % p == 0; else a[t] = a[t-p], db(t+1, p),
        // then for j = a[t-p]+1 .. sigma-1: a[t] = j, db(t+1, t).
        int[] stackJ = new int[k + 2];
        stackJ[0] = -1;
        while (top > 0) {
            int f = top - 1;
            int t = stackT[f];
            int p = stackP[f];
            if (t > k) {
                if (k % p == 0) {
                    for (int i = 1; i <= p; i++) {
                        out.add(a[i] + 1);
                    }
                }
                top--;
                continue;
            }
            int j = stackJ[f];
            if (j == -1) {
                a[t] = a[t - p];
                stackJ[f] = a[t - p] + 1;
                stackT[top] = t + 1;
                stackP[top] = p;
                stackJ[top] = -1;
                top++;
            } else if (j < sigma) {
                a[t] = j;
                stackJ[f] = j + 1;
                stackT[top] = t + 1;
                stackP[top] = t;
                stackJ[top] = -1;
                top++;
            } else {
                top--;
            }
        }
        int[] word = new int[out.size()];
        for (int i = 0; i < word.length; i++) {
            word[i] = out.get(i);
        }
        return word;
    }
}
