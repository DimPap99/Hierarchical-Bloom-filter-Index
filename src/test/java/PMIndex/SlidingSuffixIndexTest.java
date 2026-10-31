package PMIndex;

import org.junit.jupiter.api.Test;
import search.Pattern;
import tree.ssws.ChildDictionary;
import tree.ssws.SuffixTreeBuilder;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

// The paper's amortized structure (SlidingSuffixIndex) must answer exactly like a brute-force scan of the window,
// and exactly like LEGACY (StreamingSlidingWindowIndex), with each occurrence reported once.
class SlidingSuffixIndexTest {

    /** Every builder with every child dictionary. */
    private record Config(SuffixTreeBuilder.Kind builder, ChildDictionary.Kind dictionary) {
        SlidingSuffixIndex index(int w) {
            return new SlidingSuffixIndex(w, builder.create(), dictionary);
        }

        @Override
        public String toString() {
            return builder + "/" + dictionary;
        }
    }

    private static final List<Config> CONFIGS = configs();

    private static List<Config> configs() {
        List<Config> out = new ArrayList<>();
        for (SuffixTreeBuilder.Kind b : SuffixTreeBuilder.Kind.values()) {
            for (ChildDictionary.Kind d : ChildDictionary.Kind.values()) {
                out.add(new Config(b, d));
            }
        }
        return out;
    }

    // BaselineCorrectnessTest's random texts.
    @Test
    void matchesBruteForceAndLegacyOnRandomText() {
        Random rnd = new Random(17);
        for (Config config : CONFIGS) {
            for (int n : new int[]{1, 3}) {
                for (String alphabet : new String[]{"ab", "acgt", "abcdefghijklmnopqrstuvwxyz"}) {
                    check(randomText(rnd, alphabet, 2500), n, 512, 211, config, rnd, "alphabet " + alphabet);
                }
            }
        }
    }

    // Periodic and highly repetitive text: deep suffix trees and many overlapping occurrences.
    @Test
    void matchesBruteForceAndLegacyOnPeriodicText() {
        Random rnd = new Random(29);
        String[] texts = {
                "ab".repeat(1500), "abc".repeat(1000), "aab".repeat(1000), "a".repeat(3000), fibonacciWord(3000)
        };
        for (Config config : CONFIGS) {
            for (String text : texts) {
                for (int n : new int[]{1, 2, 3}) {
                    for (int w : new int[]{64, 512}) {
                        check(text, n, w, 97, config, rnd, "periodic text starting " + text.substring(0, 6));
                    }
                }
            }
        }
    }

    // Every symbol (nearly) distinct: n-grams of random text, and single characters from a huge alphabet.
    @Test
    void matchesBruteForceAndLegacyWhenSymbolsAreDistinct() {
        Random rnd = new Random(31);
        StringBuilder wide = new StringBuilder();
        for (int i = 0; i < 3000; i++) {
            wide.append((char) (0x4E00 + rnd.nextInt(20_000)));
        }
        for (Config config : CONFIGS) {
            check(randomText(rnd, "abcdefghijklmnopqrstuvwxyz", 3000), 8, 512, 211, config, rnd, "8-grams");
            check(wide.toString(), 1, 512, 211, config, rnd, "20,000-symbol alphabet");
        }
    }

    // Small windows, streams much longer than the window, a query after every insert for every pattern length
    // 1..w (and w+1). Patterns are drawn from the last 2w symbols, so many also occur before the window start.
    @Test
    void everyPatternLengthAfterEveryInsert() {
        Random rnd = new Random(37);
        for (Config config : CONFIGS) {
            for (int w : new int[]{16, 32, 64}) {
                for (String alphabet : new String[]{"a", "ab", "abc"}) {
                    SlidingSuffixIndex index = config.index(w);
                    StreamingSlidingWindowIndex legacy = new StreamingSlidingWindowIndex(w);
                    StringBuilder seen = new StringBuilder();
                    for (int t = 0; t < 10 * w + 13; t++) {
                        String c = String.valueOf(alphabet.charAt(rnd.nextInt(alphabet.length())));
                        seen.append(c);
                        index.insert(c);
                        legacy.insert(c);
                        int windowStart = Math.max(0, seen.length() - w);
                        assertEquals(windowStart, index.windowStartPosition());
                        int from = Math.max(0, windowStart - w);
                        for (int m = 1; m <= Math.min(w + 1, seen.length() - from); m++) {
                            for (int q = 0; q < 2; q++) {
                                int s = from + rnd.nextInt(seen.length() - from - m + 1);
                                char[] p = seen.substring(s, s + m).toCharArray();
                                if (q == 1) {
                                    p[rnd.nextInt(m)] = alphabet.charAt(rnd.nextInt(alphabet.length()));
                                }
                                assertSameAnswers(index, legacy, seen, windowStart, new String(p), 1,
                                        config + " w=" + w + " alphabet " + alphabet + " after " + seen.length());
                            }
                        }
                    }
                }
            }
        }
    }

    // A long run of one symbol makes suffix trees as deep as they are long.
    @Test
    void longRunOfOneSymbolIsHandled() {
        int w = 1 << 16;
        for (Config config : CONFIGS) {
            SlidingSuffixIndex index = config.index(w);
            for (int i = 0; i < 3 * w; i++) {
                index.insert("a");
            }
            ArrayList<Integer> hits = index.report(new Pattern("aaaaaaaa", 1));
            TreeSet<Integer> set = new TreeSet<>(hits);
            assertEquals(w - 7, hits.size());
            assertEquals(w - 7, set.size());
            assertEquals(0, (int) set.first());
            assertEquals(w - 8, (int) set.last());
            assertEquals(2L * w, index.windowStartPosition());
        }
    }

    // Sec. 3.1: sizes are distinct powers of two, at most w, decreasing to the right; once w symbols have arrived
    // the leftmost has size w and w <= |S| <= 2w - 1, before that everything is stored.
    @Test
    void segmentsFormABinaryCounterCappedAtW() {
        for (int w : new int[]{1, 2, 4, 8, 64}) {
            SlidingSuffixIndex index = new SlidingSuffixIndex(w);
            for (int t = 1; t <= 5 * w + 3; t++) {
                index.insert("x" + (t % 3));
                int[] sizes = index.segmentSizes();
                long stored = 0;
                for (int x = 0; x < sizes.length; x++) {
                    assertEquals(1, Integer.bitCount(sizes[x]));
                    assertTrue(sizes[x] <= w);
                    assertTrue(x == 0 || sizes[x] < sizes[x - 1]);
                    stored += sizes[x];
                }
                if (t >= w) {
                    assertEquals(w, sizes[0]);
                    assertEquals(w + (t % w), stored);
                } else {
                    assertEquals(t, stored);
                }
            }
        }
    }

    @Test
    void windowMustBeAPowerOfTwo() {
        assertThrows(IllegalArgumentException.class, () -> new SlidingSuffixIndex(0));
        assertThrows(IllegalArgumentException.class, () -> new SlidingSuffixIndex(48));
        assertThrows(IllegalArgumentException.class, () -> new SlidingSuffixIndex((1 << 20) + 1));
    }

    // Streams text into both indexes as n-grams and, every `every` characters, compares sampled patterns.
    private static void check(String text, int n, int w, int every, Config config, Random rnd,
                              String label) {
        SlidingSuffixIndex index = config.index(w);
        StreamingSlidingWindowIndex legacy = new StreamingSlidingWindowIndex(w);
        StringBuilder seen = new StringBuilder();
        for (int i = 0; i < text.length(); i++) {
            seen.append(text.charAt(i));
            if (seen.length() >= n) {
                String gram = seen.substring(seen.length() - n);
                index.insert(gram);
                legacy.insert(gram);
            }
            if (seen.length() % every == 0 && seen.length() >= n) {
                int symbols = seen.length() - n + 1;
                int windowStart = Math.max(0, symbols - w);
                assertEquals(windowStart, index.windowStartPosition(), label);
                for (String pattern : patterns(seen, windowStart, n, rnd)) {
                    assertSameAnswers(index, legacy, seen, windowStart, pattern, n,
                            config + " " + label + ", n=" + n + ", w=" + w + ", after " + seen.length() + " chars");
                }
            }
        }
    }

    private static void assertSameAnswers(SlidingSuffixIndex index, StreamingSlidingWindowIndex legacy,
                                          CharSequence seen, int windowStart, String pattern, int n, String label) {
        ArrayList<Integer> got = index.report(new Pattern(pattern, n));
        TreeSet<Integer> reported = new TreeSet<>();
        for (int relative : got) {
            reported.add(windowStart + relative);
        }
        assertEquals(got.size(), reported.size(), "reported twice: " + label + ", pattern " + pattern);
        assertEquals(bruteForce(seen, windowStart, pattern), reported, label + ", pattern " + pattern);
        TreeSet<Integer> old = new TreeSet<>();
        for (int relative : legacy.report(new Pattern(pattern, n))) {
            old.add(windowStart + relative);
        }
        assertEquals(old, reported, "LEGACY differs: " + label + ", pattern " + pattern);
    }

    // Patterns from the window and from the w symbols before it, plus one mutated copy of each.
    private static List<String> patterns(CharSequence seen, int windowStart, int n, Random rnd) {
        List<String> out = new ArrayList<>();
        int from = Math.max(0, windowStart - 64);
        int available = seen.length() - from;
        for (int q = 0; q < 20; q++) {
            int r = n + rnd.nextInt(Math.min(40, available) - n + 1);
            int s = from + rnd.nextInt(available - r + 1);
            String present = seen.subSequence(s, s + r).toString();
            out.add(present);
            char[] mutated = present.toCharArray();
            mutated[rnd.nextInt(r)] ^= 1;
            out.add(new String(mutated));
        }
        return out;
    }

    private static TreeSet<Integer> bruteForce(CharSequence seen, int windowStart, String pattern) {
        String s = seen.toString();
        TreeSet<Integer> hits = new TreeSet<>();
        for (int i = s.indexOf(pattern, windowStart); i >= 0; i = s.indexOf(pattern, i + 1)) {
            hits.add(i);
        }
        return hits;
    }

    private static String randomText(Random rnd, String alphabet, int length) {
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append(alphabet.charAt(rnd.nextInt(alphabet.length())));
        }
        return sb.toString();
    }

    private static String fibonacciWord(int length) {
        String a = "a", b = "ab";
        while (b.length() < length) {
            String next = b + a;
            a = b;
            b = next;
        }
        return b.substring(0, length);
    }
}
