package PMIndex;

import org.junit.jupiter.api.Test;
import search.Pattern;
import utilities.IndexFactory;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.TreeSet;

import static org.junit.jupiter.api.Assertions.assertEquals;

// The suffix-based baselines must answer exactly like a brute-force scan, or the benchmark comparison is meaningless.
class BaselineCorrectnessTest {

    // SSWSI (Bille et al., CPM'23) reports the occurrences inside the last w symbols.
    @Test
    void slidingWindowIndexMatchesBruteForce() {
        Random rnd = new Random(17);
        for (int n : new int[]{1, 3}) {
            for (String alphabet : new String[]{"ab", "acgt", "abcdefghijklmnopqrstuvwxyz"}) {
                checkSlidingWindow(randomText(rnd, alphabet, 2500), n, 512, 211, rnd, "alphabet " + alphabet);
            }
        }
    }

    // Periodic and highly repetitive text: deep suffix trees and many overlapping occurrences.
    @Test
    void slidingWindowIndexMatchesBruteForceOnPeriodicText() {
        Random rnd = new Random(29);
        String[] texts = {
                repeat("ab", 1500), repeat("abc", 1000), repeat("aab", 1000), repeat("a", 3000), fibonacciWord(3000)
        };
        for (String text : texts) {
            for (int n : new int[]{1, 2, 3}) {
                for (int w : new int[]{64, 512}) {
                    checkSlidingWindow(text, n, w, 97, rnd, "periodic text starting " + text.substring(0, 6));
                }
            }
        }
    }

    // Every symbol (nearly) distinct: n-grams of random text, and single characters from a huge alphabet.
    @Test
    void slidingWindowIndexMatchesBruteForceWhenSymbolsAreDistinct() {
        Random rnd = new Random(31);
        checkSlidingWindow(randomText(rnd, "abcdefghijklmnopqrstuvwxyz", 3000), 8, 512, 211, rnd, "8-grams");
        StringBuilder wide = new StringBuilder();
        for (int i = 0; i < 3000; i++) {
            wide.append((char) (0x4E00 + rnd.nextInt(20_000)));
        }
        checkSlidingWindow(wide.toString(), 1, 512, 211, rnd, "20,000-symbol alphabet");
    }

    // A long run of one symbol makes suffix trees as deep as they are long. This used to overflow the stack.
    @Test
    void longRunOfOneSymbolIsHandled() {
        int w = 1 << 16;
        StreamingSlidingWindowIndex index = new StreamingSlidingWindowIndex(w);
        for (int i = 0; i < 3 * w; i++) {
            index.insert("a");
        }
        ArrayList<Integer> hits = index.report(new Pattern("aaaaaaaa", 1));
        assertEquals(w - 7, hits.size());
        assertEquals(0, (int) new TreeSet<>(hits).first());
        assertEquals(w - 8, (int) new TreeSet<>(hits).last());
        assertEquals(2L * w, index.windowStartPosition());
    }

    // The online suffix tree indexes everything it was given, as in the benchmark (one window per run).
    @Test
    void onlineSuffixTreeMatchesBruteForce() {
        Random rnd = new Random(23);
        for (int n : new int[]{1, 3}) {
            for (String alphabet : new String[]{"ab", "acgt", "abcdefghijklmnopqrstuvwxyz"}) {
                IPMIndexing index = IndexFactory.createSuffixTreeIndex(1 << 16);
                String text = randomText(rnd, alphabet, 1500);
                for (int i = 0; i + n <= text.length(); i++) {
                    index.insert(text.substring(i, i + n));
                }
                index.insert("$");
                for (String pattern : patterns(text, 0, n, rnd)) {
                    assertEquals(bruteForce(text, 0, pattern),
                            new TreeSet<>(index.report(new Pattern(pattern, n))),
                            "suffix tree n=" + n + ", pattern " + pattern);
                }
            }
        }
    }

    // Streams text into SSWSI as n-grams and, every `every` characters, compares sampled patterns with brute force.
    private static void checkSlidingWindow(String text, int n, int w, int every, Random rnd, String label) {
        StreamingSlidingWindowIndex index = new StreamingSlidingWindowIndex(w);
        StringBuilder seen = new StringBuilder();
        for (int i = 0; i < text.length(); i++) {
            seen.append(text.charAt(i));
            if (seen.length() >= n) {
                index.insert(seen.substring(seen.length() - n));
            }
            if (seen.length() % every == 0 && seen.length() >= n) {
                int symbols = seen.length() - n + 1;
                int windowStart = Math.max(0, symbols - w);
                assertEquals(windowStart, index.windowStartPosition(), label);
                for (String pattern : patterns(seen, windowStart, n, rnd)) {
                    TreeSet<Integer> reported = new TreeSet<>();
                    for (int relative : index.report(new Pattern(pattern, n))) {
                        reported.add(windowStart + relative);
                    }
                    assertEquals(bruteForce(seen, windowStart, pattern), reported,
                            "SSWSI " + label + ", n=" + n + ", w=" + w + ", after " + seen.length()
                                    + " chars, pattern " + pattern);
                }
            }
        }
    }

    private static List<String> patterns(CharSequence seen, int windowStart, int n, Random rnd) {
        List<String> out = new ArrayList<>();
        int available = seen.length() - windowStart;
        for (int q = 0; q < 20; q++) {
            int r = n + rnd.nextInt(Math.min(40, available) - n + 1);
            int s = windowStart + rnd.nextInt(available - r + 1);
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

    private static String repeat(String unit, int times) {
        return unit.repeat(times);
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
