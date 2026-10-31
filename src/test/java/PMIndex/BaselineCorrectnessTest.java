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

    // SSWSI (Bille et al., CPM'23) reports absolute positions of occurrences inside the last w symbols.
    @Test
    void slidingWindowIndexMatchesBruteForce() {
        Random rnd = new Random(17);
        for (int n : new int[]{1, 3}) {
            for (String alphabet : new String[]{"ab", "acgt", "abcdefghijklmnopqrstuvwxyz"}) {
                int w = 512;
                StreamingSlidingWindowIndex index = new StreamingSlidingWindowIndex(w);
                String text = randomText(rnd, alphabet, 2500);
                StringBuilder seen = new StringBuilder();
                for (int i = 0; i < text.length(); i++) {
                    seen.append(text.charAt(i));
                    if (seen.length() >= n) {
                        index.insert(seen.substring(seen.length() - n));
                    }
                    if (seen.length() % 211 == 0) {
                        int symbols = seen.length() - n + 1;
                        int windowStart = Math.max(0, symbols - w);
                        for (String pattern : patterns(seen, windowStart, n, rnd)) {
                            assertEquals(bruteForce(seen, windowStart, pattern),
                                    new TreeSet<>(index.report(new Pattern(pattern, n))),
                                    "SSWSI n=" + n + " after " + seen.length() + " chars, pattern " + pattern);
                        }
                    }
                }
            }
        }
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
}
