package PMIndex;

import PMIndex.HierarchicalBloomFilterIndex.ProvisioningPolicy;
import PMIndex.HierarchicalBloomFilterIndex.SearchStrategy;
import org.junit.jupiter.api.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// Differential tests: every configuration must return exactly what a brute-force scan of the window returns.
class HierarchicalBloomFilterIndexTest {

    private static final String[] ALPHABETS = {"ab", "ACGT", "abcdefghijklmnopqrstuvwxyz"};

    @Test
    void matchesBruteForceForEveryStrategyAndNgram() {
        Random rnd = new Random(7);
        for (SearchStrategy strategy : SearchStrategy.values()) {
            for (int n : new int[]{1, 2, 3, 5, 8}) {
                for (boolean strides : new boolean[]{false, true}) {
                    for (String alphabet : ALPHABETS) {
                        String config = strategy + " n=" + n + " strides=" + strides + " alphabet=" + alphabet.length();
                        HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                                .windowLength(1000)
                                .treeLength(256)
                                .nGram(n)
                                .strides(strides)
                                .searchStrategy(strategy)
                                .build();
                        checkStream(index, randomText(rnd, alphabet, 3500), rnd, config);
                    }
                }
            }
        }
    }

    @Test
    void matchesBruteForceWithProvisioningPolicies() {
        Random rnd = new Random(11);
        for (ProvisioningPolicy policy : ProvisioningPolicy.values()) {
            for (SearchStrategy strategy : SearchStrategy.values()) {
                for (int n : new int[]{1, 4, 8}) {
                    String config = policy + " " + strategy + " n=" + n;
                    HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                            .windowLength(4096)
                            .treeLength(1024)
                            .nGram(n)
                            .searchStrategy(strategy)
                            .provisioningPolicy(policy)
                            .build();
                    // A skewed stream so that the policies actually drop levels.
                    checkStream(index, skewedText(rnd, 12_000), rnd, config);
                }
            }
        }
    }

    // hbi-optimizations.md, Section 4: the levels below the shortest expected pattern are not built. Patterns of
    // every length, shorter ones included, must still match brute force.
    @Test
    void minPatternLengthKeepsAnswersExact() {
        Random rnd = new Random(13);
        for (SearchStrategy strategy : SearchStrategy.values()) {
            for (int n : new int[]{1, 3, 8}) {
                for (int minLen : new int[]{4, 10, 33}) {
                    HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                            .windowLength(1000)
                            .treeLength(256)
                            .nGram(n)
                            .searchStrategy(strategy)
                            .minPatternLength(minLen)
                            .build();
                    checkStream(index, randomText(rnd, "abcd", 3000), rnd, strategy + " n=" + n + " minLen=" + minLen);
                }
            }
        }
        assertThrows(IllegalArgumentException.class,
                () -> HierarchicalBloomFilterIndex.builder().windowLength(64).minPatternLength(-1).build());
    }

    @Test
    void lowEntropyAndPeriodicStreams() {
        Random rnd = new Random(3);
        String[] streams = {
                "a".repeat(3000),
                "ab".repeat(1500),
                "aab".repeat(1000),
                debruijnLike(3000),
        };
        for (String stream : streams) {
            for (SearchStrategy strategy : SearchStrategy.values()) {
                for (int n : new int[]{1, 3}) {
                    HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                            .windowLength(512)
                            .treeLength(128)
                            .nGram(n)
                            .searchStrategy(strategy)
                            .build();
                    checkStream(index, stream, rnd, strategy + " n=" + n + " periodic");
                }
            }
        }
    }

    @Test
    void oneTreePerWindowKeepsTheWholeWindowSearchable() {
        // With treeLength == windowLength the oldest tree must stay until it leaves the window entirely.
        Random rnd = new Random(5);
        HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                .windowLength(256)
                .treeLength(256)
                .build();
        String text = randomText(rnd, "abcd", 1000);
        StringBuilder seen = new StringBuilder();
        for (int i = 0; i < text.length(); i++) {
            index.append(text.charAt(i));
            seen.append(text.charAt(i));
            if (seen.length() > 256) {
                String oldest = seen.substring(seen.length() - 256, seen.length() - 256 + 6);
                assertArrayEquals(bruteForce(seen, 256, oldest), index.search(oldest), "position " + i);
            }
        }
    }

    @Test
    void patternsOutsideTheWindowAreNotReported() {
        HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                .windowLength(8)
                .treeLength(4)
                .build();
        index.append("xyzabcdefgh");
        assertEquals(3, index.windowStart());
        assertArrayEquals(new long[0], index.search("xyz"));
        assertArrayEquals(new long[0], index.search("zab"));
        assertArrayEquals(new long[]{3}, index.search("abc"));
        assertArrayEquals(new long[]{3}, index.search("abcdefgh"));
        assertArrayEquals(new long[0], index.search("abcdefghi"));
    }

    @Test
    void rejectsPatternsShorterThanTheNgram() {
        HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
                .windowLength(64)
                .nGram(4)
                .build();
        index.append("hello world");
        assertThrows(IllegalArgumentException.class, () -> index.search("abc"));
        assertThrows(IllegalArgumentException.class, () -> index.search(""));
        assertArrayEquals(new long[]{6}, index.search("world"));
    }

    @Test
    void emptyIndexHasNoMatches() {
        HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder().windowLength(16).build();
        assertArrayEquals(new long[0], index.search("a"));
        assertFalse(index.contains("a"));
    }

    @Test
    void builderValidatesArguments() {
        assertThrows(IllegalArgumentException.class, () -> HierarchicalBloomFilterIndex.builder().build());
        assertThrows(IllegalArgumentException.class,
                () -> HierarchicalBloomFilterIndex.builder().windowLength(64).treeLength(48).build());
        assertThrows(IllegalArgumentException.class,
                () -> HierarchicalBloomFilterIndex.builder().windowLength(4).nGram(5).build());
        assertThrows(IllegalArgumentException.class,
                () -> HierarchicalBloomFilterIndex.builder().windowLength(64).fpRate(1.0).build());
    }

    // Streams text into the index and compares query answers with brute force at several points.
    private static void checkStream(HierarchicalBloomFilterIndex index, String text, Random rnd, String config) {
        int w = index.windowLength();
        int n = index.nGram();
        StringBuilder seen = new StringBuilder(text.length());
        int nextCheck = Math.max(n, w / 3);
        for (int i = 0; i < text.length(); i++) {
            index.append(text.charAt(i));
            seen.append(text.charAt(i));
            if (seen.length() == nextCheck || i == text.length() - 1) {
                for (String pattern : queries(seen, w, n, rnd)) {
                    long[] expected = bruteForce(seen, w, pattern);
                    long[] actual = index.search(pattern);
                    assertArrayEquals(expected, actual,
                            config + " after " + seen.length() + " chars, pattern \"" + pattern + "\"");
                    assertEquals(expected.length > 0, index.contains(pattern),
                            config + " contains, after " + seen.length() + " chars, pattern \"" + pattern + "\"");
                }
                nextCheck += w / 3 + rnd.nextInt(w / 2 + 1);
            }
        }
    }

    private static List<String> queries(CharSequence seen, int w, int n, Random rnd) {
        int len = seen.length();
        int start = Math.max(0, len - w);
        int available = len - start;
        List<String> out = new ArrayList<>();
        int maxLen = Math.min(available, 48);
        if (maxLen < n) {
            return out;
        }
        for (int q = 0; q < 25; q++) {
            int r = n + rnd.nextInt(maxLen - n + 1);
            // Present in the window, including across tree boundaries.
            int s = start + rnd.nextInt(available - r + 1);
            String present = seen.subSequence(s, s + r).toString();
            out.add(present);
            // A mutated copy, usually absent.
            char[] mutated = present.toCharArray();
            mutated[rnd.nextInt(r)] ^= 1;
            out.add(new String(mutated));
            // Text that just slid out of the window, when there is some.
            if (start >= r) {
                int old = Math.max(0, start - r + 1 + rnd.nextInt(r));
                out.add(seen.subSequence(old, old + r).toString());
            }
        }
        // The whole window and the window plus one character.
        out.add(seen.subSequence(start, len).toString());
        if (available >= 1) {
            out.add(seen.subSequence(start, len) + "x");
        }
        return out;
    }

    // Offsets s with s >= windowStart and seen[s, s + |pattern|) == pattern.
    private static long[] bruteForce(CharSequence seen, int w, String pattern) {
        String s = seen.toString();
        int start = Math.max(0, s.length() - w);
        List<Long> hits = new ArrayList<>();
        for (int i = s.indexOf(pattern, start); i >= 0; i = s.indexOf(pattern, i + 1)) {
            hits.add((long) i);
        }
        return hits.stream().mapToLong(Long::longValue).toArray();
    }

    private static String randomText(Random rnd, String alphabet, int length) {
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append(alphabet.charAt(rnd.nextInt(alphabet.length())));
        }
        return sb.toString();
    }

    // Mostly 'e' and 't', with a long tail of rare letters.
    private static String skewedText(Random rnd, int length) {
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            double u = rnd.nextDouble();
            if (u < 0.5) sb.append('e');
            else if (u < 0.8) sb.append('t');
            else sb.append((char) ('a' + rnd.nextInt(26)));
        }
        return sb.toString();
    }

    // Repeats every 3-letter word over {a,b,c} back to back, like the adversarial stream of Sec. 5.6.3.
    private static String debruijnLike(int length) {
        StringBuilder cycle = new StringBuilder();
        for (char x = 'a'; x <= 'c'; x++)
            for (char y = 'a'; y <= 'c'; y++)
                for (char z = 'a'; z <= 'c'; z++)
                    cycle.append(x).append(y).append(z);
        StringBuilder sb = new StringBuilder(length);
        while (sb.length() < length) sb.append(cycle);
        return sb.substring(0, length);
    }
}
