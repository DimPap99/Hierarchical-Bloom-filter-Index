package PMIndex;

import estimators.IntCountSketch;
import membership.BloomFilter;
import membership.Membership;
import org.junit.jupiter.api.Test;
import search.Pattern;
import tree.ImplicitTree;
import tree.TreeLayout;
import utilities.IndexFactory;
import utilities.MathUtils;
import utilities.MultiBenchmarkOptions;
import utilities.RingBuffer;
import utilities.Utils;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// Regression tests for the code-review findings, exercised through the benchmark code path
// (IndexFactory + string tokens + Pattern), which the public facade does not use.
class HbiRegressionTest {

    @Test
    void log2IsBaseTwo() {
        assertEquals(10.0, MathUtils.log2(1024.0), 1e-12);
    }

    @Test
    void pruningLevelStaysInsideTheTree() {
        ImplicitTree<Membership> tree = tree(1 << 12, 3);
        int lmax = tree.maxDepth() - 1;
        for (double p : new double[]{0.0, 1e-9, 1e-4, 0.01, 0.2, 0.9, 1.0}) {
            int bloom = MathUtils.pruningLevelBloom(tree, 0.99, p, 0.2);
            int plain = MathUtils.pruningLevel(tree, 0.99, p);
            assertTrue(bloom >= 0 && bloom <= lmax, "pruningLevelBloom(" + p + ") = " + bloom);
            assertTrue(plain >= 0 && plain <= lmax, "pruningLevel(" + p + ") = " + plain);
        }
    }

    @Test
    void dropFiltersUsesAbsoluteLevels() {
        ImplicitTree<Membership> tree = tree(1 << 10, 1);
        tree.dropFiltersUpToLp(3);
        assertEquals(3, tree.effectiveRoot());
        tree.dropFiltersUpToLp(3);  // the predictive policy calls it again when the tree fills
        assertEquals(3, tree.effectiveRoot());
        tree.dropFiltersUpToLp(2);
        assertEquals(3, tree.effectiveRoot());
        tree.dropFiltersUpToLp(100);
        assertEquals(tree.maxDepth() - 1, tree.effectiveRoot(), "the deepest level is always kept");
        assertDoesNotThrow(() -> tree.getMembershipFpRate(0));
        assertThrows(IllegalArgumentException.class, () -> tree.getMembershipFpRate(tree.maxDepth()));
    }

    @Test
    void countSketchMedianWithEvenDepthIsTheMeanOfTheMiddleRows() {
        IntCountSketch sketch = new IntCountSketch(1 << 16, 8, 42L);
        for (int i = 0; i < 100; i++) sketch.add(7L);
        assertEquals(100L, sketch.estimate(7L));
    }

    @Test
    void segmentNgramsDoNotCollide() {
        RingBuffer<String> window = new RingBuffer<>(2);
        window.append("80");
        window.append("443");
        Pattern other = new Pattern(new String[]{"804", "43"}, 2);
        Pattern same = new Pattern(new String[]{"80", "443"}, 2);
        assertNotEquals(window.join(Pattern.SEGMENT_JOINER), other.nGramArr[0]);
        assertEquals(window.join(Pattern.SEGMENT_JOINER), same.nGramArr[0]);
    }

    @Test
    void windowTokenDeterminesTheWindowLength() {
        assertEquals(1 << 18, MultiBenchmarkOptions.deriveWindowLength("w18", 1 << 21));
        assertEquals(1 << 21, MultiBenchmarkOptions.deriveWindowLength("wzipf_21_e13_128", 1 << 20));
        assertEquals(1 << 19, MultiBenchmarkOptions.deriveWindowLength("adversarial/blocks/blocks_w19", 1 << 21));
        assertEquals(1 << 21, MultiBenchmarkOptions.deriveWindowLength("1023", 1 << 21));
    }

    @Test
    void configurationDefaultsAreNotExperimentMode() {
        HbiConfiguration config = HbiConfiguration.builder()
                .searchAlgorithm(new search.BlockSearch())
                .windowLength(64).treeLength(64).fpRate(0.2).alphabetSize(16)
                .estimatorSupplier(() -> new estimators.CSEstimator(64, 3, 1))
                .membershipSupplier(BloomFilter::new)
                .pruningPlanSupplier(() -> new search.MostFreqPruning(0.99, 0.2))
                .verifier(new search.VerifierLinearLeafProbe())
                .confidence(0.99)
                .build();
        assertFalse(config.experimentMode());
        assertEquals(1, config.nGram());
        assertEquals(Utils.MemPolicy.NONE, config.memPolicy());
    }

    // Findings #1 and #2: PREDICTIVE pruned twice and crashed on insert; CF + REACTIVE crashed on query.
    @Test
    void benchmarkIndexMatchesBruteForceUnderPolicies() {
        Random rnd = new Random(1);
        for (Utils.MemPolicy policy : Utils.MemPolicy.values()) {
            for (String algorithm : new String[]{"bs", "cp", "cgbs", "cf"}) {
                for (int n : new int[]{1, 3}) {
                    HBI hbi = IndexFactory.createHbi(2048, 512, 1 << 16, 0.2, 0.99, policy, n, algorithm, 0.05, 4096);
                    String config = policy + " " + algorithm + " n=" + n;
                    String text = skewed(rnd, 9000);
                    StringBuilder seen = new StringBuilder();
                    for (int i = 0; i < text.length(); i++) {
                        seen.append(text.charAt(i));
                        if (seen.length() >= n) {
                            hbi.insert(seen.substring(seen.length() - n));
                        }
                        if (seen.length() % 1500 == 0) {
                            if (policy != Utils.MemPolicy.NONE) {
                                hbi.forceApplyMemoryPolicy();
                            }
                            int windowStart = Math.max(0, seen.length() - n + 1 - 2048);
                            for (int q = 0; q < 20; q++) {
                                int r = n + rnd.nextInt(30);
                                int s = windowStart + rnd.nextInt(seen.length() - r - windowStart + 1);
                                String pattern = seen.substring(s, s + r);
                                List<Integer> want = expected(seen, 2048, n, pattern);
                                String where = config + " after " + seen.length() + " chars, pattern " + pattern;
                                assertEquals(want, hbi.report(new Pattern(pattern, n)), where);
                                // exists stops at the first hit that the predicate takes.
                                assertEquals(!want.isEmpty(), hbi.exists(pattern), where);
                                int first = want.isEmpty() ? -1 : want.get(0);
                                assertEquals(want.size() > 1, hbi.exists(new Pattern(pattern, n), pos -> pos > first), where);
                                assertFalse(hbi.exists(pattern + "#"), where);
                            }
                        }
                    }
                }
            }
        }
    }

    // Filter seeds come from the index's seed, its trees and levels: the same seed probes exactly the same way.
    @Test
    void sameSeedGivesTheSameProbes() {
        Random rnd = new Random(8);
        String text = skewed(rnd, 6000);
        long[] totals = new long[3];
        long[] seeds = {1L, 1L, 2L};
        List<List<Integer>> answers = new ArrayList<>();
        for (int run = 0; run < 3; run++) {
            HBI hbi = IndexFactory.createHbi(2048, 512, 1 << 16, 0.3, 0.99, Utils.MemPolicy.NONE, 2, "bs", 0.05, 0,
                    new IndexFactory.HbiTuning(0, membership.FilterKind.BLOOM, seeds[run], search.BlockSearch.Options.DEFAULT, false));
            for (int i = 2; i <= text.length(); i++) hbi.insert(text.substring(i - 2, i));
            Random queries = new Random(9);
            hbi.getAllprobes();
            for (int q = 0; q < 200; q++) {
                int r = 2 + queries.nextInt(20);
                int s = text.length() - 2000 + queries.nextInt(2000 - r);
                List<Integer> hits = hbi.report(new Pattern(text.substring(s, s + r), 2));
                if (run == 0) answers.add(hits); else assertEquals(answers.get(q), hits);
            }
            totals[run] = hbi.getAllprobes();
        }
        assertEquals(totals[0], totals[1], "same seed");
        assertNotEquals(totals[0], totals[2], "another seed");
    }

    // Block Search's optional probe strategies (query plan, sketch order, move-to-front, neighbour check) change
    // only which filters are probed, never the answers: every combination, with and without strides, BS and CP-BS.
    @Test
    void probeStrategiesMatchBruteForce() {
        Random rnd = new Random(12);
        for (int bits = 0; bits < 16; bits++) {
            search.BlockSearch.Options o = new search.BlockSearch.Options(
                    (bits & 1) != 0, (bits & 2) != 0, (bits & 4) != 0, (bits & 8) != 0);
            for (String algorithm : new String[]{"bs", "cp"}) {
                for (int n : new int[]{1, 3}) {
                    for (boolean strides : new boolean[]{false, true}) {
                        HBI hbi = IndexFactory.createHbi(2048, 512, 1 << 16, 0.2, 0.99, Utils.MemPolicy.NONE, n,
                                algorithm, 0.05, 0, new IndexFactory.HbiTuning(0, membership.FilterKind.BLOOM, 1L, o, false));
                        hbi.strides = strides;
                        String config = o + " " + algorithm + " n=" + n + " strides=" + strides;
                        String text = skewed(rnd, 4000);
                        for (int i = n; i <= text.length(); i++) hbi.insert(text.substring(i - n, i));
                        int windowStart = Math.max(0, text.length() - n + 1 - 2048);
                        for (int q = 0; q < 40; q++) {
                            int r = n + rnd.nextInt(40);
                            int s = windowStart + rnd.nextInt(text.length() - r - windowStart + 1);
                            String pattern = text.substring(s, s + r);
                            if (q % 3 == 0) pattern = pattern.substring(0, r - 1) + "#";
                            assertEquals(expected(text, 2048, n, pattern), hbi.report(new Pattern(pattern, n)),
                                    config + ", pattern " + pattern);
                        }
                    }
                }
            }
        }
    }

    // The KMP verifier (Section 7) reports exactly what the linear one does, also on repetitive streams where
    // matches overlap and ranges run back to back, across tree boundaries, with and without strides.
    @Test
    void kmpVerifierMatchesBruteForce() {
        Random rnd = new Random(21);
        String[] streams = {skewed(rnd, 5000), "a".repeat(5000), "ab".repeat(2500), "aab".repeat(1700)};
        for (String text : streams) {
            for (int n : new int[]{1, 2, 3}) {
                for (boolean strides : new boolean[]{false, true}) {
                    HBI hbi = new HBI(HbiConfiguration.builder()
                            .searchAlgorithm(new search.BlockSearch())
                            .windowLength(2048).treeLength(256).fpRate(0.2).alphabetSize(1 << 16)
                            .estimatorSupplier(() -> new estimators.CSEstimator(64, 3, 1))
                            .membershipSupplier(BloomFilter::new)
                            .pruningPlanSupplier(() -> new search.MostFreqPruning(0.99, 0.2))
                            .verifier(new search.VerifierKmp())
                            .costFunction(new estimators.CostFunctionDefaultRoot())
                            .confidence(0.99).nGram(n).build());
                    hbi.strides = strides;
                    for (int i = n; i <= text.length(); i++) hbi.insert(text.substring(i - n, i));
                    int windowStart = Math.max(0, text.length() - n + 1 - 2048);
                    for (int q = 0; q < 40; q++) {
                        int r = n + rnd.nextInt(40);
                        int s = windowStart + rnd.nextInt(text.length() - r - windowStart + 1);
                        String pattern = text.substring(s, s + r);
                        if (q % 4 == 0) pattern = pattern.substring(0, r - 1) + "#";
                        assertEquals(expected(text, 2048, n, pattern), hbi.report(new Pattern(pattern, n)),
                                "kmp n=" + n + " strides=" + strides + " stream " + text.substring(0, 6) + ", pattern " + pattern);
                    }
                }
            }
        }
        assertArrayEquals(new int[]{0, 1, 2, 0, 1, 2, 3}, search.VerifierKmp.failureTable(new long[]{1, 1, 1, 2, 1, 1, 1}));
    }

    // The benchmark path with blocked filters: every policy (they read getFpRate) and strategy stays exact.
    @Test
    void blockedFiltersMatchBruteForceUnderPolicies() {
        Random rnd = new Random(2);
        for (Utils.MemPolicy policy : Utils.MemPolicy.values()) {
            for (String algorithm : new String[]{"bs", "cp", "cgbs", "cf"}) {
                HBI hbi = IndexFactory.createHbi(2048, 512, 1 << 16, 0.2, 0.99, policy, 2, algorithm, 0.05, 4096,
                        new IndexFactory.HbiTuning(0, membership.FilterKind.BLOCKED, 1L, search.BlockSearch.Options.DEFAULT, false));
                String config = policy + " " + algorithm + " blocked";
                String text = skewed(rnd, 6000);
                StringBuilder seen = new StringBuilder();
                for (int i = 0; i < text.length(); i++) {
                    seen.append(text.charAt(i));
                    if (seen.length() >= 2) hbi.insert(seen.substring(seen.length() - 2));
                    if (seen.length() % 1500 == 0) {
                        if (policy != Utils.MemPolicy.NONE) hbi.forceApplyMemoryPolicy();
                        int windowStart = Math.max(0, seen.length() - 1 - 2048);
                        for (int q = 0; q < 20; q++) {
                            int r = 2 + rnd.nextInt(30);
                            int s = windowStart + rnd.nextInt(seen.length() - r - windowStart + 1);
                            String pattern = seen.substring(s, s + r);
                            assertEquals(expected(seen, 2048, 2, pattern), hbi.report(new Pattern(pattern, 2)),
                                    config + " after " + seen.length() + " chars, pattern " + pattern);
                        }
                    }
                }
            }
        }
    }

    // Finding #9: a partial probe on an interval shorter than the pattern must not push the candidate
    // start back over the frontier (or before the interval, whose earlier starts belong to its left sibling).
    @Test
    void partialProbeNeverStartsBeforeTheFrontier() {
        ImplicitTree<Membership> tree = new ImplicitTree<>(new TreeLayout(7, 1), level -> new membership.MockMembership(),
                new membership.KeyPackingService(7, 16), null);
        for (int i = 0; i < 64; i++) tree.append('x', i);
        Pattern pattern = new Pattern("x".repeat(16) + "y", 1);
        long[] tokens = new long[17];
        java.util.Arrays.fill(tokens, 'x');
        tokens[16] = 'y';
        pattern.setTokens(tokens);
        pattern.buildStrideTokens();

        search.Frame frame = new search.Frame(3, 2);   // interval [16, 23]
        int frontier = 18;
        for (search.SearchAlgorithm algorithm : new search.SearchAlgorithm[]{new search.BlockSearch(), new search.BlockSearchCharSet()}) {
            search.CandidateRange range = algorithm.search(frame, pattern, tree, new java.util.ArrayDeque<>(), frontier);
            assertNotNull(range, algorithm.getClass().getSimpleName());
            assertTrue(range.startPos() >= frontier, algorithm.getClass().getSimpleName() + " started at " + range.startPos());
            assertEquals(23, range.endPos());
        }
    }

    // hbi-optimizations.md, Section 4: with patterns of at least minPatternLength, the levels it removes are never
    // probed, so the cut layout makes exactly the same probes and answers. Exact filters (no false positives) make
    // the probe counts comparable. Shorter patterns stay exact.
    @Test
    void minPatternLengthRemovesOnlyLevelsTheSearchNeverEnters() {
        Random rnd = new Random(4);
        int w = 2048, tree = 512;
        for (int n : new int[]{1, 2, 3, 5}) {
            for (int minLen : new int[]{6, 10, 16, 40}) {
                for (boolean strides : new boolean[]{false, true}) {
                    HBI full = exactHbi(w, tree, n, 0);
                    HBI cut = exactHbi(w, tree, n, minLen);
                    full.strides = strides;
                    cut.strides = strides;
                    String config = "n=" + n + " minLen=" + minLen + " strides=" + strides;
                    String text = skewed(rnd, 5000);
                    StringBuilder seen = new StringBuilder();
                    for (int i = 0; i < text.length(); i++) {
                        seen.append(text.charAt(i));
                        if (seen.length() >= n) {
                            String symbol = seen.substring(seen.length() - n);
                            full.insert(symbol);
                            cut.insert(symbol);
                        }
                    }
                    int leaf = Integer.highestOneBit(Math.max(n, minLen) - 1) << 1;   // smallest power of two >= it
                    assertEquals(tree / leaf, 1 << (cut.trees.getLast().maxDepth() - 1), config + ": leaf span");
                    assertTrue(cut.trees.getLast().maxDepth() <= full.trees.getLast().maxDepth(), config);

                    int windowStart = Math.max(0, seen.length() - n + 1 - w);
                    for (int q = 0; q < 60; q++) {
                        int r = n + rnd.nextInt(50);
                        int s = windowStart + rnd.nextInt(seen.length() - r - windowStart + 1);
                        String pattern = seen.substring(s, s + r);
                        if (q % 3 == 0) {
                            pattern = pattern.substring(0, r - 1) + "#";   // usually absent
                        }
                        String where = config + ", pattern " + pattern;
                        full.getAllprobes();
                        cut.getAllprobes();
                        List<Integer> want = full.report(new Pattern(pattern, n));
                        int fullProbes = full.getAllprobes();
                        assertEquals(expected(seen, w, n, pattern), want, where);
                        assertEquals(want, cut.report(new Pattern(pattern, n)), where);
                        if (r >= minLen) {
                            assertEquals(fullProbes, cut.getAllprobes(), where + ": probes");
                        }
                    }
                }
            }
        }
        assertThrows(IllegalArgumentException.class, () -> exactHbi(w, tree, 1, -1));
    }

    private static HBI exactHbi(int w, int tree, int n, int minPatternLength) {
        return new HBI(HbiConfiguration.builder()
                .searchAlgorithm(new search.BlockSearch())
                .windowLength(w).treeLength(tree).fpRate(0.2).alphabetSize(1 << 16)
                .estimatorSupplier(() -> new estimators.CSEstimator(64, 3, 1))
                .membershipSupplier(membership.MockMembership::new)
                .pruningPlanSupplier(() -> new search.MostFreqPruning(0.99, 0.2))
                .verifier(new search.VerifierLinearLeafProbe())
                .costFunction(new estimators.CostFunctionDefaultRoot())
                .confidence(0.99)
                .nGram(n)
                .minPatternLength(minPatternLength)
                .build());
    }

    @Test
    void reportRejectsPatternsTheIndexCannotAnswer() {
        HBI hbi = IndexFactory.createHbi(256, 256, 1 << 16, 0.2, 0.99, Utils.MemPolicy.NONE, 3, "cp", 0.05, 0);
        for (char c : "abcdefgh".toCharArray()) hbi.insert(String.valueOf(c).repeat(3));
        assertThrows(IllegalArgumentException.class, () -> hbi.report(new Pattern("ab", 3)));
        assertThrows(IllegalArgumentException.class, () -> hbi.report(new Pattern("", 3)));
        assertThrows(IllegalArgumentException.class, () -> hbi.report(new Pattern("abcd", 2)));
    }

    // Window-relative start positions of pattern among the last w n-gram symbols.
    private static List<Integer> expected(CharSequence seen, int w, int n, String pattern) {
        String s = seen.toString();
        int symbols = s.length() - n + 1;
        int start = Math.max(0, symbols - w);
        List<Integer> hits = new ArrayList<>();
        for (int i = s.indexOf(pattern, start); i >= 0; i = s.indexOf(pattern, i + 1)) {
            hits.add(i - start);
        }
        return hits;
    }

    private static String skewed(Random rnd, int length) {
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            double u = rnd.nextDouble();
            sb.append(u < 0.6 ? 'e' : u < 0.85 ? 't' : (char) ('a' + rnd.nextInt(26)));
        }
        return sb.toString();
    }

    private static ImplicitTree<Membership> tree(int length, int nGram) {
        int levels = 1;
        int span = length;
        while (span % 2 == 0 && span / 2 >= nGram) {
            span /= 2;
            levels++;
        }
        return new ImplicitTree<>(new TreeLayout(levels, span), level -> {
            BloomFilter bf = new BloomFilter();
            bf.init(Math.max(1, length >> level), 0.2);
            return bf;
        }, new membership.KeyPackingService(levels, 16), null);
    }
}
