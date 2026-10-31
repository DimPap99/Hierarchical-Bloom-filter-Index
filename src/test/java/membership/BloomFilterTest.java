package membership;

import PMIndex.HBI;
import org.junit.jupiter.api.Test;
import search.Pattern;
import utilities.IndexFactory;
import utilities.Utils;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// The filter's work, not just its answers. The two-word hash used to add the token without
// multiplying it, so the tokens of one node got the same bit positions and HBI probed far more than
// it should (docs/hbi-optimizations.md, Section 0). Answers stayed exact, so only tests like these
// see it.
class BloomFilterTest {

    @Test
    void tokensOfOneNodeGetTheirOwnBits() {
        BloomFilter filter = new BloomFilter(1 << 20, 0.2);
        filter.insert(5L, 101L);
        assertTrue(filter.contains(5L, 101L));
        int positives = 0;
        for (long token = 0; token < 10_000; token++) {
            if (token != 101L && filter.contains(5L, token)) positives++;
        }
        assertTrue(positives <= 2, positives + " of 9,999 other tokens of the node test positive");
    }

    @Test
    void lastChunkOfALongKeyMovesItsBits() {
        BloomFilter filter = new BloomFilter(1 << 20, 0.2);
        filter.insert(new long[]{7L, 8L, 9L});
        assertTrue(filter.contains(new long[]{7L, 8L, 9L}));
        int positives = 0;
        for (long last = 0; last < 10_000; last++) {
            if (last != 9L && filter.contains(new long[]{7L, 8L, last})) positives++;
        }
        assertTrue(positives <= 2, positives + " of 9,999 keys that differ in the last chunk test positive");
    }

    // Keys shaped like HBI's in the library API at n = 1: consecutive nodes, small tokens. A linear
    // hash leaves them on a lattice: a*k0 + c*k1 + b gave up to 1.8x the design rate for some seeds.
    @Test
    void falsePositiveRateOfAGridOfKeysIsTheDesignRate() {
        int nodes = 256, tokens = 256;
        for (double p : new double[]{0.05, 0.2, 0.3}) {
            for (int seed = 0; seed < 3; seed++) {
                BloomFilter filter = new BloomFilter(nodes * tokens, p);
                for (long node = 0; node < nodes; node++) {
                    for (long token = 0; token < tokens; token++) filter.insert(node, token);
                }
                long positives = 0;
                for (long node = 0; node < nodes; node++) {
                    for (long token = 0; token < tokens; token++) {
                        assertTrue(filter.contains(node, token), "false negative");
                        if (filter.contains(node, token + tokens)) positives++;
                    }
                }
                assertRate(positives, nodes * tokens, filter, "grid, p = " + p);
            }
        }
    }

    // Plain Carter-Wegman on consecutive keys gave up to 1.8x the design rate.
    @Test
    void falsePositiveRateOfConsecutiveKeysIsTheDesignRate() {
        int n = 1 << 16;
        for (double p : new double[]{0.05, 0.2, 0.3}) {
            for (int seed = 0; seed < 3; seed++) {
                BloomFilter filter = new BloomFilter(n, p);
                for (long key = 0; key < n; key++) filter.insert(key);
                long positives = 0;
                for (long key = 0; key < n; key++) {
                    assertTrue(filter.contains(key), "false negative");
                    if (filter.contains(key + n)) positives++;
                }
                assertRate(positives, n, filter, "consecutive keys, p = " + p);
            }
        }
    }

    // 65,536 tests: one standard deviation is under 2% of the rate at p = 0.05.
    private static void assertRate(long positives, long tested, BloomFilter filter, String what) {
        double rate = positives / (double) tested;
        assertTrue(rate < 1.15 * filter.getDesignFpRate(),
                what + ": false positive rate " + rate + ", designed for " + filter.getDesignFpRate());
    }

    // Murmur ids from the benchmark path, as they are, against the same ids put through a bijective
    // mixer before they reach the index: the hits must be identical and the probe counts close.
    // With the old hash, n = 1 probed 4x and n = 2 35-64x more without the mixer on Gutenberg.
    @Test
    void premixedIdsDoNotChangeTheWork() {
        Random random = new Random(11);
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < (1 << 14); i++) sb.append((char) ('a' + random.nextInt(26)));
        String text = sb.toString();
        int w = 1 << 14;
        for (int n : new int[]{1, 2}) {
            HBI plain = IndexFactory.createHbi(w, w, Math.min((int) Math.pow(26, n), w), 0.2, 0.99,
                    Utils.MemPolicy.NONE, n, "bs");
            HBI mixed = IndexFactory.createHbi(w, w, Math.min((int) Math.pow(26, n), w), 0.2, 0.99,
                    Utils.MemPolicy.NONE, n, "bs");
            for (int i = 0; i + n <= text.length(); i++) {
                String gram = text.substring(i, i + n);
                plain.insert(gram);
                mixed.insertToken(mix64(mixed.keyMapper.mapToLong(gram)));
            }
            long plainProbes = 0, mixedProbes = 0, plainHits = 0, mixedHits = 0;
            plain.getAllprobes();
            mixed.getAllprobes();
            for (int q = 0; q < 200; q++) {
                String query;
                if (q % 2 == 0) {
                    int start = random.nextInt(text.length() - 16);
                    query = text.substring(start, start + 16);
                } else {
                    StringBuilder absent = new StringBuilder();
                    for (int i = 0; i < 16; i++) absent.append((char) ('a' + random.nextInt(26)));
                    query = absent.toString();
                }
                plainHits += plain.report(new Pattern(query, n)).size();
                plainProbes += plain.getAllprobes();
                Pattern p = new Pattern(query, n);
                long[] ids = new long[p.nGramToLong.length];
                for (int i = 0; i < ids.length; i++) ids[i] = mix64(mixed.keyMapper.mapToLong(p.nGramArr[i]));
                p.setTokens(ids);
                mixedHits += mixed.report(p).size();
                mixedProbes += mixed.getAllprobes();
            }
            assertEquals(plainHits, mixedHits, "n = " + n);
            double ratio = plainProbes / (double) mixedProbes;
            assertTrue(ratio < 1.5 && ratio > 1 / 1.5,
                    "n = " + n + ": " + plainProbes + " probes with Murmur ids, " + mixedProbes + " with mixed ids");
        }
    }

    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }
}
