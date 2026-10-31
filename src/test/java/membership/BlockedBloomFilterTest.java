package membership;

import org.junit.jupiter.api.Test;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// BlockedBloomFilter (docs/hbi-optimizations.md, Sections 2 and 3): no false negatives, and sized so that its false
// positive rate is BloomFilter's design rate for the same (n, p), on structured and on random keys.
class BlockedBloomFilterTest {

    private static final double[] RATES = {0.4, 0.3, 0.2, 0.1, 0.05, 0.01};

    @Test
    void blockSizeFollowsK() {
        for (double p : RATES) {
            BlockedBloomFilter f = filter(1 << 16, p);
            assertEquals(f.hashCount() <= 2 ? 64 : 512, f.blockBits(), "p = " + p);
        }
        assertEquals(64, filter(1000, 0.3).blockBits());
        assertEquals(64, filter(1000, 0.2).blockBits());
        assertEquals(512, filter(1000, 0.1).blockBits());
    }

    // The sizing solves expectedRate(bits per key) = design rate, and at HBI's rates it costs at most a few percent
    // more bits than BloomFilter (0.5% at 0.3, 1.2% at 0.2, 0.5% at 0.1).
    @Test
    void sizedForTheDesignRateWithFewExtraBits() {
        int n = 1 << 18;
        for (double p : new double[]{0.3, 0.2, 0.1, 0.05}) {
            BlockedBloomFilter f = filter(n, p);
            BloomFilter standard = new BloomFilter(n, p);
            double perKey = f.bitCount() / (double) n;
            assertEquals(f.getDesignFpRate(), BlockedBloomFilter.expectedRate(perKey, f.hashCount(), f.blockBits()),
                    0.002, "p = " + p);
            double extra = f.bitCount() / (double) standard.m - 1.0;
            assertTrue(extra >= -0.01 && extra < 0.02, "p = " + p + ": " + (100 * extra) + "% more bits");
            assertEquals(standard.getDesignFpRate(), f.getDesignFpRate(), 0.002, "p = " + p);
        }
    }

    // HBI's keys: consecutive nodes with small tokens (the library API's packed characters at n = 1).
    @Test
    void gridKeysHaveNoFalseNegativesAndTheDesignRate() {
        int nodes = 256, tokens = 256;
        for (double p : RATES) {
            for (int seed = 0; seed < 3; seed++) {
                BlockedBloomFilter f = filter(nodes * tokens, p);
                for (long node = 0; node < nodes; node++)
                    for (long token = 0; token < tokens; token++) f.insert(node, token);
                long positives = 0;
                for (long node = 0; node < nodes; node++) {
                    for (long token = 0; token < tokens; token++) {
                        assertTrue(f.contains(node, token), "false negative");
                        if (f.contains(node, token + tokens)) positives++;
                    }
                }
                assertRate(positives, (long) nodes * tokens, f, "grid, p = " + p);
            }
        }
    }

    @Test
    void consecutiveAndRandomKeysHaveTheDesignRate() {
        int n = 1 << 16;
        Random random = new Random(3);
        for (double p : RATES) {
            BlockedBloomFilter f = filter(n, p);
            for (long key = 0; key < n; key++) f.insert(key);
            long positives = 0;
            for (long key = 0; key < n; key++) {
                assertTrue(f.contains(key), "false negative");
                if (f.contains(key + n)) positives++;
            }
            assertRate(positives, n, f, "consecutive keys, p = " + p);

            BlockedBloomFilter g = filter(n, p);
            long[][] keys = new long[n][];
            for (int i = 0; i < n; i++) {
                keys[i] = new long[]{random.nextLong(), random.nextLong(), random.nextLong()};
                g.insert(keys[i]);
            }
            positives = 0;
            for (int i = 0; i < n; i++) {
                assertTrue(g.contains(keys[i]), "false negative");
                if (g.contains(new long[]{random.nextLong(), random.nextLong(), random.nextLong()})) positives++;
            }
            assertRate(positives, n, g, "random three-word keys, p = " + p);
        }
    }

    // getFpRate reads the fill of every block; for a filter at its design load it is close to the measured rate.
    @Test
    void getFpRateTracksTheMeasuredRate() {
        int n = 1 << 16;
        for (double p : new double[]{0.3, 0.2, 0.1, 0.05}) {
            BlockedBloomFilter f = filter(n, p);
            assertEquals(0.0, f.getFpRate());
            for (long key = 0; key < n; key++) f.insert(7L, key);
            long positives = 0;
            for (long key = n; key < 2L * n; key++) if (f.contains(7L, key)) positives++;
            double measured = positives / (double) n;
            assertEquals(measured, f.getFpRate(), 0.1 * measured + 0.003, "p = " + p);
            assertEquals(n, f.estimateDistinct(), 0.05 * n, "p = " + p);
        }
    }

    @Test
    void rejectsBadArguments() {
        assertThrows(IllegalArgumentException.class, () -> filter(0, 0.2));
        assertThrows(IllegalArgumentException.class, () -> filter(10, 1.0));
        assertThrows(IllegalArgumentException.class, () -> new BlockedBloomFilter().insert(new long[0]));
        assertEquals(FilterKind.BLOCKED, FilterKind.parse("Blocked"));
        assertEquals(FilterKind.BLOOM, FilterKind.parse("bloom"));
        assertEquals(FilterKind.BLOCKED, FilterKind.DEFAULT);
        assertEquals(FilterKind.DEFAULT, FilterKind.parse(null));
        assertInstanceOf(BlockedBloomFilter.class, FilterKind.BLOCKED.supplier().get());
        assertThrows(IllegalArgumentException.class, () -> FilterKind.parse("cuckoo"));
    }

    private static BlockedBloomFilter filter(int n, double p) {
        BlockedBloomFilter f = new BlockedBloomFilter();
        f.init(n, p);
        return f;
    }

    // 65,536 tests: one standard deviation is about 4% of the rate at p = 0.01, under 2% at 0.05.
    private static void assertRate(long positives, long tested, BlockedBloomFilter f, String what) {
        double rate = positives / (double) tested;
        double design = f.getDesignFpRate();
        double slack = design < 0.02 ? 1.2 : 1.1;
        assertTrue(rate < slack * design, what + ": false positive rate " + rate + ", designed for " + design);
    }
}
