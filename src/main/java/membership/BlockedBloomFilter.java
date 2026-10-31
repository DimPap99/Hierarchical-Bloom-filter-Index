package membership;

import java.security.SecureRandom;
import java.util.SplittableRandom;
import java.util.concurrent.ConcurrentHashMap;
import java.util.random.RandomGenerator;

// Blocked Bloom filter (hbi-optimizations.md, Sections 2 and 3): all k bits of a key lie in one block, so a lookup
// touches one block instead of k random positions of the whole array.
// - k <= 2 (design rate >= ~0.15): 64-bit blocks, one word; a lookup is one load and a mask compare.
// - k >= 3: 512-bit blocks (8 words), which lose much less accuracy when k is larger.
// Blocks receive uneven numbers of keys, so for the same bits the false positive rate is a little higher than a
// standard filter's. init sizes the bits so that the rate equals BloomFilter's design rate for the same (n, p):
// at most 0.5% more bits at p = 0.3, 1.2% at 0.2, 0.5% at 0.1.
// Hashing: the key is prehashed into one word exactly as in BloomFilter, then hashed by one Carter-Wegman instance
// mod 2^61 - 1; its top bits pick the block and its low bits (or, for 512-bit blocks, a mix of it) the bits inside.
// Like any Bloom filter it has no false negatives. Fixed capacity, no deletion.
public class BlockedBloomFilter implements Membership {

    private long[] words;
    private int blocks;
    private int blockBits;        // 64 or 512
    private int k;
    private int n;
    private long m;               // bits
    private double designRate;    // BloomFilter's design rate for (n, p), which this filter matches

    // Bits per key for each design rate p; it depends on p only, and HBI creates many filters with the same p.
    private static final ConcurrentHashMap<Double, Double> BITS_PER_KEY = new ConcurrentHashMap<>();

    private static final long P61 = (1L << 61) - 1;
    private static final long GOLDEN = 0x9E3779B97F4A7C15L;
    private long a, b;

    public BlockedBloomFilter() {}

    @Override
    public void init(int n, double p) {
        init(n, p, new SecureRandom());
    }

    @Override
    public void init(int n, double p, long seed) {
        init(n, p, new SplittableRandom(seed));
    }

    private void init(int n, double p, RandomGenerator rng) {
        if (n <= 0) throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1) throw new IllegalArgumentException("p must be in (0,1)");
        // The same k and design rate as BloomFilter.init (before its rounding of m to whole bits).
        double ln2 = Math.log(2);
        double standardPerKey = -Math.log(p) / (ln2 * ln2);
        this.k = (int) Math.max(1, Math.round(standardPerKey * ln2));
        this.designRate = Math.pow(1 - Math.exp(-k / standardPerKey), k);
        this.n = n;
        this.blockBits = (k <= 2) ? 64 : 512;
        final int kk = k, bb = blockBits;
        final double target = designRate;
        double perKey = BITS_PER_KEY.computeIfAbsent(p, q -> bitsPerKeyFor(target, kk, bb, standardPerKey));
        long bits = (long) Math.ceil(n * perKey);
        long nBlocks = Math.max(1, (bits + blockBits - 1) / blockBits);
        if (nBlocks * (blockBits / 64) > Integer.MAX_VALUE - 8) {
            throw new IllegalArgumentException("filter too large: " + nBlocks + " blocks");
        }
        this.blocks = (int) nBlocks;
        this.m = nBlocks * blockBits;
        this.words = new long[(int) (nBlocks * (blockBits / 64))];
        do { a = rng.nextLong() & P61; } while (a == 0L);
        b = rng.nextLong() & P61;
    }

    // Bits per key at which a filter with blocks of blockBits bits and k bits per key reaches the target rate.
    static double bitsPerKeyFor(double target, int k, int blockBits, double standardPerKey) {
        double lo = standardPerKey, hi = standardPerKey;
        while (expectedRate(hi, k, blockBits) > target) {
            lo = hi;
            hi *= 1.25;
        }
        for (int i = 0; i < 50; i++) {
            double mid = 0.5 * (lo + hi);
            if (expectedRate(mid, k, blockBits) > target) lo = mid; else hi = mid;
        }
        return hi;
    }

    // Expected false positive rate with bitsPerKey bits per key: a block's load is Poisson with mean
    // blockBits / bitsPerKey, and each key sets k bits of its block drawn with replacement (the 64-bit layout;
    // the 512-bit layout draws them the same way).
    static double expectedRate(double bitsPerKey, int k, int blockBits) {
        double lambda = blockBits / bitsPerKey;
        double qk = Math.pow(1.0 - 1.0 / blockBits, k);   // a bit stays clear after one key's k draws
        double sum = 0.0;
        double logLambda = Math.log(lambda);
        double logP = -lambda;             // log of the Poisson probability of the load
        double clear = 1.0;                // qk^load
        int max = (int) (lambda + 12 * Math.sqrt(lambda) + 30);
        for (int load = 0; load <= max; load++) {
            if (load > 0) {
                logP += logLambda - Math.log(load);
                clear *= qk;
            }
            double fill = 1.0 - clear;
            double pass = fill;
            for (int i = 1; i < k; i++) pass *= fill;
            sum += Math.exp(logP) * pass;
        }
        return sum;
    }

    // Prehash, as BloomFilter: SplitMix64's finalizer, a bijection.
    private static long mix64(long z) {
        z += GOLDEN;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    private static long prehash(long k0, long k1) {
        return mix64(k1 + GOLDEN * k0);
    }

    private static long prehash(long[] chunks) {
        long x = chunks[0];
        for (int i = 1; i < chunks.length; i++) {
            x = mix64(chunks[i] + GOLDEN * x);
        }
        return x;
    }

    // a * x + b mod 2^61 - 1 over the prehashed key.
    private long hash(long x) {
        long h = (x & P61) + (x >>> 61);
        if (h >= P61) h -= P61;
        long lo = h * a;
        long hi = Math.multiplyHigh(h, a);
        long s = (lo & P61) + (lo >>> 61) + ((hi << 3) & P61) + (hi >>> 58);
        s = (s & P61) + (s >>> 61);
        if (s >= P61) s -= P61;
        long y = s + b;
        return y >= P61 ? y - P61 : y;
    }

    // Block of a hash value: its top bits, scaled to the number of blocks.
    private int block(long y) {
        return (int) Math.unsignedMultiplyHigh(y << 3, blocks);
    }

    // Bits of a 64-bit block: k <= 2 six-bit fields from the low bits of y (the block uses the top bits).
    private long mask64(long y) {
        long mask = 1L << (y & 63);
        if (k == 2) mask |= 1L << ((y >>> 6) & 63);
        return mask;
    }

    private void set(long x) {
        long y = hash(x);
        int blk = block(y);
        if (blockBits == 64) {
            words[blk] |= mask64(y);
            return;
        }
        int base = blk << 3;
        long z = mix64(y);
        int avail = 64;
        for (int i = 0; i < k; i++) {
            if (avail < 9) { z = mix64(z); avail = 64; }
            int pos = (int) (z & 511);
            z >>>= 9;
            avail -= 9;
            words[base + (pos >>> 6)] |= 1L << pos;
        }
    }

    private boolean test(long x) {
        long y = hash(x);
        int blk = block(y);
        if (blockBits == 64) {
            long mask = mask64(y);
            return (words[blk] & mask) == mask;
        }
        int base = blk << 3;
        long z = mix64(y);
        int avail = 64;
        for (int i = 0; i < k; i++) {
            if (avail < 9) { z = mix64(z); avail = 64; }
            int pos = (int) (z & 511);
            z >>>= 9;
            avail -= 9;
            if ((words[base + (pos >>> 6)] & (1L << pos)) == 0) return false;
        }
        return true;
    }

    @Override public void insert(long key) { set(mix64(key)); }
    @Override public boolean contains(long key) { return test(mix64(key)); }
    @Override public void insert(long k0, long k1) { set(prehash(k0, k1)); }
    @Override public boolean contains(long k0, long k1) { return test(prehash(k0, k1)); }

    @Override
    public void insert(long[] keyChunks) {
        if (keyChunks == null || keyChunks.length == 0)
            throw new IllegalArgumentException("keyChunks must be non-empty");
        if (keyChunks.length == 1) { insert(keyChunks[0]); return; }
        if (keyChunks.length == 2) { insert(keyChunks[0], keyChunks[1]); return; }
        set(prehash(keyChunks));
    }

    @Override
    public boolean contains(long[] keyChunks) {
        if (keyChunks == null || keyChunks.length == 0)
            throw new IllegalArgumentException("keyChunks must be non-empty");
        if (keyChunks.length == 1) return contains(keyChunks[0]);
        if (keyChunks.length == 2) return contains(keyChunks[0], keyChunks[1]);
        return test(prehash(keyChunks));
    }

    // Realized false positive rate from the fill of each block: a key absent from the filter passes when its
    // k bits, drawn in its block, are all set. The standard formula (global fill)^k would underestimate it.
    @Override
    public double getFpRate() {
        int wordsPerBlock = blockBits / 64;
        double sum = 0.0;
        for (int blk = 0; blk < blocks; blk++) {
            int ones = 0;
            for (int w = blk * wordsPerBlock, end = w + wordsPerBlock; w < end; w++) {
                ones += Long.bitCount(words[w]);
            }
            if (ones > 0) sum += Math.pow(ones / (double) blockBits, k);
        }
        return sum / blocks;
    }

    // BloomFilter's design rate for the same (n, p); this filter is sized to match it.
    public double getDesignFpRate() { return designRate; }

    // Distinct keys from the global fill, as BloomFilter; approximate here, since blocks fill unevenly.
    public long estimateDistinct() {
        long ones = 0;
        for (long w : words) ones += Long.bitCount(w);
        double rho = ones / (double) m;
        if (rho <= 0) return 0;
        if (rho >= 1) return Long.MAX_VALUE;
        return Math.round(-(m / (double) k) * Math.log(1.0 - rho));
    }

    public int blockBits() { return blockBits; }
    public int hashCount() { return k; }
    public long bitCount() { return m; }

    @Override public String toString() {
        return String.format("n=%d, m=%d, k=%d, block=%d, design rate=%.4g", n, m, k, blockBits, designRate);
    }
}
