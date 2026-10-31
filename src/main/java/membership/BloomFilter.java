package membership;

import java.security.SecureRandom;
import java.util.BitSet;

public class BloomFilter implements Membership {

    BitSet filter;
    int[] seeds;                                  // preserved (unused)
    private final byte[] leByteArr = new byte[8]; // preserved (unused)
    private final int btLength = leByteArr.length;// preserved (unused)
    int n;
    double p;
    int m;
    int k;

    @Override public String toString() {
        return String.format("n=%d, m=%d, k=%d, p(actual)=%.4g", n, m, k, p);
    }

    public BloomFilter(){}

    public void init(int n, double p){
        if (n <= 0)           throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1) throw new IllegalArgumentException("p must be in (0,1)");

        double ln2 = Math.log(2);
        this.n = n;
        this.m = (int) Math.ceil(-(n * Math.log(p)) / (ln2 * ln2));
        this.k = (int) Math.max(1, Math.round((m / (double) n) * ln2));
        this.p = Math.pow(1 - Math.exp(-k / (double) m * n), k);

        this.filter = new BitSet(m);
        initCWSeeds();
    }

    public BloomFilter(int n, double p) {
        if (n <= 0)           throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1) throw new IllegalArgumentException("p must be in (0,1)");

        double ln2 = Math.log(2);
        this.n = n;
        this.m = (int) Math.ceil(-(n * Math.log(p)) / (ln2 * ln2));
        this.k = (int) Math.max(1, Math.round((m / (double) n) * ln2));
        this.p = Math.pow(1 - Math.exp(-k / (double) m * n), k);

        this.filter = new BitSet(m);
        initCWSeeds();
    }

    @Override
    public double getFpRate() {
        double rho = (double) filter.cardinality() / m;
        if (rho <= 0) return 0.0;
        if (rho >= 1) return 1.0;
        return Math.pow(rho, k);
    }

    public double getDesignFpRate() { return p; }

    public long estimateDistinct() {
        double rho = (double) filter.cardinality() / m;
        if (rho <= 0) return 0;
        if (rho >= 1) return Long.MAX_VALUE;
        return Math.round(-(m / (double) k) * Math.log(1.0 - rho));
    }

    // Carter-Wegman modulo p = 2^61 - 1
    private static final long P61  = (1L << 61) - 1;
    private static final long MASK = P61;

    private long a1, b1, a2, b2;

    private void initCWSeeds() {
        SecureRandom rng = new SecureRandom();
        do { a1 = rng.nextLong() & MASK; } while (a1 == 0L);
        b1 = rng.nextLong() & MASK;
        do { a2 = rng.nextLong() & MASK; } while (a2 == 0L);
        b2 = rng.nextLong() & MASK;
    }

    private static long fold61(long x) {
        long t = (x & MASK) + (x >>> 61);
        return t >= P61 ? t - P61 : t;
    }

    private static long reduce128(long hi, long lo) {
        long lo_lo = lo & MASK;   // low 61 bits
        long lo_hi = lo >>> 61;   // high 3 bits
        long hi8   = hi << 3;     // * 2^3
        long hi8_lo = hi8 & MASK;
        long hi8_hi = hi >>> 58;  // >> 61 then * 2^3
        long s = lo_lo + lo_hi + hi8_lo + hi8_hi;
        s = (s & MASK) + (s >>> 61);
        return s >= P61 ? s - P61 : s;
    }

    private static long mulAdd61(long h, long a, long x) {
        long lo = h * a;
        long hi = Math.multiplyHigh(h, a);
        long r  = reduce128(hi, lo);
        long y  = r + x;
        return y >= P61 ? y - P61 : y;
    }

    // h1 and h2 are independent Carter-Wegman instances a*x + b mod p over a prehashed 64-bit key x.
    // The bit positions come from the top bits of the hash (toIndex), and structured keys, such as
    // consecutive nodes with small tokens, stay on a lattice under a linear hash, so every key is
    // mixed first. The two-word hash used to be b*a^2 + a*k0 + k1, with the token k1 only added:
    // all tokens of one node then got nearly the same bit positions.
    private long h1(long x) { return mulAdd61(fold61(x), a1, b1); }
    private long h2(long x) { return mulAdd61(fold61(x), a2, b2); }

    private static final long GOLDEN = 0x9E3779B97F4A7C15L;

    // SplitMix64's finalizer, a bijection.
    private static long mix64(long z) {
        z += GOLDEN;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    private static long prehash(long k0, long k1) {
        return mix64(k1 + GOLDEN * k0);
    }

    // Agrees with prehash(k0, k1) on two chunks.
    private static long prehash(long[] chunks) {
        long x = chunks[0];
        for (int i = 1; i < chunks.length; i++) {
            x = mix64(chunks[i] + GOLDEN * x);
        }
        return x;
    }

    // Uniform index and stride.
    private static long umulh(long x, long y) {
        long x0 = x & 0xFFFFFFFFL, x1 = x >>> 32;
        long y0 = y & 0xFFFFFFFFL, y1 = y >>> 32;
        long p11 = x1 * y1;
        long p01 = x0 * y1;
        long p10 = x1 * y0;
        long p00 = x0 * y0;
        long middle = (p10 & 0xFFFFFFFFL) + (p01 & 0xFFFFFFFFL) + (p00 >>> 32);
        return p11 + (p10 >>> 32) + (p01 >>> 32) + (middle >>> 32);
    }

    private int toIndex(long y) {
        if ((m & (m - 1)) == 0) return (int) (y & (m - 1L));
        long yShift = y << 3; // proportional to y / 2^61
        long hi = umulh(yShift, Integer.toUnsignedLong(m));
        return (int) hi;
    }

    private int stride(long y) {
        if (m <= 1) return 0;
        if ((m & (m - 1)) == 0) return (int) ((y & (m - 1L)) | 1L); // odd stride
        long yShift = y << 3;
        long hi = umulh(yShift, Integer.toUnsignedLong(m - 1));
        return 1 + (int) hi; // [1, m-1]
    }

    // Public API.

    // The k positions of a prehashed key, by double hashing.
    private void set(long x) {
        int idx = toIndex(h1(x)), step = stride(h2(x));
        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += step; if (idx >= m) idx -= m;
        }
    }

    private boolean test(long x) {
        int idx = toIndex(h1(x)), step = stride(h2(x));
        for (int i = 0; i < k; i++) {
            if (!filter.get(idx)) return false;
            idx += step; if (idx >= m) idx -= m;
        }
        return true;
    }

    // Single 64-bit key.
    public void insert(long key) { set(mix64(key)); }

    // Single 64-bit key.
    public boolean contains(long key) { return test(mix64(key)); }

    // Two-word key without array allocation.
    public void insert(long k0, long k1) { set(prehash(k0, k1)); }

    // Two-word key without array allocation.
    public boolean contains(long k0, long k1) { return test(prehash(k0, k1)); }

    // Variable-length key as 64-bit chunks.
    // Chunks are hashed in order.
    public void insert(long[] keyChunks) {
        if (keyChunks == null || keyChunks.length == 0)
            throw new IllegalArgumentException("keyChunks must be non-empty");
        if (keyChunks.length == 1) { insert(keyChunks[0]); return; }
        if (keyChunks.length == 2) { insert(keyChunks[0], keyChunks[1]); return; }
        set(prehash(keyChunks));
    }

    public boolean contains(long[] keyChunks) {
        if (keyChunks == null || keyChunks.length == 0)
            throw new IllegalArgumentException("keyChunks must be non-empty");
        if (keyChunks.length == 1) return contains(keyChunks[0]);
        if (keyChunks.length == 2) return contains(keyChunks[0], keyChunks[1]);
        return test(prehash(keyChunks));
    }
}
