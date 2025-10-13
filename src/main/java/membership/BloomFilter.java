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

    // =========================================================
    // === Carter–Wegman modulo p = 2^61 - 1 (unchanged core) ===
    // =========================================================
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

    /** Single-word h1/h2 (fast path you already had). */
    private long h1(long key) {
        long x = fold61(key);
        long lo = a1 * x, hi = Math.multiplyHigh(a1, x);
        long ax = reduce128(hi, lo);
        long y  = ax + b1;
        return y >= P61 ? y - P61 : y;
    }
    private long h2(long key) {
        long x = fold61(key);
        long lo = a2 * x, hi = Math.multiplyHigh(a2, x);
        long ax = reduce128(hi, lo);
        long y  = ax + b2;
        return y >= P61 ? y - P61 : y;
    }

    /** NEW: two-word polynomial (no array allocation). */
    private long h1(long k0, long k1) {
        long h = mulAdd61(b1, a1, fold61(k0));
        return mulAdd61(h,  a1, fold61(k1));
    }
    private long h2(long k0, long k1) {
        long h = mulAdd61(b2, a2, fold61(k0));
        return mulAdd61(h,  a2, fold61(k1));
    }

    /** NEW: variable-length polynomial over 64-bit chunks. */
    private long h1(long[] chunks) {
        long h = b1;
        for (int i = 0; i < chunks.length; i++) {
            h = mulAdd61(h, a1, fold61(chunks[i]));
        }
        return h;
    }
    private long h2(long[] chunks) {
        long h = b2;
        for (int i = 0; i < chunks.length; i++) {
            h = mulAdd61(h, a2, fold61(chunks[i]));
        }
        return h;
    }

    // -------- uniform index + stride (same as before) --------
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

    // ======================
    // === PUBLIC API =======
    // ======================

    /** Fast path: single 64-bit key (unchanged). */
    public void insert(long key) {
        long y1 = h1(key);
        long y2 = h2(key);
        int idx = toIndex(y1), step = stride(y2);
        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += step; if (idx >= m) idx -= m;
        }
    }

    /** Fast path: single 64-bit key (unchanged). */
    public boolean contains(long key) {
        long y1 = h1(key);
        long y2 = h2(key);
        int idx = toIndex(y1), step = stride(y2);
        for (int i = 0; i < k; i++) {
            if (!filter.get(idx)) return false;
            idx += step; if (idx >= m) idx -= m;
        }
        return true;
    }

    /** NEW: two-chunk key without allocating an array. */
    public void insert(long k0, long k1) {
        long y1 = h1(k0, k1);
        long y2 = h2(k0, k1);
        int idx = toIndex(y1), step = stride(y2);
        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += step; if (idx >= m) idx -= m;
        }
    }

    /** NEW: two-chunk key without allocating an array. */
    public boolean contains(long k0, long k1) {
        long y1 = h1(k0, k1);
        long y2 = h2(k0, k1);
        int idx = toIndex(y1), step = stride(y2);
        for (int i = 0; i < k; i++) {
            if (!filter.get(idx)) return false;
            idx += step; if (idx >= m) idx -= m;
        }
        return true;
    }

    /**
     * NEW: variable-length key as 64-bit chunks.
     * Order matters: we hash chunks[0], then chunks[1], … using a CW polynomial.
     */
    public void insert(long[] keyChunks) {
        if (keyChunks == null || keyChunks.length == 0)
            throw new IllegalArgumentException("keyChunks must be non-empty");
        if (keyChunks.length == 1) { insert(keyChunks[0]); return; }
        if (keyChunks.length == 2) { insert(keyChunks[0], keyChunks[1]); return; }

        long y1 = h1(keyChunks);
        long y2 = h2(keyChunks);
        int idx = toIndex(y1), step = stride(y2);
        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += step; if (idx >= m) idx -= m;
        }
    }

    public boolean contains(long[] keyChunks) {
        if (keyChunks == null || keyChunks.length == 0)
            throw new IllegalArgumentException("keyChunks must be non-empty");
        if (keyChunks.length == 1) return contains(keyChunks[0]);
        if (keyChunks.length == 2) return contains(keyChunks[0], keyChunks[1]);

        long y1 = h1(keyChunks);
        long y2 = h2(keyChunks);
        int idx = toIndex(y1), step = stride(y2);
        for (int i = 0; i < k; i++) {
            if (!filter.get(idx)) return false;
            idx += step; if (idx >= m) idx -= m;
        }
        return true;
    }
}
