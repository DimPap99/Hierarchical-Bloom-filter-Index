package membership;

// KEEPING YOUR ORIGINAL IMPORTS, COMMENTED IF UNUSED
// import java.util.HashSet;
// import java.util.Random;
// import org.apache.commons.codec.digest.MurmurHash3;

import java.security.SecureRandom;
import java.util.BitSet;

public class BloomFilter implements Membership {

    // =========================
    // === YOUR ORIGINAL FIELDS
    // =========================
    BitSet filter;
    int[] seeds;                                  // kept for backward compatibility (unused now)
    private final byte[] leByteArr = new byte[8]; // kept for backward compatibility (unused now)
    private final int btLength = leByteArr.length; // kept for backward compatibility (unused now)
    int n;
    double p;
    int m;
    int k;

    // private static final long SALT1 = 0xC6BC279692B5CC83L; // kept (unused now)
    // private static final long SALT2 = 0xD2B74407B1CE6E93L; // kept (unused now)

    @Override public String toString() {
        return String.format("n=%d, m=%d, k=%d, p(actual)=%.4g", n, m, k, p);
    }

    public BloomFilter(){}

    public void init(int n, double p){
        if (n <= 0)                  throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1)        throw new IllegalArgumentException("p must be in (0,1)");

        double ln2 = Math.log(2);
        // distinct elements
        this.n = n;
        // num of required bits
        this.m = (int) Math.ceil(-(n * Math.log(p)) / (ln2 * ln2));
        // num of hash funcs
        this.k = (int) Math.max(1, Math.round((m / (double) n) * ln2));
        // fp rate
        this.p = Math.pow(1 - Math.exp(-k / (double) m * n), k);

        this.filter = new BitSet(m);

        // initSeeds(); // original Murmur path, preserved but unused
        initCWSeeds();  // Carter–Wegman seeding
    }

    /**
     * https://en.wikipedia.org/wiki/Bloom_filter
     */
    public BloomFilter(int n, double p) {
        if (n <= 0)                  throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1)        throw new IllegalArgumentException("p must be in (0,1)");

        double ln2 = Math.log(2);
        this.n = n;
        this.m = (int) Math.ceil(-(n * Math.log(p)) / (ln2 * ln2));
        this.k = (int) Math.max(1, Math.round((m / (double) n) * ln2));
        this.p = Math.pow(1 - Math.exp(-k / (double) m * n), k);

        this.filter = new BitSet(m);
        // this.  // original stray token kept for fidelity
        // initSeeds(); // original Murmur path, preserved but unused
        initCWSeeds();  // Carter–Wegman seeding
    }

    // ===============================
    // === YOUR ORIGINAL UTIL METHODS
    // ===============================

    // Original Murmur seed init, preserved but unused
    /*
    private void initSeeds() {
        // Use only 2 hf for Kirsch and Mitzenmacher Optimization -> instead of hashes use a running sum
        this.seeds = new int[2];
        HashSet<Integer> seedHashset = new HashSet<Integer>();
        int seedsFound = 0;
        Random rand = new Random();
        while(seedsFound < seeds.length){
            int seed = rand.nextInt();
            if(!seedHashset.contains(seed)){
                seedHashset.add(seed);
                this.seeds[seedsFound] = seed;
                seedsFound++;
            }
        }
    }
    */

    // int → little-endian bytes (kept, unused)
    private void intToLE(int v, byte[] buf) {
        buf[0] = (byte)  v;
        buf[1] = (byte) (v >>> 8);
        buf[2] = (byte) (v >>> 16);
        buf[3] = (byte) (v >>> 24);
    }

    // long → little-endian bytes (kept, unused)
    private void longToLE(long v, byte[] buf) {
        buf[0] = (byte)  v;
        buf[1] = (byte) (v >>> 8);
        buf[2] = (byte) (v >>> 16);
        buf[3] = (byte) (v >>> 24);
        buf[4] = (byte) (v >>> 32);
        buf[5] = (byte) (v >>> 40);
        buf[6] = (byte) (v >>> 48);
        buf[7] = (byte) (v >>> 56);
    }

    // SplitMix64 mixer (kept, unused)
    /*
    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }
    */

    @Override
    public double getFpRate() {
        double rho = (double) filter.cardinality() / m;
        if (rho <= 0) return 0.0;
        if (rho >= 1) return 1.0;
        return Math.pow(rho, k);
    }

    // keep the design target around for logging
    public double getDesignFpRate() {
        return p;
    }

    // estimate distinct inserts so far
    public long estimateDistinct() {
        double rho = (double) filter.cardinality() / m;
        if (rho <= 0) return 0;
        if (rho >= 1) return Long.MAX_VALUE;
        return Math.round(-(m / (double) k) * Math.log(1.0 - rho));
    }

    // ============================================
    // === YOUR ORIGINAL INSERT/CONTAINS VERSIONS
    // === PRESERVED BUT COMMENTED OUT
    // ============================================

    /*
    public void insert(long key) {
        longToLE(key, leByteArr);   // fill the 8-byte buffer

        int h1 = Math.floorMod(
                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[0]), m);
        int h2 = Math.floorMod(
                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[1]), m);

        int idx = h1;
        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += h2;
            if (idx >= m) idx -= m;
        }
    }
    */

    /*
    public boolean contains(long key) {
        longToLE(key, leByteArr);   // reuse same buffer

        int h1 = Math.floorMod(
                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[0]), m);
        int h2 = Math.floorMod(
                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[1]), m);

        int idx = h1;
        for (int i = 0; i < k; i++) {
            if (!filter.get(idx))
                return false;
            idx += h2;
            if (idx >= m) idx -= m;
        }
        return true;                // could still be FP
    }
    */

    /*
    // SplitMix64 + SALT double hashing version, preserved but unused
    public void insert(long key) {
        long h1 = mix64(key ^ SALT1);
        long h2 = mix64(key ^ SALT2);

        int idx = (int) Math.floorMod(h1, m);
        int inc = (int) Math.floorMod(h2, m);

        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += inc;
            if (idx >= m) idx -= m;
        }
    }

    public boolean contains(long key) {
        long h1 = mix64(key ^ SALT1);
        long h2 = mix64(key ^ SALT2);

        int idx = (int) Math.floorMod(h1, m);
        int inc = (int) Math.floorMod(h2, m);

        for (int i = 0; i < k; i++) {
            if (!filter.get(idx)) return false;
            idx += inc;
            if (idx >= m) idx -= m;
        }
        return true; // maybe FP
    }
    */

    // =========================================================
    // === INTEGRATED CARTER–WEGMAN (2^61 - 1) HASHING (ACTIVE)
    // =========================================================

    // ---- Carter–Wegman modulo the Mersenne prime p = 2^61 - 1 ----
    private static final long P61  = (1L << 61) - 1;
    private static final long MASK = P61;

    // Two independent seed pairs (a, b) for h1 and h2
    private long a1, b1, a2, b2;

    /** Draw seeds for the two independent Carter–Wegman functions h1 and h2. */
    private void initCWSeeds() {
        SecureRandom rng = new SecureRandom();

        do { a1 = rng.nextLong() & MASK; } while (a1 == 0L);
        b1 = rng.nextLong() & MASK;

        do { a2 = rng.nextLong() & MASK; } while (a2 == 0L);
        b2 = rng.nextLong() & MASK;
    }

    /** One step fold of any 64-bit pattern to [0, P61]. */
    private static long fold61(long x) {
        long t = (x & MASK) + (x >>> 61);
        return t >= P61 ? t - P61 : t;
    }

    /**
     * Reduce a 128-bit value given as (hi:lo) modulo p using two folds.
     * Uses 2^64 = 2^3 * 2^61 and 2^61 ≡ 1 modulo p.
     */
    private static long reduce128(long hi, long lo) {
        long lo_lo = lo & MASK;   // lower 61 bits of lo
        long lo_hi = lo >>> 61;   // upper 3 bits of lo
        long hi8   = hi << 3;     // hi * 2^3
        long hi8_lo = hi8 & MASK; // lower 61 bits of hi*8
        long hi8_hi = hi >>> 58;  // upper bits of hi*8 beyond 61

        long s = lo_lo + lo_hi + hi8_lo + hi8_hi;
        s = (s & MASK) + (s >>> 61);
        return s >= P61 ? s - P61 : s;
    }

    /** h1(key) = (a1 * key + b1) mod p, in [0, p). */
    private long h1(long key) {
        long x = fold61(key);
        long lo = a1 * x;
        long hi = Math.multiplyHigh(a1, x);
        long ax = reduce128(hi, lo);
        long y  = ax + b1;
        return y >= P61 ? y - P61 : y;
    }

    /** h2(key) = (a2 * key + b2) mod p, in [0, p). */
    private long h2(long key) {
        long x = fold61(key);
        long lo = a2 * x;
        long hi = Math.multiplyHigh(a2, x);
        long ax = reduce128(hi, lo);
        long y  = ax + b2;
        return y >= P61 ? y - P61 : y;
    }

    // -------- unsigned multiply-high helper --------
    // Returns floor((x * y) / 2^64) treating x and y as unsigned 64-bit integers
    private static long umulh(long x, long y) {
        long x0 = x & 0xFFFFFFFFL;
        long x1 = x >>> 32;
        long y0 = y & 0xFFFFFFFFL;
        long y1 = y >>> 32;

        long p11 = x1 * y1;
        long p01 = x0 * y1;
        long p10 = x1 * y0;
        long p00 = x0 * y0;

        long middle = (p10 & 0xFFFFFFFFL) + (p01 & 0xFFFFFFFFL) + (p00 >>> 32);
        return p11 + (p10 >>> 32) + (p01 >>> 32) + (middle >>> 32);
    }

    /**
     * Map y in [0, p) to an exactly uniform index in [0, m).
     * Computes floor(y * m / 2^61) by evaluating high((y << 3) * m) with unsigned semantics.
     * Power-of-two m uses a direct mask.
     */
    private int toIndex(long y) {
        if ((m & (m - 1)) == 0) {
            return (int) (y & (m - 1L));
        }
        long yShift = y << 3; // now proportional to y / 2^61 in 64-bit fixed point
        long hi = umulh(yShift, Integer.toUnsignedLong(m));
        return (int) hi; // in [0, m)
    }

    /**
     * Derive a nonzero stride in [1, m-1] from y in [0, p).
     * Power-of-two m uses an odd stride to guarantee a full cycle.
     */
    private int stride(long y) {
        if (m <= 1) return 0;
        if ((m & (m - 1)) == 0) {
            return (int) ((y & (m - 1L)) | 1L); // odd stride
        }
        long yShift = y << 3;
        long hi = umulh(yShift, Integer.toUnsignedLong(m - 1));
        return 1 + (int) hi; // in [1, m-1]
    }

    // ======================
    // === ACTIVE METHODS ===
    // ======================

    /** Insert a 64-bit key using Carter–Wegman double hashing. */
    public void insert(long key) {
        long y1 = h1(key);
        long y2 = h2(key);

        int idx  = toIndex(y1);
        int step = stride(y2); // nonzero by construction

        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += step;
            if (idx >= m) idx -= m;
        }
    }

    /** Query a 64-bit key using Carter–Wegman double hashing. */
    public boolean contains(long key) {
        long y1 = h1(key);
        long y2 = h2(key);

        int idx  = toIndex(y1);
        int step = stride(y2);

        for (int i = 0; i < k; i++) {
            if (!filter.get(idx)) return false;
            idx += step;
            if (idx >= m) idx -= m;
        }
        return true; // may be a false positive
    }
}
