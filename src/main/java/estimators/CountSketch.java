package estimators;

import java.util.Arrays;
import java.util.Objects;
import java.util.Random;

public final class CountSketch {

    private final int width;
    private final int depth;
    private final long[][] table;      // [depth][width], long counters
    private final long seed;           // global seed
    private final long[] rowSeeds;     // per-row independent seeds
    private final boolean widthIsPow2;
    private final int widthMask;       // valid only if width is power of two


    /**
     * Create with explicit width (columns), depth (rows), and random seed.
     * Prefer power-of-two width for speed and slightly cleaner hashing.
     */
    public CountSketch(int width, int depth, long seed) {
        if (width <= 0 || depth <= 0) throw new IllegalArgumentException("width and depth must be > 0");
        this.width = width;
        this.depth = depth;
        this.seed = seed;
        this.table = new long[depth][width];
        this.widthIsPow2 = isPowerOfTwo(width);
        this.widthMask = widthIsPow2 ? (width - 1) : 0;

        // Derive per-row seeds deterministically from the global seed
        this.rowSeeds = new long[depth];
        long s = seed ^ 0x9E3779B97F4A7C15L;
        for (int i = 0; i < depth; i++) {
            s = mix64(s + i * 0xD1342543DE82EF95L);
            rowSeeds[i] = s;
        }
    }

    /**
     * Create with explicit width/depth and a random seed from java.util.Random.
     */
    public CountSketch(int width, int depth) {
        this(width, depth, new Random().nextLong());
    }

    /**
     * Factory based on target additive-error epsilon (ε) and failure probability delta (δ).
     * width  ≈ ceil_to_pow2( 3 / ε^2 )
     * depth  ≈ ceil( ln(1/δ) )
     */
    public static CountSketch fromEpsDelta(double epsilon, double delta, long seed) {
        if (!(epsilon > 0.0) || !(delta > 0.0) || delta >= 1.0) {
            throw new IllegalArgumentException("epsilon must be > 0 and delta in (0,1)");
        }
        int w = (int)Math.ceil(3.0 / (epsilon * epsilon));
        w = ceilPow2(Math.max(2, w)); // at least 2
        int d = (int)Math.ceil(Math.log(1.0 / delta));
        d = Math.max(1, d);
        return new CountSketch(w, d, seed);
        // Note: constants (3 and 1) are practical; you can scale for tighter bounds.
    }

    public static CountSketch fromEpsDelta(double epsilon, double delta) {
        return fromEpsDelta(epsilon, delta, new Random().nextLong());
    }

    // ---------- Public API ----------

    /**
     * Update the sketch by 'delta' occurrences of 'item'.
     * Overloads provided for int, long, String, and byte[] keys.
     */
    public void update(int item, long delta) {
        updateHash(hashInt(item), delta);
    }

    public void update(long item, long delta) {
        updateHash(hashLong(item), delta);
    }

    public void update(String item, long delta) {
        Objects.requireNonNull(item, "item");
        updateHash(hashString(item), delta);
    }

    public void update(byte[] item, long delta) {
        Objects.requireNonNull(item, "item");
        updateHash(hashBytes(item), delta);
    }

    /**
     * Convenience: increment by +1.
     */
    public void add(int item)    { update(item, 1L); }
    public void add(long item)   { update(item, 1L); }
    public void add(String item) { update(item, 1L); }
    public void add(byte[] item) { update(item, 1L); }

    /**
     * Point query estimate for a key. Unbiased. Returns a long (median of signed row estimates).
     */
    public long estimate(int item)    { return estimateHash(hashInt(item)); }
    public long estimate(long item)   { return estimateHash(hashLong(item)); }
    public long estimate(String item) { return estimateHash(hashString(item)); }
    public long estimate(byte[] item) { return estimateHash(hashBytes(item)); }

    /**
     * Merge another CountSketch in-place (must have same width, depth, and seeds).
     * Useful for parallel aggregation.
     */
    public void mergeInPlace(CountSketch other) {
        requireSameShape(other);
        for (int i = 0; i < depth; i++) {
            long[] a = this.table[i];
            long[] b = other.table[i];
            for (int j = 0; j < width; j++) {
                a[j] += b[j];
            }
        }
    }

    /** Reset all counters to zero. */
    public void clear() {
        for (int i = 0; i < depth; i++) {
            Arrays.fill(table[i], 0L);
        }
    }

    /** @return number of rows (depth). */
    public int depth() { return depth; }

    /** @return number of columns (width). */
    public int width() { return width; }

    /** @return read-only snapshot of internal table (defensive copy). */
    public long[][] snapshot() {
        long[][] cp = new long[depth][];
        for (int i = 0; i < depth; i++) {
            cp[i] = Arrays.copyOf(table[i], width);
        }
        return cp;
    }

    @Override
    public String toString() {
        return "CountSketch{width=" + width + ", depth=" + depth + ", seed=" + seed + "}";
    }

    // update/query on a precomputed item hash

    private void updateHash(long h, long delta) {
        for (int row = 0; row < depth; row++) {
            int bucket = bucket(h, row);
            int sign = sign(h, row);
            table[row][bucket] += (sign * delta);
        }
    }

    private long estimateHash(long h) {
        // Collect per-row signed bucket values
        long[] vals = new long[depth];
        for (int row = 0; row < depth; row++) {
            int bucket = bucket(h, row);
            int sign = sign(h, row);
            vals[row] = sign * table[row][bucket];
        }
        // Median across rows (robust to outliers)
        Arrays.sort(vals);
        int mid = depth >>> 1; // floor(depth/2)
        if ((depth & 1) == 1) {
            return vals[mid];
        } else {
            // even: choose upper median to keep bias-free
            return vals[mid];
        }
    }

    //  Hashing: bucket and sign per row

    private int bucket(long itemHash, int row) {
        long h = mix64(itemHash ^ rowSeeds[row]);
        if (widthIsPow2) {
            return (int)(h & widthMask);
        } else {
            // keep non-negative and modulo width
            long u = h & 0x7fffffffffffffffL;
            return (int)(u % width);
        }
    }

    private int sign(long itemHash, int row) {
        // independent hash for sign: use a disjoint tweak
        long h = mix64((itemHash + 0x9E3779B97F4A7C15L) ^ ~rowSeeds[row]);
        // use highest bit to choose ±1
        return ((h >>> 63) == 0) ? +1 : -1;
    }

    // ---------- Item hashing to 64-bit base hash ----------

    private long hashInt(int x) {
        long z = (long)x ^ seed;
        return mix64(z);
    }

    private long hashLong(long x) {
        long z = x ^ seed;
        return mix64(z);
    }

    private long hashString(String s) {
        // Simple deterministic mixing over UTF-16 chars; sufficient with SplitMix finalization.
        long z = seed ^ 0xD6E8FEB86659FD93L;
        for (int i = 0; i < s.length(); i++) {
            z = mix64(z ^ s.charAt(i));
        }
        return z;
    }

    private long hashBytes(byte[] a) {
        long z = seed ^ 0xA0761D6478BD642FL;
        for (byte b : a) {
            z = mix64(z ^ (long)b);
        }
        return z;
    }


    /**
     * SplitMix64 finalizer (public domain). Good 64-bit avalanche.
     * See: Steele et al., "Fast Splittable Pseudorandom Number Generators" (2014).
     */
    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    private static boolean isPowerOfTwo(int x) {
        return (x & (x - 1)) == 0;
    }

    private static int ceilPow2(int x) {
        int v = x - 1;
        v |= v >>> 1;
        v |= v >>> 2;
        v |= v >>> 4;
        v |= v >>> 8;
        v |= v >>> 16;
        return v + 1;
    }

    private void requireSameShape(CountSketch other) {
        if (this.width != other.width || this.depth != other.depth || this.seed != other.seed) {
            throw new IllegalArgumentException("Shapes/seeds differ; cannot merge.");
        }
        for (int i = 0; i < depth; i++) {
            if (this.rowSeeds[i] != other.rowSeeds[i]) {
                throw new IllegalArgumentException("Row seeds differ; cannot merge.");
            }
        }
    }

}
