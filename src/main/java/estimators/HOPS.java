package estimators;

import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import java.io.IOException;
import java.util.*;


public final class HOPS {

    /** Number of hash buckets (partitions). */
    private final int bucketNum;

    /** Optional: user’s current estimate of distinct keys (not required to operate). */
    private final int effectiveSigma;

    /** Seeds for bucket hash and priority hash (independent). */
    private final long bucketSeed;
    private final long prioritySeed;

    /** Fast path when bucketNum is a power of two. */
    private final boolean bucketsIsPow2;
    private final int bucketsMask; // valid only if power of two

    /** Per-bucket state: representative key and its priority. */
    private final long[] repKey;    // representative key per bucket
    private final long[] repPrio;    // representative priority per bucket (lower is better)
    private final boolean[] filled;  // whether a bucket has any representative

    /** Count of non-empty buckets (sample size). */
    private int nonEmpty;

    // Constructors

    /**
     * Create with a given number of buckets. Seeds are random.
     */
    public HOPS(int bucketNum) {
        this(bucketNum, 0 , new Random());
    }

    /**
     * Create with a given number of buckets and RNG for seeds.
     */
    public HOPS(int bucketNum, int effectiveSigma, Random rng) {
        this(bucketNum, effectiveSigma, rng.nextLong(), rng.nextLong());
    }

    /**
     * Create with explicit seeds (fully deterministic).
     */
    public HOPS(int bucketNum, int effectiveSigma, long bucketSeed, long prioritySeed) {
        if (bucketNum <= 0) throw new IllegalArgumentException("bucketNum must be > 0");
        this.bucketNum = bucketNum;
        this.effectiveSigma = Math.max(0, effectiveSigma);
        this.bucketSeed = bucketSeed;
        this.prioritySeed = prioritySeed;

        this.bucketsIsPow2 = (bucketNum & (bucketNum - 1)) == 0;
        this.bucketsMask = bucketsIsPow2 ? (bucketNum - 1) : 0;

        this.repKey = new long[bucketNum];
        this.repPrio = new long[bucketNum];
        this.filled = new boolean[bucketNum];
        Arrays.fill(this.repPrio, Long.MAX_VALUE);
        this.nonEmpty = 0;
    }

    // Streaming update

    /**
     * Offer one key from the stream. O(1) time.
     * Keeps the minimum-priority representative per bucket.
     */
    public void insert(long key) {
        final int b = bucketIndex(key);
        final long pr = priorityOf(key);
        if (!filled[b]) {
            filled[b] = true;
            repKey[b] = key;
            repPrio[b] = pr;
            nonEmpty++;
        } else if (unsignedLessThan(pr, repPrio[b])) {
            // Smaller (unsigned) priority wins
            repKey[b] = key;
            repPrio[b] = pr;
        }
        // else: do nothing
    }

    /** @return number of non-empty buckets (size of the current sample). */
    public int sampleSize() {
        return nonEmpty;
    }

    /** @return number of buckets (B). */
    public int buckets() {
        return bucketNum;
    }

    /** @return copy of seeds (bucket, priority) for reproducibility/merging. */
    public long[] seeds() {
        return new long[]{bucketSeed, prioritySeed};
    }

    /** Get a compact array of current representative keys (order arbitrary). */
    public long[] getRepresentatives() {
        long[] out = new long[nonEmpty];
        int j = 0;
        for (int i = 0; i < bucketNum; i++) {
            if (filled[i]) {
                out[j++] = repKey[i];
            }
        }
        return out;
    }

    /**
     * Get representatives and their (unsigned) priorities.
     * Useful for debugging or for stable tie-breaking outside.
     */
    public List<Rep> getRepresentativesWithPriorities() {
        List<Rep> res = new ArrayList<>(nonEmpty);
        for (int i = 0; i < bucketNum; i++) {
            if (filled[i]) {
                res.add(new Rep(repKey[i], repPrio[i], i));
            }
        }
        return res;
    }

    /** Clear the sample (keep same seeds and bucket count). */
    public void clear() {
        Arrays.fill(filled, false);
        Arrays.fill(repPrio, Long.MAX_VALUE);
        // repKey can be left as-is; 'filled' gates visibility
        nonEmpty = 0;
    }

    // Mergeability

    /**
     * Merge another HOPS into this one IN-PLACE.
     * Requires the SAME bucket count and SAME seeds (bucket & priority).
     * For each bucket, keep the representative with the smaller priority.
     */
    public void mergeFrom(HOPS other) {
        Objects.requireNonNull(other, "other");
        if (this.bucketNum != other.bucketNum) {
            throw new IllegalArgumentException("bucketNum mismatch");
        }
        if (this.bucketSeed != other.bucketSeed || this.prioritySeed != other.prioritySeed) {
            throw new IllegalArgumentException("seed mismatch; cannot merge");
        }
        for (int i = 0; i < bucketNum; i++) {
            if (!this.filled[i]) {
                if (other.filled[i]) {
                    this.filled[i] = true;
                    this.repKey[i] = other.repKey[i];
                    this.repPrio[i] = other.repPrio[i];
                    this.nonEmpty++;
                }
            } else if (other.filled[i] && unsignedLessThan(other.repPrio[i], this.repPrio[i])) {
                this.repKey[i] = other.repKey[i];
                this.repPrio[i] = other.repPrio[i];
            }
        }
    }

    //  hashing and comparisons

    /**
     * Bucket index h1(key) ∈ {0..B-1}.
     * Uses SplitMix64 finalizer for mixing; maps to bucket via mask or unsigned mod.
     */
    private int bucketIndex(long key) {
        long z = key ^ bucketSeed;
        long h = mix64(z);
        if (bucketsIsPow2) {
            return (int) (h & bucketsMask);
        } else {
            long u = h & 0x7fffffffffffffffL; // keep non-negative
            return (int) (u % bucketNum);
        }
    }

    /**
     * Priority h2(key): 64-bit pseudo-uniform; lower is better. Deterministic for a given key+seed.
     */
    private long priorityOf(long key) {
        long z = key ^ prioritySeed;
        return mix64(z);
    }

    /**
     * Unsigned comparison for 64-bit priorities.
     * Returns true iff a < b in unsigned 64-bit ordering.
     */
    private static boolean unsignedLessThan(long a, long b) {
        return (a ^ Long.MIN_VALUE) < (b ^ Long.MIN_VALUE);
    }

    /**
     * SplitMix64 finalizer (public domain).
     * Source: Steele et al., "Fast Splittable Pseudorandom Number Generators", 2014.
     */
    public static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }



    /** Representative tuple for a bucket (key, priority, bucketIndex). */
    public static final class Rep {
        public final long key;
        public final long priority;
        public final int bucket;

        public Rep(long key, long priority, int bucket) {
            this.key = key;
            this.priority = priority;
            this.bucket = bucket;
        }

        @Override public String toString() {
            // Show unsigned priority as hex for readability
            return "Rep{bucket=" + bucket +
                    ", key=" + key +
                    ", prio=0x" + Long.toUnsignedString(priority, 16) + "}";
        }
    }




}
