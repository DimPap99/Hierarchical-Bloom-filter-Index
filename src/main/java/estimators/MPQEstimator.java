package estimators;

import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import java.io.IOException;
import java.util.*;

/**
 * MPQEstimator = Partitioned Min-Hash Quantile sampler.
 *
 * Maintains B buckets. For each bucket, stores the key with the smallest
 * independent "priority" (a deterministic 64-bit pseudo-uniform value).
 * This class ONLY maintains the sample of keys; it does NOT count frequencies.
 *
 * Typical use:
 *   MPQEstimator mpq = new MPQEstimator(256);  // 256 buckets
 *   for (int key : stream) {
 *       mpq.offer(key); // O(1)
 *   }
 *   int[] reps = mpq.getRepresentatives();  // keys to query in CountSketch
 *
 * Rationale:
 * - offer(key): compute bucket = h1(key), priority = h2(key).
 *   If bucket empty or priority < currentMin, replace representative.
 * - Priorities come from SplitMix64 finalizer (Steele et al. 2014), which
 *   is a fast 64-bit mixer with good avalanche properties.
 *
 * Notes:
 * - Smaller priority is "better". Priorities are longs treated as unsigned.
 * - If the same key appears again, its priority is identical; no churn.
 * - Mergeability: per-bucket minima are trivially mergeable (see mergeFrom).
 */
public final class MPQEstimator {

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
    private final int[] repKey;     // representative key per bucket
    private final long[] repPrio;    // representative priority per bucket (lower is better)
    private final boolean[] filled;  // whether a bucket has any representative

    /** Count of non-empty buckets (sample size). */
    private int nonEmpty;

    // ---------------------------------------------------------------------
    // Constructors
    // ---------------------------------------------------------------------

    /**
     * Create with a given number of buckets. Seeds are random.
     */
    public MPQEstimator(int bucketNum) {
        this(bucketNum, 0 /*effectiveSigma unused*/, new Random());
    }

    /**
     * Create with a given number of buckets and RNG for seeds.
     */
    public MPQEstimator(int bucketNum, int effectiveSigma, Random rng) {
        this(bucketNum, effectiveSigma, rng.nextLong(), rng.nextLong());
    }

    /**
     * Create with explicit seeds (fully deterministic).
     */
    public MPQEstimator(int bucketNum, int effectiveSigma, long bucketSeed, long prioritySeed) {
        if (bucketNum <= 0) throw new IllegalArgumentException("bucketNum must be > 0");
        this.bucketNum = bucketNum;
        this.effectiveSigma = Math.max(0, effectiveSigma);
        this.bucketSeed = bucketSeed;
        this.prioritySeed = prioritySeed;

        this.bucketsIsPow2 = (bucketNum & (bucketNum - 1)) == 0;
        this.bucketsMask = bucketsIsPow2 ? (bucketNum - 1) : 0;

        this.repKey = new int[bucketNum];
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
    public void insert(int key) {
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
    public int[] getRepresentatives() {
        int[] out = new int[nonEmpty];
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
     * Merge another MPQEstimator into this one IN-PLACE.
     * Requires the SAME bucket count and SAME seeds (bucket & priority).
     * For each bucket, keep the representative with the smaller priority.
     */
    public void mergeFrom(MPQEstimator other) {
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
    private int bucketIndex(int key) {
        long z = Integer.toUnsignedLong(key) ^ bucketSeed;
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
    private long priorityOf(int key) {
        long z = Integer.toUnsignedLong(key) ^ prioritySeed;
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
    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    // ---------------------------------------------------------------------
    // Helper record for debugging/export
    // ---------------------------------------------------------------------

    /** Representative tuple for a bucket (key, priority, bucketIndex). */
    public static final class Rep {
        public final int key;
        public final long priority;
        public final int bucket;

        public Rep(int key, long priority, int bucket) {
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


    // ---------------------- Experiment helpers (static) ----------------------

    private enum Dist { UNIFORM, ZIPF }

    /** Run a single experiment and print a compact report. */
    private static void runExperiment(
            long streamLen,            // total items in stream
            int alphabetSize,          // keys from {0..alphabetSize-1}
            int buckets,               // MPQ buckets (B)
            Dist dist,                 // UNIFORM or ZIPF
            double zipfS,              // Zipf exponent (only used for ZIPF)
            long seed,                 // RNG seed for reproducibility
            double delta               // confidence for DKW, e.g., 0.05
    ) {
        System.out.println("=== MPQ vs. exact min experiment ===");
        System.out.printf("N=%d, alphabet=%d, B=%d, dist=%s, s=%.3f, seed=%d, delta=%.3f%n",
                streamLen, alphabetSize, buckets, dist, zipfS, seed, delta);

        Random rng = new Random(seed);

        // Build MPQ (deterministic seeds derived from 'seed')
        long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
        long prioSeed   = mix64(seed ^ 0xf00dbabe42aa1357L);
        MPQEstimator mpq = new MPQEstimator(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

        // Exact frequencies
        HashMap<Integer, Integer> freq = new HashMap<>(Math.max(16, alphabetSize * 2));

        // Optional: Zipf CDF precompute for discrete sampling
        double[] zipfCdf = null;
        if (dist == Dist.ZIPF) {
            zipfCdf = buildZipfCdf(alphabetSize, zipfS);
        }
        ZipfDistribution z = makeZipfIntDist(alphabetSize, zipfS, seed);

        // Stream
        for (long t = 0; t < streamLen; t++) {
            int key =  (dist == Dist.UNIFORM)
                    ? rng.nextInt(alphabetSize)
                    : z.sample();



            // exact counts
            freq.merge(key, 1, Integer::sum);

            // MPQ
            mpq.insert(key);
        }

        // Ground truth stats
        int sActive = freq.size();
        int trueMin = Integer.MAX_VALUE;
        ArrayList<Integer> trueMinKeys = new ArrayList<>();
        for (Map.Entry<Integer,Integer> e : freq.entrySet()) {
            int c = e.getValue();
            if (c < trueMin) {
                trueMin = c;
                trueMinKeys.clear();
                trueMinKeys.add(e.getKey());
            } else if (c == trueMin) {
                trueMinKeys.add(e.getKey());
            }
        }

        // MPQ-estimated minimum (min over representatives' exact counts)
        int[] reps = mpq.getRepresentatives();
        int nSample = reps.length;
        int mpqMin = Integer.MAX_VALUE;
        for (int k : reps) {
            int c = freq.getOrDefault(k, 0); // should be >0 in this setup
            if (c < mpqMin) mpqMin = c;
        }
        if (nSample == 0) {
            System.out.println("No non-empty buckets; cannot estimate. Increase N or decrease B.");
            return;
        }

        // Percent error (value error) of the minimum
        double pctErr = percentError(mpqMin, trueMin);

        // DKW rank band for the MIN order statistic
        // The minimum of n i.i.d. samples corresponds to quantile p_hat ≈ 1/n.
        // DKW: sup_x |F_n(x) - F(x)| ≤ eps with prob ≥ 1-δ  ⇒
        // the population quantile lies in [p_hat - eps, p_hat + eps].
        double eps = dkwRankEpsilon(nSample, delta);
        double pHat = 1.0 / nSample;
        double pLo = Math.max(0.0, pHat - eps);
        double pHi = Math.min(1.0, pHat + eps);

        // Turn that rank/quantile band into a VALUE band using the ground-truth distribution
        int[] truthSorted = sortedCounts(freq);
        int loVal = valueAtQuantile(truthSorted, pLo);
        int hiVal = valueAtQuantile(truthSorted, pHi);
        boolean withinDKW = (mpqMin >= loVal && mpqMin <= hiVal);

        double lambda = (double) sActive / buckets;

        // Report
        System.out.printf("Active distinct S=%d, sample size n=%d, lambda=S/B=%.3f%n",
                sActive, nSample, lambda);
        System.out.printf("True min = %d  (|argmin|=%d)%n", trueMin, trueMinKeys.size());
        System.out.printf("MPQ  min = %d%n", mpqMin);
        System.out.printf("Percent error (value): %.3f%%%n", 100.0 * pctErr);
        System.out.printf("DKW rank epsilon (n=%d, delta=%.3f): eps=%.4f ; p_hat=1/n=%.4f ; band=[%.4f, %.4f]%n",
                nSample, delta, eps, pHat, pLo, pHi);
        System.out.printf("DKW -> value band from truth: [%d, %d] ; MPQ min is %s the band%n",
                loVal, hiVal, withinDKW ? "INSIDE" : "OUTSIDE");
    }







    /** Run one MPQ-vs-truth quantile experiment and print a compact report. */
    private static void runExperiment(
            long streamLen,            // total items
            int alphabetSize,          // keys in {0..alphabetSize-1}
            int buckets,               // MPQ buckets B (also used as k for Bottom-k)
            Dist dist,                 // UNIFORM or ZIPF
            double zipfS,              // Zipf exponent if ZIPF
            double p,                  // target quantile in (0,1], e.g., 0.01
            long seed,                 // RNG seed
            double delta               // DKW confidence (1-delta)
    ) {
        System.out.println("=== MPQ & BottomK p-quantile experiment ===");
        System.out.printf("N=%d, alphabet=%d, B=%d, dist=%s, s=%.3f, p=%.5f, seed=%d, delta=%.3f%n",
                streamLen, alphabetSize, buckets, dist, zipfS, p, seed, delta);

        if (!(p > 0.0 && p <= 1.0)) {
            throw new IllegalArgumentException("p must be in (0,1]");
        }

        // ---------------------------
        // Build stream once
        // ---------------------------
        Random rng = new Random(seed);
        double[] zipfCdf = (dist == Dist.ZIPF) ? buildZipfCdf(alphabetSize, zipfS) : null;
        ArrayList<Integer> stream = new ArrayList<>((int)Math.min(streamLen, Integer.MAX_VALUE));
        for (long t = 0; t < streamLen; t++) {
            int key = (dist == Dist.UNIFORM) ? rng.nextInt(alphabetSize)
                    : sampleZipf(rng, zipfCdf); // 0..alphabetSize-1
            stream.add(key);
        }

        // ---------------------------
        // Exact frequencies (truth)
        // ---------------------------
        long t0 = System.currentTimeMillis();
        HashMap<Integer, Integer> freq = new HashMap<>(Math.max(16, alphabetSize * 2));
        for (int key : stream) {
            freq.merge(key, 1, Integer::sum);
        }
        int[] truthSorted = sortedCounts(freq);
        long tTruthMs = System.currentTimeMillis() - t0;

        int Sactive = freq.size();                // distinct keys that actually appeared
        double lambda = (double) Sactive / buckets; // occupancy ratio for MPQ

        // True p-quantile (value)
        int xTrue = valueAtQuantile(truthSorted, p);

        // ---------------------------
        // MPQ sampler
        // ---------------------------
        long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
        long prioSeed   = mix64(seed ^ 0xf00dbabe42aa1357L);
        MPQEstimator mpq = new MPQEstimator(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

        long t1 = System.currentTimeMillis();
        for (int key : stream) {
            mpq.insert(key);  // O(1)
        }
        int[] repsMPQ = mpq.getRepresentatives();
        int nbMPQ = repsMPQ.length;               // sample size for MPQ
        int[] sampleCountsMPQ = new int[nbMPQ];
        for (int i = 0; i < nbMPQ; i++) {
            sampleCountsMPQ[i] = freq.getOrDefault(repsMPQ[i], 0); // exact count
        }
        Arrays.sort(sampleCountsMPQ);
        long tMPQMs = System.currentTimeMillis() - t1;

        // MPQ quantile estimate and DKW band
        int xMPQ = valueAtQuantile(sampleCountsMPQ, p);
        double epsMPQ = dkwRankEpsilon(nbMPQ, delta);
        double pLoMPQ = Math.max(0.0, p - epsMPQ);
        double pHiMPQ = Math.min(1.0, p + epsMPQ);
        int xLoMPQ = valueAtQuantile(truthSorted, pLoMPQ);
        int xHiMPQ = valueAtQuantile(truthSorted, pHiMPQ);
        boolean mpqInside = (xMPQ >= xLoMPQ && xMPQ <= xHiMPQ);

        // ---------------------------
        // Bottom-k sampler (k = buckets)
        // ---------------------------
        // Using deduplicate=false to match "explicit representation" cost model (every occurrence tests),
        // but duplicates do not affect membership because priority is per distinct key.
        BottomKSampler bk = new BottomKSampler(/*k=*/buckets, /*seed=*/seed ^ 0x5eed5eedL, /*deduplicate=*/false);

        long t2 = System.currentTimeMillis();
        for (int key : stream) {
            bk.offer(key);  // O(1) test; O(log k) only on replacement
        }
        int[] repsBK = bk.sampleKeys();
        int nbBK = repsBK.length;                 // actual bottom-k size (<= k, e.g., if Sactive < k)
        int[] sampleCountsBK = new int[nbBK];
        for (int i = 0; i < nbBK; i++) {
            sampleCountsBK[i] = freq.getOrDefault(repsBK[i], 0);
        }
        Arrays.sort(sampleCountsBK);
        long tBKMs = System.currentTimeMillis() - t2;

        // Bottom-k quantile estimate and DKW band
        int xBK = valueAtQuantile(sampleCountsBK, p);
        double epsBK = dkwRankEpsilon(nbBK, delta);
        double pLoBK = Math.max(0.0, p - epsBK);
        double pHiBK = Math.min(1.0, p + epsBK);
        int xLoBK = valueAtQuantile(truthSorted, pLoBK);
        int xHiBK = valueAtQuantile(truthSorted, pHiBK);
        boolean bkInside = (xBK >= xLoBK && xBK <= xHiBK);

        // ---------------------------
        // Print report
        // ---------------------------
        System.out.printf("Active distinct S=%d, lambda=S/B=%.3f%n", Sactive, lambda);
        System.out.printf("Truth build time: %d ms   MPQ time: %d ms   Bottom-k time: %d ms%n",
                tTruthMs, tMPQMs, tBKMs);

        // Truth
        System.out.printf("True  p-quantile x_p      = %d%n", xTrue);

        // MPQ block
        System.out.printf("MPQ   sample size nb      = %d%n", nbMPQ);
        System.out.printf("MPQ   p-quantile \\hat{x}_p = %d   (value error = %.3f%%)%n",
                xMPQ, 100.0 * percentError(xMPQ, xTrue));
        System.out.printf("MPQ   DKW eps(nb=%d, delta=%.3f) = %.5f  ⇒ rank band [%.5f, %.5f]%n",
                nbMPQ, delta, epsMPQ, pLoMPQ, pHiMPQ);
        System.out.printf("MPQ   DKW → value band (truth)   = [%d, %d] ; %s%n",
                xLoMPQ, xHiMPQ, mpqInside ? "INSIDE" : "OUTSIDE");

        // Bottom-k block
        System.out.printf("BK    sample size nb      = %d%n", nbBK);
        System.out.printf("BK    p-quantile \\hat{x}_p = %d   (value error = %.3f%%)%n",
                xBK, 100.0 * percentError(xBK, xTrue));
        System.out.printf("BK    DKW eps(nb=%d, delta=%.3f) = %.5f  ⇒ rank band [%.5f, %.5f]%n",
                nbBK, delta, epsBK, pLoBK, pHiBK);
        System.out.printf("BK    DKW → value band (truth)   = [%d, %d] ; %s%n",
                xLoBK, xHiBK, bkInside ? "INSIDE" : "OUTSIDE");

        // Optional sanity: minimum case for both
        int xTrueMin   = truthSorted[0];
        int xMPQMin    = (nbMPQ > 0) ? sampleCountsMPQ[0] : Integer.MAX_VALUE;
        int xBKMin     = (nbBK  > 0) ? sampleCountsBK[0]  : Integer.MAX_VALUE;

        double pMinHatMPQ = (nbMPQ > 0) ? 1.0 / nbMPQ : 1.0;
        double epsMinMPQ  = dkwRankEpsilon(nbMPQ, delta);
        int xMinLoMPQ     = valueAtQuantile(truthSorted, Math.max(0.0, pMinHatMPQ - epsMinMPQ));
        int xMinHiMPQ     = valueAtQuantile(truthSorted, Math.min(1.0, pMinHatMPQ + epsMinMPQ));
        boolean mpqMinInside = (xMPQMin >= xMinLoMPQ && xMPQMin <= xMinHiMPQ);

        double pMinHatBK = (nbBK > 0) ? 1.0 / nbBK : 1.0;
        double epsMinBK  = dkwRankEpsilon(nbBK, delta);
        int xMinLoBK     = valueAtQuantile(truthSorted, Math.max(0.0, pMinHatBK - epsMinBK));
        int xMinHiBK     = valueAtQuantile(truthSorted, Math.min(1.0, pMinHatBK + epsMinBK));
        boolean bkMinInside = (xBKMin >= xMinLoBK && xBKMin <= xMinHiBK);

        System.out.printf("Sanity (min): true=%d | MPQ=%d in [%d,%d] → %s | BK=%d in [%d,%d] → %s%n",
                xTrueMin,
                xMPQMin, xMinLoMPQ, xMinHiMPQ, mpqMinInside ? "INSIDE" : "OUTSIDE",
                xBKMin,  xMinLoBK,  xMinHiBK,  bkMinInside  ? "INSIDE" : "OUTSIDE");
    }

    /** Percent error: |est - truth| / max(1, truth). */
    private static double percentError(int est, int truth) {
        int denom = Math.max(1, truth);
        return Math.abs(est - truth) / (double) denom;
    }

    /** DKW rank epsilon for confidence 1-delta over n samples. */
    private static double dkwRankEpsilon(int n, double delta) {
        if (n <= 0) return 1.0;
        if (delta <= 0.0 || delta >= 1.0) throw new IllegalArgumentException("delta must be in (0,1)");
        return Math.sqrt(Math.log(2.0 / delta) / (2.0 * n));
    }

    /** Sorted array (ascending) of exact per-key counts. */
    private static int[] sortedCounts(HashMap<Integer,Integer> freq) {
        int[] a = new int[freq.size()];
        int i = 0;
        for (int v : freq.values()) a[i++] = v;
        Arrays.sort(a);
        return a;
    }

    /**
     * Value at quantile p from an ascending array 'a' using the left-continuous definition:
     * return the smallest x such that empirical CDF ≥ p.
     * p=0 returns a[0]; p in (0,1] maps to index ceil(p*m)-1 with m=a.length.
     */
    private static int valueAtQuantile(int[] a, double p) {
        if (a.length == 0) return 0;
        if (p <= 0.0) return a[0];
        if (p >= 1.0) return a[a.length - 1];
        int rank = (int) Math.ceil(p * a.length) - 1;
        if (rank < 0) rank = 0;
        if (rank >= a.length) rank = a.length - 1;
        return a[rank];
    }

    /** Build a discrete Zipf CDF over {0..alphabetSize-1} with exponent s>0. */
    private static double[] buildZipfCdf(int alphabetSize, double s) {
        if (alphabetSize <= 0) throw new IllegalArgumentException("alphabetSize must be > 0");
        if (s <= 0.0) throw new IllegalArgumentException("Zipf exponent s must be > 0");
        double[] w = new double[alphabetSize];
        double H = 0.0;
        for (int i = 1; i <= alphabetSize; i++) {
            double wi = 1.0 / Math.pow(i, s);
            w[i - 1] = wi;
            H += wi;
        }
        double[] cdf = new double[alphabetSize];
        double acc = 0.0;
        for (int i = 0; i < alphabetSize; i++) {
            acc += w[i] / H;
            cdf[i] = acc;
        }
        cdf[alphabetSize - 1] = 1.0;
        return cdf;
    }

    /** Sample from the precomputed Zipf CDF (binary search). Returns 0..alphabet-1. */
    private static int sampleZipf(Random rng, double[] cdf) {
        double u = rng.nextDouble();
        int lo = 0, hi = cdf.length - 1;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (u <= cdf[mid]) hi = mid;
            else lo = mid + 1;
        }
        return lo;
    }

    // Build once, reuse for many samples (cheap to sample)
    static ZipfDistribution makeZipfIntDist(int alphabetSize, double exponent, long seed) {
        if (alphabetSize <= 0) throw new IllegalArgumentException("alphabetSize must be > 0");
        if (exponent <= 0.0)    throw new IllegalArgumentException("exponent must be > 0");
        // ZipfDistribution is defined on {1..N}; we will shift to {0..N-1}.
        RandomGenerator rg = new Well19937c(seed); // deterministic
        return new ZipfDistribution(rg, alphabetSize, exponent);
    }

    // Sample one int in 0..alphabetSize-1
    static int sampleZipfInt(ZipfDistribution dist) {
        return dist.sample() - 1;
    }



    public static void main(String[] args) throws IOException {
        // Defaults (override via args)
        long   N      = 1_000_000L; // stream length
        int    A      = 5000;     // alphabet size (effective Σ)
        int    B      = 512;        // MPQ buckets
        Dist   dist   = Dist.UNIFORM; // or Dist.ZIPF
        double s      = 1.0;        // Zipf exponent if ZIPF
        double p      = 0.01;       // target quantile (e.g., 1%)
        long   seed   = 123456789L; // RNG seed
        double delta  = 0.05;       // DKW (1 - delta) confidence

        // CLI: N A B dist s p seed delta
        // Example: 1000000 10000 256 UNIFORM 1.0 0.01 42 0.05
        if (args.length >= 1) N     = Long.parseLong(args[0]);
        if (args.length >= 2) A     = Integer.parseInt(args[1]);
        if (args.length >= 3) B     = Integer.parseInt(args[2]);
        if (args.length >= 4) dist  = Dist.valueOf(args[3].toUpperCase());
        if (args.length >= 5) s     = Double.parseDouble(args[4]);
        if (args.length >= 6) p     = Double.parseDouble(args[5]);
        if (args.length >= 7) seed  = Long.parseLong(args[6]);
        if (args.length >= 8) delta = Double.parseDouble(args[7]);

        runExperiment(N, A, B, dist, s, p, seed, delta);
    }

}
