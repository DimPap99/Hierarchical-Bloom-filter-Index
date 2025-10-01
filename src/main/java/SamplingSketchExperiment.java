import estimators.BottomKSampler;
import estimators.HOPS;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import java.io.IOException;
import java.util.*;

import static estimators.HOPS.mix64;

public class SamplingSketchExperiment {


    //  Experiment helpers

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
        HOPS mpq = new HOPS(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

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
    /** Run one MPQ-vs-BottomK-vs-Truth quantile experiment and (optionally) print a report. */
    private static TrialResult runExperiment(
            long streamLen,            // total items
            int alphabetSize,          // keys in {0..alphabetSize-1}
            int buckets,               // MPQ buckets B (also used as k for Bottom-k)
            Dist dist,                 // UNIFORM or ZIPF
            double zipfS,              // Zipf exponent if ZIPF
            double p,                  // target quantile in (0,1], e.g., 0.01
            long seed,                 // RNG seed
            double delta,              // DKW confidence (1-delta)
            boolean verbose            // whether to print per-trial details
    ) {
        if (!(p > 0.0 && p <= 1.0)) {
            throw new IllegalArgumentException("p must be in (0,1]");
        }

        if (verbose) {
            System.out.println("=== MPQ & BottomK p-quantile experiment ===");
            System.out.printf("N=%d, alphabet=%d, B=%d, dist=%s, s=%.3f, p=%.5f, seed=%d, delta=%.3f%n",
                    streamLen, alphabetSize, buckets, dist, zipfS, p, seed, delta);
        }

        //
        // Build stream once (shared)
        //
        Random rng = new Random(seed);
        double[] zipfCdf = (dist == Dist.ZIPF) ? buildZipfCdf(alphabetSize, zipfS) : null;
        ArrayList<Integer> stream = new ArrayList<>((int)Math.min(streamLen, Integer.MAX_VALUE));
        for (long t = 0; t < streamLen; t++) {
            int key = (dist == Dist.UNIFORM) ? rng.nextInt(alphabetSize)
                    : sampleZipf(rng, zipfCdf); // 0..alphabetSize-1
            stream.add(key);
        }

        //
        // Exact frequencies (truth)
        //
        long t0 = System.currentTimeMillis();
        HashMap<Integer, Integer> freq = new HashMap<>(Math.max(16, alphabetSize * 2));
        for (int key : stream) {
            freq.merge(key, 1, Integer::sum);
        }
        int[] truthSorted = sortedCounts(freq);
        long tTruthMs = System.currentTimeMillis() - t0;

        int Sactive = freq.size();                 // distinct keys that actually appeared
        double lambda = (double) Sactive / buckets;

        // True p-quantile (value)
        int xTrue = valueAtQuantile(truthSorted, p);

        //
        // MPQ sampler
        //
        long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
        long prioSeed   = mix64(seed ^ 0xf00dbabe42aa1357L);
        HOPS mpq = new HOPS(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

        long t1 = System.currentTimeMillis();
        for (int key : stream) {
            mpq.insert(key);  // O(1)
        }
        int[] repsMPQ = mpq.getRepresentatives();
        int nbMPQ = repsMPQ.length;
        int[] sampleCountsMPQ = new int[nbMPQ];
        for (int i = 0; i < nbMPQ; i++) {
            sampleCountsMPQ[i] = freq.getOrDefault(repsMPQ[i], 0);
        }
        Arrays.sort(sampleCountsMPQ);
        long tMPQMs = System.currentTimeMillis() - t1;

        int xMPQ = valueAtQuantile(sampleCountsMPQ, p);
        double epsMPQ = dkwRankEpsilon(nbMPQ, delta);
        double pLoMPQ = Math.max(0.0, p - epsMPQ);
        double pHiMPQ = Math.min(1.0, p + epsMPQ);
        int xLoMPQ = valueAtQuantile(truthSorted, pLoMPQ);
        int xHiMPQ = valueAtQuantile(truthSorted, pHiMPQ);
        boolean mpqInside = (xMPQ >= xLoMPQ && xMPQ <= xHiMPQ);

        //
        // Bottom-k sampler (k = buckets)
        //
        BottomKSampler bk = new BottomKSampler(/*k=*/buckets, /*seed=*/seed ^ 0x5eed5eedL);

        long t2 = System.currentTimeMillis();
        for (int key : stream) {
            bk.offer(key);  // O(1) test; O(log k) on replacement
        }
        int[] repsBK = bk.sampleKeys();
        int nbBK = repsBK.length;
        int[] sampleCountsBK = new int[nbBK];
        for (int i = 0; i < nbBK; i++) {
            sampleCountsBK[i] = freq.getOrDefault(repsBK[i], 0);
        }
        Arrays.sort(sampleCountsBK);
        long tBKMs = System.currentTimeMillis() - t2;

        int xBK = valueAtQuantile(sampleCountsBK, p);
        double epsBK = dkwRankEpsilon(nbBK, delta);
        double pLoBK = Math.max(0.0, p - epsBK);
        double pHiBK = Math.min(1.0, p + epsBK);
        int xLoBK = valueAtQuantile(truthSorted, pLoBK);
        int xHiBK = valueAtQuantile(truthSorted, pHiBK);
        boolean bkInside = (xBK >= xLoBK && xBK <= xHiBK);

        // Optional per-trial prints
        if (verbose) {
            System.out.printf("Active distinct S=%d, lambda=S/B=%.3f%n", Sactive, lambda);
            System.out.printf("Truth build time: %d ms   MPQ time: %d ms   Bottom-k time: %d ms%n",
                    tTruthMs, tMPQMs, tBKMs);

            System.out.printf("True  p-quantile x_p      = %d%n", xTrue);

            System.out.printf("MPQ   sample size nb      = %d%n", nbMPQ);
            System.out.printf("MPQ   p-quantile \\hat{x}_p = %d   (value err = %.3f%%)%n",
                    xMPQ, 100.0 * percentError(xMPQ, xTrue));
            System.out.printf("MPQ   DKW eps(nb=%d, delta=%.3f) = %.5f  ⇒ rank band [%.5f, %.5f]%n",
                    nbMPQ, delta, epsMPQ, pLoMPQ, pHiMPQ);
            System.out.printf("MPQ   DKW → value band (truth)   = [%d, %d] ; %s%n",
                    xLoMPQ, xHiMPQ, mpqInside ? "INSIDE" : "OUTSIDE");

            System.out.printf("BK    sample size nb      = %d%n", nbBK);
            System.out.printf("BK    p-quantile \\hat{x}_p = %d   (value err = %.3f%%)%n",
                    xBK, 100.0 * percentError(xBK, xTrue));
            System.out.printf("BK    DKW eps(nb=%d, delta=%.3f) = %.5f  ⇒ rank band [%.5f, %.5f]%n",
                    nbBK, delta, epsBK, pLoBK, pHiBK);
            System.out.printf("BK    DKW → value band (truth)   = [%d, %d] ; %s%n",
                    xLoBK, xHiBK, bkInside ? "INSIDE" : "OUTSIDE");
        }

        return new TrialResult(
                Sactive, nbMPQ, nbBK, lambda,
                xTrue, xMPQ, xBK,
                percentError(xMPQ, xTrue), percentError(xBK, xTrue),
                epsMPQ, pLoMPQ, pHiMPQ, xLoMPQ, xHiMPQ, mpqInside,
                epsBK, pLoBK, pHiBK, xLoBK, xHiBK, bkInside,
                tTruthMs, tMPQMs, tBKMs
        );
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

    /** Immutable result for one trial. */
    private static final class TrialResult {
        // Size & occupancy
        final int Sactive;
        final int nbMPQ;
        final int nbBK;
        final double lambda;

        // Truth and estimates
        final int xTrue;
        final int xMPQ;
        final int xBK;

        // Percent value errors (|est - truth| / max(1, truth))
        final double pctErrMPQ;
        final double pctErrBK;

        // DKW parameters and value bands
        final double epsMPQ, pLoMPQ, pHiMPQ;
        final int xLoMPQ, xHiMPQ;
        final boolean mpqInside;

        final double epsBK, pLoBK, pHiBK;
        final int xLoBK, xHiBK;
        final boolean bkInside;

        // Timings (milliseconds)
        final long tTruthMs, tMPQMs, tBKMs;

        TrialResult(
                int Sactive, int nbMPQ, int nbBK, double lambda,
                int xTrue, int xMPQ, int xBK,
                double pctErrMPQ, double pctErrBK,
                double epsMPQ, double pLoMPQ, double pHiMPQ, int xLoMPQ, int xHiMPQ, boolean mpqInside,
                double epsBK, double pLoBK, double pHiBK, int xLoBK, int xHiBK, boolean bkInside,
                long tTruthMs, long tMPQMs, long tBKMs) {
            this.Sactive = Sactive;
            this.nbMPQ = nbMPQ;
            this.nbBK = nbBK;
            this.lambda = lambda;
            this.xTrue = xTrue;
            this.xMPQ = xMPQ;
            this.xBK = xBK;
            this.pctErrMPQ = pctErrMPQ;
            this.pctErrBK = pctErrBK;
            this.epsMPQ = epsMPQ;
            this.pLoMPQ = pLoMPQ;
            this.pHiMPQ = pHiMPQ;
            this.xLoMPQ = xLoMPQ;
            this.xHiMPQ = xHiMPQ;
            this.mpqInside = mpqInside;
            this.epsBK = epsBK;
            this.pLoBK = pLoBK;
            this.pHiBK = pHiBK;
            this.xLoBK = xLoBK;
            this.xHiBK = xHiBK;
            this.bkInside = bkInside;
            this.tTruthMs = tTruthMs;
            this.tMPQMs = tMPQMs;
            this.tBKMs = tBKMs;
        }
    }
    /** Print averages over many trials. */
    private static void summarize(List<TrialResult> rs, double p, double delta) {
        int T = rs.size();
        if (T == 0) return;

        double avgS = 0, avgLam = 0;
        double avgNbMPQ = 0, avgNbBK = 0;
        double avgErrMPQ = 0, avgErrBK = 0;
        double insideMPQ = 0, insideBK = 0;
        double avgTruthMs = 0, avgMPQMs = 0, avgBKMs = 0;

        for (TrialResult r : rs) {
            avgS += r.Sactive;
            avgLam += r.lambda;
            avgNbMPQ += r.nbMPQ;
            avgNbBK += r.nbBK;
            avgErrMPQ += r.pctErrMPQ;
            avgErrBK += r.pctErrBK;
            insideMPQ += r.mpqInside ? 1 : 0;
            insideBK += r.bkInside ? 1 : 0;
            avgTruthMs += r.tTruthMs;
            avgMPQMs += r.tMPQMs;
            avgBKMs += r.tBKMs;
        }

        avgS /= T; avgLam /= T;
        avgNbMPQ /= T; avgNbBK /= T;
        avgErrMPQ = 100.0 * (avgErrMPQ / T);
        avgErrBK  = 100.0 * (avgErrBK  / T);
        double rateMPQ = 100.0 * insideMPQ / T;
        double rateBK  = 100.0 * insideBK  / T;
        avgTruthMs /= T; avgMPQMs /= T; avgBKMs /= T;

        System.out.println("\n=== Aggregate over " + T + " trials ===");
        System.out.printf("Avg Sactive = %.1f, Avg lambda=S/B = %.3f%n", avgS, avgLam);
        System.out.printf("Avg sample sizes: MPQ nb=%.1f, BK nb=%.1f%n", avgNbMPQ, avgNbBK);
        System.out.printf("Avg value error: MPQ=%.3f%%, BK=%.3f%%%n", avgErrMPQ, avgErrBK);
        System.out.printf("Fraction inside DKW band (p=%.5f, delta=%.3f): MPQ=%.1f%%, BK=%.1f%%%n",
                p, delta, rateMPQ, rateBK);
        System.out.printf("Avg times (ms): Truth=%.1f, MPQ=%.1f, Bottom-k=%.1f%n",
                avgTruthMs, avgMPQMs, avgBKMs);
    }


    public static void main(String[] args) throws IOException {
        // Defaults (override via args)
        long N = 20_000_000L;  // stream length
        int A = 20000;         // alphabet size
        int B = 16;           // MPQ buckets (also k for Bottom-k)
        Dist dist = Dist.ZIPF; // or Dist.ZIPF
        double s = 0.5;          // Zipf exponent if ZIPF
        double p = 0.01;         // target quantile (e.g., 1%)
        long seed = 123456789L;   // base RNG seed
        double delta = 0.05;         // DKW (1 - delta) confidence
        int trials = 10;            // number of independent trials to average
        boolean verbose = false;      // per-trial prints

        // CLI: N A B dist s p seed delta trials verbose
        // Example: 1000000 10000 256 UNIFORM 1.0 0.01 42 0.05 10 true
        if (args.length >= 1) N = Long.parseLong(args[0]);
        if (args.length >= 2) A = Integer.parseInt(args[1]);
        if (args.length >= 3) B = Integer.parseInt(args[2]);
        if (args.length >= 4) dist = Dist.valueOf(args[3].toUpperCase());
        if (args.length >= 5) s = Double.parseDouble(args[4]);
        if (args.length >= 6) p = Double.parseDouble(args[5]);
        if (args.length >= 7) seed = Long.parseLong(args[6]);
        if (args.length >= 8) delta = Double.parseDouble(args[7]);
        if (args.length >= 9) trials = Integer.parseInt(args[8]);
        if (args.length >= 10) verbose = Boolean.parseBoolean(args[9]);

        for (int alphabet = 5000; alphabet <= A; alphabet *= 2) {
            for (int buckets = B; buckets <= 256; buckets *= 2) {
                List<TrialResult> results = new ArrayList<>(trials);

                trials = 10;

                for (int i = 0; i < trials; i++) {
                    long trialSeed = mix64(seed + i * 0x9E3779B97F4A7C15L); // good spacing
                    TrialResult r = runExperiment(N, alphabet, buckets, dist, s, p, trialSeed, delta, false);
                    results.add(r);
                }
                summarize(results, p, delta);

            }
        }
    }
}
