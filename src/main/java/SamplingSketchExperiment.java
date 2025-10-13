import estimators.BottomKSampler;
import estimators.HOPS;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

import java.io.IOException;
import java.util.*;

import static estimators.HOPS.mix64;

public class SamplingSketchExperiment {

    // =========================
    //  Experiment configuration
    // =========================

    private enum Dist { UNIFORM, ZIPF }

    // =========================
    //  Core utilities
    // =========================

    /** Percent error on values using max(1, truth) to avoid division by zero. */
    private static double percentValueError(int est, int truth) {
        int denom = Math.max(1, truth);
        return Math.abs(est - truth) / (double) denom;
    }

    /** Empirical CDF F_m(x) from an ascending array 'a': proportion of entries <= x. */
    private static double empiricalCdfLeft(int[] a, int x) {
        int idx = upperBound(a, x) - 1; // last index with a[i] <= x
        if (idx < 0) return 0.0;
        return (idx + 1) / (double) a.length;
    }

    /** Upper bound: first index with a[i] > x in ascending 'a'. */
    private static int upperBound(int[] a, int x) {
        int lo = 0, hi = a.length;
        while (lo < hi) {
            int mid = (lo + hi) >>> 1;
            if (a[mid] <= x) lo = mid + 1;
            else hi = mid;
        }
        return lo;
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

    // =========================
    //  Occupancy math (exact)
    // =========================

    /** Expected non-empty buckets for D distinct keys into B buckets. */
    private static double occupancyExpectation(int D, int B) {
        if (B <= 0) return 0.0;
        return B * (1.0 - Math.pow(1.0 - 1.0 / B, D));
    }

    /** Variance of non-empty buckets (exact balls-into-bins formula). */
    private static double occupancyVariance(int D, int B) {
        if (B <= 0) return 0.0;
        double t1 = Math.pow(1.0 - 1.0 / B, D);
        double t2 = Math.pow(1.0 - 2.0 / B, D);
        double q = 1.0 - t1;
        return B * q * (1.0 - q) + B * (B - 1.0) * (1.0 - 2.0 * t1 + t2 - q * q);
    }

    /**
     * Chebyshev lower bound: with probability at least 1 - delta
     * n_b ≥ E - sqrt(Var/delta).
     */
    private static int occupancyLowerBoundChebyshev(int D, int B, double delta_samp) {
        double mu = occupancyExpectation(D, B);
        double var = occupancyVariance(D, B);
        double nstar = mu - Math.sqrt(var / Math.max(1e-12, delta_samp));
        return (int) Math.floor(Math.max(0.0, nstar));
    }

    /** Required sample size for DKW rank error eps at confidence 1-delta. */
    private static int requiredSampleSizeForDKW(double eps, double delta_q) {
        return (int) Math.ceil(Math.log(2.0 / delta_q) / (2.0 * eps * eps));
    }

    // =========================
    //  Auto-design for B (Chebyshev only, as requested)
    // =========================

    /** Small immutable record for design outputs. */
    private static final class DesignResult {
        final int Bsuggested;    // smallest B meeting the target (by Chebyshev LB)
        final int nReq;          // required sample size for DKW
        final int nLB;           // Chebyshev lower bound at Bsuggested
        final double mu;         // E[n_b] at Bsuggested
        final double var;        // Var[n_b] at Bsuggested
        final boolean impossible; // true if Dhat < nReq

        DesignResult(int Bsuggested, int nReq, int nLB, double mu, double var, boolean impossible) {
            this.Bsuggested = Bsuggested;
            this.nReq = nReq;
            this.nLB = nLB;
            this.mu = mu;
            this.var = var;
            this.impossible = impossible;
        }
    }

    /**
     * Design the minimal number of buckets B such that, with distinct-count guess Dhat,
     * Chebyshev’s lower bound n_b ≥ n_req holds with probability at least 1 - delta_samp.
     * We use doubling + binary search to find the smallest B that works.
     */
    private static DesignResult designBucketsForRankTargetChebyshev(
            int Dhat, double epsTarget, double delta_q, double delta_samp
    ) {
        if (Dhat <= 0) throw new IllegalArgumentException("Dhat must be > 0");
        if (!(epsTarget > 0 && epsTarget < 1)) throw new IllegalArgumentException("epsTarget must be in (0,1)");
        if (!(delta_q > 0 && delta_q < 1)) throw new IllegalArgumentException("delta_q must be in (0,1)");
        if (!(delta_samp > 0 && delta_samp < 1)) throw new IllegalArgumentException("delta_samp must be in (0,1)");

        int nReq = requiredSampleSizeForDKW(epsTarget, delta_q);
        boolean impossible = Dhat < nReq;

        int lo = 1, hi = 1;
        while (occupancyLowerBoundChebyshev(Dhat, hi, delta_samp) < nReq && hi < (1 << 29)) {
            hi <<= 1;
        }
        int best = hi;
        while (lo <= hi) {
            int mid = lo + ((hi - lo) >>> 1);
            int nLB = occupancyLowerBoundChebyshev(Dhat, mid, delta_samp);
            if (nLB >= nReq) {
                best = mid;
                hi = mid - 1;
            } else {
                lo = mid + 1;
            }
        }
        int nLBbest = occupancyLowerBoundChebyshev(Dhat, best, delta_samp);
        double mu = occupancyExpectation(Dhat, best);
        double var = occupancyVariance(Dhat, best);
        return new DesignResult(best, nReq, nLBbest, mu, var, impossible);
    }

    // =========================
    //  One trial (truth + samplers)
    // =========================

    /**
     * Run one MPQ-vs-Bottom-k vs Truth quantile experiment and return detailed results.
     * Reports both value errors and rank errors, plus guarantee checks.
     */
    private static TrialResult runExperiment(
            long streamLen,
            int alphabetSize,
            int buckets,                 // used for both MPQ and Bottom-k unless autoB overrides
            Dist dist,
            double zipfS,
            double p,                    // target quantile in (0,1], for example 0.01
            long seed,
            double delta,                // DKW delta for rank band
            boolean verbose,
            // Auto-design knobs
            boolean autoDesignB,         // if true, compute B from rank target using Chebyshev
            Integer distinctGuess,       // if non-null, use this Dhat for design, else use Sactive
            Double rankEpsTarget,        // absolute rank error target
            Double deltaSamp             // sampling failure probability
    ) {
        if (!(p > 0.0 && p <= 1.0)) {
            throw new IllegalArgumentException("p must be in (0,1]");
        }

        Random rng = new Random(seed);

        // Generate stream once
        ArrayList<Integer> stream = new ArrayList<>((int) Math.min(streamLen, Integer.MAX_VALUE));
        ZipfDistribution z = makeZipfIntDist(alphabetSize, zipfS, seed);
        for (long t = 0; t < streamLen; t++) {
            int key = (dist == Dist.UNIFORM) ? rng.nextInt(alphabetSize) : sampleZipfInt(z);
            stream.add(key);
        }

        // Exact truth
        long tTruth0 = System.currentTimeMillis();
        HashMap<Integer, Integer> freq = new HashMap<>(Math.max(16, alphabetSize * 2));
        for (int key : stream) freq.merge(key, 1, Integer::sum);
        int[] truthSorted = sortedCounts(freq);
        long tTruthMs = System.currentTimeMillis() - tTruth0;

        int Sactive = freq.size();
        int xTrue = valueAtQuantile(truthSorted, p);

        // Auto-design B if requested (Chebyshev)
        int Bused = buckets;
        DesignResult design = null;
        double epsT = (rankEpsTarget != null) ? rankEpsTarget : 0.01;
        double deltaS = (deltaSamp != null) ? deltaSamp : 0.05;
        if (autoDesignB) {
            int Dhat = (distinctGuess != null && distinctGuess > 0) ? distinctGuess : Math.max(1, Sactive);
            design = designBucketsForRankTargetChebyshev(Dhat, epsT, /*delta_q=*/delta, /*delta_samp=*/deltaS);
            Bused = design.Bsuggested;

            if (verbose) {
                System.out.printf("Auto-design Chebyshev: Dhat=%d, epsTarget=%.6f, delta_q=%.3f, delta_samp=%.3f%n",
                        Dhat, epsT, delta, deltaS);
                System.out.printf("Suggested B=%d, n_req=%d, E[n_b]=%.2f, Var=%.2f, LB(n_b)=%d%s%n",
                        design.Bsuggested, design.nReq, design.mu, design.var, design.nLB,
                        design.impossible ? "  (WARNING: Dhat < n_req ⇒ target impossible)" : "");
            }
        }

        // --- MPQ sampler (HOPS) ---
        long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
        long prioSeed   = mix64(seed ^ 0xf00dbabe42aa1357L);
        HOPS mpq = new HOPS(Bused, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

        long tMPQ0 = System.currentTimeMillis();
        for (int key : stream) mpq.insert(key);
        int[] repsMPQ = mpq.getRepresentatives();
        int nbMPQ = repsMPQ.length;

        int[] sampleCountsMPQ = new int[nbMPQ];
        for (int i = 0; i < nbMPQ; i++) sampleCountsMPQ[i] = freq.getOrDefault(repsMPQ[i], 0);
        Arrays.sort(sampleCountsMPQ);
        long tMPQMs = System.currentTimeMillis() - tMPQ0;

        int xMPQ = valueAtQuantile(sampleCountsMPQ, p);
        double epsMPQ = dkwRankEpsilon(nbMPQ, delta);
        double pLoMPQ = Math.max(0.0, p - epsMPQ);
        double pHiMPQ = Math.min(1.0, p + epsMPQ);
        int xLoMPQ = valueAtQuantile(truthSorted, pLoMPQ);
        int xHiMPQ = valueAtQuantile(truthSorted, pHiMPQ);
        boolean mpqInDKWValueBand = (xMPQ >= xLoMPQ && xMPQ <= xHiMPQ);

        double pAchMPQ = empiricalCdfLeft(truthSorted, xMPQ);
        double rankErrMPQ = Math.abs(pAchMPQ - p);

        // --- Bottom-k sampler (k=Bused) ---
        BottomKSampler bk = new BottomKSampler(/*k=*/Bused, /*seed=*/seed ^ 0x5eed5eedL);

        long tBK0 = System.currentTimeMillis();
        for (int key : stream) bk.offer(key);
        int[] repsBK = bk.sampleKeys();
        int nbBK = repsBK.length;
        int[] sampleCountsBK = new int[nbBK];
        for (int i = 0; i < nbBK; i++) sampleCountsBK[i] = freq.getOrDefault(repsBK[i], 0);
        Arrays.sort(sampleCountsBK);
        long tBKMs = System.currentTimeMillis() - tBK0;

        int xBK = valueAtQuantile(sampleCountsBK, p);
        double epsBK = dkwRankEpsilon(nbBK, delta);
        double pLoBK = Math.max(0.0, p - epsBK);
        double pHiBK = Math.min(1.0, p + epsBK);
        int xLoBK = valueAtQuantile(truthSorted, pLoBK);
        int xHiBK = valueAtQuantile(truthSorted, pHiBK);
        boolean bkInDKWValueBand = (xBK >= xLoBK && xBK <= xHiBK);

        double pAchBK = empiricalCdfLeft(truthSorted, xBK);
        double rankErrBK = Math.abs(pAchBK - p);

        // Design checks
        int nLB = -1;
        boolean occLBMet = true;
        if (design != null) {
            nLB = design.nLB;                // Chebyshev LB at suggested B
            occLBMet = nbMPQ >= nLB;
        }

        if (verbose) {
            System.out.printf("Active distinct S=%d, Bused=%d%n", Sactive, Bused);
            System.out.printf("Truth time: %d ms   MPQ time: %d ms   Bottom-k time: %d ms%n",
                    tTruthMs, tMPQMs, tBKMs);

            System.out.printf("True  x_p = %d%n", xTrue);

            System.out.printf("MPQ   nb=%d, x̂_p=%d, value error=%.3f%%, achieved percentile=%.5f, rank error=%.5f%n",
                    nbMPQ, xMPQ, 100.0 * percentValueError(xMPQ, xTrue), pAchMPQ, rankErrMPQ);
            System.out.printf("MPQ   DKW eps=%.5f ⇒ rank band [%.5f, %.5f], value band [%d, %d] ⇒ %s%n",
                    epsMPQ, pLoMPQ, pHiMPQ, xLoMPQ, xHiMPQ, mpqInDKWValueBand ? "INSIDE" : "OUTSIDE");
            if (design != null) {
                System.out.printf("MPQ   occupancy LB check: observed nb=%d, LB(n_b)=%d ⇒ %s%n",
                        nbMPQ, nLB, occLBMet ? "OK" : "LOW");
            }

            System.out.printf("BK    nb=%d, x̂_p=%d, value error=%.3f%%, achieved percentile=%.5f, rank error=%.5f%n",
                    nbBK, xBK, 100.0 * percentValueError(xBK, xTrue), pAchBK, rankErrBK);
            System.out.printf("BK    DKW eps=%.5f ⇒ rank band [%.5f, %.5f], value band [%d, %d] ⇒ %s%n",
                    epsBK, pLoBK, pHiBK, xLoBK, xHiBK, bkInDKWValueBand ? "INSIDE" : "OUTSIDE");
        }

        return new TrialResult(
                // sizes and occupancy
                Sactive, Bused, nbMPQ, nbBK,
                // truth and estimates
                xTrue, xMPQ, xBK,
                // value errors (percentage)
                100.0 * percentValueError(xMPQ, xTrue),
                100.0 * percentValueError(xBK, xTrue),
                // rank metrics
                p, pAchMPQ, pAchBK, rankErrMPQ, rankErrBK,
                // DKW parameters and value bands and checks
                delta, dkwRankEpsilon(nbMPQ, delta), dkwRankEpsilon(nbBK, delta),
                pLoMPQ, pHiMPQ, xLoMPQ, xHiMPQ, mpqInDKWValueBand,
                pLoBK, pHiBK, xLoBK, xHiBK, bkInDKWValueBand,
                // occupancy LB check
                nLB, occLBMet,
                // timings
                tTruthMs, tMPQMs, tBKMs,
                // design info
                design
        );
    }

    // =========================
    //  Trial result and summary
    // =========================

    /** Immutable result for one trial with both value and rank diagnostics and guarantee checks. */
    private static final class TrialResult {
        // Sizes
        final int Sactive;
        final int Bused;
        final int nbMPQ;
        final int nbBK;

        // Truth and estimates
        final int xTrue;
        final int xMPQ;
        final int xBK;

        // Value errors in percent
        final double pctErrMPQ;
        final double pctErrBK;

        // Rank targets and achieved percentiles
        final double pTarget;
        final double pAchMPQ;
        final double pAchBK;
        final double rankErrMPQ;
        final double rankErrBK;

        // DKW configuration and bands
        final double delta;
        final double epsMPQ;
        final double epsBK;

        final double pLoMPQ, pHiMPQ;
        final int xLoMPQ, xHiMPQ;
        final boolean mpqInDKWValueBand;

        final double pLoBK, pHiBK;
        final int xLoBK, xHiBK;
        final boolean bkInDKWValueBand;

        // Occupancy lower bound check
        final int nLB;            // Chebyshev LB at Bused (or -1 if design disabled)
        final boolean occLBMet;

        // Timings in milliseconds
        final long tTruthMs, tMPQMs, tBKMs;

        // Optional design info
        final DesignResult design;

        TrialResult(
                int Sactive, int Bused, int nbMPQ, int nbBK,
                int xTrue, int xMPQ, int xBK,
                double pctErrMPQ, double pctErrBK,
                double pTarget, double pAchMPQ, double pAchBK, double rankErrMPQ, double rankErrBK,
                double delta, double epsMPQ, double epsBK,
                double pLoMPQ, double pHiMPQ, int xLoMPQ, int xHiMPQ, boolean mpqInDKWValueBand,
                double pLoBK, double pHiBK, int xLoBK, int xHiBK, boolean bkInDKWValueBand,
                int nLB, boolean occLBMet,
                long tTruthMs, long tMPQMs, long tBKMs,
                DesignResult design
        ) {
            this.Sactive = Sactive;
            this.Bused = Bused;
            this.nbMPQ = nbMPQ;
            this.nbBK = nbBK;
            this.xTrue = xTrue;
            this.xMPQ = xMPQ;
            this.xBK = xBK;
            this.pctErrMPQ = pctErrMPQ;
            this.pctErrBK = pctErrBK;
            this.pTarget = pTarget;
            this.pAchMPQ = pAchMPQ;
            this.pAchBK = pAchBK;
            this.rankErrMPQ = rankErrMPQ;
            this.rankErrBK = rankErrBK;
            this.delta = delta;
            this.epsMPQ = epsMPQ;
            this.epsBK = epsBK;
            this.pLoMPQ = pLoMPQ;
            this.pHiMPQ = pHiMPQ;
            this.xLoMPQ = xLoMPQ;
            this.xHiMPQ = xHiMPQ;
            this.mpqInDKWValueBand = mpqInDKWValueBand;
            this.pLoBK = pLoBK;
            this.pHiBK = pHiBK;
            this.xLoBK = xLoBK;
            this.xHiBK = xHiBK;
            this.bkInDKWValueBand = bkInDKWValueBand;
            this.nLB = nLB;
            this.occLBMet = occLBMet;
            this.tTruthMs = tTruthMs;
            this.tMPQMs = tMPQMs;
            this.tBKMs = tBKMs;
            this.design = design;
        }
    }

    /** Print averages and guarantee satisfaction rates over many trials. */
    private static void summarize(List<TrialResult> rs, double p, double delta) {
        int T = rs.size();
        if (T == 0) return;

        double avgS = 0, avgB = 0, avgNbMPQ = 0, avgNbBK = 0;
        double avgValErrMPQ = 0, avgValErrBK = 0;
        double avgRankErrMPQ = 0, avgRankErrBK = 0;
        double insideMPQ = 0, insideBK = 0;
        double occOK = 0, occTotal = 0;
        double avgTruthMs = 0, avgMPQMs = 0, avgBKMs = 0;

        for (TrialResult r : rs) {
            avgS += r.Sactive;
            avgB += r.Bused;
            avgNbMPQ += r.nbMPQ;
            avgNbBK += r.nbBK;

            avgValErrMPQ += r.pctErrMPQ;
            avgValErrBK  += r.pctErrBK;

            avgRankErrMPQ += r.rankErrMPQ;
            avgRankErrBK  += r.rankErrBK;

            insideMPQ += r.mpqInDKWValueBand ? 1 : 0;
            insideBK  += r.bkInDKWValueBand ? 1 : 0;

            if (r.nLB >= 0) {
                occTotal += 1;
                occOK += r.occLBMet ? 1 : 0;
            }

            avgTruthMs += r.tTruthMs;
            avgMPQMs   += r.tMPQMs;
            avgBKMs    += r.tBKMs;
        }

        avgS /= T; avgB /= T; avgNbMPQ /= T; avgNbBK /= T;
        avgValErrMPQ /= T; avgValErrBK /= T;
        avgRankErrMPQ /= T; avgRankErrBK /= T;
        double rateMPQ = 100.0 * insideMPQ / T;
        double rateBK  = 100.0 * insideBK  / T;
        double rateOcc = (occTotal > 0) ? 100.0 * occOK / occTotal : Double.NaN;
        avgTruthMs /= T; avgMPQMs /= T; avgBKMs /= T;

        System.out.println("\n=== Aggregate over " + T + " trials ===");
        System.out.printf("Avg Sactive = %.1f, Avg Bused = %.1f%n", avgS, avgB);
        System.out.printf("Avg sample sizes: MPQ nb=%.1f, BK nb=%.1f%n", avgNbMPQ, avgNbBK);

        System.out.printf("Avg VALUE error vs truth x_p:  MPQ=%.3f%%, BK=%.3f%%%n",
                avgValErrMPQ, avgValErrBK);
        System.out.printf("Avg RANK error |p̂−p|:        MPQ=%.5f,  BK=%.5f%n",
                avgRankErrMPQ, avgRankErrBK);

        System.out.printf("Inside DKW→value band (p=%.5f, delta=%.3f): MPQ=%.1f%%, BK=%.1f%%%n",
                p, delta, rateMPQ, rateBK);

        if (!Double.isNaN(rateOcc)) {
            System.out.printf("Occupancy LB satisfied (design check): %.1f%%%n", rateOcc);
        }

        System.out.printf("Avg times (ms): Truth=%.1f, MPQ=%.1f, Bottom-k=%.1f%n",
                avgTruthMs, avgMPQMs, avgBKMs);
    }

    // =========================
    //  Zipf helpers
    // =========================

    static ZipfDistribution makeZipfIntDist(int alphabetSize, double exponent, long seed) {
        if (alphabetSize <= 0) throw new IllegalArgumentException("alphabetSize must be > 0");
        if (exponent <= 0.0)    throw new IllegalArgumentException("exponent must be > 0");
        RandomGenerator rg = new Well19937c(seed); // deterministic
        return new ZipfDistribution(rg, alphabetSize, exponent);
    }

    static int sampleZipfInt(ZipfDistribution dist) {
        return dist.sample() - 1; // map {1..N} -> {0..N-1}
    }

    // =========================
    //  Main: CLI compatible + optional auto-design
    // =========================

    public static void main(String[] args) throws IOException {
        // Original defaults
        long N = 20_000_000L;    // stream length
        int A = 20_000;          // alphabet size
        int B = 16;              // buckets (also k for Bottom-k)
        Dist dist = Dist.ZIPF;   // UNIFORM or ZIPF
        double s = 0.5;          // Zipf exponent if ZIPF
        double p = 0.01;         // target quantile
        long seed = 123456789L;  // RNG seed
        double delta = 0.05;     // DKW delta
        int trials = 1;         // trials
        boolean verbose = true; // per-trial prints

        // New optional knobs with safe defaults
        boolean autoDesignB = true;       // compute B from rank target if true
        Integer distinctGuess = null;      // if null, we will use Sactive from the stream
        Double rankEpsTarget = 0.05;       // if null, default 0.01
        Double deltaSamp = null;           // if null, default 0.05

        // CLI: N A B dist s p seed delta trials verbose autoB distinctGuess rankEpsTarget deltaSamp
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
        if (args.length >= 11) autoDesignB = Boolean.parseBoolean(args[10]);
        if (args.length >= 12) distinctGuess = Integer.parseInt(args[11]);
        if (args.length >= 13) rankEpsTarget = Double.parseDouble(args[12]);
        if (args.length >= 14) deltaSamp = Double.parseDouble(args[13]);

        for (int alphabet = 5000; alphabet <= A; alphabet *= 2) {
            for (int buckets = B; buckets <= 256; buckets *= 2) {
                List<TrialResult> results = new ArrayList<>(trials);

                for (int i = 0; i < 5; i++) {
                    long trialSeed = mix64(seed + i * 0x9E3779B97F4A7C15L);
                    TrialResult r = runExperiment(
                            N, alphabet, buckets, dist, s, p, trialSeed, delta, /*verbose*/ verbose,
                            autoDesignB, distinctGuess, rankEpsTarget, deltaSamp
                    );
                    results.add(r);
                }

                System.out.printf(
                        "%n--- Setting: alphabet=%d, B=%d, dist=%s, s=%.3f, p=%.5f, delta=%.3f, autoB=%s ---\n",
                        alphabet, B, dist, s, p, delta, String.valueOf(autoDesignB)
                );
                if (autoDesignB && !results.isEmpty() && results.get(0).design != null) {
                    DesignResult d = results.get(0).design;
                    System.out.printf("Design summary (trial 1): suggested B=%d, n_req=%d, LB(n_b)=%d, E[n_b]=%.2f\n",
                            d.Bsuggested, d.nReq, d.nLB, d.mu);
                    if (d.impossible) {
                        System.out.println("WARNING: Dhat < n_req, impossible to meet the requested rank target");
                    }
                }

                summarize(results, p, delta);
            }
        }
    }
}
