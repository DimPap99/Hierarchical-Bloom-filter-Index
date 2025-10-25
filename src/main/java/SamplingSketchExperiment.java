    import estimators.BottomKSampler;
    import estimators.HOPS;
    import org.apache.commons.math3.distribution.ZipfDistribution;
    import org.apache.commons.math3.random.RandomGenerator;
    import org.apache.commons.math3.random.Well19937c;
    import utilities.CsvUtil;

    import java.nio.file.Files;
    import java.nio.file.Path;
    import java.nio.file.StandardOpenOption;
    import java.io.OutputStream;
    import java.nio.file.StandardOpenOption.*;

    import java.io.IOException;
    import java.util.*;

    import static estimators.HOPS.mix64;

    public class SamplingSketchExperiment {

        // =========================
        //  Experiment configuration
        // =========================

        private enum Dist {UNIFORM, ZIPF}
        private static double EXPONENT;
        // =========================
        //  Core utilities
        // =========================

        /**
         * Percent error on values using max(1, truth) to avoid division by zero.
         */
        private static double percentValueError(int est, int truth) {
            int denom = Math.max(1, truth);
            return Math.abs(est - truth) / (double) denom;
        }
        /** Suite-A/B CSV header in the exact order you’ve been using. */
        // --- replace suiteACsvHeader() with this ---
        private static List<?> suiteACsvHeader() {
            return List.of(
                    "eps_target","delta_q","delta_samp","delta_tot","Dhat",
                    "B_suggested","n_req","E_nb","Var_nb","LB_nb",
                    "Sactive","Bused","nb_MPQ","eps_DKW_MPQ",
                    "p_lo_MPQ","p_hi_MPQ","p_band_width_MPQ",
                    "x_lo_MPQ","x_hi_MPQ","x_band_width_MPQ",
                    "MPQ_xhat","MPQ_value_error_pct","MPQ_achieved_percentile","MPQ_rank_error",
                    "inside_band_MPQ","occupancy_ok","union_success",
                    "ns_per_insert_MPQ","ns_close_MPQ",
                    "nb_BK","eps_DKW_BK","x_lo_BK","x_hi_BK",
                    "BK_xhat","BK_value_error_pct","BK_achieved_percentile","BK_rank_error",
                    "ns_per_insert_BK","ns_close_BK",
                    "truth_x_p","mpq_time_ms","bk_time_ms","exponent"
            );
        }

        // --- replace suiteACsvRow(...) with this ---
        private static List<?> suiteACsvRow(TrialResult r, double exponent) {
            // design fields (when auto-design is ON)
            int    B_suggested = (r.design != null) ? r.design.Bsuggested : r.Bused;
            int    n_req       = (r.design != null) ? r.design.nReq       : -1;
            double E_nb        = (r.design != null) ? r.design.mu         : Double.NaN;
            double Var_nb      = (r.design != null) ? r.design.var        : Double.NaN;
            int    LB_nb       = (r.design != null) ? r.design.nLB        : -1;

            boolean unionSuccess = r.mpqInDKWValueBand && r.occLBMet;
            double  deltaTot     = r.delta + r.deltaSampUsed;
            double  pBandWidth   = r.pHiMPQ - r.pLoMPQ;
            int     xBandWidth   = r.xHiMPQ - r.xLoMPQ;

            return List.of(
                    // knobs
                    r.epsTargetUsed, r.delta, r.deltaSampUsed, deltaTot, r.DhatUsed,
                    // design
                    B_suggested, n_req, E_nb, Var_nb, LB_nb,
                    // sizes
                    r.Sactive, r.Bused, r.nbMPQ, r.epsMPQ,
                    // bands (MPQ)
                    r.pLoMPQ, r.pHiMPQ, pBandWidth,
                    r.xLoMPQ, r.xHiMPQ, xBandWidth,
                    // MPQ estimate & checks
                    r.xMPQ, r.pctErrMPQ, r.pAchMPQ, r.rankErrMPQ,
                    r.mpqInDKWValueBand ? "TRUE" : "FALSE",
                    r.occLBMet ? "TRUE" : "FALSE",
                    unionSuccess ? "TRUE" : "FALSE",
                    // costs (MPQ)
                    r.nsPerInsertMPQ, r.closeNsPerItemMPQ,
                    // BK side
                    r.nbBK, r.epsBK, r.xLoBK, r.xHiBK,
                    r.xBK, r.pctErrBK, r.pAchBK, r.rankErrBK,
                    r.nsPerInsertBK, r.closeNsPerItemBK,
                    // truth + times
                    r.xTrue, (int) r.tMPQMs, (int) r.tBKMs, exponent
            );
        }


        /** Write rows to CSV with append semantics; auto-writes header if file is new/empty. */
        private static void writeSuiteACsv(Path out, List<TrialResult> rs) throws IOException {
            if (rs == null || rs.isEmpty()) return;

            boolean needHeader = true;
            if (Files.exists(out) && Files.size(out) > 0) {
                needHeader = false;
            }

            List<List<?>> rows = new ArrayList<>(rs.size() + 1);
            if (needHeader) rows.add(suiteACsvHeader());
            for (TrialResult r : rs) rows.add(suiteACsvRow(r, EXPONENT));

            try (OutputStream os = Files.newOutputStream(out, StandardOpenOption.CREATE, StandardOpenOption.APPEND)) {
                CsvUtil.writeRows(os, rows, CsvUtil.Config.defaults());
            }
            System.out.printf("CSV: wrote %d trial row(s) to %s (append)\n", rs.size(), out.toAbsolutePath());
        }

        /**
         * Empirical CDF F_m(x) from an ascending array 'a': proportion of entries <= x.
         */
        private static double empiricalCdfLeft(int[] a, int x) {
            int idx = upperBound(a, x) - 1; // last index with a[i] <= x
            if (idx < 0) return 0.0;
            return (idx + 1) / (double) a.length;
        }

        /**
         * Upper bound: first index with a[i] > x in ascending 'a'.
         */
        private static int upperBound(int[] a, int x) {
            int lo = 0, hi = a.length;
            while (lo < hi) {
                int mid = (lo + hi) >>> 1;
                if (a[mid] <= x) lo = mid + 1;
                else hi = mid;
            }
            return lo;
        }

        /**
         * Holds per-item timing results in nanoseconds.
         */
        private static final class PerItemCost {
            final double nsPerInsertMPQ;
            final double nsPerInsertBK;
            final double closeNsPerItemMPQ;
            final double closeNsPerItemBK;

            PerItemCost(double a, double b, double c, double d) {
                this.nsPerInsertMPQ = a;
                this.nsPerInsertBK = b;
                this.closeNsPerItemMPQ = c;
                this.closeNsPerItemBK = d;
            }
        }

        /**
         * Measure per-item costs for MPQ (HOPS) and Bottom-k on the given stream.
         * We time only the insert loops for update cost, and we time the close-time
         * steps separately. We include a short warm-up to let the JIT compile.
         */
        private static PerItemCost measurePerItemCosts(
                ArrayList<Long> stream,
                int buckets,
                long seed
        ) {
            final int n = stream.size();

            // Seeds consistent with your runExperiment
            long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
            long prioSeed = mix64(seed ^ 0xf00dbabe42aa1357L);

            // -------- MPQ / HOPS update cost --------
            // Warm-up
            {
                HOPS warm = new HOPS(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);
                for (long k : stream) warm.insert(k);
            }
            // Timed
            HOPS mpq = new HOPS(buckets, /*effectiveSigma*/ 0, bucketSeed, prioSeed);
            long t0 = System.nanoTime();
            for (long k : stream) mpq.insert(k);
            long t1 = System.nanoTime();
            double nsPerInsertMPQ = (t1 - t0) / (double) n;

            // -------- Bottom-k update cost --------
            // Warm-up
            {
                BottomKSampler warm = new BottomKSampler(/*k=*/buckets, /*seed=*/seed ^ 0x5eed5eedL);
                for (long k : stream) warm.offer(k);
            }
            // Timed
            BottomKSampler bk = new BottomKSampler(/*k=*/buckets, /*seed=*/seed ^ 0x5eed5eedL);
            long t2 = System.nanoTime();
            for (long k : stream) bk.offer(k);
            long t3 = System.nanoTime();
            double nsPerInsertBK = (t3 - t2) / (double) n;

            // -------- Close-time cost (sketch-only and quantile computation) --------
            // MPQ close: get representatives and sort their frequency proxies if available later.
            long t4 = System.nanoTime();
            long[] repsMPQ = mpq.getRepresentatives();
            // Sorting representatives by their counts belongs to the *estimation* step.
            // If you later plug in your frequency sketch, include that lookup here.
            Arrays.sort(repsMPQ); // placeholder to avoid dead-code elimination
            long t5 = System.nanoTime();
            double closeNsPerItemMPQ = (t5 - t4) / (double) n;

            // Bottom-k close: extract and sort keys
            long t6 = System.nanoTime();
            long[] repsBK = bk.sampleKeys();
            Arrays.sort(repsBK);
            long t7 = System.nanoTime();
            double closeNsPerItemBK = (t7 - t6) / (double) n;

            return new PerItemCost(nsPerInsertMPQ, nsPerInsertBK, closeNsPerItemMPQ, closeNsPerItemBK);
        }

        /**
         * DKW rank epsilon for confidence 1-delta over n samples.
         */
        private static double dkwRankEpsilon(int n, double delta) {
            if (n <= 0) return 1.0;
            if (delta <= 0.0 || delta >= 1.0) throw new IllegalArgumentException("delta must be in (0,1)");
            return Math.sqrt(Math.log(2.0 / delta) / (2.0 * n));
        }

        /**
         * Sorted array (ascending) of exact per-key counts.
         */
        private static int[] sortedCounts(HashMap<Long, Integer> freq) {
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

        /**
         * Expected non-empty buckets for D distinct keys into B buckets.
         */
        private static double occupancyExpectation(int D, int B) {
            if (B <= 0) return 0.0;
            return B * (1.0 - Math.pow(1.0 - 1.0 / B, D));
        }

        /**
         * Variance of non-empty buckets (exact balls-into-bins formula).
         */
        private static double occupancyVariance(int D, int B) {
            if (B <= 0) return 0.0;
            double t1 = Math.pow(1.0 - 1.0 / B, D);
            double t2 = Math.pow(1.0 - 2.0 / B, D);
            double q = 1.0 - t1;
            double var = B * q * (1.0 - q) + B * (B - 1.0) * (1.0 - 2.0 * t1 + t2 - q * q);
            return Math.max(0.0, var); // numeric guard
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

        /**
         * Required sample size for DKW rank error eps at confidence 1-delta.
         */
        private static int requiredSampleSizeForDKW(double eps, double delta_q) {
            return (int) Math.ceil(Math.log(2.0 / delta_q) / (2.0 * eps * eps));
        }

        // =========================
        //  Auto-design for B (Chebyshev only, as requested)
        // =========================

        /**
         * Small immutable record for design outputs.
         */
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

            // --- NEW: short-circuit if target is impossible with the available distincts ---
            if (Dhat < nReq) {
                // pick a sane B to avoid collisions (e.g., 2×Dhat) but keep memory bounded
                final int Bcap = 1 << 22; // 4,194,304; tune to your heap
                int Bsafe = Math.min(Bcap, Math.max(16, 2 * Dhat));
                double mu = occupancyExpectation(Dhat, Bsafe);
                double var = Math.max(0.0, occupancyVariance(Dhat, Bsafe));
                int nLB = occupancyLowerBoundChebyshev(Dhat, Bsafe, delta_samp);
                return new DesignResult(Bsafe, nReq, nLB, mu, var, /*impossible=*/true);
            }

            // normal path: search the minimal B meeting the Chebyshev LB
            int lo = 1, hi = 1, best = -1;
            final int Bcap = 1 << 24; // cap to protect heap (16,777,216)
            while (occupancyLowerBoundChebyshev(Dhat, hi, delta_samp) < nReq && hi < Bcap) {
                hi <<= 1;
            }
            if (hi >= Bcap) {
                // failed to meet target within cap: return capped design, flagged as impossible
                double mu = occupancyExpectation(Dhat, Bcap);
                double var = Math.max(0.0, occupancyVariance(Dhat, Bcap));
                int nLB = occupancyLowerBoundChebyshev(Dhat, Bcap, delta_samp);
                return new DesignResult(Bcap, nReq, nLB, mu, var, /*impossible=*/true);
            }
            best = hi;
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
            double mu = occupancyExpectation(Dhat, best);
            double var = Math.max(0.0, occupancyVariance(Dhat, best));
            int nLB = occupancyLowerBoundChebyshev(Dhat, best, delta_samp);
            return new DesignResult(best, nReq, nLB, mu, var, /*impossible=*/false);
        }


        // =========================
        //  One trial (truth + samplers)
        // =========================

        /**
         * Run one MPQ-vs-Bottom-k vs Truth quantile experiment and return detailed results.
         * Reports both value errors and rank errors, plus guarantee checks.
         */
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
                double delta,                // DKW delta for rank band (delta_q)
                boolean verbose,
                // Auto-design knobs
                boolean autoDesignB,         // if true, compute B from rank target using Chebyshev
                Integer distinctGuess,       // if non-null, use this Dhat for design, else use Sactive
                Double rankEpsTarget,        // absolute rank error target (eps_target)
                Double deltaSamp             // sampling failure probability (delta_samp)
        ) {
            if (!(p > 0.0 && p <= 1.0)) {
                throw new IllegalArgumentException("p must be in (0,1]");
            }

            Random rng = new Random(seed);

            // Generate stream once
            ArrayList<Long> stream = new ArrayList<>((int) Math.min(streamLen, Integer.MAX_VALUE));
            ZipfDistribution z = makeZipfIntDist(alphabetSize, zipfS, seed);
            for (long t = 0; t < streamLen; t++) {
                long key = (dist == Dist.UNIFORM) ? rng.nextInt(alphabetSize) : sampleZipfInt(z);
                stream.add(key);
            }

            // Exact truth
            long tTruth0 = System.currentTimeMillis();
            HashMap<Long, Integer> freq = new HashMap<>(Math.max(16, alphabetSize * 2));
            for (long key : stream) freq.merge(key, 1, Integer::sum);
            int[] truthSorted = sortedCounts(freq);
            long tTruthMs = System.currentTimeMillis() - tTruth0;

            int Sactive = freq.size();
            int xTrue = valueAtQuantile(truthSorted, p);

            // Auto-design B if requested (Chebyshev)
            int Bused = buckets;
            DesignResult design = null;
            double epsT = (rankEpsTarget != null) ? rankEpsTarget : 0.01;
            double deltaS = (deltaSamp != null) ? deltaSamp : 0.05;
            int DhatUsed = Math.max(1, (distinctGuess != null && distinctGuess > 0) ? distinctGuess : Sactive);

            if (autoDesignB) {
                design = designBucketsForRankTargetChebyshev(DhatUsed, epsT, /*delta_q=*/delta, /*delta_samp=*/deltaS);
                Bused = design.Bsuggested;

                // Best-achievable rank error with realized sample size at design LB:
                int nAch = Math.max(1, Math.min(DhatUsed, design.nLB)); // cannot exceed Dhat
                double epsAch = Math.sqrt(Math.log(2.0 / delta) / (2.0 * nAch));
                if(verbose){
                    System.out.printf("Auto-design Chebyshev: Dhat=%d, epsTarget=%.6f, delta_q=%.3f, delta_samp=%.3f%n",
                            DhatUsed, epsT, delta, deltaS);
                    System.out.printf("Suggested B=%d, n_req=%d, E[n_b]=%.2f, Var=%.2f, LB(n_b)=%d%s%n",
                            design.Bsuggested, design.nReq, design.mu, design.var, design.nLB,
                            design.impossible ? "  (WARNING: target impossible with Dhat)" : "");
                }

                if (design.impossible) {
                    System.out.printf("NOTE: With Dhat=%d, the tightest DKW rank band at delta_q=%.3f is eps≈%.5f%n",
                            DhatUsed, delta, epsAch);
                }
            }

            // Measure per-item costs for the chosen B (uses same stream and seeds)
            PerItemCost cost = measurePerItemCosts(stream, Bused, seed);

            // --- MPQ sampler (HOPS) ---
            long bucketSeed = mix64(seed ^ 0x1234abcd5678ef00L);
            long prioSeed = mix64(seed ^ 0xf00dbabe42aa1357L);
            HOPS mpq = new HOPS(Bused, /*effectiveSigma*/ 0, bucketSeed, prioSeed);

            long tMPQ0 = System.currentTimeMillis();
            for (long key : stream) mpq.insert(key);
            long[] repsMPQ = mpq.getRepresentatives();
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
            for (long key : stream) bk.offer(key);
            long[] repsBK = bk.sampleKeys();
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
                System.out.printf("Per-item cost (ns): MPQ insert=%.1f, MPQ close=%.3f, BK insert=%.1f, BK close=%.3f%n",
                        cost.nsPerInsertMPQ, cost.closeNsPerItemMPQ, cost.nsPerInsertBK, cost.closeNsPerItemBK);

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
                    // per-item costs (ns)
                    cost.nsPerInsertMPQ, cost.closeNsPerItemMPQ, cost.nsPerInsertBK, cost.closeNsPerItemBK,
                    // deltas and design knob actually used
                    epsT, deltaS, DhatUsed,
                    // design info
                    design
            );
        }

        // =========================
        //  Trial result and summary
        // =========================

        /**
         * Immutable result for one trial with both value/rank diagnostics, guarantees, costs, and design knobs.
         */
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
            final double delta;   // delta_q
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

            // Per-item costs in nanoseconds
            final double nsPerInsertMPQ;
            final double closeNsPerItemMPQ;
            final double nsPerInsertBK;
            final double closeNsPerItemBK;

            // Design knobs actually used (for CSV)
            final double epsTargetUsed;
            final double deltaSampUsed;
            final int DhatUsed;

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
                    double nsPerInsertMPQ, double closeNsPerItemMPQ, double nsPerInsertBK, double closeNsPerItemBK,
                    double epsTargetUsed, double deltaSampUsed, int DhatUsed,
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

                this.nsPerInsertMPQ = nsPerInsertMPQ;
                this.closeNsPerItemMPQ = closeNsPerItemMPQ;
                this.nsPerInsertBK = nsPerInsertBK;
                this.closeNsPerItemBK = closeNsPerItemBK;

                this.epsTargetUsed = epsTargetUsed;
                this.deltaSampUsed = deltaSampUsed;
                this.DhatUsed = DhatUsed;

                this.design = design;
            }
        }


        /** Print averages and guarantee satisfaction rates over many trials. */
        /**
         * Print averages and guarantee satisfaction rates over many trials, including per-item costs.
         */
        /**
         * Print averages and guarantee satisfaction rates over many trials, including per-item costs.
         * NOW ALSO prints the average achieved epsilon (DKW band half-width) for MPQ and BK.
         */
        /**
         * Print averages and guarantee satisfaction rates over many trials, including per-item costs.
         * Now also prints:
         *  - avg achieved percentile p̂ (MPQ and BK)
         *  - avg rank error |p̂ - p| (MPQ and BK) with clearer labels
         *  - avg achieved epsilon from DKW (epsMPQ, epsBK)
         */
        private static void summarize(List<TrialResult> rs, double p, double delta) {
            int T = rs.size();
            if (T == 0) return;

            double avgS = 0, avgB = 0, avgNbMPQ = 0, avgNbBK = 0;

            // value error (in count space, %)
            double avgValErrMPQ = 0, avgValErrBK = 0;

            // rank stuff
            double avgRankErrMPQ = 0, avgRankErrBK = 0;
            double avgPAchMPQ    = 0, avgPAchBK    = 0; // achieved percentiles p̂
            double avgEpsMPQ     = 0, avgEpsBK     = 0; // DKW half-widths ε

            double insideMPQ = 0, insideBK = 0;
            double occOK = 0, occTotal = 0;
            double unionOK = 0, unionTotal = 0;

            double avgTruthMs = 0, avgMPQMs = 0, avgBKMs = 0;

            // per-item costs
            double avgNsInsMPQ = 0, avgNsCloseMPQ = 0, avgNsInsBK = 0, avgNsCloseBK = 0;

            for (TrialResult r : rs) {
                avgS += r.Sactive;
                avgB += r.Bused;
                avgNbMPQ += r.nbMPQ;
                avgNbBK += r.nbBK;

                avgValErrMPQ += r.pctErrMPQ;
                avgValErrBK  += r.pctErrBK;

                avgRankErrMPQ += r.rankErrMPQ;
                avgRankErrBK  += r.rankErrBK;

                avgPAchMPQ += r.pAchMPQ;
                avgPAchBK  += r.pAchBK;

                avgEpsMPQ  += r.epsMPQ;
                avgEpsBK   += r.epsBK;

                insideMPQ += r.mpqInDKWValueBand ? 1 : 0;
                insideBK  += r.bkInDKWValueBand ? 1 : 0;

                if (r.nLB >= 0) {
                    occTotal += 1;
                    occOK += r.occLBMet ? 1 : 0;

                    unionTotal += 1;
                    boolean unionSuccess = r.occLBMet && r.mpqInDKWValueBand;
                    unionOK += unionSuccess ? 1 : 0;
                }

                avgTruthMs += r.tTruthMs;
                avgMPQMs   += r.tMPQMs;
                avgBKMs    += r.tBKMs;

                avgNsInsMPQ    += r.nsPerInsertMPQ;
                avgNsCloseMPQ  += r.closeNsPerItemMPQ;
                avgNsInsBK     += r.nsPerInsertBK;
                avgNsCloseBK   += r.closeNsPerItemBK;
            }

            // divide by T to get averages
            avgS /= T;
            avgB /= T;
            avgNbMPQ /= T;
            avgNbBK /= T;

            avgValErrMPQ /= T;
            avgValErrBK  /= T;

            avgRankErrMPQ /= T;
            avgRankErrBK  /= T;

            avgPAchMPQ /= T;
            avgPAchBK  /= T;

            avgEpsMPQ /= T;
            avgEpsBK  /= T;

            double rateMPQ = 100.0 * insideMPQ / T;
            double rateBK  = 100.0 * insideBK  / T;

            double rateOcc   = (occTotal   > 0) ? 100.0 * occOK   / occTotal   : Double.NaN;
            double rateUnion = (unionTotal > 0) ? 100.0 * unionOK / unionTotal : Double.NaN;

            avgTruthMs /= T;
            avgMPQMs   /= T;
            avgBKMs    /= T;

            avgNsInsMPQ    /= T;
            avgNsCloseMPQ  /= T;
            avgNsInsBK     /= T;
            avgNsCloseBK   /= T;

            System.out.println("\n=== Aggregate over " + T + " trials ===");
            System.out.printf("Target percentile p = %.5f (delta_q = %.3f)\n", p, delta);

            System.out.printf("Avg Sactive = %.1f, Avg Bused = %.1f\n", avgS, avgB);
            System.out.printf("Avg sample sizes: MPQ nb=%.1f, BK nb=%.1f\n", avgNbMPQ, avgNbBK);

            // Value error (count space)
            System.out.printf("Avg VALUE error vs truth x_p (%%):  MPQ=%.3f%%, BK=%.3f%%\n",
                    avgValErrMPQ, avgValErrBK);

            // Achieved percentile and rank error
            System.out.printf("Avg achieved percentile p̂:        MPQ=%.5f, BK=%.5f\n",
                    avgPAchMPQ, avgPAchBK);
            System.out.printf("Avg RANK error |p̂−p|:             MPQ=%.5f, BK=%.5f\n",
                    avgRankErrMPQ, avgRankErrBK);

            // Certified epsilon from DKW given nb per trial
            System.out.printf("Avg certified DKW ε (half-width):  MPQ=%.5f, BK=%.5f\n",
                    avgEpsMPQ, avgEpsBK);

            // Band success
            System.out.printf("Inside DKW→value band:             MPQ=%.1f%% of trials, BK=%.1f%% of trials\n",
                    rateMPQ, rateBK);

            if (!Double.isNaN(rateOcc)) {
                System.out.printf("Occupancy LB satisfied:            %.1f%% of trials\n", rateOcc);
            }
            if (!Double.isNaN(rateUnion)) {
                System.out.printf("Union success (occ LB && band):    %.1f%% of trials\n", rateUnion);
            }

            System.out.printf("Avg times (ms): Truth=%.1f, MPQ=%.1f, Bottom-k=%.1f\n",
                    avgTruthMs, avgMPQMs, avgBKMs);

            System.out.printf("Avg per-item cost (ns): MPQ insert=%.1f, MPQ close=%.3f, BK insert=%.1f, BK close=%.3f\n",
                    avgNsInsMPQ, avgNsCloseMPQ, avgNsInsBK, avgNsCloseBK);
        }


        /**
         * Print Suite-A CSV header (exactly the columns you asked for).
         */
        private static void printSuiteACsvHeader() {
            System.out.println(
                    "eps_target,delta_q,delta_samp,Dhat,B_suggested,n_req,E_nb,Var_nb,LB_nb," +
                            "Sactive,Bused,nb_MPQ,eps_DKW_MPQ,p_lo_MPQ,p_hi_MPQ,x_lo_MPQ,x_hi_MPQ," +
                            "MPQ_xhat,MPQ_value_error_pct,MPQ_achieved_percentile,MPQ_rank_error," +
                            "inside_band_MPQ,occupancy_ok,ns_per_insert_MPQ,ns_close_MPQ," +
                            "nb_BK,eps_DKW_BK,x_lo_BK,x_hi_BK,BK_xhat,BK_value_error_pct," +
                            "BK_achieved_percentile,BK_rank_error,ns_per_insert_BK,ns_close_BK," +
                            "truth_x_p,mpq_time_ms,bk_time_ms"
            );
        }

        /**
         * Print one Suite-A CSV row for a single trial result.
         */
        private static void printSuiteACsvRow(TrialResult r) {
            String inside = r.mpqInDKWValueBand ? "TRUE" : "FALSE";
            String occOK = r.occLBMet ? "TRUE" : "FALSE";

            // Design fields (present when auto-design is enabled)
            int B_suggested = (r.design != null) ? r.design.Bsuggested : r.Bused;
            int n_req = (r.design != null) ? r.design.nReq : -1;
            double E_nb = (r.design != null) ? r.design.mu : Double.NaN;
            double Var_nb = (r.design != null) ? r.design.var : Double.NaN;
            int LB_nb = (r.design != null) ? r.design.nLB : -1;

            System.out.printf(
                    Locale.US,
                    "%.3f,%.2f,%.3f,%d,%d,%d,%.2f,%.2f,%d," +   // eps_target, delta_q, delta_samp, Dhat, B_suggested, n_req, E_nb, Var_nb, LB_nb
                            "%d,%d,%d,%.5f,%.5f,%.5f,%d,%d," +          // Sactive, Bused, nb_MPQ, eps_DKW_MPQ, p_lo_MPQ, p_hi_MPQ, x_lo_MPQ, x_hi_MPQ
                            "%d,%.3f,%.5f,%.5f,%s,%s,%.3f,%.3f," +      // MPQ_xhat, MPQ_value_error_pct, MPQ_achieved_percentile, MPQ_rank_error, inside_band_MPQ, occupancy_ok, ns_per_insert_MPQ, ns_close_MPQ
                            "%d,%.5f,%d,%d,%d,%.3f,%.5f,%.5f,%.3f,%.3f," + // nb_BK, eps_DKW_BK, x_lo_BK, x_hi_BK, BK_xhat, BK_value_error_pct, BK_achieved_percentile, BK_rank_error, ns_per_insert_BK, ns_close_BK
                            "%d,%d,%d%n",                                // truth_x_p, mpq_time_ms, bk_time_ms
                    r.epsTargetUsed, r.delta, r.deltaSampUsed, r.DhatUsed, B_suggested, n_req, E_nb, Var_nb, LB_nb,
                    r.Sactive, r.Bused, r.nbMPQ, r.epsMPQ, r.pLoMPQ, r.pHiMPQ, r.xLoMPQ, r.xHiMPQ,
                    r.xMPQ, r.pctErrMPQ, r.pAchMPQ, r.rankErrMPQ, inside, occOK, r.nsPerInsertMPQ, r.closeNsPerItemMPQ,
                    r.nbBK, r.epsBK, r.xLoBK, r.xHiBK, r.xBK, r.pctErrBK, r.pAchBK, r.rankErrBK, r.nsPerInsertBK, r.closeNsPerItemBK,
                    r.xTrue, r.tMPQMs, r.tBKMs
            );
        }


        // =========================
        //  Zipf helpers
        // =========================

        static ZipfDistribution makeZipfIntDist(int alphabetSize, double exponent, long seed) {
            if (alphabetSize <= 0) throw new IllegalArgumentException("alphabetSize must be > 0");
            if (exponent <= 0.0) throw new IllegalArgumentException("exponent must be > 0");
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
            // Defaults
            long N = 10_000_000L;   // stream length
            int A = 10000;          // alphabet size
            int B = 2500;             // buckets (also k for Bottom-k)
            Dist dist = Dist.ZIPF;  // UNIFORM or ZIPF
            double s = 1.5;           // Zipf exponent if ZIPF
            EXPONENT = s;
            double p = 0.05;        // target quantile
            long seed = 123456789L; // RNG seed
            int trials = 50;         // trials
            boolean verbose = false; // per-trial prints

            boolean autoDesignB = true;     // compute B from rank target if true
            Integer distinctGuess = A;   // if null, we will use Sactive from the stream

            Double rankEpsTarget = 0.05;   // eps_target
            double delta = 0.01;    // DKW delta (delta_q)
            Double deltaSamp = 0.09;        // delta_samp
//            # (δ_q, δ_samp) = (0.01, 0.09)  [sampling-heavy]
//            # (δ_q, δ_samp) = (0.05, 0.05)  [balanced]
//            # (δ_q, δ_samp) = (0.09, 0.01)  [quantile-heavy]

            //11877 2952
            //7000 4239
            //6233 2481
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

            // NEW: optional CSV out path and append flag
            String csvPath = "/home/dimpap/Desktop/util_scripts/mpq_res_auto"+EXPONENT+".csv";             // if null, no file output
            boolean csvAppend = true;          // we always append; header handled automatically

            // CLI: N A B dist s p seed delta trials verbose autoB distinctGuess rankEpsTarget deltaSamp [csvPath] [append]
            if (args.length >= 15) csvPath = args[14];
            if (args.length >= 16) csvAppend = Boolean.parseBoolean(args[15]); // kept for symmetry; we still append
//            List<Double> delta_samp    = new ArrayList<>(Arrays.asList(0.001, 0.025, 0.05, 0.075, 0.09));
//            List<Double> delta_q = new ArrayList<>(Arrays.asList(0.099, 0.075, 0.05, 0.025, 0.01));
            List<Double> delta_samp    = new ArrayList<>(Arrays.asList(0.05, 0.025, 0.09));
            List<Double> delta_q = new ArrayList<>(Arrays.asList(0.05, 0.075,  0.01));
//            List<Double> delta_q = new ArrayList<>(Arrays.asList(0.05, 0.05, 0.05, 0.05, 0.05));

            List<TrialResult> results = new ArrayList<>(trials);
            for(int z = 0; z < delta_q.size(); z++) {

                for (int i = 0; i < trials; i++) {
                    long trialSeed = mix64(
                            seed
                                    + i * 0x9E3779B97F4A7C15L   // per-trial
                                    + z * 0xD1B54A32D192ED03L   // per-split
                    );
                    TrialResult r = runExperiment(
                            N, A, B, dist, s, p, trialSeed, delta_q.get(z), /*verbose*/ verbose,
                            autoDesignB, distinctGuess, rankEpsTarget, delta_samp.get(z)
                    );
                    results.add(r);
                }
                int base = z * trials; // first row of the current split in `results`
                if (autoDesignB && results.size() > base && results.get(base).design != null) {
                    DesignResult d = results.get(base).design;
                    if(verbose) System.out.printf("Design summary (split %d): suggested B=%d, n_req=%d, LB(n_b)=%d, E[n_b]=%.2f%n",
                            z, d.Bsuggested, d.nReq, d.nLB, d.mu);

                    if (d.impossible) {
                        System.out.println("WARNING: Dhat < n_req, impossible to meet the requested rank target");
                    }
                }

            }



            // === CSV export ===
            if (csvPath != null && !csvPath.isEmpty()) {
                // write all trials (append; header if needed)
                writeSuiteACsv(Path.of(csvPath), results);
            } else if (results.size() == 1) {
                // fallback: keep your old single-trial stdout CSV
                printSuiteACsvHeader();
                printSuiteACsvRow(results.get(0));
            }

            // keep your summary prints for multi-trial runs
            summarize(results, p, delta);
        }
    }
