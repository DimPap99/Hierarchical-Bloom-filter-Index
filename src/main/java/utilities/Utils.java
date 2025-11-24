package utilities;

public class Utils {

    // Using >>> is likely faster than masking.
    public static int getIntervalIndex(int domainTotalBits, int level, int value){
//        int bitsToExtract = totalBits - level;
//        //bring one to the position of MSB for the amount of bits we want to extract
//        //everything else remains 0. -1 to flip all the other bits. This will give us as many '111...' as the
//        //total bits we want to extract. Then shift left by level in order to push them to the leftmost side
//        int mask          = ((1 << bitsToExtract) - 1) << level;
//        //Logical and between our value and mask. The mask will grab the MSB bits we want. Then shift right (divide) based
//        //on the level to find at which interval the value exists
//        //Example: Level 4, max level = 5. 2^5 = 32. So for level 4 we got 32 / 2^4 = 32/16 = 2. We got 2 16 intervals
//        //Consider value = 17. We want to find at which interval it should be placed. Since we got 2 intervals, the total
//        //bits we want to extract are 1, as it can represent 2 options.
//        //Our mask then becomes: 10000. 17: 10001. We do: 17 & 10000 = 10000 (16) and then shift as many positions as the
//        //current level, that is 4. 10000 >> 4 = 1. So it belongs in the second interval.
//
//        return (value & mask) >>> level;
        return value >>> domainTotalBits - level;
    }

    public static enum MemPolicy {
        NONE,
        REACTIVE,
        PREDICTIVE
    }

    // Returns "<interval bits><# placeholders>" for a value at a given level.
    public static String intervalWithHashes(int totalBits, int level, int value) {

        int lowerBits = totalBits - level;          // how many '#' chars
        int interval  = (level == 0) ? 0           // avoid >>> with shift = totalBits
                : value >>> lowerBits;

        // Interval bits on the left.
        String intervalBits;
        if (level > 0) {                            // <-- guard fixes %0s problem
            intervalBits = Integer.toBinaryString(interval);
            intervalBits = String.format("%" + level + "s", intervalBits)
                    .replace(' ', '0');
        } else {
            intervalBits = "";                      // level == 0, no bits
        }

        // '#' placeholders on the right.
        String hashes = "#".repeat(lowerBits);      // repeat(0) returns ""

        return intervalBits + hashes;
    }

    public static final class HopsDesignResult {
        public int suggestedBuckets = 100000;
        public final int requiredSampleSize;
        public final int occupancyLowerBound;
        public final double expectedNonEmpty;
        public final double variance;
        public final boolean impossible;

        public HopsDesignResult(int suggestedBuckets,
                                int requiredSampleSize,
                                int occupancyLowerBound,
                                double expectedNonEmpty,
                                double variance,
                                boolean impossible) {
            this.suggestedBuckets = suggestedBuckets;
            this.requiredSampleSize = requiredSampleSize;
            this.occupancyLowerBound = occupancyLowerBound;
            this.expectedNonEmpty = expectedNonEmpty;
            this.variance = variance;
            this.impossible = impossible;
        }
    }

    // Design the number of HOPS buckets for a rank error target using Chebyshev bounds.
    public static HopsDesignResult designBucketsForRankTargetChebyshev(int distinctEstimate,
                                                                       double epsTarget,
                                                                       double deltaQ,
                                                                       double deltaSample) {
        if (distinctEstimate <= 0) {
            throw new IllegalArgumentException("distinctEstimate must be > 0");
        }
        if (!(epsTarget > 0 && epsTarget < 1)) {
            throw new IllegalArgumentException("epsTarget must be in (0,1)");
        }
        if (!(deltaQ > 0 && deltaQ < 1)) {
            throw new IllegalArgumentException("deltaQ must be in (0,1)");
        }
        if (!(deltaSample > 0 && deltaSample < 1)) {
            throw new IllegalArgumentException("deltaSample must be in (0,1)");
        }

        int required = requiredSampleSizeForDKW(epsTarget, deltaQ);

        if (distinctEstimate < required) {
            final int bucketCap = 1 << 22;
            int suggested = Math.min(bucketCap, Math.max(16, 2 * distinctEstimate));
            double mu = occupancyExpectation(distinctEstimate, suggested);
            double var = occupancyVariance(distinctEstimate, suggested);
            int nLb = occupancyLowerBoundChebyshev(distinctEstimate, suggested, deltaSample);
            return new HopsDesignResult(suggested, required, nLb, mu, var, true);
        }

        final int bucketCap = 1 << 24;
        int hi = 1;
        while (occupancyLowerBoundChebyshev(distinctEstimate, hi, deltaSample) < required && hi < bucketCap) {
            hi <<= 1;
        }
        if (hi >= bucketCap) {
            double mu = occupancyExpectation(distinctEstimate, bucketCap);
            double var = occupancyVariance(distinctEstimate, bucketCap);
            int nLb = occupancyLowerBoundChebyshev(distinctEstimate, bucketCap, deltaSample);
            return new HopsDesignResult(bucketCap, required, nLb, mu, var, true);
        }

        int lo = 1;
        int best = hi;
        while (lo <= hi) {
            int mid = lo + ((hi - lo) >>> 1);
            int nLb = occupancyLowerBoundChebyshev(distinctEstimate, mid, deltaSample);
            if (nLb >= required) {
                best = mid;
                hi = mid - 1;
            } else {
                lo = mid + 1;
            }
        }

        double mu = occupancyExpectation(distinctEstimate, best);
        double var = occupancyVariance(distinctEstimate, best);
        int nLb = occupancyLowerBoundChebyshev(distinctEstimate, best, deltaSample);
        return new HopsDesignResult(best, required, nLb, mu, var, false);
    }

    private static double occupancyExpectation(int distinct, int buckets) {
        if (buckets <= 0) {
            return 0.0;
        }
        return buckets * (1.0 - Math.pow(1.0 - 1.0 / buckets, distinct));
    }

    private static double occupancyVariance(int distinct, int buckets) {
        if (buckets <= 0) {
            return 0.0;
        }
        double t1 = Math.pow(1.0 - 1.0 / buckets, distinct);
        double t2 = Math.pow(1.0 - 2.0 / buckets, distinct);
        double q = 1.0 - t1;
        double var = buckets * q * (1.0 - q) + buckets * (buckets - 1.0) * (1.0 - 2.0 * t1 + t2 - q * q);
        return Math.max(0.0, var);
    }

    private static int occupancyLowerBoundChebyshev(int distinct, int buckets, double deltaSample) {
        double mu = occupancyExpectation(distinct, buckets);
        double var = occupancyVariance(distinct, buckets);
        double nStar = mu - Math.sqrt(var / Math.max(1e-12, deltaSample));
        return (int) Math.floor(Math.max(0.0, nStar));
    }

    private static int requiredSampleSizeForDKW(double eps, double deltaQ) {
        return (int) Math.ceil(Math.log(2.0 / deltaQ) / (2.0 * eps * eps));
    }


}
