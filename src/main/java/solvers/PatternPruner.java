package solvers;

import java.util.Arrays;

public final class PatternPruner {

    /** continuous, strictly increasing objective F(b) - a */
    private static double FminusA(double b, double[] pHat, double a) {
        double logProd = 0.0;               // work in log-space to avoid underflow
        for (double p : pHat) {
            double term = 1.0 - Math.pow(1.0 - p, b);  // 1 - (1-p)^b
            if (term <= 0.0)            // numerical guard (shouldn’t happen)
                return -a;              // forces caller to increase b
            logProd += Math.log(term);
        }
        double F = Math.exp(logProd);    // back to linear domain
        return F - a;
    }

    public static double solveBbyBisection(double[] pHat,
                                           double    a,
                                           double    epsRel,
                                           double    epsAbs) {

        if (!(a > 0.0 && a < 1.0))
            throw new IllegalArgumentException("a must be in (0,1)");

        if (Arrays.stream(pHat).anyMatch(p -> p <= 0.0 || p >= 1.0))
            throw new IllegalArgumentException("all pHat must lie in (0,1)");

        // --- 1. find a finite upper bound bHigh with F(bHigh) > a ----------
        double bLow  = 0.0;              // F(0+) = 0 < a
        double bHigh = 1.0;
        while (FminusA(bHigh, pHat, a) < 0.0) {
            bHigh *= 2.0;               // exponential search
            if (bHigh > 1e12)           // safety stop
                throw new IllegalStateException("excessive bHigh (>1e12)");
        }

        // --- 2. bisection ---------------------------------------------------
        double fLow  = -a;               // F(0+) - a  < 0
        double fHigh = FminusA(bHigh, pHat, a); //  > 0 by construction

        while (true) {
            double bMid = 0.5 * (bLow + bHigh);
            double fMid = FminusA(bMid, pHat, a);

            // convergence criteria -----------------------------------------
            if (Math.abs(fMid)   <= epsAbs ||
                    (bHigh - bLow)   <= epsRel * bMid) {
                return bMid;
            }

            if (fMid < 0.0) {           // root lies to the right
                bLow  = bMid;
                fLow  = fMid;
            } else {                    // root lies to the left
                bHigh = bMid;
                fHigh = fMid;
            }
        }
    }

    // ----------------- quick demo ------------------------------------------
    public static void main(String[] args) {
        double[] pHat = {0.12, 0.033, 0.004, 0.44, 0.003, 0.014};  // pattern of length r = 3
        double   a    = 0.99;                  // desired confidence
        long start = System.currentTimeMillis();

        double   b    = solveBbyBisection(pHat, a, 1e-8, 1e-12);

        System.out.printf("b_a,P  = %.6f%n", b);

        int W = 131072;
        int Lprune = (int) Math.floor(Math.log(W / b) / Math.log(2)) + 1;
        System.out.printf("L_prune = %d%n", Lprune);
    }
}
