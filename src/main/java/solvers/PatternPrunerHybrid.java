package solvers;

import java.util.Arrays;

public final class PatternPrunerHybrid {

    // ------------------------ public API ------------------------
    public static double solveBHybrid(double[] pHat,
                                      double    a,
                                      double    epsRel,
                                      double    epsAbs) {

        validateInputs(pHat, a);

        // ---- 1. bracketing  -----------------------------------
        double bLow  = 0.0;                    // G(0+) = -a < 0
        double bHigh = 1.0;
        while (G(bHigh, pHat, a) < 0.0) {
            bHigh *= 2.0;
            if (bHigh > 1e12)
                throw new IllegalStateException("excessive bHigh (>1e12)");
        }

        // ---- 2. initial bisection stage -----------------------
        final double bisectTol = 1e-3;         // absolute width trigger
        while ((bHigh - bLow) > bisectTol) {
            double bMid = 0.5 * (bLow + bHigh);
            double gMid = G(bMid, pHat, a);
            if (gMid < 0.0)
                bLow = bMid;
            else
                bHigh = bMid;
        }

        // ---- 3. Newton stage inside the shrunken bracket -----
        double b = 0.5 * (bLow + bHigh);
        while (true) {
            double g  = G(b,  pHat, a);
            if (Math.abs(g) <= epsAbs) return b;

            double gPrime = Gprime(b, pHat);   // uses same log-space trick
            if (gPrime == 0.0) {               // fallback
                b = 0.5 * (bLow + bHigh);
            } else {
                double bNew = b - g / gPrime;

                // clamp to current bracket
                if (bNew <= bLow || bNew >= bHigh)
                    bNew = 0.5 * (bLow + bHigh);

                if (Math.abs(bNew - b) <= epsRel * b) return bNew;

                b = bNew;
            }

            // keep bracket consistent
            if (G(b, pHat, a) < 0.0)
                bLow = b;
            else
                bHigh = b;
        }
    }

    // ------------------------ helpers --------------------------
    private static void validateInputs(double[] pHat, double a) {
        if (!(a > 0.0 && a < 1.0))
            throw new IllegalArgumentException("a must be in (0,1)");
        if (Arrays.stream(pHat).anyMatch(p -> p <= 0.0 || p >= 1.0))
            throw new IllegalArgumentException("all pHat must lie in (0,1)");
    }

    /** G(b) = F(b) - a, in linear domain but computed via log-space */
    private static double G(double b, double[] pHat, double a) {
        double logProd = 0.0;
        for (double p : pHat) {
            double term = 1.0 - Math.pow(1.0 - p, b); // 1-(1-p)^b
            if (term <= 0.0) return -a;               // numeric guard
            logProd += Math.log(term);
        }
        return Math.exp(logProd) - a;
    }

    /** G'(b) = F'(b), computed safely in log-space */
    private static double Gprime(double b, double[] pHat) {
        double logProd = 0.0;
        double sumFrac = 0.0;
        for (double p : pHat) {
            double oneMinusP = 1.0 - p;
            double pow       = Math.pow(oneMinusP, b);
            double term      = 1.0 - pow;            // y_i(b)
            double termDer   = -pow * Math.log(oneMinusP); // y_i'(b)

            logProd += Math.log(term);               // for F(b)
            sumFrac += termDer / term;               // Σ y_i'/y_i
        }
        double F = Math.exp(logProd);
        return F * sumFrac;                          // F'(b)
    }

    // --------------------- demo main (optional) ----------------
    public static void main(String[] args) {
        double[] pHat = {0.12, 0.033, 0.004, 0.44, 0.003, 0.014};
        double   a    = 0.99;
        long start = System.currentTimeMillis();
        double   b    = solveBHybrid(pHat, a, 1e-10, 1e-12);
        long end = System.currentTimeMillis();
        System.out.println(end-start);
        int W = 131072;
        int Lprune = (int) Math.floor(Math.log(W / b) / Math.log(2)) + 1;

        System.out.printf("b_a,P  = %.12f%n", b);
        System.out.println("L_prune = " + Lprune);
    }
}
