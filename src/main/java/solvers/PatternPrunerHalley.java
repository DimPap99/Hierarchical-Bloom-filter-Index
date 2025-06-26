package solvers;

import java.util.Arrays;

/**
 *  Solve  ∏i [ 1 − (1 − p_i)^b ] = a   for  b > 0   via bracket-clamped Halley.
 *
 *  All computations are performed in log-space to avoid underflow.
 *
 *  Compile & run:
 *      javac solvers.PatternPrunerHalley.java
 *      java  solvers.PatternPrunerHalley
 */
public final class PatternPrunerHalley {

    /* ------------------------------------------------------------------ */
    public static double solveB(double[] pHat, double a,
                                double epsRel, double epsAbs) {

        if (!(a > 0.0 && a < 1.0))
            throw new IllegalArgumentException("a must lie in (0,1)");
        if (Arrays.stream(pHat).anyMatch(p -> p <= 0.0 || p >= 1.0))
            throw new IllegalArgumentException("each pHat in (0,1)");

        // 1. Initial bracket  [bLow, bHigh]  with opposite signs
        double bLow  = 0.0;                // G(0+) = -a < 0
        double bHigh = 1.0;
        while (G(bHigh, pHat, a) < 0.0) {  // expand until G > 0
            bHigh *= 2.0;
            if (bHigh > 1e12)
                throw new IllegalStateException("failed to bracket (bHigh > 1e12)");
        }

        // 2. Start from the midpoint
        double b = 0.5 * (bLow + bHigh);

        // 3. Halley iterations with bracket clamping
        for (int iter = 0; iter < 50; iter++) {

            double g = G(b, pHat, a);
            if (Math.abs(g) <= epsAbs) return b;      // vertical tolerance

            Derivs d   = derivs(b, pHat);             // F, S1, S2 in one pass
            double F   = d.F;
            double S1  = d.S1;
            double S2  = d.S2;
            double g1  = F * S1;                      // F′
            double g2  = F * (S1 * S1 - S2);          // F′′

            double denom = 2.0 * g1 * g1 - g * g2;
            double bNew;
            if (denom == 0.0) {                       // fallback: bisection
                bNew = 0.5 * (bLow + bHigh);
            } else {
                double delta = (2.0 * g * g1) / denom;
                bNew = b - delta;

                // clamp the step inside the bracket
                if (bNew <= bLow || bNew >= bHigh)
                    bNew = 0.5 * (bLow + bHigh);
            }

            if (Math.abs(bNew - b) <= epsRel * b) return bNew;   // relative tol
            b = bNew;

            // maintain the bracket
            if (G(b, pHat, a) < 0.0)
                bLow = b;
            else
                bHigh = b;
        }
        throw new IllegalStateException("Halley failed to converge in 50 iterations");
    }

    /* -------------------------- helper functions ----------------------- */
    /** G(b) = F(b) - a */
    private static double G(double b, double[] pHat, double a) {
        double logProd = 0.0;
        for (double p : pHat) {
            double term = 1.0 - Math.pow(1.0 - p, b);      // 1 − (1−p)^b
            if (term <= 0.0) return -a;                    // numeric guard
            logProd += Math.log(term);
        }
        return Math.exp(logProd) - a;
    }

    /** simultaneous computation of F(b),  Σ y′/y,  Σ y″/y  */
    private static class Derivs {
        final double F, S1, S2;
        Derivs(double F, double S1, double S2) {
            this.F = F; this.S1 = S1; this.S2 = S2;
        }
    }
    private static Derivs derivs(double b, double[] pHat) {
        double logProd = 0.0, S1 = 0.0, S2 = 0.0;
        for (double p : pHat) {
            double oneMinusP = 1.0 - p;
            double pow       = Math.pow(oneMinusP, b);
            double ln        = Math.log(oneMinusP);
            double y         = 1.0 - pow;           // y_i(b)
            double y1        = -pow * ln;           // y_i′
            double y2        =  pow * ln * ln;      // y_i″

            logProd += Math.log(y);
            S1      += y1 / y;
            S2      += y2 / y;
        }
        return new Derivs(Math.exp(logProd), S1, S2);
    }

    /* ------------------------------ demo ------------------------------ */
    public static void main(String[] args) {
        // Example pattern and confidence
        double[] pHat = {0.12, 0.033, 0.004};
        double   a    = 0.85;

        // Solve for b_a,P
        double b = solveB(pHat, a, 1e-10, 1e-12);
        System.out.printf("b_a,P  = %.12f%n", b);

        // Map to pruning level
        int W = 1024;                                   // window length
        int Lprune = (int) Math.floor(Math.log(W / b) / Math.log(2)) + 1;
        System.out.println("L_prune = " + Lprune);
    }
}
