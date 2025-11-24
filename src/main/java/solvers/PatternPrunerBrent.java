package solvers;

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import java.util.Arrays;

public final class PatternPrunerBrent {

    // Solve product_i [1 - (1 - p_i)^b] = a for b > 0 using a Brent solver.
    public static double solveB(double[] pHat,
                                double   a,
                                double   relTol,
                                double   absTol) {

        // Sanity checks.
        if (!(a > 0.0 && a < 1.0))
            throw new IllegalArgumentException("a must lie in (0,1)");
        if (Arrays.stream(pHat).anyMatch(p -> p <= 0.0 || p >= 1.0))
            throw new IllegalArgumentException("all pHat must be in (0,1)");

        // Objective G(b) = F(b) - a, computed in log-space to avoid underflow.
        UnivariateFunction G = b -> {
            double logProd = 0.0;
            for (double p : pHat) {
                double term = 1.0 - Math.pow(1.0 - p, b);   // 1 - (1-p)^b
                if (term <= 0.0) return -a;                 // guard: pushes root rightward
                logProd += Math.log(term);
            }
            return Math.exp(logProd) - a;
        };

        // Brent-Dekker solver (order = 5, classic Brent).
        BracketingNthOrderBrentSolver solver =
                new BracketingNthOrderBrentSolver(absTol, relTol, 5);

        // Build a bracket [lower, upper] with opposite signs.
        double lower = 0.0;          // G(0+) = -a < 0
        double upper = 1.0;          // start guess; expand until G(upper) > 0
        while (G.value(upper) < 0.0) {
            upper *= 2.0;
            if (upper > 1e12)
                throw new IllegalStateException("failed to bracket root (upper > 1e12)");
        }

        return solver.solve(100, G, lower, upper);
    }

    // Demo.
    public static void main(String[] args) {
        double[] pHat = {0.12, 0.033, 0.004, 0.44, 0.003, 0.014};  // pattern of length r = 3
        double   a    = 0.99;                  //                   // desired confidence

        double b = solveB(pHat, a, 1e-10, 1e-12);
        System.out.printf("b_a,P = %.12f%n", b);

        int W = 131072;                           // window length, say 1024
        int Lprune = (int) Math.floor(Math.log(W / b) / Math.log(2)) + 1;
        System.out.println("L_prune = " + Lprune);
    }
}
