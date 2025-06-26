package solvers;/*  solvers.PatternPrunerBrent.java
 *
 *  Solve:   ∏i [ 1 - (1 - p_i)^b ]  =  a      for  b > 0
 *  using Apache Commons Math 3's BracketingNthOrderBrentSolver
 *
 *  Maven dependency (pom.xml):
 *  <dependency>
 *      <groupId>org.apache.commons</groupId>
 *      <artifactId>commons-math3</artifactId>
 *      <version>3.6.1</version>
 *  </dependency>
 */

import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BracketingNthOrderBrentSolver;
import java.util.Arrays;

public final class PatternPrunerBrent {

    /** ------------------------------------------------------------------
     *  @param pHat    empirical probabilities 0 < p_i < 1
     *  @param a       confidence threshold 0 < a < 1
     *  @param relTol  relative tolerance on b
     *  @param absTol  absolute tolerance on F(b) - a
     *  @return        unique root b_{a,P}
     *  ------------------------------------------------------------------ */
    public static double solveB(double[] pHat,
                                double   a,
                                double   relTol,
                                double   absTol) {

        // --- sanity checks ------------------------------------------------
        if (!(a > 0.0 && a < 1.0))
            throw new IllegalArgumentException("a must lie in (0,1)");
        if (Arrays.stream(pHat).anyMatch(p -> p <= 0.0 || p >= 1.0))
            throw new IllegalArgumentException("all pHat must be in (0,1)");

        /* Objective G(b) = F(b) - a, computed in log-space so the
           product cannot underflow when b is large. */
        UnivariateFunction G = b -> {
            double logProd = 0.0;
            for (double p : pHat) {
                double term = 1.0 - Math.pow(1.0 - p, b);   // 1 − (1−p)^b
                if (term <= 0.0) return -a;                 // guard: pushes root rightward
                logProd += Math.log(term);
            }
            return Math.exp(logProd) - a;
        };

        // Brent–Dekker solver (order = 5 ⇢ classic Brent)
        BracketingNthOrderBrentSolver solver =
                new BracketingNthOrderBrentSolver(absTol, relTol, 5);

        // --- build a bracket [lower, upper] with opposite signs ----------
        double lower = 0.0;          // G(0+) = -a < 0
        double upper = 1.0;          // start guess; expand until G(upper) > 0
        while (G.value(upper) < 0.0) {
            upper *= 2.0;
            if (upper > 1e12)
                throw new IllegalStateException("failed to bracket root (upper > 1e12)");
        }

        // solve!  maxEval = 100 is plenty (Brent typically converges in <10)
        return solver.solve(100, G, lower, upper);
    }

    /* -------------------------- demo ----------------------------------- */
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
