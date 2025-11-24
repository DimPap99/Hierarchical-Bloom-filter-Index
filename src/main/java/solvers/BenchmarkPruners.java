package solvers;

import java.util.Random;

// Compare different b-a,P solvers under identical synthetic workloads.
public class BenchmarkPruners {

    // Solver functional interface.
    @FunctionalInterface
    public interface PrunerSolver {
        double solve(double[] pHat, double a);
    }

    // Run the same workload on a solver and report elapsed-time stats.
    public static void benchmark(String name,
                                 PrunerSolver solver,
                                 double[][]   patterns,
                                 double       a,
                                 int          repeats) {

        long bestNs = Long.MAX_VALUE, totalNs = 0;

        for (int k = 0; k < repeats; k++) {
            long t0 = System.currentTimeMillis();
            for (double[] patt : patterns) solver.solve(patt, a);
            long dt = System.currentTimeMillis() - t0;
            totalNs += dt;
            if (dt < bestNs) bestNs = dt;
        }
        double avgMs  = totalNs / 1e6 / repeats;
        double bestMs = bestNs   / 1e6;

        System.out.printf("%-12s : avg = %.4f ms   best = %.4f ms%n",
                name, avgMs, bestMs);
    }

    public static void main(String[] args) {

        int       numPatterns = 500;
        int       maxLen      = 15;           // pattern length r  (1..10)
        Random    rng         = new Random(42);
        double[][] patterns   = new double[numPatterns][];

        for (int i = 0; i < numPatterns; i++) {
            int r = 1 + rng.nextInt(maxLen);
            double[] p = new double[r];
            for (int j = 0; j < r; j++)
                p[j] = 0.1 + rng.nextDouble() * 0.15;   // 0.001 .. 0.151
            patterns[i] = p;
        }
        double a = 0.99;                      // same confidence for all
        int    repeats = 100;                   // run each solver 5 times

        PrunerSolver halley = (p, conf) ->
                PatternPrunerHalley.solveB(p, conf, 1e-10, 1e-12);

        PrunerSolver brent  = (p, conf) ->
                PatternPrunerBrent.solveB(p, conf, 1e-10, 1e-12);

        PrunerSolver pPrune  = (p, conf) ->
                PatternPruner.solveBbyBisection(p, conf, 1e-10, 1e-12);
        PrunerSolver phPrune  = (p, conf) ->
                PatternPrunerHybrid.solveBHybrid(p, conf, 1e-10, 1e-12);


        benchmark("Halley", halley, patterns, a, repeats);
        benchmark("Brent",  brent,  patterns, a, repeats);
        benchmark("pPrune", pPrune, patterns, a, repeats);
        benchmark("phPrune",  phPrune,  patterns, a, repeats);
    }
}
