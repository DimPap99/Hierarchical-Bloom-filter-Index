package utilities;

import tree.ImplicitTree;

public final class  MathUtils {
    private MathUtils() {
        throw new AssertionError("MathUtils must not be instantiated");
    }

    public static double h_b(int width, int lvl, double prob){
        return Math.pow( 1 - (1 - prob), width/Math.pow(2, lvl));
    }

    public static double  expectedProbesPerNode(double[] probs,
                                 double bloomFp,
                                 int     blockLen) {
        double prod = 1.0;
        double sum  = 0.0;
        for (double p : probs) {
            double h   = 1.0 - Math.pow(1.0 - p, blockLen); // hit-once prob
            double qi  = h + bloomFp * (1.0 - h);           // Bloom “yes”
            prod *= qi;                                     // q1·…·qk
            sum  += prod;                                   // accumulate Σ
        }
        return sum;  // E[# probes in this node]
    }

    public static  int pruningLevel(ImplicitTree tree, double conf, double prob){
        double b_a;

        b_a = Math.log(1 - conf) / Math.log(1 - prob);   //   bα
        double log2 = Math.log(tree.baseIntervalSize() / b_a)    //   log_e
                / Math.log(2.0);
        int    rawLp  = (int) Math.floor(log2) + 1;
        return Math.max(0, Math.min(rawLp, tree.maxDepth() - 1));
    }

}
