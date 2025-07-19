package estimators;

import search.Pattern;
import tree.ImplicitTree;

import java.util.Arrays;

import static utilities.MathUtils.*;

public class CostFunctionMaxProb implements CostFunction {

    //cost of constant operations. Typically measured in nano seconds.
    double bloomProbeCost = 3654;
    double leafSearchCost = 600;
    public double alpha;
    public CostFunctionMaxProb() {

    }

    //hit probability for a single character in an interval b of size w/2^lvl



    public int minCostLp(ImplicitTree tree,
                         double   bfFalsePosRate,   // Bloom collision per symbol (β)
                         double confInit,
                            Pattern p, double bfCost, double leafCost)       // largest p̂ in pattern
    {


        double pMax = Arrays.stream(tree.estimator.estimateALl(p)).min().getAsDouble();
        double[] pArr = new double[1];
        pArr[0] = pMax;
        double bestCost = Double.POSITIVE_INFINITY;
        int bestLp = 20;
        for (double alpha = confInit; alpha <= 1 + 1e-9; alpha += 0.05) {
            if(alpha >=1) alpha = 0.99;
            int lp   = pruningLevel(tree, alpha, pMax);
            double c = costFunc(pArr, bfFalsePosRate, tree.baseIntervalSize(), lp, p.nGramToInt.length, tree.maxDepth());

            if (c < bestCost) {
                bestCost = c;
                bestLp   = lp;
                this.alpha = alpha;
            }

            if(alpha == 0.99) {
                //System.out.println(bestLp + "niggr");
                return bestLp;
            }

        }
        System.out.println("Cyka blyat");

        return pruningLevel(tree, 0.99, pMax);          // could also return alpha or bestCost
    }


    //False Probability rate for current query. Disclaimer: Its not the same as Bloom Filter fp
    //The current index does not enforce an order in the characters of a pattern. Therefor we might be searching for
    //abc but we can be getting a positive answer for all its permutation such as bca etc without abc being inside
    //the indexed items. Along that false positive rate we also add the FP of the bloom filter as shown below.
    //This is the simple worst case cf, where we dont track an active set of characters per level and we use a universal
    //pruning level
    public double fpRate(double[] probs, double bfFalsePosRate, int width, int lvl){
        double total_probes_prob = 1;
         //assume independece between character sequences. (Basically we dont use conditional probs)
        for(double prob: probs){
            double currHB = h_b(width, lvl, prob);
            double elementFromBFFP = (1 - currHB) * bfFalsePosRate;
            //add the probabilities that a particular character is inside our interval as well as
            //the probability that the character is not inside that interval (1 - hb) but we still get a positive
            //answer due to the fp from the bloom filter.
            total_probes_prob *= currHB + elementFromBFFP;
        }
        return total_probes_prob;
    }

    double costFunc(double[] probs,
                    double   bfFalsePosRate,   // Bloom collision per symbol (β)
                    int      width,            // root-window size W = 2^d
                    int      lvl,              // current L_max
                    int      patternSize,      // r
                    int      maxDepth)         // deepest internal level (d)
    {
        // basic node counts
        double horizontalNodes   = 1L << lvl;               // 2^L
        double vertNodesPerBr    = maxDepth - lvl + 1;      // incl. start node

        // symbol-set false-positive rate
        double fpr = fpRate(probs, bfFalsePosRate, width, lvl);   // see helper

        // expected # intervals followed = true branch + FP branches
        double intervalsFollowed = 1.0 + fpr * (horizontalNodes - 1);

        // expected probes per node
        int    blockLen = width >> lvl;                      // W / 2^L
        double probesPerNode = expectedProbesPerNode(probs, bfFalsePosRate,
                blockLen);


        double horizontalCost = probesPerNode * horizontalNodes * bloomProbeCost;

        double verticalCost   = patternSize * intervalsFollowed *
                vertNodesPerBr * bloomProbeCost;

        double leafCost       = patternSize * ((intervalsFollowed * leafSearchCost) + 1) ;

        return horizontalCost + verticalCost + leafCost;
    }
//    public double costFunc(double[] probs, double bfFalsePosRate, int width, int lvl, int patternSize, int maxDepth){
//        //horizontal nodes
//        double horizontalNodes =  Math.pow(2, lvl);
//        double verticalNodes =  maxDepth - lvl;
//        double fpr = this.fpRate(probs, bfFalsePosRate, width, lvl);
//        //the cost incured by doing a the first vertical scan to determine if the pattern is inside an interval
//        double intervalsToBeChecked = 1 + fpr * (horizontalNodes - 1); //1 because one of them will be true the others are fp
//        //for the intervals that we have the true occurence as we ll as the fp, we gonna probe for all the characters of the pattern
//        //the other intervals do not have the full occurences. As such we need at least one and less that patternSize
//        //to reject them.
//        // TODO: For now assume we reject them on 1. (thats not necessarily accurate) We could even fully add them
//        //and not include the vertical cost. Need to play around with this.)
//        double horizontalCost = patternSize * intervalsToBeChecked * this.bloomProbeCost + (horizontalNodes - intervalsToBeChecked);
//
//        double verticalCost = patternSize * intervalsToBeChecked * verticalNodes * this.bloomProbeCost;
//        //this occurs only once per vertical search one we reach the leaf level
//        double leafCost = patternSize * intervalsToBeChecked * this.leafSearchCost;
//
//
//        //full cost of current plan in terms of operations multiplied with time
//        return horizontalCost + verticalCost + leafCost;
//
//
//    }
}
