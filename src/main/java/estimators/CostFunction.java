package estimators;

import tree.ImplicitTree;

import static java.lang.Math.pow;

public class CostFunction {

    //cost of constant operations. Typically measured in nano seconds.
    double bloomProbeCost;
    double leafSearchCost;

    public CostFunction() {

    }

    //hit probability for a single character in an interval b of size w/2^lvl
    public double h_b(int width, int lvl, double prob){
        return Math.pow( 1 - (1 - prob), width/Math.pow(2, lvl));
    }


    public int pruningLevel(ImplicitTree tree, double conf, double prob){
        double b_a;
        b_a = Math.log(1 - conf) / Math.log(1 - prob);   //   bα
        double log2 = Math.log(tree.baseIntervalSize() / b_a)    //   log_e
                / Math.log(2.0);
        int    rawLp  = (int) Math.floor(log2) + 1;
        return Math.max(0, Math.min(rawLp, tree.maxDepth() - 1));
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


    public double costFunc(double[] probs, double bfFalsePosRate, int width, int lvl, int patternSize, int maxDepth){
        //horizontal nodes
        double horizontalNodes =  Math.pow(2, lvl);
        double verticalNodes =  maxDepth - lvl;
        double fpr = this.fpRate(probs, bfFalsePosRate, width, lvl);
        //the cost incured by doing a the first vertical scan to determine if the pattern is inside an interval
        double intervalsToBeChecked = 1 + fpr * (horizontalNodes - 1); //1 because one of them will be true the others are fp
        //for the intervals that we have the true occurence as we ll as the fp, we gonna probe for all the characters of the pattern
        //the other intervals do not have the full occurences. As such we need at least one and less that patternSize
        //to reject them.
        // TODO: For now assume we reject them on 1. (thats not necessarily accurate) We could even fully add them
        //and not include the vertical cost. Need to play around with this.)
        double horizontalCost = patternSize * intervalsToBeChecked * this.bloomProbeCost + (horizontalNodes - intervalsToBeChecked);

        double verticalCost = patternSize * intervalsToBeChecked * verticalNodes * this.bloomProbeCost;
        //this occurs only once per vertical search one we reach the leaf level
        double leafCost = patternSize * intervalsToBeChecked * this.leafSearchCost;


        //full cost of current plan in terms of operations multiplied with time
        return horizontalCost + verticalCost + leafCost;


    }
}
