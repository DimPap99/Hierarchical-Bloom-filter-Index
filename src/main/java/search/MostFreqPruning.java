package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Deque;
import java.util.List;

public class MostFreqPruning implements  PruningPlan {

    public double conf;

    public void PruningPlan() {

    }
    public MostFreqPruning(double conf){
        this.conf = conf;

    }
    @Override
    public void prepare(Deque<Frame> stack, int lp) {
        for(int i=0; i< Math.pow(2, lp); i++){
            stack.addLast(new Frame(lp, i));
        }
    }

    @Override
    public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, int alphabetSize, double confidence ) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();

        double b_a;
        double p_max = 0;
        double p_temp;
        ArrayList<Double> probs = new ArrayList<>();
        for( int c : pattern.nGramToInt){
            p_temp = tree.estimator.estimate(c);
            probs.add(p_temp);
            if(p_temp > p_max){
                p_max = p_temp;
            }
        }

        b_a = Math.log(1 - conf) / Math.log(1 - p_max);   //   bα
        double log2 = Math.log(tree.baseIntervalSize() / b_a)    //   log_e
                / Math.log(2.0);
        int    raw  = (int) Math.floor(log2) + 1;                // ⌊·⌋ + 1
        int    lp   = Math.max(0, Math.min(raw, tree.maxDepth() - 1));

        prunedLevels.add(lp);
        return prunedLevels;
    }


}