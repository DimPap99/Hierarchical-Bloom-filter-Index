package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Deque;

public class MostFreqPruning implements  PruningPlan {



    public void PruningPlan() {

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
        for( char c : pattern.text){
            p_temp = tree.estimator.estimate(c);
            if(p_temp > p_max){
                p_max = p_temp;
            }
        }
        b_a = Math.log(1 - confidence) / Math.log(1 - p_max);
        prunedLevels.add((int) (Math.floor(Math.log(tree.baseIntervalSize()/b_a)) + 1));
        return prunedLevels;
    }


}