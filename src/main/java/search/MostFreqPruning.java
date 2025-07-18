package search;

import estimators.CostFunctionMaxProb;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Deque;

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
    public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, int alphabetSize, double confidence, CostFunctionMaxProb cf ) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();

        double b_a;
        double p_max = 0;
        double p_temp;
//        ArrayList<Double> probs = new ArrayList<>();
        for( int c : pattern.nGramToInt){
            p_temp = tree.estimator.estimate(c);
//            probs.add(p_temp);
            if(p_temp > p_max){
                p_max = p_temp;
            }
        }

        int    lp   = MathUtils.pruningLevel(tree, confidence, p_max);//Math.max(0, Math.min(raw, tree.maxDepth() - 1));

        prunedLevels.add(lp);
        return prunedLevels;
    }





}