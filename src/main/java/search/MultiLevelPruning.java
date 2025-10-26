package search;

import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Deque;

public class MultiLevelPruning implements  PruningPlan {
    public double conf;

    public void PruningPlan() {

    }
    public MultiLevelPruning(double conf){
        //this.conf = conf;

    }
    @Override
    public void prepare(Deque<Frame> stack, int lp) {}
    @Override
    public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, double confidence, boolean strides) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();
        long[] arr = strides ? pattern.effectiveNgramArr : pattern.nGramToLong;
        for(long token : arr){
            int lp = Math.min(tree.effectiveRoot(),MathUtils.pruningLevel(tree, confidence, tree.estimator.estimate(token)));
            prunedLevels.add(lp);
        }
        return prunedLevels;
    }




}
