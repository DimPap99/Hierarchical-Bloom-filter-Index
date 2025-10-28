package search;

import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Deque;

public class MultiLevelPruning implements  PruningPlan {
    public double conf;
    public double fp;
    public void PruningPlan() {

    }
    public MultiLevelPruning(double conf, double fp) {
        //this.conf = conf;
        this.fp = fp;

    }
    @Override
    public void prepare(Deque<Frame> stack, int lp) {}
    @Override
    public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, double confidence, boolean strides) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();
        long[] arr = strides ? pattern.effectiveNgramArr : pattern.nGramToLong;
        for(long token : arr){
            double p = tree.estimator.estimate(token);
            int lp = Math.max(tree.effectiveRoot(),MathUtils.pruningLevelBloom(tree, confidence, p, this.fp));
            prunedLevels.add(lp);
        }
        return prunedLevels;
    }




}
