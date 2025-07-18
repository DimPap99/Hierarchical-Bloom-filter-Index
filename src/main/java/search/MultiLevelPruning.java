package search;

import estimators.CostFunctionMaxProb;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Deque;

public class MultiLevelPruning implements  PruningPlan {
    public double conf;

    public void PruningPlan() {

    }
    public MultiLevelPruning(double conf){
        this.conf = conf;

    }
    @Override
    public void prepare(Deque<Frame> stack, int lp) {}
    @Override
    public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, int alphabetSize, double confidence, CostFunctionMaxProb cf) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();

        for(int c: pattern.nGramToInt){
            int lp = MathUtils.pruningLevel(tree, confidence, tree.estimator.estimate(c));
            prunedLevels.add(lp);
        }
        return prunedLevels;
    }




}
