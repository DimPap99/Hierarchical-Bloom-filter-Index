package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Deque;

public class MultiLevelPruning implements  PruningPlan {

    public void PruningPlan() {

    }

    @Override
    public void prepare(Deque<Frame> stack, int lp) {}
    @Override
    public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, int alphabetSize, double confidence) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();

        for(char c: pattern.text){
            int lp =tree.estimator.minPruneLevel(
                    String.valueOf(c),
                    1,
                    alphabetSize,
                    confidence);
            prunedLevels.add(lp);
        }
        return prunedLevels;
    }



}
