package search;
import estimators.CostFunctionMaxProb;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Deque;

/** Decides the lowest tree level worth exploring for a pattern. */
public interface PruningPlan {
    void prepare(Deque<Frame> stack, int lp);
    ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, int alphabetSize, double confidence, CostFunctionMaxProb cf);

}
