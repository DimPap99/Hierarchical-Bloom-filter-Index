package search;
import tree.ImplicitTree;
import estimators.Estimator;

import javax.swing.tree.TreeNode;
import java.util.ArrayList;
import java.util.Deque;

/** Decides the lowest tree level worth exploring for a pattern. */
public interface PruningPlan {
    void prepare(Deque<Frame> stack, int lp);
    ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, int alphabetSize, double confidence);

}
