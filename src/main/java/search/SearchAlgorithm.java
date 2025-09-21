package search;

import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Deque;

public interface SearchAlgorithm {


//    ArrayList<Integer> report(char[] key, ArrayList<ImplicitTree> trees, boolean existence);
    CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset);

    int getCurrentOffset();
}
