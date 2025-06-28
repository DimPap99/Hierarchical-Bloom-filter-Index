package search;

import tree.ImplicitTree;

import java.util.ArrayList;

public interface SearchAlgorithm {


    ArrayList<Integer> report(char[] key, ArrayList<ImplicitTree> trees, boolean existence);

}
