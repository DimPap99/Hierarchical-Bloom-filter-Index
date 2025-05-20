package algorithms;

import PMIndex.HBI;
import PMIndex.ImplicitTree;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;

public interface SearchAlgorithm {


    ArrayList<Integer> report(String key, ImplicitTree tree, boolean existence);

}
