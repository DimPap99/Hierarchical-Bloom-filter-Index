package algorithms;

import PMIndex.HBI;
import PMIndex.ImplicitTree;
import org.apache.commons.math3.util.Pair;

import java.util.ArrayList;

public interface SearchAlgorithm {

    boolean exists(String key, ImplicitTree tree);

    ArrayList<Integer> report(String key, ImplicitTree tree);

}
