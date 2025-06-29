package search;

import membership.Membership;
import org.apache.commons.math3.util.Pair;
import tree.ImplicitTree;

import java.util.ArrayList;

public interface Verifier {

    Pair<ArrayList<Integer>, Integer> verify(int currentTreeIdx,
                                             ArrayList<ImplicitTree<Membership>> trees,
                                             CandidateRange candidateRange,
                                             Pattern pat);
}
