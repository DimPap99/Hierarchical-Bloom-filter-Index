package search;

import tree.ImplicitTree;
import estimators.Estimator;
import org.apache.commons.math3.util.Pair;

import java.util.*;

public class BlockSearch implements SearchAlgorithm{
    public Estimator estimator;
    public int currentOffset;

    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        Probe probe = probe(tree,
                f.level(),             // unchanged helper signature
                f.intervalIdx(),
                p.nGramToInt);
        this.currentOffset = positionOffset;
        int currentIntervalSize = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval = tree.baseIntervalSize();

        if (probe.consumed() == 0) {
            this.currentOffset = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval;
            return null;
        } else {
            //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
            if (probe.complete() && childrenIntervalSize >= p.text.length) {
                //we just generate children as theres a chance that the pattern is in both of them
                tree.generateChildren(f, stack, positionOffset, tree.id);
                return null;
            } else {
                int intervalEndIdx = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval - 1;

                if(probe.complete()) return new CandidateRange(this.currentOffset, intervalEndIdx);
                else {
                    //the probes were incomplete. Meaning:
                    //Regardless of the children we are at a point where the current interval >= pattern and we have at least 1 character match from it
                    //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
                    //neighboring child (we have an overlap).
                    this.currentOffset = intervalEndIdx - probe.consumed() + 1;
                    return new CandidateRange(this.currentOffset, intervalEndIdx);
                }
            }
        }
    }


    public int getCurrentOffset(){return this.currentOffset;}
    Probe probe(ImplicitTree tree, int level, int interval, int[] pattern) {


        int matches = 0;

        long key;
        for (int i = 0; i < pattern.length; i++) {
            key =  tree.codec.pack(level, interval, pattern[i]);

            if (!tree.contains(level, key)) {
                return new Probe(matches, false);          // first mismatch at i
            }
            matches++;
            if(matches == pattern.length){
                return new Probe(pattern.length, true);
            }
        }
        return new Probe(pattern.length, true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;              // == tree.intervalSize

        /* end position of current interval inside the *global* stream */
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll
                - 1;

        return maxPos >= positionOffset;
    }
}
