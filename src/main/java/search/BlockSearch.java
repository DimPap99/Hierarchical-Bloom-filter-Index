package search;

import tree.ImplicitTree;
import estimators.Estimator;
import org.apache.commons.math3.util.Pair;

import java.util.*;

public class BlockSearch implements SearchAlgorithm{
    public Estimator estimator;

    public class GlobalInfo{
        public boolean isSplit = false;
        public int matched = 0;
        public int patternLength;
        public boolean isReset = false;
        public int positionOffset = -1;
        public GlobalInfo(int patternLength){
            this.patternLength = patternLength;
        }


    }

    record Probe(int consumed, boolean complete) {} // complete==true if no miss in this node

    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        Probe probe = probe(tree,
                f.level(),             // unchanged helper signature
                f.intervalIdx(),
                p.text,
                0,
                p.text.length);
        int currentIntervalSize = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval = tree.baseIntervalSize();

        if (probe.consumed == 0) {
            positionOffset = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval;
            return null;
        } else {
            //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
            if (probe.complete && childrenIntervalSize >= p.text.length) {
                //we just generate children as theres a chance that the pattern is in both of them
                tree.generateChildren(f, stack, positionOffset, tree.id);
                return null;
            } else { //the probes were inclomplete. Meaning:
                //Regardless of the children we are at a point where the current interval >= pattern and we have at least 1 character match from it
                //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
                //neighboring child (we have an overlap).
                int intervalEndIdx = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval - 1;
                positionOffset = intervalEndIdx - probe.consumed + 1;
                return new CandidateRange(positionOffset, intervalEndIdx);

            }
        }
    }

    Probe probe(ImplicitTree tree, int level, int interval, char[] pattern, int offset, int remainder) {


        int matches = 0;

        long key;
        for (int i = offset; i < pattern.length; i++) {
            key =  tree.codec.pack(level, interval, pattern[i]);

            if (!tree.contains(level, key)) {
                return new Probe(matches, false);          // first mismatch at i
            }
            matches++;
            if(matches == remainder){
                return new Probe(remainder, true);
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

    private static int[] prefixFunction(char[] pat) {
        int[] pi = new int[pat.length];
        int k = 0;
        for (int i = 1; i < pat.length; ++i) {
            while (k > 0 && pat[k] != pat[i]) k = pi[k - 1];
            if (pat[k] == pat[i]) ++k;
            pi[i] = k;
        }
        return pi;
    }

}
