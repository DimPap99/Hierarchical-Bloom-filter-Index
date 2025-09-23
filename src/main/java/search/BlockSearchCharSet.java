// search/BlockSearchCharSet.java
package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;

public class BlockSearchCharSet implements SearchAlgorithm {
    public Estimator estimator;
    public int currentOffset;

    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        Probe probe = probe(tree, f.level(), f.intervalIdx(), p.nGramToInt, p.charStartLp);
        this.currentOffset = positionOffset;

        int currentIntervalSize  = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval     = tree.baseIntervalSize();
        int intervalEndIdx       = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval - 1;

        // **Only** skip when we truly saw a NO at the first symbol.
        if (probe.consumed() == 0 && probe.testedFirst()) {
            this.currentOffset = intervalEndIdx + 1;
            return null;
        }

        // If we didn't test first symbol but have no NO yet → treat as “complete so far”
        boolean noDefinitiveNo = probe.complete() || !probe.testedFirst();

        if (noDefinitiveNo && childrenIntervalSize >= p.text.length) {
            tree.generateChildren(f, stack, positionOffset, tree.id);
            return null;
        }

        if (probe.complete()) {
            // Full success for what we were allowed to test → verifier over whole node
            return new CandidateRange(this.currentOffset, intervalEndIdx);
        } else {
            // Partial: verifier only where the rightmost alignment can still succeed
            this.currentOffset = intervalEndIdx - probe.consumed() + 1;
            return new CandidateRange(this.currentOffset, intervalEndIdx);
        }
    }

    public int getCurrentOffset(){ return this.currentOffset; }

    // Counts consecutive TRUEs from start
    int countStartConsecutiveMatches(boolean[] matchedArr){
        int matches = 0;
        for (boolean b : matchedArr) {
            if (!b) break;
            matches++;
        }
        return matches;
    }

    Probe probe(ImplicitTree tree, int level, int interval, int[] pattern, ArrayList<Integer> lp) {
        long key;
        boolean[] matchedArr = new boolean[pattern.length];
        Arrays.fill(matchedArr, false);

        boolean testedFirst = false;   // NEW: did we actually test i==0?

        for (int i = 0; i < pattern.length; i++) {
            if (level >= lp.get(i)) {
                // this position is enabled at this level
                key = tree.codec.pack(level, interval, pattern[i]);
                if (i == 0) testedFirst = true;  // NEW

                if (!tree.contains(level, key)) {
                    // First observed mismatch at i → ensure we know the true prefix length:
                    for (int j = 0; j < i; j++) {
                        key = tree.codec.pack(level, interval, pattern[j]);
                        matchedArr[j] = tree.contains(level, key);
                        if (j == 0) testedFirst = true;       // NEW: we tested index 0 here
                        if (!matchedArr[j]) break;             // stop at first NO in prefix
                    }
                    return new Probe(countStartConsecutiveMatches(matchedArr), false, testedFirst);
                }
                matchedArr[i] = true;
            }
        }


        return new Probe(countStartConsecutiveMatches(matchedArr), true, testedFirst);
    }

    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll - 1;
        return maxPos >= positionOffset;
    }
}
