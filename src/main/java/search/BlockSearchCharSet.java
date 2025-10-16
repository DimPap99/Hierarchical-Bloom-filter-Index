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
        Probe probe = probe(tree, f.level(), f.intervalIdx(), p.nGramToLong, p.charStartLp);
        this.currentOffset = positionOffset;

        int currentIntervalSize  = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval     = tree.baseIntervalSize();
        int intervalEndIdx       = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval - 1;

        // Only skip when we truly saw a NO at the first symbol.
        if (probe.consumed() == 0 && probe.testedFirst()) {
            this.currentOffset = intervalEndIdx + 1;
            return null;
        }

        // If we didn't test first symbol but have no NO yet -> treat as “complete so far”
        boolean noDefinitiveNo = probe.complete() || !probe.testedFirst();

        if (noDefinitiveNo && childrenIntervalSize >= p.size) {
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

    Probe probe(ImplicitTree tree, int level, int interval, long[] pattern, ArrayList<Integer> lp) {
        long key;
        boolean[] matchedArr = new boolean[pattern.length];
        Arrays.fill(matchedArr, false);

        boolean testedFirst = false;   // NEW: did we actually test i==0?
        boolean contains = false;
        for (int i = 0; i < pattern.length; i++) {
            if (level >= lp.get(i)) {
                // this position is enabled at this level
                if (tree.codec.fitsOneWord(interval, pattern[i])) {
                    long w = tree.codec.packWord(interval, pattern[i]);
                    contains = tree.contains(level, w);
                } else {
                    long hi = Integer.toUnsignedLong(interval);
                    long lo = pattern[i];
                    contains = tree.contains(level, hi, lo);
                }
                if (i == 0) testedFirst = true;  // NEW

                if (!contains) {
                    // First observed mismatch at i → ensure we know the true prefix length:
                    for (int j = 0; j < i; j++) {
//                        int prefixSymbol = Math.toIntExact(pattern[j]);
                        if (tree.codec.fitsOneWord(interval, pattern[j])) {
                            long w = tree.codec.packWord(interval, pattern[j]);
                            contains = tree.contains(level, w);
                        } else {
                            long hi = Integer.toUnsignedLong(interval);
                            long lo = pattern[j];
                            contains = tree.contains(level, hi, lo);
                        }
                        matchedArr[j] = contains;
                        if (j == 0) testedFirst = true;
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
