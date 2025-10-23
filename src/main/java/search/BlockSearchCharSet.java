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
    private boolean strides;

    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        int currentIntervalSize  = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval     = tree.baseIntervalSize();
        int intervalEndIdx       = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval - 1;

        long[] tokens = strides ? p.effectiveNgramArr : p.nGramToLong;
        Probe probe = probe(tree, f.level(), f.intervalIdx(), p, tokens, p.charStartLp, currentIntervalSize);
        this.currentOffset = positionOffset;

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

    @Override
    public void setStrides(boolean strides) {
        this.strides = strides;
    }

    @Override
    public boolean usesStrides() {
        return strides;
    }

    Probe probe(ImplicitTree tree,
                int level,
                int interval,
                Pattern pattern,
                long[] tokens,
                ArrayList<Integer> lp,
                int currentIntervalSize) {
        int limit = Math.min(tokens.length, maxTokensForInterval(currentIntervalSize, pattern));
        boolean[] matchedArr = new boolean[limit];
        Arrays.fill(matchedArr, false);

        boolean testedFirst = false;   // NEW: did we actually test i==0?
        boolean contains = false;
        for (int i = 0; i < limit; i++) {
            if(i > 0){
                int b = 1;
            }
            if (lp != null && lp.size() > i && level < lp.get(i)) {
                continue;
            }
            long token = tokens[i];
            if (tree.codec.fitsOneWord(interval, token)) {
                long w = tree.codec.packWord(interval, token);
                contains = tree.contains(level, w);
            } else {
                long hi = Integer.toUnsignedLong(interval);
                long lo = token;
                contains = tree.contains(level, hi, lo);
            }
            if (i == 0) testedFirst = true;  // NEW

            if (!contains) {
                int prefixMatches = countStartConsecutiveMatches(matchedArr);
                int consumed = charactersMatched(prefixMatches, pattern);
                return new Probe(consumed, false, testedFirst);
            }
            matchedArr[i] = true;
        }


        int prefixMatches = countStartConsecutiveMatches(matchedArr);
        int consumed = charactersMatched(prefixMatches, pattern);
        return new Probe(consumed, true, testedFirst);
    }

    private int charactersMatched(int matchedTokens, Pattern pattern) {
        if (matchedTokens <= 0) {
            return 0;
        }

        int nGramSize = Math.max(1, pattern.nGram);
        if (strides) {
            long chars = (long) matchedTokens * nGramSize;
            return (int) Math.min(pattern.originalSz, chars);
        }

        if (nGramSize == 1) {
            return Math.min(pattern.originalSz, matchedTokens);
        }

        long chars = (long) matchedTokens + nGramSize - 1L;
        return (int) Math.min(pattern.originalSz, chars);
    }

    private int maxTokensForInterval(int intervalSize, Pattern pattern) {
        if (!strides) {
            return intervalSize;
        }

        int nGramSize = Math.max(1, pattern.nGram);
        if (nGramSize <= 1) {
            return intervalSize;
        }

        int full = intervalSize / nGramSize;
        int rem = intervalSize % nGramSize;
        return full + (rem > 0 ? 1 : 0);
    }

    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll - 1;
        return maxPos >= positionOffset;
    }
}
