package search;

import tree.ImplicitTree;
import estimators.Estimator;

import java.util.*;

public class BlockSearch implements SearchAlgorithm{
    public Estimator estimator;
    public int currentOffset;
    private boolean strides;

    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        int currentIntervalSize = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval = tree.baseIntervalSize();

        long[] tokens = strides ? p.effectiveNgramArr : p.nGramToLong;
        Probe probe = probe(tree,
                f.level(),
                f.intervalIdx(),
                p,
                tokens,
                currentIntervalSize);
        this.currentOffset = positionOffset;

        if (probe.consumed() == 0) {
            this.currentOffset = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval;
            return null;
        } else {
            //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
            if (probe.complete() && childrenIntervalSize >= p.size) {
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
                int currentIntervalSize) {
        //p72xcHxRu1

        int matches = 0;
        int effectiveLookupRange = Math.min(tokens.length, maxTokensForInterval(currentIntervalSize, pattern));
        boolean contains = false;
        for (int i = 0; i < effectiveLookupRange; i+=1) {
            long token = tokens[i];
//            if (tree.codec.fitsOneWord(interval, token)) {
//                long w = tree.codec.packWord(interval, token);
//                contains = tree.contains(level, w);
//            } else {
                long hi = Integer.toUnsignedLong(interval);
                long lo = token;
                contains = tree.contains(level, hi, lo);
//            }

            if (!contains) {
                int consumed = charactersMatched(matches, pattern);
                return new Probe(consumed, false);          // first mismatch at i
            }
            matches+=1;
            if(matches == effectiveLookupRange){
                return new Probe(pattern.originalSz, true);
            }
        }
        return new Probe(pattern.originalSz, true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;              // == tree.intervalSize

        // end position of current interval inside the *global* stream
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll
                - 1;

        return maxPos >= positionOffset;
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
}
