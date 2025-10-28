package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;

/**
 * BlockSearchCharSet:
 * A block search strategy that uses Bloom-filter style membership tests
 * over node intervals, with pruning plan awareness (lp),
 * and fallback logic to ensure we actually test the first token
 * when needed.
 */
public class BlockSearchCharSet implements SearchAlgorithm {

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
                p.charStartLp,
                currentIntervalSize);
        this.currentOffset = positionOffset;

        if (probe.consumed() == 0) {
            this.currentOffset = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval;
            return null;
        } else {
            //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
            boolean canDescend = f.level() + 1 < tree.maxDepth();
            if (probe.complete() && childrenIntervalSize >= p.size && canDescend) {
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

    public int getCurrentOffset() {
        return this.currentOffset;
    }

    // Counts consecutive TRUEs starting at index 0 in matchedArr.
    int countStartConsecutiveMatches(boolean[] matchedArr) {
        int matches = 0;
        for (boolean b : matchedArr) {
            if (!b) {
                break;
            }
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
                int intervalIdx,
                Pattern pattern,
                long[] tokens,
                ArrayList<Integer> lp,
                int currentIntervalSize) {

        final int limit = Math.min(tokens.length, maxTokensForInterval(currentIntervalSize, pattern));
        if (limit == 0) {
            return new Probe(0, false, true); // degenerate, whatever
        }

        final long hi = Integer.toUnsignedLong(intervalIdx); // hoist once

        // ---- 1. Force-check token 0 ----
        long t0 = tokens[0];
        boolean firstContains = tree.contains(level, hi, t0);
        if (!firstContains) {
            return new Probe(0, false, true); // hard no at start
        }

        int confirmedPrefix = 1; // we know token0 hit
        int maxScanned = 1;      // highest token index we've *looked at* in main scan (>= confirmedPrefix)

        // ---- 2. Main scan for i = 1..limit-1 ----
        for (int i = 1; i < limit; i++) {

            // obey lp (skip early if this token isn't allowed yet at this level)
            if (lp != null && lp.size() > i && level < lp.get(i)) {
                // we didn't check it, so we can't extend confirmedPrefix
                maxScanned = i + 1;
                continue;
            }

            long ti = tokens[i];
            boolean contains = tree.contains(level, hi, ti);

            maxScanned = i + 1;

            if (!contains) {
                // miss at i -> repair gap [confirmedPrefix .. i-1]
                int repairedPrefix = confirmedPrefix;
                for (int j = confirmedPrefix; j < i; j++) {
                    // only needed if lp skipped that j
                    // We know: if j < confirmedPrefix, it's already guaranteed good
                    // So we only test j>=confirmedPrefix
                    if (lp != null && lp.size() > j && level < lp.get(j)) {
                        long tj = tokens[j];
                        if (!tree.contains(level, hi, tj)) {
                            int consumedChars = charactersMatched(repairedPrefix, pattern);
                            return new Probe(consumedChars, false, true);
                        }
                        repairedPrefix++;
                    } else {
                        // if we didn't skip j in main scan, that means confirmedPrefix
                        // should have already included j. If it didn't, we won't get here
                        // because confirmedPrefix would have caught up. So nothing to do.
                    }
                }

                int consumedChars = charactersMatched(repairedPrefix, pattern);
                return new Probe(consumedChars, false, true);
            }

            // contains == true
            // can this grow the guaranteed prefix?
            if (i == confirmedPrefix) {
                confirmedPrefix++;
            }
        }

        // ---- 3. No miss in main scan ----
        // We never saw a definitive NO anywhere we were allowed to check.
        // We still need to "repair" any skipped tail between confirmedPrefix..maxScanned-1
        // so we don't lie about prefix length.
        int repairedPrefix = confirmedPrefix;
        for (int j = confirmedPrefix; j < maxScanned; j++) {
            if (lp != null && lp.size() > j && level < lp.get(j)) {
                long tj = tokens[j];
                if (!tree.contains(level, hi, tj)) {
                    int consumedChars = charactersMatched(repairedPrefix, pattern);
                    return new Probe(consumedChars, true, true);
                }
                repairedPrefix++;
            }
        }

        int consumedChars = charactersMatched(repairedPrefix, pattern);
        return new Probe(consumedChars, true, true);
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

    public boolean isValidChild(int positionOffset,
                                int intervalIdx,
                                int level,
                                int maxLevel,
                                int workingTreeIdx) {

        int spanAll = 1 << maxLevel;
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll
                - 1;
        return maxPos >= positionOffset;
    }
}
