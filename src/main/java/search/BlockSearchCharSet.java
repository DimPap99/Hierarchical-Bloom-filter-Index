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

        int limit = Math.min(tokens.length,
                maxTokensForInterval(currentIntervalSize, pattern));

        boolean[] matchedArr = new boolean[limit];
        Arrays.fill(matchedArr, false);

        // ---- 1. Force-check token 0, ignoring lp pruning ----
        long firstTok = tokens[0];
        long hi = Integer.toUnsignedLong(intervalIdx);
        long lo0 = firstTok;
        boolean firstContains = tree.contains(level, hi, lo0);
        // we definitely tested the first token
        boolean testedFirst = true;

        if (!firstContains) {
            // Hard "no" at the very first symbol:
            // pattern cannot start at this interval boundary.
            // consumed == 0, not complete.
            return new Probe(0, /*complete*/ false, /*testedFirst*/ true);
        }
        int prefixCount = 1;
        int consumedChars;


        // MAIN SCAN (respects lp)
        for (int i = 1; i < limit; i++) {

            // pruning plan: skip checking this token i at this level if lp says "too early"
            if (lp != null && lp.size() > i && level < lp.get(i)) {
                continue;
            }

            long token = tokens[i];

            boolean contains;

            long lo = token;
            contains = tree.contains(level, hi, lo);

            if (!contains) {
                // First negative at token i in the main scan.

                // Repair: explicitly check tokens 1..i-1 ignoring lp.

                for (int j = 1; j < i; j++) {
                    if(matchedArr[j]) {
                        prefixCount +=1;
                        continue;
                    }
                    long t0 = tokens[j];

                    boolean c0;

                    lo0 = t0;
                    c0 = tree.contains(level, hi, lo0);


                    if (!c0) {
                        // prefix breaks at j
                        consumedChars = charactersMatched(prefixCount, pattern);
                        return new Probe(consumedChars, /*complete*/ false, /*testedFirst*/ true);
                    }

                    matchedArr[j] = true;
                    prefixCount++;
                }

                // All tokens 0..i-1 were present.
                consumedChars = charactersMatched(prefixCount, pattern);
                return new Probe(consumedChars, /*complete*/ false, /*testedFirst*/ true);

            }

            // still "maybe present" for this token in main scan
            matchedArr[i] = true;
            // Do not increment prefixCount here, because this might not actually be
            // contiguous from 0 if lp skipped early tokens. We only use prefixCount
            // in fallback, where we explicitly walk from 0 upward.
        }

        // If we get here, main scan saw NO Bloom filter negatives at all.
        // complete == true.
        // testedFirst may still be false if lp skipped i == 0 and we never had to fallback.

        consumedChars = charactersMatched(matchedArr.length, pattern);

        return new Probe(consumedChars, /*complete*/ true, /*testedFirst*/ testedFirst);
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