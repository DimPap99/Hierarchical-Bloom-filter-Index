package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Deque;

// Confidence-Gated Block Search (Sec. 4.2.4): block search that only tests the tokens whose gate
// level has been reached, and always tests the first token to anchor the alignment.
public class BlockSearchCharSet implements SearchAlgorithm {

    public Estimator estimator;
    public int currentOffset;
    private boolean strides;
    private boolean useTokenLengthGuard = false;

    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        int currentIntervalSize = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval = tree.baseIntervalSize();

        int intervalStartIdx = tree.id * treeBaseInterval + currentIntervalSize * f.intervalIdx();
        int intervalEndIdx = intervalStartIdx + currentIntervalSize - 1;

        // Candidate ranges never start before the frontier or before this interval.
        this.currentOffset = Math.max(positionOffset, intervalStartIdx);
        if (intervalEndIdx < positionOffset) {
            this.currentOffset = positionOffset;
            return null;
        }

        long[] tokens = strides ? p.effectiveNgramArr : p.nGramToLong;
        Probe probe = probe(tree,
                f.level(),
                f.intervalIdx(),
                p,
                tokens,
                p.charStartLp);
        int requiredSpan = useTokenLengthGuard ? tokens.length : p.size;
        requiredSpan = Math.max(1, requiredSpan);
        boolean canDescend = f.level() + 1 < tree.maxDepth();

        if (probe.consumed() == 0) {
            // The first token is absent, so no occurrence can start inside this interval.
            this.currentOffset = intervalEndIdx + 1;
            return null;
        }
        //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
        if (probe.complete() && canDescend && childrenIntervalSize >= requiredSpan) {
            //we just generate children as theres a chance that the pattern is in both of them
            tree.generateChildren(f, stack, positionOffset, tree.id);
            return null;
        }
        if (probe.complete()) {
            return new CandidateRange(this.currentOffset, intervalEndIdx);
        }
        //the probes were incomplete. Meaning:
        //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
        //neighboring interval (we have an overlap). Starts before this interval are handled by earlier intervals.
        this.currentOffset = Math.max(this.currentOffset, intervalEndIdx - probe.consumed() + 1);
        return new CandidateRange(this.currentOffset, intervalEndIdx);
    }

    public int getCurrentOffset() {
        return this.currentOffset;
    }

    @Override
    public void setStrides(boolean strides) {
        this.strides = strides;
    }

    @Override
    public boolean usesStrides() {
        return strides;
    }

    public void setUseTokenLengthGuard(boolean useTokenLengthGuard) {
        this.useTokenLengthGuard = useTokenLengthGuard;
    }

    // Algorithm 2 (TestActive). Returns the number of pattern characters covered by the longest prefix
    // of tokens that tested positive, and whether no gated-in token tested negative.
    Probe probe(ImplicitTree tree,
                int level,
                int intervalIdx,
                Pattern pattern,
                long[] tokens,
                ArrayList<Integer> lp) {

        int limit = tokens.length;
        boolean[] matchedArr = new boolean[limit];

        // Force-check token 0, ignoring lp pruning.
        long hi = Integer.toUnsignedLong(intervalIdx);
        if (!tree.contains(level, hi, tokens[0])) {
            // Hard "no" at the very first symbol:
            // pattern cannot start at this interval boundary.
            return new Probe(0, false, true);
        }
        int prefixCount = 1;

        // MAIN SCAN (respects lp)
        for (int i = 1; i < limit; i++) {

            // pruning plan: skip checking this token i at this level if lp says "too early"
            if (lp != null && lp.size() > i && level < lp.get(i)) {
                continue;
            }

            if (!tree.contains(level, hi, tokens[i])) {
                // First negative at token i in the main scan.
                // Repair: explicitly check tokens 1..i-1 ignoring lp to find where the prefix breaks.
                for (int j = 1; j < i; j++) {
                    if (!matchedArr[j] && !tree.contains(level, hi, tokens[j])) {
                        // prefix breaks at j
                        return new Probe(charactersMatched(prefixCount, pattern), false, true);
                    }
                    matchedArr[j] = true;
                    prefixCount++;
                }

                // All tokens 0..i-1 were present.
                return new Probe(charactersMatched(prefixCount, pattern), false, true);
            }

            // still "maybe present" for this token in main scan
            matchedArr[i] = true;
            // Do not increment prefixCount here, because this might not actually be
            // contiguous from 0 if lp skipped early tokens. We only use prefixCount
            // in the repair walk, which goes from 1 upward.
        }

        // The main scan saw no Bloom filter negatives at all.
        return new Probe(charactersMatched(limit, pattern), true, true);
    }

    // Upper bound on how far into the pattern the failing token reaches, in pattern positions.
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
