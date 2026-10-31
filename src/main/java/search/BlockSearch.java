package search;

import tree.ImplicitTree;
import estimators.Estimator;

import java.util.*;

public class BlockSearch implements SearchAlgorithm{
    public Estimator estimator;
    public int currentOffset;
    private boolean strides;
    private boolean useTokenLengthGuard = false;

    // Optional probe strategies (hbi-optimizations.md, Section 6). All off by default: measured on the benchmark
    // scenarios, none of them was faster everywhere.
    // - queryPlan: test each distinct token once, in QueryPlan order (smallest offset of each token).
    // - sketchOrder: with queryPlan, start the plan rarest first by the tree's frequency sketch, when it has one.
    // - moveToFront: with queryPlan, a token that tests negative is tested first from then on.
    // - neighbourCheck: when a probe leaves a tail of starts, first test the failing token in the next interval.
    public record Options(boolean queryPlan, boolean sketchOrder, boolean moveToFront, boolean neighbourCheck) {
        public static final Options DEFAULT = new Options(false, false, false, false);
    }
    private Options options = Options.DEFAULT;
    // The token whose absence ended the last probe (for neighbourCheck).
    private long failedToken;

    public void setOptions(Options options) {
        this.options = Objects.requireNonNull(options, "options");
    }

    @Override
    public CandidateRange search(Frame f,
                                 Pattern p,
                                 ImplicitTree tree,
                                 Deque<Frame> stack,
                                 int positionOffset) {

        int currentIntervalSize   = tree.intervalSize(f.level());
        int childrenIntervalSize  = currentIntervalSize / 2;
        int treeBaseInterval      = tree.baseIntervalSize();

        // Start and end of this interval among the retained trees
        int globalStart = tree.id * treeBaseInterval
                + currentIntervalSize * f.intervalIdx();

        int globalEnd   = globalStart + currentIntervalSize - 1;

        // Set currentOffset so that candidate ranges never start before frontier
        this.currentOffset = Math.max(positionOffset, globalStart);

        // If this node's entire interval lies strictly before positionOffset,
        // we can skip probing this node entirely
        if (globalEnd < positionOffset) {
            this.currentOffset = positionOffset;
            return null;
        }

        long[] tokens = strides ? p.effectiveNgramArr : p.nGramToLong;

        Probe probe = options.queryPlan()
                ? probeInPlanOrder(tree, f.level(), f.intervalIdx(),
                        p.plan(strides, options.sketchOrder() ? tree.estimator : null), p.originalSz, currentIntervalSize)
                : probe(tree, f.level(), f.intervalIdx(), p, tokens, currentIntervalSize);

        int requiredSpan = useTokenLengthGuard ? tokens.length : p.size;
        requiredSpan = Math.max(1, requiredSpan);
        boolean canDescend = f.level() + 1 < tree.maxDepth();

        if (probe.consumed() == 0) {
            // Bloom said "no match at this interval's left edge"
            // We know everything up to globalEnd is dead
            this.currentOffset = globalEnd + 1;
            return null;
        }

        if (probe.complete()
                && canDescend
                && childrenIntervalSize >= requiredSpan
                ) {

            tree.generateChildren(f, stack, positionOffset, tree.id);
            return null;
        }
        if (probe.complete()) {
            return new CandidateRange(this.currentOffset, globalEnd);
        } else {
            // Only starts near the right edge can still cross into the next interval. Starts before
            // this interval belong to earlier intervals, so the range never reaches past currentOffset.
            if (options.neighbourCheck() && absentFromNextInterval(tree, f)) {
                this.currentOffset = globalEnd + 1;
                return null;
            }
            this.currentOffset = Math.max(this.currentOffset, globalEnd - probe.consumed() + 1);
            return new CandidateRange(this.currentOffset, globalEnd);
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

    public void setUseTokenLengthGuard(boolean useTokenLengthGuard) {
        this.useTokenLengthGuard = useTokenLengthGuard;
    }

    Probe probe(ImplicitTree tree,
                int level,
                int interval,
                Pattern pattern,
                long[] tokens,
                int currentIntervalSize) {
        //p72xcHxRu1

        int matches = 0;
        int effectiveLookupRange =
                Math.min(tokens.length, maxTokensForInterval(currentIntervalSize, pattern));
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
                failedToken = token;
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


    // Option queryPlan: the plan's tokens that can lie in the interval (offset < its size), in plan order. The
    // first absent token at offset o ends the probe: an occurrence starting here must start in the last o positions.
    Probe probeInPlanOrder(ImplicitTree tree,
                           int level,
                           int interval,
                           QueryPlan plan,
                           int patternLength,
                           int currentIntervalSize) {
        long hi = Integer.toUnsignedLong(interval);
        for (int i = 0; i < plan.size(); i++) {
            int offset = plan.offset(i);
            if (offset >= currentIntervalSize) {
                continue;
            }
            if (!tree.contains(level, hi, plan.token(i))) {
                failedToken = plan.token(i);
                if (options.moveToFront()) {
                    plan.failed(i);
                }
                return new Probe(offset, false);
            }
        }
        return new Probe(patternLength, true);
    }

    // Option neighbourCheck: whether the token that ended the probe is absent from the next interval of the same
    // level in this tree. Every start in the tail would need it there. The last interval of a level is not checked.
    private boolean absentFromNextInterval(ImplicitTree tree, Frame f) {
        int next = f.intervalIdx() + 1;
        if (next >= (1 << f.level())) {
            return false;
        }
        return !tree.contains(f.level(), Integer.toUnsignedLong(next), failedToken);
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
