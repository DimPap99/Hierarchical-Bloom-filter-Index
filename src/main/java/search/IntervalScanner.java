package search;

import tree.ImplicitTree;
import java.util.*;

/** DFS over the Bloom hierarchy that yields _candidate_ intervals. */
public final class IntervalScanner implements Iterator<CandidateRange> {

    private final Deque<Frame> stack = new ArrayDeque<>();
    private final ImplicitTree<?> tree;
    private final Pattern pattern;
    public int positionOffset;
    public SearchAlgorithm searchAlgorithm;

    private int seedLevel = -1;
    private int nextSeedIdx = 0;
    private int seedLimit = 0;

    public IntervalScanner(ImplicitTree<?> tree, Pattern pattern, SearchAlgorithm searchAlgorithm, int positionOffset) {

        this.tree    = tree;
        this.pattern = pattern;
        this.searchAlgorithm = searchAlgorithm;
        this.positionOffset = positionOffset;
    }

    /** Initialise the scanner with the starting level for the DFS walk. */
    public void seedLevel(int level) {
        stack.clear();
        this.seedLevel = Math.max(0, level);
        this.nextSeedIdx = 0;
        this.seedLimit = intervalsForLevel(level);
        ensureSeedFrame();
    }

    @Override
    public boolean hasNext() {
        ensureSeedFrame();
        return !stack.isEmpty();
    }

    @Override
    public CandidateRange next() {
        while (true) {
            ensureSeedFrame();
            if (stack.isEmpty()) {
                return null;
            }
            Frame f = stack.pop();
            CandidateRange r = this.searchAlgorithm.search(f, pattern, this.tree, this.stack, this.positionOffset);
            this.positionOffset = this.searchAlgorithm.getCurrentOffset();
            if (r != null) {
                return r;
            }
        }
    }

    private void ensureSeedFrame() {
        if (!stack.isEmpty()) {
            return;
        }
        if (seedLevel < 0) {
            return;
        }
        if (nextSeedIdx < seedLimit) {
            stack.addLast(new Frame(seedLevel, nextSeedIdx++));
        }
        if (nextSeedIdx >= seedLimit) {
            seedLevel = -1;
        }
    }

    private static int intervalsForLevel(int level) {
        if (level <= 0) {
            return 1;
        }
        if (level >= 31) {
            if (level >= 63) {
                return Integer.MAX_VALUE;
            }
            long intervals = 1L << level;
            if (intervals >= Integer.MAX_VALUE) {
                return Integer.MAX_VALUE;
            }
            return (int) intervals;
        }
        return 1 << level;
    }
}
