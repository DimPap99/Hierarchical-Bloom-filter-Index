package search;

import tree.ImplicitTree;
import java.util.*;

/** DFS over the Bloom hierarchy that yields _candidate_ intervals. */
public final class IntervalScanner implements Iterator<CandidateRange> {

    private Deque<Frame> stack = new ArrayDeque<>();
    private final ImplicitTree<?> tree;
    private final Pattern pattern;
    public int positionOffset;
    public SearchAlgorithm searchAlgorithm;
    public IntervalScanner(ImplicitTree<?> tree, Pattern pattern, SearchAlgorithm searchAlgorithm, int positionOffset) {

        this.tree    = tree;
        this.pattern = pattern;
        this.searchAlgorithm = searchAlgorithm;
        this.positionOffset = positionOffset;
    }
    // Seeds the stack with initial level(s) if needed
    public void seedStack(Deque<Frame> stack) {
        this.stack = stack;
    }
    @Override
    public boolean hasNext() { return !stack.isEmpty(); }

    @Override
    public CandidateRange next() {
        while (!stack.isEmpty()) {
            Frame f = stack.pop();
            CandidateRange r = this.searchAlgorithm.search(f, pattern, this.tree, this.stack, this.positionOffset);
            this.positionOffset = this.searchAlgorithm.getCurrentOffset();
            if (r == null) continue;
            else {
//                if(pattern.nGramToInt.length > 8){
//                    int pos = this.tree.traverse(f, this.tree.maxDepth()-1, this.positionOffset, pattern.nGramToInt[0]);
//                    if(pos != -1) return new CandidateRange(pos, r.endPos());
//                }
                return r;
            }
        }
        return null;
    }


}
