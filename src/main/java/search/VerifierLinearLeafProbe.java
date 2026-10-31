package search;

import java.util.ArrayList;

import membership.Membership;
import org.apache.commons.math3.util.Pair;
import tree.ImplicitTree;
import tree.StreamBuffer;

// Checks every start in a candidate range directly against the buffered symbols (the raw data level).
// Positions are relative to the retained trees: tree t holds positions [t.id * span, (t.id + 1) * span).
public class VerifierLinearLeafProbe implements Verifier {

    public VerifierLinearLeafProbe() {}
    int leafProbesCounter = 0;

    @Override
    public Pair<ArrayList<Integer>, Integer> verify(int currentTreeIdx,
                                                    ArrayList<ImplicitTree<Membership>> trees,
                                                    CandidateRange candidateRange,
                                                    Pattern pat)
    {
        ArrayList<Integer> matches = new ArrayList<>();
        int start = Math.max(0, candidateRange.startPos());
        int stopPos = candidateRange.endPos();
        int next = Math.max(start, stopPos + 1);
        long[] tokens = pat.nGramToLong;
        if (trees.isEmpty() || tokens.length == 0) {
            return new Pair<>(matches, next);
        }

        final int span = trees.get(0).baseIntervalSize();
        final int baseTreeId = trees.get(0).id;

        // Last start whose occurrence would still end inside the indexed data.
        ImplicitTree<Membership> lastTree = trees.get(trees.size() - 1);
        int lastPos = lastTree.id * span + lastTree.buffer.lastIndex();
        int lastStart = Math.min(stopPos, lastPos - tokens.length + 1);

        for (int s = start; s <= lastStart; s++) {
            if (matchesAt(trees, span, baseTreeId, s, tokens)) {
                matches.add(s);
            }
        }
        return new Pair<>(matches, next);
    }

    @Override
    public int getLeafProbes() {
        return this.leafProbesCounter;
    }

    @Override
    public void reset() {
        this.leafProbesCounter = 0;
    }

    // Compare the pattern tokens with the buffered tokens starting at pos, crossing into later trees
    // when needed. Every tree but the newest is full, so a tree ends exactly at its buffer length.
    private boolean matchesAt(ArrayList<ImplicitTree<Membership>> trees,
                              int span,
                              int baseTreeId,
                              int pos,
                              long[] tokens) {
        int treeIdx = pos / span - baseTreeId;
        if (treeIdx < 0 || treeIdx >= trees.size()) {
            return false;
        }
        int local = pos - (treeIdx + baseTreeId) * span;
        StreamBuffer buffer = trees.get(treeIdx).buffer;

        for (long token : tokens) {
            if (local >= buffer.length()) {
                treeIdx++;
                if (treeIdx >= trees.size()) {
                    return false;
                }
                buffer = trees.get(treeIdx).buffer;
                local = 0;
            }
            this.leafProbesCounter++;
            if (buffer.get(local) != token) {
                return false;
            }
            local++;
        }
        return true;
    }
}
