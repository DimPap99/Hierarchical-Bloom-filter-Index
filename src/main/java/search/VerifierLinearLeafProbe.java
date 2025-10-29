package search;

import java.util.ArrayList;

import membership.Membership;
import org.apache.commons.math3.util.Pair;
import tree.ImplicitTree;
import tree.StreamBuffer;

public class VerifierLinearLeafProbe implements Verifier {

    public VerifierLinearLeafProbe() {}
    int leafProbesCounter = 0;

    @Override
    public Pair<ArrayList<Integer>, Integer> verify(int currentTreeIdx,
                                                    ArrayList<ImplicitTree<Membership>> trees,
                                                    CandidateRange candidateRange,
                                                    Pattern pat)
    {
        int leafStartIdx = candidateRange.startPos();
        int stopPos      = candidateRange.endPos();

        // span is the number of symbols covered by one ImplicitTree block
        // This is usually 2^(maxDepth-1) or treeLength
        final int span      = trees.get(0).baseIntervalSize();
        final int spanShift = Integer.numberOfTrailingZeros(span);

        // The first active tree in the sliding window might no longer have id == 0
        // We take that into account when mapping absolute positions back to an index in 'trees'
        final int baseTreeId = trees.get(0).id;

        if (leafStartIdx < 0) {
            leafStartIdx = 0;
        }

        ArrayList<Integer> matches = new ArrayList<>();
        MatchResult result;

        // We do not cache workingTree/buffer up here because we may have to jump
        // to the correct relative index if the caller passed a stale currentTreeIdx
        while (leafStartIdx <= stopPos) {

            // Figure out which absolute tree id contains this starting position
            int absoluteTreeId = leafStartIdx >>> spanShift; // == floor(leafStartIdx / span)
            int relIdx         = absoluteTreeId - baseTreeId;

            // If relIdx is outside [0, trees.size()) then we ran past the active window
            if (relIdx < 0 || relIdx >= trees.size()) {
                return new Pair<>(matches, leafStartIdx);
            }

            // Run the naive verifier starting exactly at leafStartIdx
            result = verifyAtLeavesNaive(relIdx, trees, leafStartIdx, pat.nGramToLong);

            if (result.matched) {
                matches.add(result.pos);
                // jump past this match, classic exact string search semantics
                leafStartIdx = result.pos + 1;
            } else {
                // advance one symbol in global coordinates
                leafStartIdx++;
            }
        }

        return new Pair<>(matches, leafStartIdx);
    }

    @Override
    public int getLeafProbes() {
        return this.leafProbesCounter;
    }

    @Override
    public void reset() {
        this.leafProbesCounter = 0;
    }

    public MatchResult verifyAtLeavesNaive(
            int startTreeRelIdx,
            ArrayList<ImplicitTree<Membership>> trees,
            int leafStartIdx,
            long[] pat)
    {
        final int span      = trees.get(0).baseIntervalSize();
        final int spanShift = Integer.numberOfTrailingZeros(span);
        final int m         = pat.length;

        // The first active tree may have id > 0
        final int baseTreeId = trees.get(0).id;

        // Compute the last valid absolute position in the active window
        int globalLastPos = -1;
        if (!trees.isEmpty()) {
            ImplicitTree<Membership> lastTree = trees.get(trees.size() - 1);
            // lastTree.id is the absolute tree id
            // lastTree.buffer.lastIndex() is the last local offset inside that block
            globalLastPos = lastTree.id * span + lastTree.buffer.lastIndex();
        }

        // If the pattern would run past the active window, we can fail fast
        if (leafStartIdx + m - 1 > globalLastPos) {
            return new MatchResult(false, globalLastPos);
        }

        // Initialise to the correct working tree
        int workingTreeRelIdx = startTreeRelIdx;
        ImplicitTree<?> tree  = trees.get(workingTreeRelIdx);
        StreamBuffer buffer   = tree.buffer;
        int bufferLength      = buffer.length();

        for (int i = 0; i < m; i++) {
            int leafIdx = leafStartIdx + i;

            // Compute which absolute tree id holds position leafIdx
            int absTreeId   = leafIdx >>> spanShift;
            int relTreeIdx  = absTreeId - baseTreeId;

            // Crossed into a new ImplicitTree block
            if (relTreeIdx != workingTreeRelIdx) {
                if (relTreeIdx < 0 || relTreeIdx >= trees.size()) {
                    // we left the active window
                    return new MatchResult(false, leafIdx - 1);
                }
                workingTreeRelIdx = relTreeIdx;
                tree              = trees.get(workingTreeRelIdx);
                buffer            = tree.buffer;
                bufferLength      = buffer.length();
            }

            // Map leafIdx to local index inside this tree block
            // Because span is a power of two, leafIdx & (span - 1) is leafIdx % span
            int localPos = leafIdx & (span - 1);

            // Bounds check against the last tree which may be only partially full
            if (localPos >= bufferLength) {
                return new MatchResult(false, leafIdx - 1);
            }

            // Actual comparison
            this.leafProbesCounter++;
            if (buffer.get(localPos) != pat[i]) {
                int failPos = (i == 0) ? leafStartIdx : leafIdx - 1;
                return new MatchResult(false, failPos);
            }
        }

        // All characters matched
        return new MatchResult(true, leafStartIdx);
    }
}
