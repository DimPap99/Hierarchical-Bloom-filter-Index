package search;

import java.util.ArrayList;

import membership.Membership;
import org.apache.commons.math3.util.Pair;
import tree.ImplicitTree;

/** Verifier that performs "leaf" probing exhaustively over a range. Testing all possible alignment of the pattern
    until all possible */
public class VerifierLinearLeafProbe implements Verifier {

    public VerifierLinearLeafProbe() {}   // static utility


    public Pair<ArrayList<Integer>, Integer> verify(int currentTreeIdx,
                                                    ArrayList<ImplicitTree<Membership>> trees,
                                                    CandidateRange candidateRange,
                                                    Pattern pat){
        int leafStartIdx = candidateRange.startPos();
        int stopPos = candidateRange.endPos();
        int span       = trees.get(0).baseIntervalSize();              // same for every tree
        int spanShift  = Integer.numberOfTrailingZeros(span);  // log₂(span)

        if (leafStartIdx < 0) leafStartIdx = 0;

        ArrayList<Integer> matches = new ArrayList<>();
        MatchResult result;

        /* current working tree – will be reloaded whenever we cross a boundary */
        ImplicitTree workingTree = trees.get(currentTreeIdx);

        while (leafStartIdx <= stopPos) {

            /* ---------- 1. determine which tree the leaf belongs to -------- */
            int treeIdx = leafStartIdx >>> spanShift;      // == leafStartIdx / SPAN

            if (treeIdx != currentTreeIdx) {                // crossed ≥ 1 tree(s)
                if (treeIdx >= trees.size()) {              // hit stream end
                    return new Pair<>(matches, leafStartIdx);
                }
                currentTreeIdx = treeIdx;
                workingTree    = trees.get(currentTreeIdx);
            }

            /* ---------- 2. run the naive verifier ------------------------- */
            result = verifyAtLeavesNaive(currentTreeIdx, trees, leafStartIdx, pat.nGramToInt);

            if (result.matched) {                           // full hit
                matches.add(result.pos);
                leafStartIdx = result.pos + 1;              // jump past hit
            } else {
                leafStartIdx++;                             // slide window
            }
        }
        return new Pair<>(matches, leafStartIdx);
    }

    public MatchResult verifyAtLeavesNaive(
            int currentTreeIdx,
            ArrayList<ImplicitTree<Membership>> trees,
            int leafStartIdx,
            int[] pat)
    {
        /* ---------- constants (same for every tree) ------------------- */
        final int span      = trees.get(0).baseIntervalSize();               // 2^maxDepth
        final int spanShift = Integer.numberOfTrailingZeros(span);   // log₂(span)
        final int m         = pat.length;

        /* ---------- quick global bound check -------------------------- */
        /* Highest usable absolute position in the *whole* stream */
        int globalLastPos = (trees.size() - 1) * span
                + trees.get(trees.size() - 1).buffer.data.size() - 1;

        if (leafStartIdx + m - 1 > globalLastPos)
            return new MatchResult(false, globalLastPos);

        /* ---------- initialise current tree --------------------------- */
        int workingTreeIdx  = currentTreeIdx;
        ImplicitTree tree   = trees.get(workingTreeIdx);
        ArrayList<Integer> sb    = tree.buffer.data;                 // convenience ref

        /* ---------- main loop ----------------------------------------- */
        for (int i = 0; i < m; i++) {
            int leafIdx  = leafStartIdx + i;
            int treeIdx  = leafIdx >>> spanShift;          // == leafIdx / span

            /* --- did we cross into a new implicit-tree block? --------- */
            if (treeIdx != workingTreeIdx) {
                if (treeIdx >= trees.size())               // safety net
                    return new MatchResult(false, leafIdx - 1);

                workingTreeIdx = treeIdx;
                tree           = trees.get(treeIdx);
                sb             = tree.buffer.data;              // update ref
            }

            /* --- local offset inside this tree’s stream --------------- */
            int localPos = leafIdx & (span - 1);           // faster than % span

            /* --- bounds-check (last tree may be shorter than span) ---- */
            if (localPos >= sb.size())
                return new MatchResult(false, leafIdx - 1);

            /* --- actual character comparison -------------------------- */
            if (sb.get(localPos) != pat[i]) {
                int failPos = (i == 0) ? leafStartIdx : leafIdx - 1;
                return new MatchResult(false, failPos);
            }
        }

        return new MatchResult(true, leafStartIdx);
    }

}
