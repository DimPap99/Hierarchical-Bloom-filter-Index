package search;

import membership.Membership;
import org.apache.commons.math3.util.Pair;
import tree.ImplicitTree;
import tree.StreamBuffer;

import java.util.ArrayList;

// Verifies candidate ranges with the Knuth-Morris-Pratt automaton of the pattern's tokens (hbi-optimizations.md,
// Section 7): one left-to-right pass over the buffered tokens, O(range + pattern) per range, where
// VerifierLinearLeafProbe restarts the comparison at every start, O(range * pattern) in the worst case (a run of
// one symbol). A query's ranges arrive in increasing order, so when a range starts where the last scan ended, the
// automaton carries on instead of rescanning. Positions are relative to the retained trees, as in
// VerifierLinearLeafProbe.
public class VerifierKmp implements Verifier {

    private int leafProbes = 0;
    // Per query, built at the first range: the pattern's failure table, the last position fed to the automaton,
    // and the number of tokens it had matched there.
    private int[] failure;
    private int scannedTo;
    private int state;

    @Override
    public void beginQuery(Pattern pattern) {
        failure = null;
    }

    @Override
    public Pair<ArrayList<Integer>, Integer> verify(int currentTreeIdx,
                                                    ArrayList<ImplicitTree<Membership>> trees,
                                                    CandidateRange candidateRange,
                                                    Pattern pat) {
        ArrayList<Integer> matches = new ArrayList<>();
        int start = Math.max(0, candidateRange.startPos());
        int stopPos = candidateRange.endPos();
        int next = Math.max(start, stopPos + 1);
        long[] tokens = pat.nGramToLong;
        int r = tokens.length;
        if (trees.isEmpty() || r == 0) {
            return new Pair<>(matches, next);
        }

        final int span = trees.get(0).baseIntervalSize();
        final int baseTreeId = trees.get(0).id;
        ImplicitTree<Membership> lastTree = trees.get(trees.size() - 1);
        int lastPos = lastTree.id * span + lastTree.buffer.lastIndex();
        int lastStart = Math.min(stopPos, lastPos - r + 1);
        if (lastStart < start) {
            return new Pair<>(matches, next);
        }

        if (failure == null) {
            failure = failureTable(tokens);
            scannedTo = Integer.MIN_VALUE;
        }
        // Continue from the last scan if it reached the position just before this range, else start afresh.
        int from;
        if (scannedTo >= start - 1) {
            from = scannedTo + 1;
        } else {
            from = start;
            state = 0;
        }
        int end = lastStart + r - 1;

        int treeIdx = from / span - baseTreeId;
        int local = from - (treeIdx + baseTreeId) * span;
        StreamBuffer buffer = trees.get(treeIdx).buffer;
        for (int p = from; p <= end; p++) {
            if (local >= buffer.length()) {   // every tree but the newest is full, so it ends at its buffer length
                buffer = trees.get(++treeIdx).buffer;
                local = 0;
            }
            long token = buffer.get(local++);
            leafProbes++;
            while (state > 0 && tokens[state] != token) {
                state = failure[state - 1];
            }
            if (tokens[state] == token) {
                state++;
            }
            if (state == r) {
                int s = p - r + 1;
                if (s >= start) {   // a continued scan can complete matches that start before this range
                    matches.add(s);
                }
                state = failure[r - 1];
            }
        }
        scannedTo = end;
        return new Pair<>(matches, next);
    }

    // failure[i]: length of the longest proper prefix of tokens[0..i] that is also a suffix of it.
    public static int[] failureTable(long[] tokens) {
        int[] failure = new int[tokens.length];
        int k = 0;
        for (int i = 1; i < tokens.length; i++) {
            while (k > 0 && tokens[i] != tokens[k]) {
                k = failure[k - 1];
            }
            if (tokens[i] == tokens[k]) {
                k++;
            }
            failure[i] = k;
        }
        return failure;
    }

    @Override
    public int getLeafProbes() {
        return leafProbes;
    }

    @Override
    public void reset() {
        leafProbes = 0;
    }
}
