package algorithms;

import PMIndex.ImplicitTree;
import org.apache.commons.math3.util.Pair;
import utilities.Utils;

import java.lang.reflect.Array;
import java.util.*;

public class BlockSearch implements SearchAlgorithm{

    public class GlobalInfo{
        public boolean isSplit = false;
        public int matched = 0;
        public int patternLength;
        public boolean isReset = false;
        public int positionOffset = -1;
        public GlobalInfo(int patternLength){
            this.patternLength = patternLength;
        }


    }
    static final class Frame {
        final int level;          // depth
        final int intervalIdx;       // index within that depth
//        final int charsToMatch;
//        final boolean isRemainder;
        Frame(int level, int interval){//, int charsToMatch, boolean isRemainder ) {
            this.level = level;
            this.intervalIdx = interval;
        }
    }

    record Probe(int consumed, boolean complete) {} // complete==true if no miss in this node


    public class MatchResult{
        public boolean matched;
        public int charsMatched;
        public int currentTreeIdx;
        public int pos;
        public boolean reset = false;
        MatchResult(boolean matched, int pos, int charsMatched){
            this.matched = matched;
            this.pos = pos;
            this.charsMatched = charsMatched;
        }
    }



    Probe probe(ImplicitTree tree, int level, int interval, char[] pattern, int offset, int remainder) {


        int matches = 0;

        long key;
        for (int i = offset; i < pattern.length; i++) {
            key =  tree.keyService.pack(level, interval, pattern[i]);

            if (!tree.membership.contains(key)) {
                return new Probe(matches, false);          // first mismatch at i
            }
            matches++;
            if(matches == remainder){
                return new Probe(remainder, true);
            }
        }
        return new Probe(pattern.length, true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1) + workingTreeIdx * 8 - 1;
        if(maxPos >= positionOffset) return true;
        else return false;
    }

    private static int[] prefixFunction(char[] pat) {
        int[] pi = new int[pat.length];
        int k = 0;
        for (int i = 1; i < pat.length; ++i) {
            while (k > 0 && pat[k] != pat[i]) k = pi[k - 1];
            if (pat[k] == pat[i]) ++k;
            pi[i] = k;
        }
        return pi;
    }

//    public Pair<ArrayList<Integer>, Integer> exhaustChild(int currentTreeIdx, ArrayList<ImplicitTree> trees, int intervalSize, int intervalIdx, int leafStartIdx, char[] pat, int[] pi, int positionChecks) {
//        int workingTreeIdx = currentTreeIdx;
//        ImplicitTree workingTree = trees.get(currentTreeIdx);
//        int stopPos = leafStartIdx + positionChecks;
//        ArrayList<Integer> matches = new ArrayList<>();
//        MatchResult result; //= verifyAtLeaves(t, leafStartIdx, pat, pi);
//        if(leafStartIdx == -1){ leafStartIdx = 0; }
//
//        while(leafStartIdx <= stopPos){
//            int idxx = (int)(leafStartIdx / (workingTree.maxIdx + 1));
//            if((int)(leafStartIdx / (workingTree.maxIdx + 1)) != workingTreeIdx) {
//                workingTreeIdx++;
//                if (workingTreeIdx < trees.size()) {
//                    workingTree = trees.get(workingTreeIdx);
//                } else {
//                    //No more stream characters to check || We are the end of the stream
//                    return new Pair<ArrayList<Integer>, Integer>(matches, leafStartIdx);
//
//                }
//            }
//            result = verifyAtLeavesNaive(currentTreeIdx, trees, leafStartIdx, pat);
//            if(result.matched){
//                matches.add(result.pos + 1 - pat.length);
//                leafStartIdx = result.pos + 1;
//            }else{
//                leafStartIdx++;
//            }
//        }
//        return new Pair<ArrayList<Integer>, Integer>(matches, leafStartIdx);
//    }
//
//    public MatchResult verifyAtLeavesNaive(int currentTreeIdx, ArrayList<ImplicitTree> trees, int leafStartIdx,
//                                           char[] pat) {
//        int m = pat.length;
//        int workingTreeIdx = currentTreeIdx;
//        int leafIdx;
//        ImplicitTree workingTree = trees.get(currentTreeIdx);;
//        for (int i = 0; i < m; i++) {
//            // build the composite key for leaf index = leafStartIdx + i,
//            // and the single pattern character pat[i]
//            leafIdx = leafStartIdx + i;
//            if((int)(leafIdx / (workingTree.maxIdx + 1)) != workingTreeIdx){
//                workingTreeIdx++;
//                if(workingTreeIdx < trees.size()){
//                    workingTree = trees.get(workingTreeIdx);
//                }else{
//                    //No more stream characters to check || We are the end of the stream
//                    return new MatchResult(false,
//                            leafStartIdx + i - 1,
//                            i);
//                }
//
//            }
//            String key = workingTree.createCompositeKey(workingTree.maxDepth,
//                    leafIdx,
//                    pat[i]);
//            if (!workingTree.membership.contains(key)) {
//                // mismatch at the i‐th character of pat:
//                //   - pos = absolute leaf where it failed = leafStartIdx + i
//                //   - charsMatched = i (how many matched before failing)
//                if(i == 0){
//                    return new MatchResult(false,
//                            leafStartIdx,
//                            i);
//                }
//                else{
//                    return new MatchResult(false,
//                            leafStartIdx + i - 1,
//                            i);
//                }
//
//            }
//        }
//        // if we reach here, all m characters matched in a row
//        int lastPos = leafStartIdx + (m - 1);
//        return new MatchResult(true,
//                lastPos,
//                m);
//    }

    /* If you don’t have a global constant, compute them once in the ctor:
     *   SPAN       = trees.get(0).maxIdx + 1;
     *   SPAN_SHIFT = Integer.numberOfTrailingZeros(SPAN);
     */

    /* ================================================================== */
    /*  1.  exhaustChild – division-free                                   */
    /* ================================================================== */
    public Pair<ArrayList<Integer>, Integer> exhaustChild(
            int currentTreeIdx,
            ArrayList<ImplicitTree> trees,
            int intervalSize,           //  not touched
            int intervalIdx,            //  not touched
            int leafStartIdx,
            char[] pat,
            int[] pi,
            int positionChecks) {
        int span       = trees.get(0).maxIdx + 1;              // same for every tree
        int spanShift  = Integer.numberOfTrailingZeros(span);  // log₂(span)

        if (leafStartIdx < 0) leafStartIdx = 0;
        int stopPos = leafStartIdx + positionChecks;

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
            result = verifyAtLeavesNaive(currentTreeIdx, trees, leafStartIdx, pat);

            if (result.matched) {                           // full hit
                matches.add(result.pos);
                leafStartIdx = result.pos + 1;              // jump past hit
            } else {
                leafStartIdx++;                             // slide window
            }
        }
        return new Pair<>(matches, leafStartIdx);
    }

    /* ================================================================== */
    /*  2.  verifyAtLeavesNaive – division-free                            */
    /* ================================================================== */
    /**
     * Character-by-character verifier, but reading directly from
     * workingTree.stream at the deepest level instead of using the Bloom filter.
     *
     * @return MatchResult  – identical contract as before
     */
    public MatchResult verifyAtLeavesNaive(
            int currentTreeIdx,
            ArrayList<ImplicitTree> trees,
            int leafStartIdx,
            char[] pat)
    {
        /* ---------- constants (same for every tree) ------------------- */
        final int span      = trees.get(0).maxIdx + 1;               // 2^maxDepth
        final int spanShift = Integer.numberOfTrailingZeros(span);   // log₂(span)
        final int m         = pat.length;

        /* ---------- quick global bound check -------------------------- */
        /* Highest usable absolute position in the *whole* stream */
        int globalLastPos = (trees.size() - 1) * span
                + trees.get(trees.size() - 1).stream.length() - 1;

        if (leafStartIdx + m - 1 > globalLastPos)
            return new MatchResult(false, globalLastPos, 0);

        /* ---------- initialise current tree --------------------------- */
        int workingTreeIdx  = currentTreeIdx;
        ImplicitTree tree   = trees.get(workingTreeIdx);
        StringBuilder sb    = tree.stream;                 // convenience ref

        /* ---------- main loop ----------------------------------------- */
        for (int i = 0; i < m; i++) {
            int leafIdx  = leafStartIdx + i;
            int treeIdx  = leafIdx >>> spanShift;          // == leafIdx / span

            /* --- did we cross into a new implicit-tree block? --------- */
            if (treeIdx != workingTreeIdx) {
                if (treeIdx >= trees.size())               // safety net
                    return new MatchResult(false, leafIdx - 1, i);

                workingTreeIdx = treeIdx;
                tree           = trees.get(treeIdx);
                sb             = tree.stream;              // update ref
            }

            /* --- local offset inside this tree’s stream --------------- */
            int localPos = leafIdx & (span - 1);           // faster than % span

            /* --- bounds-check (last tree may be shorter than span) ---- */
            if (localPos >= sb.length())
                return new MatchResult(false, leafIdx - 1, i);

            /* --- actual character comparison -------------------------- */
            if (sb.charAt(localPos) != pat[i]) {
                int failPos = (i == 0) ? leafStartIdx : leafIdx - 1;
                return new MatchResult(false, failPos, i);
            }
        }

        return new MatchResult(true, leafStartIdx, m);
    }


    //    public MatchResult verifyAtLeaves(int currentTreeIdx, ArrayList<ImplicitTree> trees, int leafStartIdx,
//                                              char[] pat,
//                                              int[] pi) {
//        int i = 0;               // index in text (leaf offset)
//        int j = 0;               // index in pattern
//        int m = pat.length*2 - 1;
//
//        ImplicitTree t = trees.get(currentTreeIdx);
//        while (i < m) {
//            String key = t.createCompositeKey(t.maxDepth, leafStartIdx + i, String.valueOf(pat[j]));
//            if (t.membership.contains(key)) {
//                ++i; ++j;
//                if (j == m) {   // full match
//                    int lastPos = leafStartIdx + m - 1;
//                    return new MatchResult(true, lastPos, m);
//                }
//            } else {
//                if (j == 0) {   // miss on first char
//                    return new MatchResult(false, leafStartIdx + i, 0);
//                }
//                j = pi[j - 1];  // fallback;
//            }
//        }
//        // should not reach, but return mismatch safeguard
//        return new MatchResult(false, leafStartIdx + i - 1, j);
//    }
    //generate children for the current frame. Essentially get the right/left split of the current interval
    void generateChildren(Frame currentFrame, Deque<Frame> framesStack, int positionOffset, ImplicitTree tree, int workingTreeIdx){
        int rightChild = tree.getRightChild(currentFrame.intervalIdx);
        int leftChild = tree.getLeftChild(currentFrame.intervalIdx);

        Frame leftFrame = new Frame(currentFrame.level + 1, leftChild);
        Frame rightFrame = new Frame(currentFrame.level + 1, rightChild);

        //add right child
        if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth, workingTreeIdx)) {
            framesStack.push(rightFrame);
        }
        //add left child
        if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth, workingTreeIdx)) {
            framesStack.push(leftFrame);
        }

    }


    public ArrayList<Integer> report(char[] key, ArrayList<ImplicitTree> trees, boolean existence) {
        ArrayList<Integer> results = new ArrayList<>();
        char[] pat = key;
        int[] pi   = prefixFunction(pat);
        Probe bfProbe;
        int currentTreeIdx = 0;
        GlobalInfo info = new GlobalInfo(key.length);
        for(ImplicitTree tree : trees){
            int iterationsUntilPossible =0;
            int childrenIntervalSize;
            int currentIntervalSize;
            int belongingTree = (int)(info.positionOffset  / tree.intervalSize);
            Deque<Frame> stack = new ArrayDeque<>();
            //check root
            Frame rootFrame = new Frame(0, 0);
            bfProbe = probe(tree,  rootFrame.level, rootFrame.intervalIdx, key, 0, key.length - info.matched);
            if(bfProbe.complete){
                generateChildren(rootFrame, stack, info.positionOffset, tree, currentTreeIdx);
            }
            while(!stack.isEmpty() && belongingTree == currentTreeIdx) {
                Frame currentFrame = stack.pop();

                currentIntervalSize = tree.getIntervalSize(currentFrame.level);
                childrenIntervalSize = currentIntervalSize/2;
                bfProbe = probe(tree,  currentFrame.level, currentFrame.intervalIdx, key, info.matched, key.length - info.matched);
//                System.out.println("Checking " + Utils.intervalWithHashes(4, currentFrame.level, currentFrame.intervalIdx*currentIntervalSize)
//                        + " Matched: " + info.matched + " Probe: " + bfProbe.consumed);
                //skip entiner interval - the pattern is not inside it
                if(bfProbe.consumed == 0){
                    info.positionOffset = currentIntervalSize * (currentFrame.intervalIdx + 1) + (currentTreeIdx)*8;
                    info.matched = 0;
                    iterationsUntilPossible++;
                }else{
                    //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
                    if(bfProbe.complete && childrenIntervalSize >= key.length) {
                        //we just generate children as theres a chance that the pattern is in both of them
                        generateChildren(currentFrame, stack, info.positionOffset, tree, currentTreeIdx);
                        iterationsUntilPossible++;
                    }
                    else if(bfProbe.complete && childrenIntervalSize < key.length ){
                        //In this case we know that the probe matched all characters in the current interval.
                        //But also that the children intervals do not completely fit the pattern. If the pattern exists
                        //It will overlap between Left Interval and Right.

                        Pair<ArrayList<Integer>, Integer> res =  exhaustChild(currentTreeIdx, trees, currentIntervalSize, currentFrame.intervalIdx, info.positionOffset, pat, pi, currentIntervalSize);

                        for(int r : res.getFirst()){
                            results.add(r);
                        }
                        info.positionOffset = res.getSecond();
                        belongingTree = (int)(info.positionOffset  / tree.intervalSize);
                        iterationsUntilPossible = 0;
                    }
                    else{ //the probes were inclomplete. Meaning:
                        //Regardless of the children we are at a point where the current interval >= pattern and we have at least 1 character match from it
                        //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
                        //neighboring child (we have an overlap).
                        int intervalEndIdx = currentIntervalSize * (currentFrame.intervalIdx + 1) + (currentTreeIdx)*8 - 1;
                        int checkpos = bfProbe.consumed;
                        if(info.positionOffset == -1){
                            checkpos = intervalEndIdx;
                        }
                        info.positionOffset = intervalEndIdx - bfProbe.consumed + 1;
                        Pair<ArrayList<Integer>, Integer> res =  exhaustChild(currentTreeIdx, trees, currentIntervalSize, currentFrame.intervalIdx, info.positionOffset, pat, pi, checkpos);

                        for(int r : res.getFirst()){
                            results.add(r);
                        }
                        info.positionOffset = res.getSecond();
                        belongingTree = (int)(info.positionOffset  / tree.intervalSize);
                        iterationsUntilPossible = 0;

                    }
                }
            }
            currentTreeIdx++;
        }

        return results;
    }
}
