package algorithms;

import PMIndex.ImplicitTree;
import org.apache.commons.math3.util.Pair;
import utilities.Utils;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
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
        for (int i = offset; i < pattern.length; i++) {
            String key = tree.createCompositeKey(level, interval,
                    pattern[i]);
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
        int maxPos = (1 << maxLevel - level) * (intervalIdx + 1) + (workingTreeIdx)*8 - 1;
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

    public Pair<ArrayList<Integer>, Integer> exhaustChild(int currentTreeIdx, ArrayList<ImplicitTree> trees, int intervalSize, int intervalIdx, int leafStartIdx, char[] pat, int[] pi, int positionChecks) {
        int workingTreeIdx = currentTreeIdx;
        ImplicitTree workingTree = trees.get(currentTreeIdx);
        int stopPos = leafStartIdx + positionChecks;
        ArrayList<Integer> matches = new ArrayList<>();
        MatchResult result; //= verifyAtLeaves(t, leafStartIdx, pat, pi);
        if(leafStartIdx == -1){ leafStartIdx = 0; }

        while(leafStartIdx <= stopPos){
            int idxx = (int)(leafStartIdx / (workingTree.maxIdx + 1));
            if((int)(leafStartIdx / (workingTree.maxIdx + 1)) != workingTreeIdx) {
                workingTreeIdx++;
                if (workingTreeIdx < trees.size()) {
                    workingTree = trees.get(workingTreeIdx);
                } else {
                    //No more stream characters to check || We are the end of the stream
                    return new Pair<ArrayList<Integer>, Integer>(matches, leafStartIdx);

                }
            }
            result = verifyAtLeavesNaive(currentTreeIdx, trees, leafStartIdx, pat);
            if(result.matched){
            matches.add(result.pos + 1 - pat.length);
            leafStartIdx = result.pos + 1;
            }else{
                leafStartIdx++;
            }
        }
        return new Pair<ArrayList<Integer>, Integer>(matches, leafStartIdx);
    }

    public MatchResult verifyAtLeavesNaive(int currentTreeIdx, ArrayList<ImplicitTree> trees, int leafStartIdx,
                                           char[] pat) {
        int m = pat.length;
        int workingTreeIdx = currentTreeIdx;
        int leafIdx;
        ImplicitTree workingTree = trees.get(currentTreeIdx);;
        for (int i = 0; i < m; i++) {
            // build the composite key for leaf index = leafStartIdx + i,
            // and the single pattern character pat[i]
            leafIdx = leafStartIdx + i;
            if((int)(leafIdx / (workingTree.maxIdx + 1)) != workingTreeIdx){
                workingTreeIdx++;
                if(workingTreeIdx < trees.size()){
                    workingTree = trees.get(workingTreeIdx);
                }else{
                    //No more stream characters to check || We are the end of the stream
                    return new MatchResult(false,
                            leafStartIdx + i - 1,
                            i);
                }

            }
            String key = workingTree.createCompositeKey(workingTree.maxDepth,
                    leafIdx,
                    pat[i]);
            if (!workingTree.membership.contains(key)) {
                // mismatch at the i‐th character of pat:
                //   - pos = absolute leaf where it failed = leafStartIdx + i
                //   - charsMatched = i (how many matched before failing)
                if(i == 0){
                    return new MatchResult(false,
                            leafStartIdx,
                            i);
                }
                else{
                    return new MatchResult(false,
                            leafStartIdx + i - 1,
                            i);
                }

            }
        }
        // if we reach here, all m characters matched in a row
        int lastPos = leafStartIdx + (m - 1);
        return new MatchResult(true,
                lastPos,
                m);
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
    void generateChildren(Frame currentFrame, Stack<Frame> framesStack, int positionOffset, ImplicitTree tree, int workingTreeIdx){
        int rightChild = tree.getRightChild(currentFrame.intervalIdx);
        int leftChild = tree.getLeftChild(currentFrame.intervalIdx);

        Frame leftFrame = new Frame(currentFrame.level + 1, leftChild);
        Frame rightFrame = new Frame(currentFrame.level + 1, rightChild);

        //add right child
        if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth, workingTreeIdx)) {
            framesStack.add(rightFrame);
        }
        //add left child
        if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth, workingTreeIdx)) {
            framesStack.add(leftFrame);
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

            int childrenIntervalSize;
            int currentIntervalSize;
            int belongingTree = (int)(info.positionOffset  / tree.intervalSize);
            Stack<Frame> stack = new Stack<>();
            //check root
            Frame rootFrame = new Frame(0, 0);
            bfProbe = probe(tree,  rootFrame.level, rootFrame.intervalIdx, key, 0, key.length - info.matched);
            if(bfProbe.complete){
                generateChildren(rootFrame, stack, info.positionOffset, tree, currentTreeIdx);
            }
            while(!stack.isEmpty() && belongingTree == currentTreeIdx) {
                Frame currentFrame = stack.pop();

                currentIntervalSize = tree.getIntervalSize(currentFrame.level);
                childrenIntervalSize = tree.getIntervalSize(currentFrame.level + 1);
                bfProbe = probe(tree,  currentFrame.level, currentFrame.intervalIdx, key, info.matched, key.length - info.matched);
//                System.out.println("Checking " + Utils.intervalWithHashes(4, currentFrame.level, currentFrame.intervalIdx*currentIntervalSize)
//                        + " Matched: " + info.matched + " Probe: " + bfProbe.consumed);
                //skip entiner interval - the pattern is not inside it
                if(bfProbe.consumed == 0){
                    info.positionOffset = currentIntervalSize * (currentFrame.intervalIdx + 1) + (currentTreeIdx)*8;
                    info.matched = 0;
                }else{
                    //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
                    if(bfProbe.complete & childrenIntervalSize >= key.length) {
                        //we just generate children as theres a chance that the pattern is in both of them
                        generateChildren(currentFrame, stack, info.positionOffset, tree, currentTreeIdx);
                    }
                    else if(bfProbe.complete & childrenIntervalSize < key.length ){
                        //In this case we know that the probe matched all characters in the current interval.
                        //But also that the children intervals do not completely fit the pattern. If the pattern exists
                        //It will overlap between Left Interval and Right.

                        Pair<ArrayList<Integer>, Integer> res =  exhaustChild(currentTreeIdx, trees, currentIntervalSize, currentFrame.intervalIdx, info.positionOffset, pat, pi, currentIntervalSize);

                        for(int r : res.getFirst()){
                            results.add(r);
                        }
                        info.positionOffset = res.getSecond();
                        belongingTree = (int)(info.positionOffset  / tree.intervalSize);
                    }
                    else{ //the probes were inclomplete. Meaning:
                        //Regardless of the children we are at a point where the current interval >= pattern and we have at least 1 character match from it
                        //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
                        //neighboring child (we have an overlap).
                        int intervalEndIdx = currentIntervalSize * (currentFrame.intervalIdx + 1) + (currentTreeIdx)*8 - 1;
                        int checkpos = intervalEndIdx - info.positionOffset;
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
                    }
                }
            }
            currentTreeIdx++;
        }

        return results;
    }
}
