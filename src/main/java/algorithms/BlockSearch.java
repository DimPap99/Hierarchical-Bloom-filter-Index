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
        public int pos;
        public boolean reset = false;
        MatchResult(boolean matched, int pos, int charsMatched){
            this.matched = matched;
            this.pos = pos;
            this.charsMatched = charsMatched;
        }
    }



    Probe probe(ImplicitTree tree, int level, int interval, String pattern, int offset, int remainder) {


        int matches = 0;
        for (int i = offset; i < pattern.length(); i++) {
            String key = tree.createCompositeKey(level, interval,
                    pattern.substring(i, i + 1));
            if (!tree.membership.contains(key)) {
                return new Probe(i, false);          // first mismatch at i
            }
            matches++;
            if(matches == remainder){
                return new Probe(remainder, true);

            }
        }
        return new Probe(pattern.length(), true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel){
        if((1 << maxLevel - level) * (intervalIdx + 1) >= positionOffset) return true;
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

    public Pair<ArrayList<Integer>, Integer> exhaustChild(ImplicitTree t, int intervalSize, int intervalIdx, int leafStartIdx, char[] pat, int[] pi) {
        int intervalEnd = (intervalIdx + 1) * intervalSize - 1;
        ArrayList<Integer> matches = new ArrayList<>();
        MatchResult result; //= verifyAtLeaves(t, leafStartIdx, pat, pi);
        if(leafStartIdx == -1){ leafStartIdx = 0; }
        while(leafStartIdx <= intervalEnd){
            result = verifyAtLeavesNaive(t, leafStartIdx, pat);
            if(result.matched) matches.add(result.pos + 1 - pat.length);
            leafStartIdx = result.pos + 1;
        }
        return new Pair<ArrayList<Integer>, Integer>(matches, leafStartIdx);
    }

    public MatchResult verifyAtLeavesNaive(ImplicitTree t,
                                           int leafStartIdx,
                                           char[] pat) {
        int m = pat.length;
        for (int i = 0; i < m; i++) {
            // build the composite key for leaf index = leafStartIdx + i,
            // and the single pattern character pat[i]
            String key = t.createCompositeKey(t.maxDepth,
                    leafStartIdx + i,
                    String.valueOf(pat[i]));
            if (!t.membership.contains(key)) {
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


    public MatchResult verifyAtLeaves(ImplicitTree t,
                                              int leafStartIdx,
                                              char[] pat,
                                              int[] pi) {
        int i = 0;               // index in text (leaf offset)
        int j = 0;               // index in pattern
        int m = pat.length;
        while (i < m) {
            String key = t.createCompositeKey(t.maxDepth, leafStartIdx + i, String.valueOf(pat[j]));
            if (t.membership.contains(key)) {
                ++i; ++j;
                if (j == m) {   // full match
                    int lastPos = leafStartIdx + m - 1;
                    return new MatchResult(true, lastPos, m);
                }
            } else {
                if (j == 0) {   // miss on first char
                    return new MatchResult(false, leafStartIdx + i, 0);
                }
                j = pi[j - 1];  // fallback;
            }
        }
        // should not reach, but return mismatch safeguard
        return new MatchResult(false, leafStartIdx + i - 1, j);
    }

    void generateChildren(Frame currentFrame, Stack<Frame> framesStack, int positionOffset, ImplicitTree tree){
        int rightChild = tree.getRightChild(currentFrame.intervalIdx);
        int leftChild = tree.getLeftChild(currentFrame.intervalIdx);

        Frame leftFrame = new Frame(currentFrame.level + 1, leftChild);
        Frame rightFrame = new Frame(currentFrame.level + 1, rightChild);

        //add right child
        if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth)) {
            framesStack.add(rightFrame);
        }
        //add left child
        if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth)) {
            framesStack.add(leftFrame);
        }

    }

    ArrayList<Frame> getChildren(Frame currentFrame, int positionOffset, ImplicitTree tree){
        ArrayList<Frame> frames = new ArrayList<>();
        int rightChild = tree.getRightChild(currentFrame.intervalIdx);
        int leftChild = tree.getLeftChild(currentFrame.intervalIdx);

        Frame leftFrame = new Frame(currentFrame.level + 1, leftChild );
        Frame rightFrame = new Frame(currentFrame.level + 1, rightChild);

        //add right child
        if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth)) {
            frames.add(rightFrame);
        }else frames.add(null);
        //add left child
        if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth)) {
            frames.add(leftFrame);
        }else frames.add(null);
        return frames;
    }

    @Override
    public ArrayList<Integer> report(String key, ImplicitTree tree, boolean existence) {
        ArrayList<Integer> results = new ArrayList<>();
        char[] pat = key.toCharArray();
        int[] pi   = prefixFunction(pat);
        GlobalInfo info = new GlobalInfo(key.length());

        int childrenIntervalSize;
        int currentIntervalSize;
//        childrenIntervalSize = tree.getIntervalSize(currentFrame.level + 1);
//        int intervalFinalIdx = currentIntervalSize * (currentFrame.intervalIdx + 1) - 1;
//        int possibleMatches = intervalFinalIdx - positionOffset;
        Probe bfProbe;
        Stack<Frame> stack = new Stack<>();
        //check root
        Frame rootFrame = new Frame(0, 0);
        bfProbe = probe(tree,  rootFrame.level, rootFrame.intervalIdx, key, 0, key.length() - info.matched);
        if(bfProbe.complete){
            generateChildren(rootFrame, stack, info.positionOffset, tree);
        }
        while(!stack.isEmpty()) {
            Frame currentFrame = stack.pop();

            currentIntervalSize = tree.getIntervalSize(currentFrame.level);
            childrenIntervalSize = tree.getIntervalSize(currentFrame.level + 1);


            bfProbe = probe(tree,  currentFrame.level, currentFrame.intervalIdx, key, info.matched, key.length() - info.matched);
            System.out.println("Checking " + Utils.intervalWithHashes(4, currentFrame.level, currentFrame.intervalIdx*currentIntervalSize)
                    + " Matched: " + info.matched + " Probe: " + bfProbe.consumed);

            if(bfProbe.consumed == 0){
                info.positionOffset = currentIntervalSize * (currentFrame.intervalIdx + 1) - 1;
                info.matched = 0;
            }else{
                //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
                if(bfProbe.complete & childrenIntervalSize >= key.length()) {
                    //we just generate children as theres a chance that the pattern is in both of them
                    generateChildren(currentFrame, stack, info.positionOffset, tree);
                }
                else if(bfProbe.complete & childrenIntervalSize < key.length() ){
                    //In this case we know that the probe matched all characters in the current interval.
                    //But also that the children intervals do not completely fit the pattern. If the pattern exists
                    //It will overlap between Left Interval and Right.
                    Pair<ArrayList<Integer>, Integer> res =  exhaustChild(tree, currentIntervalSize, currentFrame.intervalIdx, info.positionOffset, pat, pi);

                    for(int r : res.getFirst()){
                        results.add(r);
                    }
                    info.positionOffset = res.getSecond();
                }
                else{ //the probes were inclomplete. Meaning:
                    //Regardless of the children we are at a point where the current interval >= pattern and we have at least 1 character match from it
                    //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
                    //neighboring child (we have an overlap).
                    int intervalEndIdx = currentIntervalSize * (currentFrame.intervalIdx + 1) - 1;
                    info.positionOffset = intervalEndIdx - bfProbe.consumed + 1;
                    Pair<ArrayList<Integer>, Integer> res =  exhaustChild(tree, currentIntervalSize, currentFrame.intervalIdx, info.positionOffset, pat, pi);

                    for(int r : res.getFirst()){
                        results.add(r);
                    }
                    info.positionOffset = res.getSecond();
                }

            }




        }
        return results;
    }
}
