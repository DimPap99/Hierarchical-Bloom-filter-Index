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
        if((1 << maxLevel - level) * (intervalIdx + 1) > positionOffset) return true;
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

        while(leafStartIdx <= intervalEnd){
            result = verifyAtLeaves(t, leafStartIdx, pat, pi);
            if(result.matched) matches.add(result.pos + 1 - pat.length);
            leafStartIdx = result.pos + 1;
        }
        return new Pair<ArrayList<Integer>, Integer>(matches, leafStartIdx);
    }

    /**
     * Starting at the given leafStartIdx, exhaustively scan for all
     * occurrences of `pat` whose first character lies within the
     * leaf‐level interval [intervalIdx*intervalSize .. (intervalIdx+1)*intervalSize−1].
     *
     * @param t             the ImplicitTree (the HBI block containing these leaves)
     * @param intervalSize  size of each leaf‐level interval (a power of two)
     * @param intervalIdx   index of the interval (0-based) that we are scanning
     * @param leafStartIdx  the absolute leaf index at which to begin scanning;
     *                      must satisfy: intervalStart ≤ leafStartIdx ≤ intervalEnd
     * @param pat           the pattern as a char[]
     * @param pi            the KMP prefix‐function array for pat
     * @return a list of MatchResult for each (probable) match whose first
     *         character falls between intervalStart and intervalEnd
     */
    public ArrayList<MatchResult> exhaustChild2(
            ImplicitTree t,
            int intervalSize,
            int intervalIdx,
            int leafStartIdx,
            char[] pat,
            int[] pi
    ) {
        ArrayList<MatchResult> matches = new ArrayList<>();

        // Compute absolute boundaries of this leaf interval:
        int intervalStart = intervalIdx * intervalSize;
        int intervalEnd   = intervalStart + intervalSize - 1;
        int m             = pat.length;

        // We only attempt a start if the full pattern can fit inside [intervalStart..intervalEnd].
        // So the largest valid start is intervalEnd - (m-1).
        int maxStart = intervalEnd - (m - 1);

        // If the caller gave us a start beyond maxStart, there's nothing to do.
        if (leafStartIdx > maxStart) {
            return matches;
        }

        // Iterate over all candidate leafStart positions (with KMP skipping).
        while (leafStartIdx <= maxStart) {
            // 1) Verify as far as possible at leaf level:
            MatchResult res = verifyAtLeaves(t, leafStartIdx, pat, pi);

            if (res.matched) {
                // Full pattern match. Record it.
                matches.add(res);

                // Advance by 1 to allow overlapping matches inside the same interval.
                // (If you prefer non-overlapping matches, use: leafStartIdx += m;)
                leafStartIdx = leafStartIdx + 1;
            } else {
                // Mismatch occurred. Compute how many chars matched:
                int j   = res.charsMatched; // 0..m-1
                int pos = res.pos;          // absolute position of the failed leaf

                // KMP‐style shift: if j==0, skip 1; else skip j - pi[j-1].
                int shift = (j == 0 ? 1 : j - pi[j - 1]);

                // The original start was leafStartIdx, so:
                int origStart = leafStartIdx;

                // The next “safe” start is:
                leafStartIdx = origStart + shift;
            }
        }

        return matches;
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
                j = pi[j - 1];  // fallback; note: no Bloom probe here
            }
        }
        // should not reach, but return mismatch safeguard
        return new MatchResult(false, leafStartIdx + i - 1, j);
    }

    public MatchResult checkImmediately(Frame frame, ImplicitTree tree, String pattern, int possibleMatches, int offset, int posAligned){

        int intSz = tree.getIntervalSize(frame.level);
        //starting index of current interval
//        int pos = intSz * frame.intervalIdx;
        int charsConsumed = 0;
        int passes = 0;
        //boolean canReset = false;
        int initPos;

        while (passes < possibleMatches) {
            passes++;
            posAligned++;
            String key = tree.createCompositeKey(tree.maxDepth, posAligned, pattern.substring(charsConsumed+offset, charsConsumed+1 + offset));
            if(!tree.membership.contains(key)){
                //if we only have 1 match from the prev interval and it fails its possible that the first char of the new interval
                //is correct but it is used at the second position in the pattern. E.g. ***b|bdca , Match started at last b
                //of 1st interval.
                if(offset == 1 && intSz - possibleMatches == offset){
                    possibleMatches = key.length();
                    //bring back 1 pos
                    posAligned--;

                }else{
                    possibleMatches = Math.max(0, possibleMatches-passes);
                }
                offset = 0;
                charsConsumed = 0;
                //reset the possible matches
                passes = 0;


            }else{
                charsConsumed++;
            }

            if(charsConsumed == pattern.length()) return new MatchResult(true, posAligned, charsConsumed);
            if(charsConsumed == possibleMatches) return new MatchResult(false, posAligned, charsConsumed);

        }
        posAligned--;

        return null;
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

            if(currentIntervalSize < key.length()){

            }else{
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

                    }

                }

            }


        }
        return results;
    }
}
