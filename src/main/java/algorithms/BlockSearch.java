package algorithms;

import PMIndex.ImplicitTree;
import utilities.Utils;

import java.util.ArrayList;
import java.util.Stack;
public class BlockSearch implements SearchAlgorithm{

    public class TraversalInfo{
        public int level;
        public int intervalIdx;
        public String str;

        public TraversalInfo(int level, int intervalIdx, String str){
            this.level = level;
            this.intervalIdx = intervalIdx;
            this.str = str;
        }


    }
    static final class Frame {
        final int level;          // depth
        final int intervalIdx;       // index within that depth
        final int charsToMatch;
        Frame(int level, int interval, int charsToMatch ) {
            this.level = level;
            this.intervalIdx = interval;
            this.charsToMatch = charsToMatch;
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
                int b = 2;
                return new Probe(remainder, true);

            }
        }
        return new Probe(pattern.length(), true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel){
        if((1 << maxLevel - level) * (intervalIdx + 1) > positionOffset) return true;
        else return false;
    }

    public MatchResult checkImmediately(Frame frame, ImplicitTree tree, String pattern, int possibleMatches, int offset, int posAligned){

        int intSz = tree.getIntervalSize(frame.level);
        //starting index of current interval
//        int pos = intSz * frame.intervalIdx;
        int charsConsumed = 0;
        int passes = 0;
        //boolean canReset = false;


        while (passes < possibleMatches) {
            passes++;
            posAligned++;
            String key = tree.createCompositeKey(tree.maxDepth, posAligned, pattern.substring(charsConsumed+offset, charsConsumed+1 + offset));
            if(!tree.membership.contains(key)){
                charsConsumed = 0;
                //reset the possible matches
                possibleMatches = Math.max(0, possibleMatches-passes);
                passes = 0;
            }else{
                charsConsumed++;
            }

            if(charsConsumed == pattern.length()) return new MatchResult(true, posAligned, charsConsumed);
            if(charsConsumed == possibleMatches) return new MatchResult(false, posAligned, charsConsumed);

        }
        posAligned--;
//        for(int i=offset; i< intSz; i++){
//
//            String key = tree.createCompositeKey(tree.maxDepth, pos, pattern.substring(charsConsumed+offset, charsConsumed+1 + offset));
//            //reset
//            if(!tree.membership.contains(key)){
//                charsConsumed = 0;
//                //reset the possible matches
//
//            }else{
//                charsConsumed++;
//            }
//
//            if(charsConsumed == pattern.length()) return new MatchResult(true, pos, charsConsumed);
//            if(charsConsumed == possibleMatches) return new MatchResult(true, pos, charsConsumed);
//            pos++;
//            passes++;
//        }
//        pos--;
        if(charsConsumed != 0) return new MatchResult(false, posAligned, charsConsumed);


        return null;
    }

    void generateChildren(Frame currentFrame, Stack<Frame> framesStack, int remainder, int positionOffset, ImplicitTree tree, String key, boolean fits){
        int rightChild = tree.getRightChild(currentFrame.intervalIdx);
        int leftChild = tree.getLeftChild(currentFrame.intervalIdx);

        Frame leftFrame = new Frame(currentFrame.level + 1, leftChild, key.length() );
        Frame rightFrame = new Frame(currentFrame.level + 1, rightChild, key.length());
        if(fits){
            //add right child
            if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth)) {
                framesStack.add(rightFrame);
            }
            //add left child
            if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth)) {
                framesStack.add(leftFrame);
            }
        }
    }
    @Override
    public ArrayList<Integer> report(String key, ImplicitTree tree, boolean existence) {
        ArrayList<Integer> results = new ArrayList<>();
        int positionOffset = -1;
        //Shows how many chars we have not matched yet from our total pattern
        //because search starts from left, the remainder for the left child will be 0.
        int matched = 0;
        int remainder = key.length();
        int rightChild;
        int leftChild;
        int childrenIntervalSize;
        int currentIntervalSize;
        Probe bfProbe;
        Stack<Frame> stack = new Stack<>();
//        Frame leftFrame = null;
//        Frame rightFrame = null;
        int leftSideChars = 0;
        int rightSideChars = 0;
        //check root
        Frame rootFrame = new Frame(0, 0, key.length());
        bfProbe = probe(tree,  rootFrame.level, rootFrame.intervalIdx, key, matched, remainder);
        if(bfProbe.complete){
            generateChildren(rootFrame, stack, remainder, positionOffset, tree, key, true);
        }
        while(!stack.isEmpty()){
            Frame currentFrame = stack.pop();

            currentIntervalSize = tree.getIntervalSize(currentFrame.level);
            //if for any reason we got in the stack frames that are before the position offset drop them
           // if(tree.getIntervalSize(currentFrame.level) * (currentFrame.intervalIdx + 1) > positionOffset && positionOffset != -1) continue;
            //probe the filter.
            bfProbe = probe(tree,  currentFrame.level, currentFrame.intervalIdx, key, matched, remainder);
            System.out.println("Checking " + Utils.intervalWithHashes(4, currentFrame.level, currentFrame.intervalIdx*currentIntervalSize) + "\nRemainder: " + remainder
                    + " Matched: " + matched + " Probe: " + bfProbe.consumed);
            //No matches at all
            if(bfProbe.consumed == 0){
                positionOffset += currentIntervalSize;
                matched = 0;
                remainder = key.length();
                continue;
            }
            else{
                childrenIntervalSize = tree.getIntervalSize(currentFrame.level + 1);
                //2 cases
                //Full match
                if(bfProbe.complete){

                    if(childrenIntervalSize >= key.length()){
                        generateChildren(currentFrame, stack, remainder, positionOffset, tree, key, true);
                    }
                    else{
                        //look immediately for the pattern. Find how the characters are distributed along the children
                        //and if we have any remaining chars update the remainder
                        System.out.println("Immediate Check on: " + Utils.intervalWithHashes(4, currentFrame.level, currentFrame.intervalIdx*currentIntervalSize));
                        MatchResult result = checkImmediately(currentFrame, tree, key, remainder, matched, positionOffset);
                        if(result.matched || remainder - result.charsMatched == 0){
                            results.add(result.pos-key.length());
                            positionOffset = result.pos;
                            matched = 0;
                            remainder = key.length();
                            //At this point check immediately if there can be the start of another pattern from this interval
                            if(positionOffset < (currentIntervalSize * currentFrame.level) - 1){
                                int possibleMatches = (currentIntervalSize * currentFrame.level) - 1 - positionOffset;
                                result = checkImmediately(currentFrame, tree, key, possibleMatches, matched, positionOffset );
                                positionOffset = result.pos;
                                remainder -= result.charsMatched;
                                matched = result.charsMatched;
                            }
                        }else{
                            positionOffset = result.pos;

                            if(result.charsMatched == 0) {

                                //reset
                                remainder = key.length();
                                matched = 0;
                                if(positionOffset < (currentIntervalSize * currentFrame.level) - 1){
                                    int possibleMatches = (currentIntervalSize * currentFrame.level) - 1 - positionOffset;
                                    result = checkImmediately(currentFrame, tree, key, possibleMatches, matched, positionOffset );
                                    positionOffset = result.pos;
                                    remainder -= result.charsMatched;
                                    matched = result.charsMatched;
                                }

                            }
                            else {
                                remainder -= result.charsMatched;
                                matched = result.charsMatched;
                            }
                            continue;
                        }
                    }
                }
                else {
                    //MatchResult result = checkImmediately(currentFrame, tree, key, remainder, matched);
                    int b = 2;
                    }
                }
            }





        return results;
    }
}
