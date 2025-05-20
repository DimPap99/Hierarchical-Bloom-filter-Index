package algorithms;

import PMIndex.ImplicitTree;

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
        final int offset;         // how many chars of P already aligned to the left edge
        final int charsToMatch;
        Frame(int level, int interval, int offset, int charsToMatch) {
            this.level = level;
            this.intervalIdx = interval;
            this.offset  = offset;
            this.charsToMatch = charsToMatch;
        }
    }

    record Probe(int consumed, boolean complete) {} // complete==true if no miss in this node


    public class MatchResult{
        public ArrayList<String> matches;
        public int matchedChars;

        MatchResult(){
            this.matches = new ArrayList<String>();
            this.matchedChars = 0;
        }
    }

    Probe probe(ImplicitTree tree, int width, int level, int interval, String pattern, int charsToMatch, int offset) {



        for (int i = offset; i < charsToMatch; i++) {
            String key = tree.createCompositeKey(level, interval,
                    pattern.substring(i, i + 1));
            if (!tree.membership.contains(key)) {
                return new Probe(i, false);          // first mismatch at i
            }
        }
        return new Probe(pattern.length(), true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel){
        if((1 << maxLevel - level) * (intervalIdx + 1) > positionOffset) return true;
        else return false;
    }

    @Override
    public ArrayList<Integer> report(String key, ImplicitTree tree, boolean existence) {
        ArrayList<Integer> results = new ArrayList<>();
        int positionOffset = -1;
        //Shows how many chars we have not matched yet from our total pattern
        //because search starts from left, the remainder for the left child will be 0.
        int matched = 0;
        int rightChild;
        int leftChild;
        int childrenIntervalSize;
        int currentIntervalSize;
        Probe bfProbe;
        Stack<Frame> stack = new Stack<>();
        stack.add(new Frame(0, 0, 0, key.length()));
        Frame leftFrame = null;
        Frame rightFrame = null;
        int leftSideChars = 0;
        while(!stack.isEmpty()){
            Frame currentFrame = stack.pop();

            currentIntervalSize = tree.getIntervalSize(currentFrame.level);
            //if for any reason we got in the stack frames that are before the position offset drop them
            if(tree.getIntervalSize(currentFrame.level) * (currentFrame.intervalIdx + 1) > positionOffset && positionOffset != -1) continue;
            //probe the filter.
            bfProbe = probe(tree, currentIntervalSize,  currentFrame.level, currentFrame.intervalIdx, key, currentFrame.charsToMatch, currentFrame.offset);
            //No matches at all
            if(bfProbe.consumed == 0){
                positionOffset += tree.getIntervalSize(currentFrame.level);
                matched = 0;
                continue;
            }
            else{
                childrenIntervalSize = tree.getIntervalSize(currentFrame.level + 1);
                rightChild = tree.getRightChild(currentFrame.intervalIdx);
                leftChild = tree.getLeftChild(currentFrame.intervalIdx);
                //2 cases
                //Full match -> No remainder
                if(bfProbe.complete){
                    matched = 0;
                    leftFrame = new Frame(currentFrame.level + 1, leftChild, matched, bfProbe.consumed);
                    rightFrame = new Frame(currentFrame.level + 1, rightChild, matched, bfProbe.consumed);

                    //add right child
                    if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth)) {
                        stack.add(rightFrame);
                    }
                    //add left child
                    if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth)) {
                        stack.add(leftFrame);
                    }
                }else {
                    //We try to fit the consumed characters to the right most sides of the children intervals
                    leftSideChars = Math.max(0, currentIntervalSize - bfProbe.consumed);

                    if (leftSideChars > 0) {
                        matched += leftSideChars;
                        //create leftFrame
                        leftFrame = new Frame(currentFrame.level + 1, leftChild, matched, leftSideChars);
                    }

                    //create rightFrame
                    rightFrame = new Frame(currentFrame.level + 1, rightChild, matched, bfProbe.consumed - leftSideChars);
                    matched += bfProbe.consumed - leftSideChars;

                    //add right child
                    if (isValidChild(positionOffset, rightChild, currentFrame.level + 1, tree.maxDepth)) {
                        stack.add(rightFrame);
                    }
                    //add left child
                    if (isValidChild(positionOffset, leftChild, currentFrame.level + 1, tree.maxDepth) && leftSideChars > 0) {
                        stack.add(leftFrame);
                    }
                }
                }
            }





        return results;
    }
}
