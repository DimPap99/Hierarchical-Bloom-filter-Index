package algorithms;

import PMIndex.HBI;
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

    public class MatchResult{
        public ArrayList<String> matches;
        public int matchedChars;

        MatchResult(){
            this.matches = new ArrayList<String>();
            this.matchedChars = 0;
        }
    }
    /// Perform bloom probes on blocks (intervals) of the tree. Separate matches based on potential consecutive matches.
    /// E.g. For support that for the current interval up to 8 chars fit in a block, maxIntervalChars = 8. If our key contains 10 chars and all of them
    /// match when probing the bloom filter, then we gonna have 8 consecutive characters as one result entry and
    ///  2 consecutive characters as a second result entry. If 9 characters match and then 10th doesnt. We will
    /// still have 2 intervals 8chars, 1 char. If 7 chars match and then 3 last dont, we have 1 interval with the 7 chars.
    /// etc.
    public MatchResult getPotentialConsecutiveMatches(TraversalInfo info, ImplicitTree tree){
        int maxIntervalChars = 1 << (tree.maxDepth - info.level);
        int run = 0;
        MatchResult res = new MatchResult();

        for (int i = 0; i < info.str.length(); i++) {
            String ch = info.str.substring(i, i + 1);
            String key = tree.createCompositeKey(info.level, info.intervalIdx, ch);

            if (run < maxIntervalChars && tree.membership.contains(key)) {
                ++run;
                ++res.matchedChars;
                if (run == maxIntervalChars) {                // full block
                    int start = i - run + 1;
                    res.matches.add(info.str.substring(start, i + 1));
                    run = 0;
                }
            } else {                                          // first mismatch
                if (run > 0) {
                    int start = i - run;
                    res.matches.add(info.str.substring(start, i));
                }
                return res;                                   // stop at first miss
            }
        }
        if (run > 0) {                                       // tail fragment
            int start = info.str.length() - run;
            res.matches.add(info.str.substring(start));
        }
        return res;
    }
    public int findInstance(int level, int intervalIdx, String pattern, ImplicitTree tree){
        int position = -1;
        Stack<TraversalInfo> stack = new Stack<>();
        stack.add(new TraversalInfo(level, intervalIdx, pattern));
        int nextLevel;
        int maxElements;
        while(true){
            TraversalInfo currentInfo = stack.pop();
            //if the matches we have for that interval are consecutive and equal to the length of the key
            //we check the left and right child
            MatchResult res = getPotentialConsecutiveMatches(currentInfo, tree);
            if(res.matchedChars == 0 && stack.size() == 0) break;
            else{
                int leftChild = tree.getLeftChild(currentInfo.intervalIdx);
                int rightChild = tree.getRightChild(currentInfo.intervalIdx);
                nextLevel = currentInfo.level + 1;
                maxElements = 1 << (tree.maxDepth - nextLevel);
                //check if the intervals we have are
                if(res.matches.size() == 1 || res.matches.size() == 2){

                }


            }
            //left


            //right

            break;
        }
        return position;
    }
    @Override
    public boolean exists(String key, ImplicitTree tree) {
        return findInstance(0, 0,key) == -1 ? false : true;
    }

    @Override
    public ArrayList<Integer> report(String key, ImplicitTree tree) {
        return null;
    }
}
