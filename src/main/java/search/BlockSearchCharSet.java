package search;

import estimators.Estimator;
import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;

public class BlockSearchCharSet implements SearchAlgorithm {
    public Estimator estimator;
    public int currentOffset;


    @Override
    public CandidateRange search(Frame f, Pattern p, ImplicitTree tree, Deque<Frame> stack, int positionOffset) {
        Probe probe = probe(tree,
                f.level(),             // unchanged helper signature
                f.intervalIdx(),
                p.nGramToInt, p.charStartLp);
        this.currentOffset = positionOffset;
        int currentIntervalSize = tree.intervalSize(f.level());
        int childrenIntervalSize = currentIntervalSize / 2;
        int treeBaseInterval = tree.baseIntervalSize();

        if (probe.consumed() == 0) {
            this.currentOffset = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval;
            return null;
        } else {
            //the intervals of the children are bigger or equal to the pattern and the probe matched all characters
            if (probe.complete() && childrenIntervalSize >= p.text.length) {
                //we just generate children as theres a chance that the pattern is in both of them
                tree.generateChildren(f, stack, positionOffset, tree.id);
                return null;
            } else {
                int intervalEndIdx = currentIntervalSize * (f.intervalIdx() + 1) + tree.id * treeBaseInterval - 1;

                if(probe.complete()) return new CandidateRange(this.currentOffset, intervalEndIdx);
                else {
                    //the probes were inclomplete. Meaning:
                    //Regardless of the children we are at a point where the current interval >= pattern and we have at least 1 character match from it
                    //For a pattern to truly exist it must reside in the right most positions of the current interval. The remaining characters will be at the
                    //neighboring child (we have an overlap).
                    this.currentOffset = intervalEndIdx - probe.consumed() + 1;
                    return new CandidateRange(this.currentOffset, intervalEndIdx);
                }
            }
        }
    }


    public int getCurrentOffset(){return this.currentOffset;}

    //Counts the consecutive matches from the start of the pattern
    int countStartConsecutiveMatches(boolean[] matchedArr){
        int matches = 0;
        for(int i = 0; i < matchedArr.length; i++){
            if(!matchedArr[i]){
               return matches;
            }else{
                matches++;
            }
        }
        return matches;
    }
    Probe probe(ImplicitTree tree, int level, int interval, int[] pattern, ArrayList<Integer> lp) {

        long key;

        boolean[] matchedArr = new boolean[pattern.length];
        Arrays.fill(matchedArr, Boolean.FALSE);

        for (int i = 0; i < pattern.length; i++) {

            if ( level >= lp.get(i)) {
                key =  tree.codec.pack(level, interval, pattern[i]);

                if (!tree.contains(level, key)) {
                    //in case of a mismatch we need to know how many consecutive occurences we have from the start of the pattern
                    //for anything else we dont care. Because there is a chance we havent checked characters from the beggining
                    //of the pattern at all, we perform probes until the character we are currently at.
                    for(int j =0; j < i; j++){
                        key =  tree.codec.pack(level, interval, pattern[j]);
                        matchedArr[j] = tree.contains(level, key);
                        //on first mismatch from the starting char we break
                        if(!matchedArr[j])break;
                    }
                    return new Probe(countStartConsecutiveMatches(matchedArr), false);          // first mismatch at i
                }
                matchedArr[i] =  true;


            }

        }
        return new Probe(countStartConsecutiveMatches(matchedArr), true);                 // checked len chars, all good
    }


    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;              // == tree.intervalSize

        /* end position of current interval inside the *global* stream */
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll
                - 1;

        return maxPos >= positionOffset;
    }

}
