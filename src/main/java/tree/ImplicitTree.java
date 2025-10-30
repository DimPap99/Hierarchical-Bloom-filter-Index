package tree;

import estimators.Estimator;
import membership.*;
import search.BlockSearch;
import search.Frame;
import search.PruningPlan;

import java.util.ArrayList;
import java.util.Deque;
import java.util.function.IntFunction;



/**
 * An implicit binary tree of Membership data structures filters that indexes
 * a contiguous fragment of the text stream.
 *
 * @param <M> concrete {@link Membership} implementation
 */
public final class ImplicitTree< M extends Membership> {

    public int containCounter = 0;
    private final TreeLayout           layout;
    public int indexedItemsCounter=-1;
    private final int[] visitedPerLevel = new int[64];

    private final LevelDirectory< M> levels;
    public final StreamBuffer         buffer;
    public final KeyPackingService codec;
    public boolean alwaysFits = false;
    public Estimator estimator;
    int treeSize;
    // Maximum number of characters the root interval spans.
    private final int capacity;
    public PruningPlan pruningPlan;
//    public int treeIdx;
    // Global position of the last appended character.  Starts at –1.
    private long endPos = -1;
    public int id = 0;

    public ImplicitTree(TreeLayout layout,
                        IntFunction<M> filterFactory,
                        KeyPackingService codec, Estimator est) {
        this.layout   = layout;
        this.levels   = new LevelDirectory<>(layout, filterFactory);
        this.capacity = layout.intervalSize(0);
        this.buffer   = new StreamBuffer(true, this.capacity);
        this.codec    = codec;
        this.estimator = est;
        containCounter = 0;
//        this.treeSize = size;
//        this.tr
    }


    /**
     * Insert a character that occurs at {@code globalPos} in the full stream.
     * The character is packed into a key for each level and inserted into
     * the corresponding membership structure.
     */
    public void append(long symbol, long globalPos) {
        buffer.append(symbol);
        endPos = globalPos;
        indexedItemsCounter++;


//        int span       = layout.intervalSize(layout.getEffectiveRootLevel());
        for (int level = layout.getEffectiveRootLevel(); level < layout.getEffectiveLeafLevel(); level++) {
            int span       = layout.intervalSize(level);
            int intervalId = indexedItemsCounter / span;      // node id
//            if (codec.fitsOneWord(intervalId, symbol)) {
//                long w = codec.packWord(intervalId, symbol);
//                levels.insert(level, w);           // calls BloomFilter.insert(long)
//            } else {
            long hi = Integer.toUnsignedLong(intervalId);
            long lo = symbol;
            levels.insert(level, hi, lo);      // calls BloomFilter.insert(long,long)
            //}

        }
    }

//    public void append(long c, long globalPos) {
//        buffer.append(c);
//        endPos = globalPos;
//        indexedItemsCounter++;
//
//        for (int level = layout.getEffectiveRootLevel(); level < layout.getEffectiveLeafLevel(); level++) {
//            int span       = layout.intervalSize(level);
//            int intervalId = (int) (indexedItemsCounter / span);      // LOCAL id
//            long key       = codec.pack(level, intervalId, c);
//            levels.insert(level, key);
//        }
//    }

    /** Returns {@code true} once the tree’s root interval is full. */
    public boolean isFull() {
        return buffer.length() >= capacity;
    }

    public int baseIntervalSize(){ return this.capacity; }
    /** Inclusive global position of the last character held. */
    public long endPos() { return endPos; }

    //  Geometry helpers forwarded from TreeLayout
    public int intervalSize(int level)      { return layout.intervalSize(level); }
    public int leftChild(int idx) { return layout.leftChild(idx); }
    public int rightChild(int idx){ return layout.rightChild(idx); }
    public int maxDepth(){return layout.levels();}
    public int effectiveDepth(){
        return layout.getEffectiveLeafLevel();
    }
    //  Expose per-level membership and level-count (read-only)
    public M filterAtLevel(int level) { return levels.filter(level); }

    public int totalFilters(){return  levels.depth();}
    public void dropFiltersUpToLp(int lp){

        int effectiveRoot = this.effectiveRoot();
        while(lp > 0){
            this.levels.dropLevel(effectiveRoot);
            effectiveRoot++;
            lp--;
        }
        this.layout.setEffectiveRootLevel(effectiveRoot);

    }
    public int levelCount() { return levels.depth(); }

    public int effectiveRoot(){
        return layout.getEffectiveRootLevel();
    }
    public void generateChildren(Frame currentFrame,
                                 Deque<Frame> framesStack,
                                 int positionOffset,
                                 int workingTreeIdx) {

        int nextLevel = currentFrame.level() + 1;
        int maxDepth  = this.maxDepth();

        int rightChildIdx = this.rightChild(currentFrame.intervalIdx());
        int leftChildIdx  = this.leftChild(currentFrame.intervalIdx());

        // Push right, then left, but only if that subtree can still
        // contain a position >= positionOffset.
//        if (this.isValidChild(positionOffset, rightChildIdx, nextLevel, maxDepth, workingTreeIdx)) {
            framesStack.push(new Frame(nextLevel, rightChildIdx));
//        }

//        if (this.isValidChild(positionOffset, leftChildIdx, nextLevel, maxDepth, workingTreeIdx)) {
            framesStack.push(new Frame(nextLevel, leftChildIdx));
//        }
    }



    public ArrayList<Frame> generateChildren(Frame currentFrame, int positionOffset, int workingTreeIdx){
        int rightChild = this.rightChild(currentFrame.intervalIdx());
        int leftChild = this.leftChild(currentFrame.intervalIdx());
        ArrayList<Frame> children = new  ArrayList<Frame>();
        int nextLevel = currentFrame.level() + 1;
        Frame leftFrame = new Frame(nextLevel, leftChild);
        Frame rightFrame = new Frame(nextLevel, rightChild);
        int maxDepth = this.maxDepth();
        //add right child
//        if (this.isValidChild(positionOffset, rightChild, nextLevel, maxDepth, workingTreeIdx)) {
            children.add(rightFrame);
        //}
        //add left child
        //if (isValidChild(positionOffset, leftChild, nextLevel, maxDepth, workingTreeIdx)) {
            children.add(leftFrame);
        //}
        return children;
    }

//    public int traverse(Frame f,
//                 int maxLevel,
//                 int posOffset,
//                 long firstSymbol) {
//
//        // Leaf: check the actual buffer position
//        final int packedSymbol = Math.toIntExact(firstSymbol);
//
//        if (f.level() == maxLevel) {
//            return (f.intervalIdx() >= posOffset &&
//                    this.buffer.get(f.intervalIdx()) == firstSymbol)
//                    ? f.intervalIdx()
//                    : -1;
//        }
//
//        // Internal node: depth-first search of children
//        for (Frame child : this.generateChildren(f, posOffset, this.id)) {
//
//            long key = this.codec.pack(child.level(), child.intervalIdx(), packedSymbol);
//
//            if (!this.contains(child.level(), key)) {
//                // This child cannot possibly contain firstChar – skip it
//                continue;
//            }
//
//            int res = traverse(child, maxLevel, posOffset, firstChar);
//            if (res != -1) {          // found it in this subtree
//                return res;
//            }
//        }
//
//        // Not found in any child
//        return -1;
//    }

    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;              // == tree.intervalSize

        // end position of current interval inside the global stream
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll
                - 1;

        return maxPos >= positionOffset;
    }


    //  Membership lookup
    public boolean contains(int level, long key) {

        containCounter++;
        if (level >= 0 && level < visitedPerLevel.length) {
            visitedPerLevel[level]++;
        }
        return levels.contains(level, key);
    }
    public int[] snapshotVisitedPerLevel() {
        return java.util.Arrays.copyOf(visitedPerLevel, visitedPerLevel.length);
    }
    /** Clear both the per-level visits and the aggregate Bloom counter */
    public void clearVisitCounters() {
        java.util.Arrays.fill(visitedPerLevel, 0);
        containCounter = 0;
    }

    public boolean contains(int level, long hi, long lo) {
        containCounter++;
        if (level >= 0 && level < visitedPerLevel.length) {
            visitedPerLevel[level]++;
        }
        return levels.contains(level, hi, lo);
    }

    public double getMembershipFpRate(int level) {
        if(level <= maxDepth()) return this.levels.getLevelFpRate(level);
        else return -1;
    }
}
