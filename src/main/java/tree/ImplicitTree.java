package tree;

import estimators.Estimator;
import membership.Key;
import membership.Key64;
import membership.LongKey;
import membership.Membership;
import search.BlockSearch;
import search.Frame;
import search.PruningPlan;

import java.util.ArrayList;
import java.util.Deque;
import java.util.function.IntFunction;



/**
 * An implicit binary tree of Bloom/Cuckoo filters that indexes
 * a contiguous fragment of the text stream.
 *
 * @param <M> concrete {@link Membership} implementation
 */
public final class ImplicitTree< M extends Membership> {

    public int containCounter = 0;
    private final TreeLayout           layout;
    public int indexedItemsCounter=-1;

    private final LevelDirectory< M> levels;
    public final StreamBuffer         buffer;
    public final LongKey codec;
    public Estimator estimator;
    /** Maximum number of characters the root interval spans. */
    private final int capacity;
    public PruningPlan pruningPlan;
//    public int treeIdx;
    /** Global position of the last appended character.  Starts at –1. */
    private long endPos = -1;
    public int id = 0;
    public ImplicitTree(TreeLayout layout,
                        IntFunction<M> filterFactory,
                        LongKey codec, Estimator est) {
        this.layout   = layout;
        this.levels   = new LevelDirectory<>(layout.levels(), filterFactory);
        this.buffer   = new StreamBuffer();
        this.codec    = codec;
        this.capacity = layout.intervalSize(0);
        this.estimator = est;
        containCounter = 0;
//        this.tr
    }

    /**
     * Insert a character that occurs at {@code globalPos} in the full stream.
     * The character is packed into a key for each level and inserted into
     * the corresponding membership structure.
     */
    public void append(int c, long globalPos) {
        buffer.append(c);
        endPos = globalPos;
        indexedItemsCounter++;

        for (int level = 0; level < layout.levels(); level++) {
            int span       = layout.intervalSize(level);
            int intervalId = (int) (indexedItemsCounter / span);      // LOCAL id
            long key       = codec.pack(level, intervalId, c);
            levels.insert(level, key);
        }
    }

    /** Returns {@code true} once the tree’s root interval is full. */
    public boolean isFull() {
        return buffer.length() >= capacity;
    }

    public int baseIntervalSize(){ return this.capacity; }
    /** Inclusive global position of the last character held. */
    public long endPos() { return endPos; }

    // ---- Geometry helpers forwarded from TreeLayout ----
    public int intervalSize(int level)      { return layout.intervalSize(level); }
    public int leftChild(int idx) { return layout.leftChild(idx); }
    public int rightChild(int idx){ return layout.rightChild(idx); }
    public int maxDepth(){return layout.levels();}
    public void generateChildren(Frame currentFrame, Deque<Frame> framesStack, int positionOffset, int workingTreeIdx){
        int rightChild = this.rightChild(currentFrame.intervalIdx());
        int leftChild = this.leftChild(currentFrame.intervalIdx());

        int nextLevel = currentFrame.level() + 1;
        Frame leftFrame = new Frame(nextLevel, leftChild);
        Frame rightFrame = new Frame(nextLevel, rightChild);
        int maxDepth = this.maxDepth();
        //add right child
        //if (this.isValidChild(positionOffset, rightChild, nextLevel, maxDepth, workingTreeIdx)) {
            framesStack.push(rightFrame);
        //}
        //add left child
        //if (isValidChild(positionOffset, leftChild, nextLevel, maxDepth, workingTreeIdx)) {
            framesStack.push(leftFrame);
        //}
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
        //if (this.isValidChild(positionOffset, rightChild, nextLevel, maxDepth, workingTreeIdx)) {
            children.add(rightFrame);
        //}
        //add left child
        //if (isValidChild(positionOffset, leftChild, nextLevel, maxDepth, workingTreeIdx)) {
            children.add(leftFrame);
        //}
        return children;
    }

    public int traverse(Frame f,
                 int maxLevel,
                 int posOffset,
                 int firstChar) {

        // Leaf: check the actual buffer position
        if (f.level() == maxLevel) {
            return (f.intervalIdx() >= posOffset &&
                    this.buffer.data.get(f.intervalIdx()) == firstChar)
                    ? f.intervalIdx()
                    : -1;
        }

        // Internal node: depth-first search of children
        for (Frame child : this.generateChildren(f, posOffset, this.id)) {

            long key = this.codec.pack(child.level(), child.intervalIdx(), firstChar);

            if (!this.contains(child.level(), key)) {
                // This child cannot possibly contain firstChar – skip it
                continue;
            }

            int res = traverse(child, maxLevel, posOffset, firstChar);
            if (res != -1) {          // found it in this subtree
                return res;
            }
        }

        // Not found in any child
        return -1;
    }

    public boolean isValidChild(int positionOffset, int intervalIdx, int level, int maxLevel, int workingTreeIdx){
        int spanAll = 1 << maxLevel;              // == tree.intervalSize

        // end position of current interval inside the global stream
        int maxPos = (1 << (maxLevel - level)) * (intervalIdx + 1)
                + workingTreeIdx * spanAll
                - 1;

        return maxPos >= positionOffset;
    }


    // ---- Membership lookup ------------------------------------
    public boolean contains(int level, long key) {
        containCounter++;
        return levels.contains(level, key);
    }

    public double getMembershipFpRate(int level) {
        if(level <= maxDepth()) return this.levels.getLevelFpRate(level);
        else return -1;
    }
}
