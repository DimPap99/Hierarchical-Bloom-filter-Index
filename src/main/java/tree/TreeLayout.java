package tree;

//JIT friendly helper for the layout of the tree
public record TreeLayout(int levels,
                         int baseIntervalSize) {

    /** Size (in characters) of one interval on {@code level}. Essentially gives the size of a block when we divide
     the window by W/2^level intervals*/
    public int intervalSize(int level) {
        if (level < 0 || level >= levels)
            throw new IllegalArgumentException("level " + level + " out of range");
        // Each level halves the span relative to its parent.
        return baseIntervalSize << (levels - 1 - level);
    }

    /** Index of the left child for (level, idx). */
    public int leftChild(int idx) { return idx << 1; }

    /** Index of the right child for (level, idx). */
    public int rightChild(int idx) { return (idx << 1) | 1; }

    /** Number of intervals that exist on {@code level}. */
    public int intervalsOnLevel(int level) { return 1 << level; }


    public int levels() { return levels; }


}
