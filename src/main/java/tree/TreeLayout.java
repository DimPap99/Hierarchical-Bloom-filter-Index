package tree;

// Helper for the layout of the tree.
public class TreeLayout{

    private int levels;
    public int baseIntervalSize;

    // When using an estimator we can change the effective root and leaf levels.
    private int effectiveRootLevel;
    private int effectiveLeafLevel;

    public TreeLayout(int levels, int leafSpan) {
        this.levels = levels;
        this.effectiveRootLevel = 0;
        this.effectiveLeafLevel = levels;
        this.baseIntervalSize = leafSpan;
    }

        // Size (in characters) of one interval on level.
        public int intervalSize(int level) {
            if (level < 0 || level >= levels)
                throw new IllegalArgumentException("level " + level + " out of range");
            // Each level halves the span relative to its parent.
            return baseIntervalSize << (levels - 1 - level);
        }

        //Index of the left child for (level, idx).
        public int leftChild(int idx) { return idx << 1; }

        //Index of the right child for (level, idx).
        public int rightChild(int idx) { return (idx << 1) | 1; }

        //Number of intervals that exist on {@code level}.
        public int intervalsOnLevel(int level) { return 1 << level; }


        public int levels() { return levels; }


        public void setEffectiveRootLevel(int effRoot){
            this.effectiveRootLevel = effRoot;
        }

        public void setEffectiveLeafLevel(int effLeaf){
            this.effectiveLeafLevel = effLeaf;
        }

        public int getEffectiveRootLevel(){ return effectiveRootLevel; }
        public int getEffectiveLeafLevel(){ return effectiveLeafLevel; }

    }
