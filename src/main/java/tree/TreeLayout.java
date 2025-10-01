package tree;

//JIT friendly helper for the layout of the tree
public class TreeLayout{

    private int levels;
    public int baseIntervalSize;

    //when using an estimator we can change the effective root level as well as the effective leaf level
    //we can do the pre emptively, assuming no distribution shift to discard known suboptimal levels close to the roots
    //as well as close to the leafs in case we have ngrams. We can also do this as a post processing step after we finish inserting in the tree. Either
    //way we can achieve memory savings as well as insertion speed up.
    private int effectiveRootLevel;
    private int effectiveLeafLevel;

    public TreeLayout(int levels, int leafSpan) {
        this.levels = levels;
        this.effectiveRootLevel = 0;
        this.effectiveLeafLevel = levels;
        this.baseIntervalSize = leafSpan;
    }



        //Size (in characters) of one interval on {@code level}. Essentially gives the size of a block when we divide
        // the window by W/2^level intervals*/
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

        public void setEffectiveLeafLevel(int effRoot){
            this.effectiveLeafLevel = effectiveLeafLevel;
        }

        public int getEffectiveRootLevel(){ return effectiveRootLevel; }
        public int getEffectiveLeafLevel(){ return effectiveLeafLevel; }

    }

