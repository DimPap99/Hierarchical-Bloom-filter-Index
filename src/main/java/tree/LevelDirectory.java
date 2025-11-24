package tree;


import membership.Membership;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntFunction;

// Owns one Membership filter per tree level.
public final class LevelDirectory<M extends Membership> {

    private final List<M> levelFilters;
    private int firstLevel;

    public LevelDirectory(TreeLayout layout,
                          IntFunction<M> factory) {
        this.firstLevel = layout.getEffectiveRootLevel();
        int totalLevels = layout.getEffectiveLeafLevel() - this.firstLevel;
        this.levelFilters = new ArrayList<>(totalLevels);
        for (int l = 0; l < totalLevels; l++) {
            levelFilters.add(factory.apply(this.firstLevel + l));
        }
    }

    private int toLocalIndex(int level) {
        int idx = level - firstLevel;
        if (idx < 0 || idx >= levelFilters.size()) {
            throw new IllegalArgumentException("level " + level + " out of range");
        }
        return idx;
    }

    public void insert(int level, long key) {
        levelFilters.get(toLocalIndex(level)).insert(key);
    }

    public void insert(int level, long hi, long lo) {
        levelFilters.get(toLocalIndex(level)).insert(hi, lo);
    }


    public boolean contains(int level, long key) {
        return levelFilters.get(toLocalIndex(level)).contains(key);
    }

    public boolean contains(int level, long hi, long lo) {
        return levelFilters.get(toLocalIndex(level)).contains(hi, lo);
    }


    public M filter(int level) { return levelFilters.get(toLocalIndex(level)); }

    public void dropLevel(int level){
        int idx = toLocalIndex(level);
        levelFilters.remove(idx);
        if (idx == 0) {
            firstLevel++;
        }
    }

    public int depth() { return levelFilters.size(); }


    public double getLevelFpRate(int level) {
        return levelFilters.get(toLocalIndex(level)).getFpRate();
    }
}
