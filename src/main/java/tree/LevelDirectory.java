package tree;


import membership.Membership;

import java.util.ArrayList;
import java.util.List;
import java.util.function.IntFunction;


/**
 * Owns one {@link Membership} structure per tree level.
 * Generic type {@code K} is the packed-key representation
 */
public final class LevelDirectory<M extends Membership> {

    private final List<M> levelFilters;

    public LevelDirectory(TreeLayout layout,
                          IntFunction<M> factory) {
        int totalLevels = layout.getEffectiveLeafLevel() - layout.getEffectiveRootLevel();
        this.levelFilters = new ArrayList<>(totalLevels);
        for (int l = layout.getEffectiveRootLevel(); l < layout.getEffectiveLeafLevel(); l++) {
            levelFilters.add(factory.apply(l));
        }
    }

    public void insert(int level, long key) {
        levelFilters.get(level).insert(key);
    }

    public boolean contains(int level, long key) {
        return levelFilters.get(level).contains(key);
    }

    public M filter(int level) { return levelFilters.get(level); }

    public void dropFilter(int idx){
        levelFilters.remove(idx);
    }

    public int depth() { return levelFilters.size(); }


    public double getLevelFpRate(int level) {
        return levelFilters.get(level).getFpRate();
    }
}
