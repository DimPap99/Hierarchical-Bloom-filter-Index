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

    public LevelDirectory(int depth,
                          IntFunction<M> factory) {
        this.levelFilters = new ArrayList<>(depth);
        for (int l = 0; l < depth; l++) {
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

    public int depth() { return levelFilters.size(); }
}
