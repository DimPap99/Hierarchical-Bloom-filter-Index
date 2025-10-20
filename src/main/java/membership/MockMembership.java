package membership;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class MockMembership implements Membership {

    private final HashSet<Long> singleKeys;
    private final HashMap<Long, HashSet<Long>> pairKeys;
    private final HashSet<List<Long>> sequenceKeys;

    public MockMembership() {
        this.singleKeys = new HashSet<>();
        this.pairKeys = new HashMap<>();
        this.sequenceKeys = new HashSet<>();
    }

    @Override
    public void init(int n, double p) {
        singleKeys.clear();
        pairKeys.clear();
        sequenceKeys.clear();
    }

    @Override
    public void insert(long key) {
        singleKeys.add(key);
    }

    @Override
    public boolean contains(long key) {
        return singleKeys.contains(key);
    }

    @Override
    public void insert(long hi, long lo) {
        pairKeys.computeIfAbsent(hi, k -> new HashSet<>()).add(lo);
    }

    @Override
    public boolean contains(long hi, long lo) {
        HashSet<Long> bucket = pairKeys.get(hi);
        return bucket != null && bucket.contains(lo);
    }

    @Override
    public void insert(long[] key) {
        if (key == null || key.length == 0) {
            return;
        }
        if (key.length == 1) {
            insert(key[0]);
            return;
        }
        if (key.length == 2) {
            insert(key[0], key[1]);
            return;
        }
        ArrayList<Long> seq = new ArrayList<>(key.length);
        for (long v : key) {
            seq.add(v);
        }
        sequenceKeys.add(seq);
    }

    @Override
    public boolean contains(long[] key) {
        if (key == null || key.length == 0) {
            return false;
        }
        if (key.length == 1) {
            return contains(key[0]);
        }
        if (key.length == 2) {
            return contains(key[0], key[1]);
        }
        ArrayList<Long> seq = new ArrayList<>(key.length);
        for (long v : key) {
            seq.add(v);
        }
        return sequenceKeys.contains(seq);
    }

    @Override
    public double getFpRate() {
        return 0.0;
    }
}
