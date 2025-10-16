package membership;

import java.util.HashSet;

public class MockMembership implements Membership{

    HashSet<Long> hashSet;

    public MockMembership() {
        this.hashSet = new HashSet<>();
    }


    public void init(int n, double p) {
        this.hashSet = new HashSet<>();
    }

    public void insert(long key) {
        this.hashSet.add(key);
    }

    public boolean contains(long key) {
        return this.hashSet.contains(key);
    }

    @Override
    public void insert(long[] key) {

    }

    @Override
    public void insert(long hi, long lo) {

    }

    @Override
    public boolean contains(long[] key) {
        return false;
    }

    @Override
    public boolean contains(long hi, long lo) {
        return false;
    }

    @Override
    public double getFpRate() {
        return 0;
    }
}
