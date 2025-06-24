package membership;

import java.util.HashSet;

public class MockMembership implements Membership{

    HashSet<Long> hashSet;

    public MockMembership() {
        this.hashSet = new HashSet<>();
    }


    @Override
    public void init(int n, double p) {
        this.hashSet = new HashSet<>();
    }

    @Override
    public void insert(long key) {
        this.hashSet.add(key);
    }

    @Override
    public boolean contains(long key) {
        return this.hashSet.contains(key);
    }
}
