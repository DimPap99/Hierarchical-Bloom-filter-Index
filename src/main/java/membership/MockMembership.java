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
    public double getFpRate() {
        return 0;
    }
}
