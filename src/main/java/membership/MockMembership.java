package membership;

import java.util.HashSet;

public class MockMembership implements Membership{

    HashSet<Integer> hashSet;

    public MockMembership() {
        this.hashSet = new HashSet<>();
    }


    @Override
    public void init(int n, double p) {
        this.hashSet = new HashSet<>();
    }

    @Override
    public void insert(int key) {
        this.hashSet.add(key);
    }

    @Override
    public boolean contains(int key) {
        return this.hashSet.contains(key);
    }
}
