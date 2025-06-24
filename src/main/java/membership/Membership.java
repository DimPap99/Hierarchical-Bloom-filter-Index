package membership;

public interface Membership {
    void init(int n, double p);
    void insert(long key);
    boolean contains(long key);
}
