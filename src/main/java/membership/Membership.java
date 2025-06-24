package membership;

public interface Membership {
    void init(int n, double p);
    void insert(int key);
    boolean contains(int key);
}
