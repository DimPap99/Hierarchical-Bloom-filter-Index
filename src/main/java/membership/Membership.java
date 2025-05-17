package membership;

public interface Membership {
    void init(int n, double p);
    void insert(String key);
    boolean contains(String key);
}
