package membership;

public interface Membership {
    void init(int n, double p);

    // As init(n, p), with the hash functions drawn from seed, so the same seed gives the same filter.
    default void init(int n, double p, long seed) {
        init(n, p);
    }
    void insert(long key);
    boolean contains(long key);

    void insert(long[] key);
    void insert(long hi, long lo);

    boolean contains(long[] key);
    boolean contains(long hi, long lo);

    double getFpRate();
}
