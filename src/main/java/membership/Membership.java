package membership;

public interface Membership {
    void init(int n, double p);
    void insert(long key);
    boolean contains(long key);

    void insert(long[] key);
    void insert(long hi, long lo);

    boolean contains(long[] key);
    boolean contains(long hi, long lo);

    double getFpRate();
}
