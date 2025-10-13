package membership;

public class PackedKey {

    public final long[] key;

    public PackedKey(long[] key) {
        this.key = java.util.Objects.requireNonNull(key);
    }
    public PackedKey(int chunks) {
        if (chunks <= 0) throw new IllegalArgumentException("chunks must be > 0");
        this.key = new long[chunks];
    }

    public static PackedKey of2(long lo, long hi) {
        PackedKey pk = new PackedKey(2);
        pk.key[0] = lo;
        pk.key[1] = hi;
        return pk;
    }
}
