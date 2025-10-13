package utilities;

public final class PackedKeyScratch {
    private static final ThreadLocal<long[]> TL2 = ThreadLocal.withInitial(() -> new long[2]);

    public static long[] two() {
        return TL2.get(); // reuse per thread
    }
}
