package membership;

/** Specialised for 64-bit primitive keys. */
public interface LongKey extends Key {
    long pack(int level, int intervalIdx, char ch);
    int  level(long key);
    int  interval(long key);
    char ch(long key);
}



