package membership;

/** Specialised for 32-bit primitive keys (avoids boxing). */
public interface IntKey extends Key {
    int  pack(int level, int intervalIdx, char ch);
    int  level(int key);
    int  interval(int key);
    char ch(int key);
}


