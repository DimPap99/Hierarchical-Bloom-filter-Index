package membership;

/**
 * Dynamically sized 32-bit composite key.
 *
 *
 *  │  level   │   interval     │   char     │
 *
 *     lvBits        ivBits        chBits
 *
 *  Bits are allocated as:
 *      lvBits  = bits needed to store maxLevel (ceil(log2(maxLevel+1)))
 *      ivBits  = maxLevel          (because deepest level has 2^d intervals)
 *      chBits  = remaining bits    (32 − lvBits − ivBits)
 */
public final class Key32 implements IntKey {

    private final int LV_BITS, IV_BITS, CH_BITS;
    private final int IV_MASK, CH_MASK;
    private final int IV_SHIFT, LV_SHIFT;

    public Key32(int maxLevel, int minCharBits) {
        if (maxLevel < 0)
            throw new IllegalArgumentException("maxLevel must be ≥0");
        if (minCharBits <= 0 || minCharBits >= 32)
            throw new IllegalArgumentException("minCharBits out of range");

        LV_BITS = Integer.SIZE - Integer.numberOfLeadingZeros(maxLevel);
        IV_BITS = maxLevel;                    // 2^d intervals at depth d
        CH_BITS = 32 - LV_BITS - IV_BITS;

        if (CH_BITS < minCharBits) {
            throw new IllegalArgumentException(
                    "32-bit key cannot host: level=" + maxLevel +
                            ", need ≥" + minCharBits + " char-bits (have only " + CH_BITS + ')');
        }

        CH_MASK  = (1 << CH_BITS) - 1;
        IV_MASK  = (1 << IV_BITS) - 1;

        IV_SHIFT = CH_BITS;
        LV_SHIFT = CH_BITS + IV_BITS;
    }
    @Override
    public int pack(int level, int intervalIdx, char ch) {
        return (level       << LV_SHIFT) |
                (intervalIdx << IV_SHIFT) |
                (ch & CH_MASK);
    }

    @Override public int  level(int key)    { return  key >>> LV_SHIFT; }
    @Override public int  interval(int key) { return (key >>> IV_SHIFT) & IV_MASK; }
    @Override public char ch(int key)       { return (char) (key & CH_MASK); }

    public int charBits() { return CH_BITS; }
    public int intervalBits() { return IV_BITS; }
    public int levelBits() { return LV_BITS; }
}
