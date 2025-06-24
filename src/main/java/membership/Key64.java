package membership;

/**
 * Dynamically sized 64-bit composite key.
 *
 * Layout (high → low):   [ level |    interval    |   char   ]
 *                        63 … LSHIFT … IVSHIFT … 0
 *
 *   lvBits = bitsNeeded(maxLevel)      // fits depth number d
 *   ivBits = maxLevel                  // 2^d intervals at deepest level
 *   chBits = 64 - lvBits - ivBits      // whatever is left
 *
 * Throws if chBits < minCharBits or if 64 bits are not enough.
 */
public final class Key64 implements LongKey {

    /* -------------------- immutable layout -------------------------- */
    int LV_BITS, IV_BITS, CH_BITS;
    long IV_MASK, CH_MASK;
    int  IV_SHIFT, LV_SHIFT;

    /* helper: bits required to encode a positive int v (0 → 0 bits) */
    private static int bitsNeeded(int v) {
        return (v == 0) ? 0 : Integer.SIZE - Integer.numberOfLeadingZeros(v);
    }

    /**
     * @param maxLevel     deepest level d you must support (≥0)
     * @param minCharBits  minimum bits you need for the character (typically 8 or 16)
     */
    public Key64(int maxLevel, int minCharBits) {
        if (maxLevel < 0)                       throw new IllegalArgumentException("maxLevel < 0");
        if (minCharBits <= 0 || minCharBits >= 64)
            throw new IllegalArgumentException("minCharBits out of range");

        LV_BITS = bitsNeeded(maxLevel);         // store the depth number
        IV_BITS = maxLevel;                     // need up to 2^d intervals
        CH_BITS = 64 - LV_BITS - IV_BITS;       // remaining bits

        if (CH_BITS < minCharBits) {
            throw new IllegalArgumentException(
                    "64-bit key cannot host: level=" + maxLevel +
                            " with ≥" + minCharBits + " char bits (only " + CH_BITS + " left)");
        }

        CH_MASK  = (1L << CH_BITS) - 1L;
        IV_MASK  = (IV_BITS == 64) ? -1L : (1L << IV_BITS) - 1L;

        IV_SHIFT = CH_BITS;
        LV_SHIFT = CH_BITS + IV_BITS;
    }

    /* ==================== LongKey implementation ==================== */

    @Override
    public long pack(int level, int intervalIdx, char ch) {
        return  ((long) level       << LV_SHIFT) |
                ((long) intervalIdx << IV_SHIFT) |
                ((long) ch & CH_MASK);
    }

    @Override
    public int level(long key) {
        return (int) (key >>> LV_SHIFT);
    }

    @Override
    public int interval(long key) {
        return (int) ((key >>> IV_SHIFT) & IV_MASK);
    }

    @Override
    public char ch(long key) {
        return (char) (key & CH_MASK);
    }

    /* ---------------- optional introspection ------------------------ */
    public int charBits()     { return CH_BITS; }
    public int intervalBits() { return IV_BITS; }
    public int levelBits()    { return LV_BITS; }
}
