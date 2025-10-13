package membership;

/**
 * Dynamically sized 64-bit composite key without level.
 *
 * Layout (high → low):   [     interval     |      char      ]
 *                        63 … IVSHIFT ……………………………… 0
 *
 *   IV_BITS = maxLevel - 1            // deepest effective level has 2^(maxLevel-1) intervals
 *   CH_BITS = 64 - IV_BITS
 *
 * Throws if CH_BITS cannot accommodate the requested alphabet size.
 *
 * This variant does not encode the level:
 *   - pack(level, intervalIdx, ch) ignores the level
 *   - level(key) always returns 0
 *   - levelBits() returns 0
 */
public final class Key64 implements LongKey {

    // immutable layout
    private final int  IV_BITS;
    private final int  CH_BITS;
    private final long IV_MASK;
    private final long CH_MASK;
    private final int  IV_SHIFT;

    /** Build a mask with the lowest n bits set. Works for n in [0, 64]. */
    private static long mask64(int n) {
        if (n <= 0)  return 0L;
        if (n >= 64) return -1L;               // all ones
        return -1L >>> (64 - n);
    }

    /** Bits required to encode values in [0, alphabetSize-1]. */
    private static int bitsNeededForAlphabet(int alphabetSize) {
        if (alphabetSize <= 1) return 1;       // at least one bit
        return 32 - Integer.numberOfLeadingZeros(alphabetSize - 1);
    }

    /**
     * @param maxLevel  deepest tree level you plan to support (≥ 0)
     *                  this class uses IV_BITS = maxLevel - 1 to match your existing convention
     * @param alphabet  number of distinct token values that must fit in the char field
     */
    public Key64(int maxLevel, int alphabet) {
        if (maxLevel < 0) {
            throw new IllegalArgumentException("maxLevel must be ≥ 0");
        }
        if (alphabet < 0) {
            throw new IllegalArgumentException("alphabet must be ≥ 0");
        }

        // no-level design: only interval and char are packed
        this.IV_BITS = Math.max(0, maxLevel - 1);  // deepest effective level has 2^(maxLevel-1) intervals
        if (IV_BITS > 63) {
            throw new IllegalArgumentException("interval field would need " + IV_BITS + " bits which does not fit in 64 bits");
        }

        this.CH_BITS = 64 - IV_BITS;
        if (CH_BITS <= 0) {
            throw new IllegalArgumentException("no bits left for char field. reduce maxLevel");
        }

        int requiredCharBits = bitsNeededForAlphabet(alphabet);
        if (CH_BITS < requiredCharBits) {
            throw new IllegalArgumentException(
                    "64-bit key cannot host alphabet size " + alphabet +
                            " with maxLevel=" + maxLevel +
                            " because only " + CH_BITS + " char bits are available");
        }

        this.IV_MASK  = mask64(IV_BITS);
        this.CH_MASK  = mask64(CH_BITS);
        this.IV_SHIFT = CH_BITS;
    }

    /** Pack interval index and char. The level parameter is ignored by design. */
    @Override
    public long pack(int level, int intervalIdx, int ch) {
        long ivPart = (Integer.toUnsignedLong(intervalIdx) & IV_MASK) << IV_SHIFT;
        long chPart = (Integer.toUnsignedLong(ch)         & CH_MASK);
        return ivPart | chPart;
    }

    /** Extract level which is not encoded in this variant and therefore always zero. */
    @Override
    public int level(long key) {
        return 0;
    }

    /** Extract the interval index from the packed key. */
    @Override
    public int interval(long key) {
        long iv = (key >>> IV_SHIFT) & IV_MASK;
        return (int) iv;
    }

    /** Extract the char or token identifier from the packed key. */
    @Override
    public char ch(long key) {
        return (char) (key & CH_MASK);
    }

    /** Number of bits reserved for the char field. */
    public int charBits() {
        return CH_BITS;
    }

    /** Number of bits reserved for the interval field. */
    public int intervalBits() {
        return IV_BITS;
    }

    /** Always zero in this no-level design. */
    public int levelBits() {
        return 0;
    }

    /**
     * Maximum number of distinct char values that can be represented for a given maxLevel.
     * Clamped to Integer.MAX_VALUE to keep the original return type.
     */
    public static int getMaxPossibleChar(int maxLevel) {
        int ivBits = Math.max(0, maxLevel - 1);
        int chBits = 64 - ivBits;
        if (chBits >= 31) return Integer.MAX_VALUE;  // avoid overflow of 1<<chBits in int
        return 1 << chBits;
    }
}
