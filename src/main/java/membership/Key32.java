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

    private final int IV_BITS;
    private final int CH_BITS;

    private final int IV_MASK;   // mask with IV_BITS low bits set
    private final int CH_MASK;   // mask with CH_BITS low bits set

    private final int IV_SHIFT;  // number of bits the interval field is shifted left

    /**
     * @param maxLevel     maximum tree depth used to size the interval field.
     *                     the deepest level has 2^maxLevel intervals, so we need exactly maxLevel bits
     * @param minCharBits  minimum number of bits you require for the character or token field
     */
    public Key32(int maxLevel, int minCharBits) {
        if (maxLevel < 0) {
            throw new IllegalArgumentException("maxLevel must be greater than or equal to zero");
        }
        if (minCharBits <= 0 || minCharBits > 32) {
            throw new IllegalArgumentException("minCharBits must be in the range [1, 32]");
        }

        // allocate bits
        //in case where we dont use 1 bloom filter for all the levels then we can ommits those bits from the key cause the assignment happens manually
        this.IV_BITS = maxLevel;                 // needs exactly maxLevel bits to count up to 2^maxLevel - 1
        if (IV_BITS > 31) {
            // shifting an int by thirty two is undefined in Java, so we disallow IV_BITS == 32
            throw new IllegalArgumentException("interval field would require " + IV_BITS + " bits, which does not fit into thirty two bit packing");
        }
        this.CH_BITS = 32 - IV_BITS;

        if (CH_BITS < minCharBits) {
            throw new IllegalArgumentException(
                    "thirty two bit key cannot satisfy requirements. have CH_BITS=" + CH_BITS
                            + " after reserving IV_BITS=" + IV_BITS
                            + " but need at least minCharBits=" + minCharBits);
        }

        // safe masks for any bit width in [0, 32]
        this.IV_MASK = maskNBits(IV_BITS);
        this.CH_MASK = maskNBits(CH_BITS);

        // layout: [ interval bits | character bits ]
        this.IV_SHIFT = CH_BITS;
    }

    /** Pack interval index and character bits. The level parameter is ignored by design. */
    @Override
    public int pack(int level, int intervalIdx, char ch) {
        // keep only the configured number of character bits
        int chPart = ch & CH_MASK;
        int ivPart = (intervalIdx & IV_MASK) << IV_SHIFT;
        return ivPart | chPart;
    }

    /** Convenience overload that makes the no level design explicit. */
    public int pack(int intervalIdx, int tokenBits) {
        int chPart = tokenBits & CH_MASK;
        int ivPart = (intervalIdx & IV_MASK) << IV_SHIFT;
        return ivPart | chPart;
    }

    /** Returns zero, since level is not encoded in this variant. */
    @Override
    public int level(int key) {
        return 0;
    }

    /** Extract the interval index from the packed key. */
    @Override
    public int interval(int key) {
        return (key >>> IV_SHIFT) & IV_MASK;
    }

    /** Extract the character or token identifier from the packed key. */
    @Override
    public char ch(int key) {
        return (char) (key & CH_MASK);
    }

    /** Number of bits reserved for the character or token field. */
    public int charBits() { return CH_BITS; }

    /** Number of bits reserved for the interval index. */
    public int intervalBits() { return IV_BITS; }

    /** Always zero in this design. Included for symmetry with the original class. */
    public int levelBits() { return 0; }

    /** Build an integer mask with the lowest {@code n} bits set. Works for n in [0, 32]. */
    private static int maskNBits(int n) {
        if (n <= 0) return 0;
        if (n >= 32) return -1;                  // all ones
        // use unsigned shift to avoid overflow edge cases
        return -1 >>> (32 - n);
    }
}
