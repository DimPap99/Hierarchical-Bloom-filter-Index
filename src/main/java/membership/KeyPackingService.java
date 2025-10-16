package membership;

/**
 * KeyPackingService
 *
 * Bit layout of one 64-bit word:
 *   [ interval field (IV_BITS) in the high bits | symbol field (CH_BITS) in the low bits ]
 *
 * Configure once, validate against alphabet size, then pack or unpack.
 */
public final class KeyPackingService {

    // Configuration
    private final int  IV_BITS;
    private final int  CH_BITS;
    private final long IV_MASK;
    private final long CH_MASK;
    private final int  IV_SHIFT;

    // Public container, thin wrapper over long[]
    public static final class Key {
        private final long[] words;

        public Key(long[] words) {
            if (words == null) throw new NullPointerException("words must not be null");
            if (words.length != 1 && words.length != 2) {
                throw new IllegalArgumentException("Key must contain exactly 1 or 2 words");
            }
            this.words = words; // do not clone for performance; treat as read-only
        }

        public long[] words() { return words; }
        public int size()     { return words.length; }
        public long word0()   { return words[0]; }
        public long word1()   { if (words.length != 2) throw new IllegalStateException("not a two-word key"); return words[1]; }

        @Override public String toString() {
            return (words.length == 1)
                    ? "Key[" + Long.toUnsignedString(words[0]) + "]"
                    : "Key[" + Long.toUnsignedString(words[0]) + ", " + Long.toUnsignedString(words[1]) + "]";
        }
    }

    /** Build from maximum level and alphabet size; deepest effective level is maxLevel-1 with 2^{maxLevel-1} intervals */
    public KeyPackingService(int maxLevel, int alphabetSize) {
        if (maxLevel < 0)  throw new IllegalArgumentException("maxLevel must be nonnegative");
        if (alphabetSize < 0) throw new IllegalArgumentException("alphabetSize must be nonnegative");

        int ivBits = Math.max(0, maxLevel - 1);
        if (ivBits > 63)  throw new IllegalArgumentException("interval field would need " + ivBits + " bits");
        int chBits = 64 - ivBits;
        if (chBits <= 0)  throw new IllegalArgumentException("no bits left for symbol field");
        int need = bitsNeededForAlphabet(alphabetSize);
        if (chBits < need) {
            throw new IllegalArgumentException("alphabet needs " + need + " bits but only " + chBits + " available");
        }

        this.IV_BITS  = ivBits;
        this.CH_BITS  = chBits;
        this.IV_MASK  = mask64(IV_BITS);
        this.CH_MASK  = mask64(CH_BITS);
        this.IV_SHIFT = CH_BITS;
    }

    /** Direct constructor when you already know how many high bits go to the interval field */
    public KeyPackingService(int intervalBits, int alphabetSize, boolean directBits) {
        if (!directBits) throw new IllegalArgumentException("use directBits=true for this constructor");
        if (intervalBits < 0 || intervalBits > 63) throw new IllegalArgumentException("intervalBits must be in [0, 63]");
        int chBits = 64 - intervalBits;
        if (chBits <= 0) throw new IllegalArgumentException("no bits left for symbol field");
        int need = bitsNeededForAlphabet(alphabetSize);
        if (chBits < need) throw new IllegalArgumentException("alphabet needs " + need + " bits but only " + chBits + " available");

        this.IV_BITS  = intervalBits;
        this.CH_BITS  = chBits;
        this.IV_MASK  = mask64(IV_BITS);
        this.CH_MASK  = mask64(CH_BITS);
        this.IV_SHIFT = CH_BITS;
    }


    /**
     * always accepts one int and one long and returns a Key
     * If both fields fit the configured split it returns a one-word Key
     * Otherwise it returns a two-word Key {hi, lo} with hi = interval as unsigned 64 bit and lo = symbol as 64 bit
     */
    public Key packAuto(int intervalIdx, long symbolCode) {
        if (fitsOneWord(intervalIdx, symbolCode)) {
            long word = packWord(intervalIdx, symbolCode);
            return new Key(new long[]{ word });
        } else {
            long hi = Integer.toUnsignedLong(intervalIdx);
            long lo = symbolCode;
            return new Key(new long[]{ hi, lo }); // order {hi, lo} as requested
        }
    }


    /** Pack into one word without allocation; caller can use this to fill their own arrays */
    public long packWord(int intervalIdx, long symbolCode) {
        long ivPart = (Integer.toUnsignedLong(intervalIdx) & IV_MASK) << IV_SHIFT;
        long chPart = symbolCode & CH_MASK;
        return ivPart | chPart;
    }

    /** Validate against level then pack into one word Key */
    public Key packAtLevel(int level, int intervalIdx, long symbolCode) {
        validateLevelIndex(level, intervalIdx);
        return new Key(new long[]{ packWord(intervalIdx, symbolCode) });
    }

    // -------- unpacking for one-word keys --------

    public int  interval(Key key) { ensureOneWord(key); return intervalFromWord(key.word0()); }
    public long symbol  (Key key) { ensureOneWord(key); return symbolFromWord  (key.word0()); }

    public int  intervalFromWord(long word) { return (int) ((word >>> IV_SHIFT) & IV_MASK); }
    public long symbolFromWord  (long word) { return (word & CH_MASK); }

    // -------- helpers --------

    public boolean fitsOneWord(int intervalIdx, long symbolCode) {
        // interval fits if its unsigned value uses at most IV_BITS bits
        boolean intervalFits = (Integer.toUnsignedLong(intervalIdx) >>> IV_BITS) == 0;
        // symbol fits if its unsigned value uses at most CH_BITS bits
        boolean symbolFits   = (symbolCode >>> CH_BITS) == 0;
        return intervalFits && symbolFits;
    }

    private static long mask64(int n) {
        if (n <= 0)  return 0L;
        if (n >= 64) return -1L;
        return -1L >>> (64 - n);
    }

    private static int bitsNeededForAlphabet(int alphabetSize) {
        if (alphabetSize <= 1) return 1;
        return 32 - Integer.numberOfLeadingZeros(alphabetSize - 1);
    }

    private void validateLevelIndex(int level, int intervalIdx) {
        if (level < 0) throw new IllegalArgumentException("level must be nonnegative");
        if (level > IV_BITS) {
            throw new IllegalArgumentException("level " + level + " needs " + level + " bits but only " + IV_BITS + " are allocated");
        }
        if (level == 0) {
            if (intervalIdx != 0) throw new IllegalArgumentException("intervalIdx must be zero at level zero");
            return;
        }
        int perLevelMask = (level == 31) ? 0x7FFF_FFFF : ((1 << level) - 1);
        if ((intervalIdx & ~perLevelMask) != 0) {
            throw new IllegalArgumentException("intervalIdx=" + intervalIdx + " does not fit level " + level);
        }
    }

    private static void ensureOneWord(Key key) {
        if (key.size() != 1) throw new IllegalArgumentException("operation requires a one-word key");
    }

    // Introspection, if you need them in tests or logs
    public int  intervalBits()  { return IV_BITS; }
    public int  symbolBits()    { return CH_BITS; }
    public int  intervalShift() { return IV_SHIFT; }
    public long intervalMaskLowBits() { return IV_MASK; }
    public long symbolMaskLowBits()   { return CH_MASK; }
}
