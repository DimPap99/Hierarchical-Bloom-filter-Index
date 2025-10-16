package utilities;

import org.apache.commons.codec.digest.MurmurHash3;

import java.nio.charset.StandardCharsets;
import java.security.SecureRandom;

/**
 * StringKeyMapper
 *
 * Chooses the hash width from a target collision probability epsilon and the expected
 * number of distinct tokens n using the birthday bound
 *
 *   M  >=  n(n - 1) / (2 * -ln(1 - epsilon))
 *   bits = ceil(log2(M))
 *
 * If bits <= 31 it maps to int, otherwise to long.
 * You can choose Apache Commons Codec MurmurHash3 or Carter–Wegman universal hashing.
 *
 * NOTE: If requiredBits > 31, mapToInt(String) will still return a 31-bit value (best-effort)
 * and print a one-time warning with the implied collision probability for 31 bits.
 */
public final class StringKeyMapper {

    public enum Algo { MURMUR3_APACHE, CARTER_WEGMAN }

    private final Algo algo;
    private final long expectedDistinct;
    private final double epsilon;
    private final int requiredBits;     // computed from epsilon and n
    private final boolean use32bit;     // true when requiredBits <= 31
    private final long mask64;          // (1<<requiredBits)-1 or all ones for 64 bits
    private final int mask32;           // when use32bit: (1<<requiredBits)-1; otherwise unused
    private final int mask32BestEffort; // always 0x7FFFFFFF for safe 31-bit ints

    private boolean warnedOnce = false;

    // ---- Carter–Wegman state (only used when algo == CARTER_WEGMAN) ----
    private static final long P61  = (1L << 61) - 1;
    private static final long MASK = P61;

    // two independent affine parameters produce two independent 64-bit hashes
    private final long cwA1, cwB1;
    private final long cwA2, cwB2;

    /** Construct with expected distinct tokens n and target collision probability epsilon. Uses Murmur by default. */
    public StringKeyMapper(long expectedDistinct, double epsilon) {
        this(expectedDistinct, epsilon, Algo.MURMUR3_APACHE);

    }

    public static StringKeyMapper.Algo chooseAlgo(int avgLen, boolean haveBytes) {
        if (haveBytes) return StringKeyMapper.Algo.MURMUR3_APACHE;
        if (avgLen >= 24) return StringKeyMapper.Algo.MURMUR3_APACHE;
        return StringKeyMapper.Algo.CARTER_WEGMAN;
    }

    /** Same as above but lets you pick the hashing algorithm. */
    public StringKeyMapper(long expectedDistinct, double epsilon, Algo algo) {
        if (expectedDistinct < 0) throw new IllegalArgumentException("expectedDistinct must be nonnegative");
        if (!(epsilon > 0.0 && epsilon < 1.0)) throw new IllegalArgumentException("epsilon must be in (0, 1)");

        this.expectedDistinct = expectedDistinct;
        this.epsilon = epsilon;
        this.requiredBits = bitsForEpsilon(expectedDistinct, epsilon);
        this.use32bit = requiredBits <= 31;

        // Build masks once
        if (requiredBits >= 64) {
            this.mask64 = -1L; // all bits
        } else {
            this.mask64 = (1L << requiredBits) - 1L;
        }
        if (use32bit) {
            this.mask32 = (requiredBits == 31) ? 0x7FFFFFFF : ((1 << requiredBits) - 1);
        } else {
            this.mask32 = 0; // not used when !use32bit
        }
        this.mask32BestEffort = 0x7FFFFFFF; // always safe 31-bit positive mask

        this.algo = algo;

        // Seed Carter–Wegman parameters only if requested
        if (algo == Algo.CARTER_WEGMAN) {
            SecureRandom rng = new SecureRandom();

            long a1, b1, a2, b2;
            do { a1 = rng.nextLong() & MASK; } while (a1 == 0L);
            b1 = rng.nextLong() & MASK;
            do { a2 = rng.nextLong() & MASK; } while (a2 == 0L);
            b2 = rng.nextLong() & MASK;

            this.cwA1 = a1; this.cwB1 = b1;
            this.cwA2 = a2; this.cwB2 = b2;
        } else {
            this.cwA1 = this.cwB1 = this.cwA2 = this.cwB2 = 0L;
        }
    }

    /** Return true if this mapper will produce int identifiers while honoring epsilon (i.e., requiredBits <= 31). */
    public boolean is32bit() { return use32bit; }

    /** The number of hash bits chosen from epsilon and n. */
    public int requiredBits() { return requiredBits; }

    /**
     * Map a string to an int. If requiredBits <= 31, masks exactly to requiredBits.
     * If requiredBits > 31, prints a one-time warning and returns a best-effort 31-bit value.
     */
    public int mapToInt(String s) {
        long h64 = hash64(s);

        if (use32bit) {
            long v = h64 & (long) mask32;   // exact mask to requiredBits
            return (int) v;
        } else {
            if (!warnedOnce) {
                double err31 = birthdayCollisionProbability(this.expectedDistinct, 31);
                System.out.println(
                        "Mapper requires " + requiredBits +
                                " bits to meet epsilon=" + epsilon +
                                " for n=" + expectedDistinct +
                                ". Returning 31-bit int (best-effort). Expected collision prob at 31 bits: " + err31
                );
                warnedOnce = true;
            }
            // best-effort: always non-negative 31-bit id
            return (int) (h64 & (long) mask32BestEffort);
        }
    }

    /** Map a string to a long in the range [0, 2^{requiredBits}-1]. Always valid and honors epsilon. */
    public long mapToLong(String s) {
        if (!warnedOnce) {
            double err31 = birthdayCollisionProbability(this.expectedDistinct, 63);
            System.out.println(
                    "Mapper requires " + requiredBits +
                            " bits to meet epsilon=" + epsilon +
                            " for n=" + expectedDistinct +
                            ". Returning 31-bit int (best-effort). Expected collision prob at 63 bits: " + err31
            );
            warnedOnce = true;
        }
        long h64 = hash64(s);
        return h64 & mask64;
    }

    /**
     * Produce two independent 64-bit hashes for the same string.
     * Useful for double hashing in Bloom filters. Unmasked (caller decides use).
     */
    public long[] mapToTwoLongs(String s) {
        if (algo == Algo.MURMUR3_APACHE) {
            byte[] bytes = s.getBytes(StandardCharsets.UTF_8);
            return MurmurHash3.hash128x64(bytes); // {h1, h2}
        } else {
            long h1 = cwHashUtf8To61(s, cwA1, cwB1);
            long h2 = cwHashUtf8To61(s, cwA2, cwB2);
            return new long[]{h1, h2};
        }
    }

    /** Convenience: estimated collision probability for the configured (epsilon-derived) bit width. */
    public double collisionProbability() {
        return birthdayCollisionProbability(expectedDistinct, requiredBits);
    }

    /* ====================== hashing backends ====================== */

    private long hash64(String s) {
        if (algo == Algo.MURMUR3_APACHE) {
            // hash128 then fold to one 64-bit word
            byte[] bytes = s.getBytes(StandardCharsets.UTF_8);
            long[] hh = MurmurHash3.hash128x64(bytes);
            return hh[0] ^ hh[1];
        } else {
            return cwHashUtf8To61(s, cwA1, cwB1);
        }
    }

    // Carter–Wegman polynomial over UTF-8 then affine a*x + b modulo 2^61 - 1
    private static long cwHashUtf8To61(String s, long a, long b) {
        long h = b;
        for (int i = 0, len = s.length(); i < len; i++) {
            int ch = s.charAt(i);
            if (ch < 0x80) {
                h = cwMulAdd(h, a, ch);
            } else if (ch < 0x800) {
                h = cwMulAdd(h, a, 0xC0 | (ch >>> 6));
                h = cwMulAdd(h, a, 0x80 | (ch & 0x3F));
            } else if (ch < 0xD800 || ch > 0xDFFF) {
                h = cwMulAdd(h, a, 0xE0 | (ch >>> 12));
                h = cwMulAdd(h, a, 0x80 | ((ch >>> 6) & 0x3F));
                h = cwMulAdd(h, a, 0x80 | (ch & 0x3F));
            } else {
                if (i + 1 < len) {
                    int low = s.charAt(i + 1);
                    if ((ch & 0xFC00) == 0xD800 && (low & 0xFC00) == 0xDC00) {
                        int cp = 0x10000 + ((ch & 0x3FF) << 10) + (low & 0x3FF);
                        i++;
                        h = cwMulAdd(h, a, 0xF0 | (cp >>> 18));
                        h = cwMulAdd(h, a, 0x80 | ((cp >>> 12) & 0x3F));
                        h = cwMulAdd(h, a, 0x80 | ((cp >>> 6)  & 0x3F));
                        h = cwMulAdd(h, a, 0x80 | (cp & 0x3F));
                        continue;
                    }
                }
                // replacement character U+FFFD encoded in UTF-8
                h = cwMulAdd(h, a, 0xEF);
                h = cwMulAdd(h, a, 0xBF);
                h = cwMulAdd(h, a, 0xBD);
            }
        }
        return h;
    }

    private static long cwMulAdd(long h, long a, int x) {
        long lo = h * a;
        long hi = Math.multiplyHigh(h, a);
        long r  = reduce128mod61(hi, lo);
        long y  = r + x;
        return y >= P61 ? y - P61 : y;
    }

    private static long reduce128mod61(long hi, long lo) {
        long lo_lo = lo & MASK;  // low 61 bits
        long lo_hi = lo >>> 61;  // high 3 bits
        long hi8   = hi << 3;    // * 2^3
        long hi8_lo = hi8 & MASK;
        long hi8_hi = hi >>> 58; // >> 61 then * 2^3
        long s = lo_lo + lo_hi + hi8_lo + hi8_hi;
        s = (s & MASK) + (s >>> 61);
        return s >= P61 ? s - P61 : s;
    }

    /* ====================== birthday bound helpers ====================== */

    /** Collision probability for n distinct items in a space of 2^{bits}. */
    public static double birthdayCollisionProbability(long n, int bits) {
        if (n <= 1) return 0.0;
        if (bits <= 0) return 1.0;

        // M = 2^bits as a double; safe up to 2^1023 in double exponent
        final double space = Math.scalb(1.0, bits); // 2^{bits}
        // Use double for numerator; this is ~ n^2/2 and can exceed 2^53 for large n,
        // but we only need it as a floating ratio below.
        final double nn = (double) n;
        final double num = 0.5d * nn * (nn - 1.0d); // n(n-1)/2

        // If num/space is huge, probability is essentially 1.
        // Use a threshold where exp(-x) underflows or is numerically negligible.
        final double x = num / space; // nonnegative
        if (x >= 40.0) { // exp(-40) ~ 4.2e-18, far below double precision concerns
            return 1.0;
        }

        // P = 1 - exp(-x); use -expm1(-x) for accuracy when x is small
        return -Math.expm1(-x);
    }


    /** Bits required so that collision probability <= epsilon for n items. */
    public static int bitsForEpsilon(long n, double epsilon) {
        if (n <= 1) return 1;
        if (!(epsilon > 0.0 && epsilon < 1.0)) {
            throw new IllegalArgumentException("epsilon must be in (0, 1)");
        }

        // When epsilon is very small, -ln(1 - eps) ~ eps; using log1p improves accuracy.
        final double denom = -Math.log1p(-epsilon); // = -ln(1 - eps)
        final double nn = (double) n;
        final double requiredStates = (nn * (nn - 1.0d)) / (2.0d * denom); // M

        // bits = ceil(log2(M))
        final double bitsExact = Math.log(requiredStates) / Math.log(2.0d);
        int bits = (int) Math.ceil(bitsExact);

        if (bits < 1)  bits = 1;
        if (bits > 64) bits = 64; // clamp to 64 since you only mask to 64 in this class
        return bits;
    }

}
