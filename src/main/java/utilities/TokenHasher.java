package utilities;

import org.apache.commons.codec.digest.MurmurHash3;

import java.nio.charset.StandardCharsets;

/**
 * Stateless hashing utility that maps arbitrary strings to stable positive integers/longs.
 * We expose both 63-bit and 31-bit flavours so callers that must interoperate with legacy
 * int-based APIs can do so, while the streaming suffix structures can keep the full 64-bit
 * entropy and rely on per-tree remapping to obtain dense alphabets.
 */
public final class TokenHasher {

    private TokenHasher() {
        // Utility class; do not instantiate.
    }

    /**
     * Hash a token to a strictly positive 63-bit value. Returns a number in [1, 2^63 - 1].
     */
    public static long hashToPositiveLong(String token) {
        if (token == null) {
            throw new IllegalArgumentException("token must not be null");
        }
        byte[] bytes = token.getBytes(StandardCharsets.UTF_8);
        long[] h = MurmurHash3.hash128x64(bytes);
        long mixed = mix64(h[0] ^ h[1]);
        long value = mixed & Long.MAX_VALUE;
        return (value == 0L) ? 1L : value;
    }

    /**
     * Hash a token to a positive 31-bit value for compatibility with components that still
     * expect ints. The value is derived from the 63-bit hash; collisions are therefore as
     * rare as they would be with a simple xor-fold.
     */
    public static int hashToPositiveInt(String token) {
        long value = hashToPositiveLong(token);
        int folded = (int) (value ^ (value >>> 32));
        folded &= 0x7FFF_FFFF;
        return (folded == 0) ? 1 : folded;
    }

    private static long mix64(long z) {
        z ^= (z >>> 33);
        z *= 0xff51afd7ed558ccdL;
        z ^= (z >>> 33);
        z *= 0xc4ceb9fe1a85ec53L;
        z ^= (z >>> 33);
        return z;
    }
}
