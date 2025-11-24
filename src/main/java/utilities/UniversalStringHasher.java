package utilities;

import java.security.SecureRandom;
import java.util.Objects;
import java.util.SplittableRandom;

// Deterministic 64-bit string hasher with randomly-sampled seeds (not cryptographic).
public final class UniversalStringHasher {

    private static final long FNV_OFFSET_BASIS = 0xcbf29ce484222325L;
    private static final long FNV_PRIME = 0x100000001b3L;

    private final long offsetBasis;
    private final long xorSeed;

    private UniversalStringHasher(long offsetBasis, long xorSeed) {
        this.offsetBasis = offsetBasis;
        this.xorSeed = xorSeed;
    }


    public static UniversalStringHasher create() {
        SecureRandom random = new SecureRandom();
        long seed = random.nextLong();
        long tweak = random.nextLong();
        long offset = mix(FNV_OFFSET_BASIS ^ seed);
        long xor = mix(tweak);
        if (xor == 0) {
            xor = 0x9e3779b97f4a7c15L; // golden ratio constant avoids degenerate xor
        }
        return new UniversalStringHasher(offset, xor);
    }

    // Creates a hasher with deterministic seeds derived from seed (useful for tests).
    public static UniversalStringHasher withSeed(long seed) {
        SplittableRandom random = new SplittableRandom(seed);
        long offset = mix(FNV_OFFSET_BASIS ^ random.nextLong());
        long xor = mix(random.nextLong());
        if (xor == 0) {
            xor = 0x9e3779b97f4a7c15L;
        }
        return new UniversalStringHasher(offset, xor);
    }

    public long hash(String value) {
        Objects.requireNonNull(value, "value");
        long hash = offsetBasis;
        for (int i = 0; i < value.length(); i++) {
            hash ^= value.charAt(i);
            hash *= FNV_PRIME;
        }
        hash ^= xorSeed;
        return mix(hash) & Long.MAX_VALUE; // ensure non-negative domain
    }

    private static long mix(long z) {
        z ^= (z >>> 33);
        z *= 0xff51afd7ed558ccdL;
        z ^= (z >>> 33);
        z *= 0xc4ceb9fe1a85ec53L;
        z ^= (z >>> 33);
        return z;
    }
}
