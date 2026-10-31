package membership;

import java.util.Locale;
import java.util.function.Supplier;

// Which Bloom filter HBI's levels use (--filter in the harness and MemoryBenchmark).
public enum FilterKind {
    // membership.BloomFilter: k positions over the whole bit array.
    BLOOM,
    // membership.BlockedBloomFilter: all k bits of a key in one 64- or 512-bit block, sized for the same rate.
    BLOCKED;

    public Supplier<Membership> supplier() {
        return switch (this) {
            case BLOOM -> BloomFilter::new;
            case BLOCKED -> BlockedBloomFilter::new;
        };
    }

    public String label() {
        return name().toLowerCase(Locale.ROOT);
    }

    public static FilterKind parse(String value) {
        if (value == null) return BLOOM;
        return switch (value.trim().toLowerCase(Locale.ROOT)) {
            case "bloom", "standard" -> BLOOM;
            case "blocked" -> BLOCKED;
            default -> throw new IllegalArgumentException("--filter must be bloom or blocked, not " + value);
        };
    }
}
