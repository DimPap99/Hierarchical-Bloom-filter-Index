package tree.ssws;

import utilities.UniversalStringHasher;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Objects;
import java.util.SplittableRandom;

/**
 * Maps the (possibly large) token alphabet present in a window substring into a dense integer range
 * with high probability of being collision free. The implementation follows the high-probability
 * universal hashing scheme outlined in the SSWSI paper (Section 5). A random affine hash in the
 * 64-bit domain is sampled; collisions are checked explicitly and resampling occurs when needed.
 */
final class TokenRanker {

    private static final int MAX_ATTEMPTS = 32;

    private TokenRanker() {
        // utility class
    }

    static MappingResult rank(String[] text, int windowSize) {
        if (text == null) {
            throw new IllegalArgumentException("text cannot be null");
        }
        if (text.length == 0) {
            return new MappingResult(new int[0], Translator.fromDirect(new HashMap<>()));
        }

        String[] distinct = Arrays.copyOf(text, text.length);
        Arrays.sort(distinct);
        int distinctCount = deduplicate(distinct);
        distinct = Arrays.copyOf(distinct, distinctCount);

        SplittableRandom rng = new SplittableRandom();

        for (int attempt = 0; attempt < MAX_ATTEMPTS; attempt++) {
            HashFunction hash = HashFunction.random(rng);
            if (isInjective(hash, distinct)) {
                LongIntHashTable table = new LongIntHashTable(distinctCount);
                for (int i = 0; i < distinctCount; i++) {
                    long hashed = hash.hash(distinct[i]);
                    table.put(hashed, i);
                }

                int[] mapped = new int[text.length];
                for (int i = 0; i < text.length; i++) {
                    long hashed = hash.hash(text[i]);
                    int rank = table.get(hashed);
                    if (rank < 0) {
                        throw new IllegalStateException("Missing rank for token: " + text[i]);
                    }
                    mapped[i] = rank;
                }
                return new MappingResult(mapped, Translator.fromHash(hash, table));
            }
        }

        // Extremely unlikely fallback: rank by sorted order deterministically.
        HashMap<String, Integer> direct = new HashMap<>(distinctCount * 2);
        for (int i = 0; i < distinctCount; i++) {
            direct.put(distinct[i], i);
        }
        int[] mapped = new int[text.length];
        for (int i = 0; i < text.length; i++) {
            Integer rank = direct.get(text[i]);
            if (rank == null) {
                throw new IllegalStateException("Token not found during deterministic ranking");
            }
            mapped[i] = rank;
        }
        return new MappingResult(mapped, Translator.fromDirect(direct));
    }

    private static int deduplicate(String[] values) {
        if (values.length == 0) {
            return 0;
        }
        int write = 1;
        String last = values[0];
        for (int i = 1; i < values.length; i++) {
            String v = values[i];
            if (!Objects.equals(v, last)) {
                values[write++] = v;
                last = v;
            }
        }
        return write;
    }

    private static boolean isInjective(HashFunction hash, String[] distinct) {
        LongHashSet seen = new LongHashSet(distinct.length * 2);
        for (String token : distinct) {
            long hashed = hash.hash(token);
            if (!seen.add(hashed)) {
                return false;
            }
        }
        return true;
    }

    static final class MappingResult {
        final int[] mappedText;
        final Translator translator;

        MappingResult(int[] mappedText, Translator translator) {
            this.mappedText = mappedText;
            this.translator = translator;
        }
    }

    static final class Translator {
        private final HashFunction hash;
        private final LongIntHashTable hashedTable;
        private final HashMap<String, Integer> directMap;

        private Translator(HashFunction hash, LongIntHashTable hashedTable, HashMap<String, Integer> directMap) {
            this.hash = hash;
            this.hashedTable = hashedTable;
            this.directMap = directMap;
        }

        static Translator fromHash(HashFunction hash, LongIntHashTable table) {
            return new Translator(hash, table, null);
        }

        static Translator fromDirect(HashMap<String, Integer> direct) {
            return new Translator(null, null, direct);
        }

        int[] mapPattern(String[] pattern) {
            int[] mapped = new int[pattern.length];
            if (directMap != null) {
                for (int i = 0; i < pattern.length; i++) {
                    Integer rank = directMap.get(pattern[i]);
                    if (rank == null) {
                        return null;
                    }
                    mapped[i] = rank;
                }
            } else {
                for (int i = 0; i < pattern.length; i++) {
                    long key = hash.hash(pattern[i]);
                    int rank = hashedTable.get(key);
                    if (rank < 0) {
                        return null;
                    }
                    mapped[i] = rank;
                }
            }
            return mapped;
        }
    }

    private static final class HashFunction {
        private static final long GOLDEN_GAMMA = 0x9e3779b97f4a7c15L;

        private final long mul;
        private final long add;
        private final UniversalStringHasher stringHasher;

        HashFunction(long mul, long add, UniversalStringHasher stringHasher) {
            if ((mul & 1L) == 0) {
                mul |= 1L; // enforce odd multiplier for better mixing
            }
            this.mul = mul;
            this.add = add;
            this.stringHasher = stringHasher;
        }

        static HashFunction random(SplittableRandom random) {
            long m = random.nextLong() | 1L;
            long a = random.nextLong();
            UniversalStringHasher hasher = UniversalStringHasher.withSeed(random.nextLong());
            return new HashFunction(m, a, hasher);
        }

        long hash(String value) {
            long base = stringHasher.hash(value);
            long z = base + add;
            z ^= mul;
            z = mix64(z);
            return z & Long.MAX_VALUE;
        }

        private static long mix64(long z) {
            z = (z ^ (z >>> 33)) * 0xff51afd7ed558ccdL;
            z = (z ^ (z >>> 33)) * 0xc4ceb9fe1a85ec53L;
            z ^= z >>> 33;
            z += GOLDEN_GAMMA;
            return z;
        }
    }

    private static final class LongHashSet {
        private static final long EMPTY = Long.MIN_VALUE;

        private final long[] keys;
        private final int mask;

        LongHashSet(int expected) {
            int capacity = 1;
            while (capacity < expected) {
                capacity <<= 1;
            }
            keys = new long[capacity];
            Arrays.fill(keys, EMPTY);
            mask = capacity - 1;
        }

        boolean add(long key) {
            int idx = mix(key) & mask;
            while (true) {
                long current = keys[idx];
                if (current == EMPTY) {
                    keys[idx] = key;
                    return true;
                }
                if (current == key) {
                    return false;
                }
                idx = (idx + 1) & mask;
            }
        }

        private static int mix(long key) {
            long z = key ^ (key >>> 33);
            z *= 0xff51afd7ed558ccdL;
            z ^= z >>> 33;
            return (int) z;
        }
    }

    private static final class LongIntHashTable {
        private static final long EMPTY = Long.MIN_VALUE;

        private final long[] keys;
        private final int[] values;
        private final int mask;

        LongIntHashTable(int expected) {
            int capacity = 1;
            while (capacity < expected * 2) {
                capacity <<= 1;
            }
            keys = new long[capacity];
            Arrays.fill(keys, EMPTY);
            values = new int[capacity];
            mask = capacity - 1;
        }

        void put(long key, int value) {
            int idx = mix(key) & mask;
            while (true) {
                long current = keys[idx];
                if (current == EMPTY) {
                    keys[idx] = key;
                    values[idx] = value;
                    return;
                }
                if (current == key) {
                    values[idx] = value;
                    return;
                }
                idx = (idx + 1) & mask;
            }
        }

        int get(long key) {
            int idx = mix(key) & mask;
            while (true) {
                long current = keys[idx];
                if (current == EMPTY) {
                    return -1;
                }
                if (current == key) {
                    return values[idx];
                }
                idx = (idx + 1) & mask;
            }
        }

        private static int mix(long key) {
            long z = key ^ (key >>> 33);
            z *= 0xff51afd7ed558ccdL;
            z ^= z >>> 33;
            return (int) z;
        }
    }
}
