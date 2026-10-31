package tree.ssws;

import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

class RankReducerTest {

    // Names must be 1..sigma in order of first occurrence, followed by the terminator 0, on every path:
    // direct comparison, 8-, 11- and 16-bit digits, and skipped passes.
    @Test
    void namesFollowFirstOccurrence() {
        Random rnd = new Random(5);
        int[] sizes = {0, 1, 2, 3, 31, 32, 33, 100,
                RankReducer.MEDIUM - 1, RankReducer.MEDIUM, RankReducer.MEDIUM + 1, 4096,
                RankReducer.LARGE - 1, RankReducer.LARGE, RankReducer.LARGE + 1, 200_000};
        for (int n : sizes) {
            for (int kind = 0; kind < 6; kind++) {
                long[] pool = new long[(kind == 2) ? 3 : 100];
                for (int i = 0; i < pool.length; i++) {
                    pool[i] = rnd.nextLong() & Long.MAX_VALUE;
                }
                long[] text = new long[n];
                for (int i = 0; i < n; i++) {
                    text[i] = switch (kind) {
                        case 0 -> rnd.nextLong() & Long.MAX_VALUE;          // all distinct, 63-bit hashes
                        case 1 -> 42L;                                      // one symbol
                        case 2, 3 -> pool[rnd.nextInt(pool.length)];        // few distinct
                        case 4 -> rnd.nextInt(1000);                        // only the low digits vary
                        default -> rnd.nextLong();                          // negative values too
                    };
                }
                assertArrayEquals(reference(text), RankReducer.rankWithTerminator(text),
                        "n=" + n + " kind=" + kind);
            }
        }
    }

    private static int[] reference(long[] text) {
        HashMap<Long, Integer> name = new HashMap<>();
        int[] out = new int[text.length + 1];
        for (int i = 0; i < text.length; i++) {
            out[i] = name.computeIfAbsent(text[i], v -> name.size() + 1);
        }
        return out;
    }
}
