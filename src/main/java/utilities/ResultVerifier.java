package utilities;

import PMIndex.IPMIndexing;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// Cross-checks the indexes of a benchmark run (--verify true): every index is exact, so for the same
// query all of them must report the same occurrences, and any difference is a bug.
//
// The first index to answer a query records a fingerprint of its sorted positions; every later answer
// to that query, from any index or run, is compared with it. Positions are window-relative for HBI and
// SSWSI, and the online suffix tree indexes exactly one window, so they are directly comparable.
// The query loops keep their answers until the timer stops and only then call check(), so verifying
// adds no work inside the timed region.
public final class ResultVerifier {

    private static final int MAX_EXAMPLES = 10;

    private static volatile boolean enabled = false;
    private static final Map<String, Expected> EXPECTED = new HashMap<>();
    private static final List<String> EXAMPLES = new ArrayList<>();
    private static long checked = 0;
    private static long mismatches = 0;

    private ResultVerifier() {
    }

    public static synchronized void reset(boolean on) {
        enabled = on;
        EXPECTED.clear();
        EXAMPLES.clear();
        checked = 0;
        mismatches = 0;
    }

    public static boolean enabled() {
        return enabled;
    }

    public static synchronized void check(IPMIndexing index, int nGram, Path queryFile, List<? extends List<Integer>> answers) {
        String who = index.getClass().getSimpleName() + " n=" + nGram;
        for (int q = 0; q < answers.size(); q++) {
            long[] positions = answers.get(q).stream().mapToLong(Integer::longValue).sorted().toArray();
            long fingerprint = positions.length * 0x9E3779B97F4A7C15L ^ Arrays.hashCode(positions);
            String key = queryFile + "#" + q;
            Expected first = EXPECTED.get(key);
            if (first == null) {
                EXPECTED.put(key, new Expected(fingerprint, positions.length, who));
                continue;
            }
            checked++;
            if (first.fingerprint != fingerprint) {
                mismatches++;
                if (EXAMPLES.size() < MAX_EXAMPLES) {
                    EXAMPLES.add(String.format("%s query %d: %s found %d occurrences, %s found %d",
                            queryFile, q + 1, first.who, first.count, who, positions.length));
                }
            }
        }
    }

    public static synchronized String summary() {
        if (!enabled) {
            return "Result verification off (--verify true turns it on).";
        }
        StringBuilder sb = new StringBuilder();
        sb.append(String.format("Result verification: %d answers compared across indexes and runs, %d mismatches.",
                checked, mismatches));
        if (checked == 0) {
            sb.append(" Nothing to compare: run more than one index, or more than one run.");
        }
        for (String example : EXAMPLES) {
            sb.append(System.lineSeparator()).append("  MISMATCH ").append(example);
        }
        if (mismatches > EXAMPLES.size()) {
            sb.append(System.lineSeparator()).append("  ... and ").append(mismatches - EXAMPLES.size()).append(" more");
        }
        return sb.toString();
    }

    public static synchronized long mismatches() {
        return mismatches;
    }

    private record Expected(long fingerprint, int count, String who) {
    }
}
