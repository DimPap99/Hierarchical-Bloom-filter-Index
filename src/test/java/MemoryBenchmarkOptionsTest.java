import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import utilities.SswsiSettings;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

import static org.junit.jupiter.api.Assertions.*;

// MemoryBenchmark takes the same three SSWSI flags as the harness and names its CSV after the setting.
class MemoryBenchmarkOptionsTest {

    @Test
    void sswsiFlags(@TempDir Path dir) throws IOException {
        Files.createDirectories(dir.resolve("w10"));
        String[] base = {"--data-root", dir.toString(), "--window", "w10", "--window-power", "10", "--tree-power", "10",
                "--policy", "NONE"};

        MemoryBenchmark.MemoryOptions defaults = MemoryBenchmark.MemoryOptions.parse(base);
        assertEquals(SswsiSettings.DEFAULT, defaults.sswsi());
        assertTrue(defaults.csvOutputFile().toString().endsWith("_sswsi-amortized-dc3-fks.csv"),
                defaults.csvOutputFile().toString());

        MemoryBenchmark.MemoryOptions fc = MemoryBenchmark.MemoryOptions.parse(concat(base,
                "--sswsi-builder", "fc", "--sswsi-dict", "fastutil", "--suffix-match-ngram", "true"));
        assertEquals("amortized-fc-fastutil", fc.sswsi().label());
        assertTrue(fc.csvOutputFile().toString().endsWith("_suffix-match_sswsi-amortized-fc-fastutil.csv"),
                fc.csvOutputFile().toString());

        MemoryBenchmark.MemoryOptions legacy = MemoryBenchmark.MemoryOptions.parse(concat(base, "--sswsi-mode", "legacy"));
        assertEquals(SswsiSettings.LEGACY, legacy.sswsi());
        assertTrue(legacy.csvOutputFile().toString().endsWith("_sswsi-legacy.csv"), legacy.csvOutputFile().toString());

        assertThrows(IllegalArgumentException.class,
                () -> MemoryBenchmark.MemoryOptions.parse(concat(base, "--sswsi-mode", "legacy", "--sswsi-builder", "dc3")));
        assertThrows(IllegalArgumentException.class, () -> MemoryBenchmark.MemoryOptions.parse(concat(base,
                "--window-length", "1000", "--tree-length", "1000")));
        assertDoesNotThrow(() -> MemoryBenchmark.MemoryOptions.parse(concat(base,
                "--window-length", "1000", "--tree-length", "1000", "--sswsi-mode", "legacy")));
    }

    @Test
    void minPatternLengthFlag(@TempDir Path dir) throws IOException {
        Files.createDirectories(dir.resolve("w10"));
        String[] base = {"--data-root", dir.toString(), "--window", "w10", "--window-power", "10", "--tree-power", "10",
                "--policy", "NONE"};
        assertEquals(0, MemoryBenchmark.MemoryOptions.parse(base).minPatternLength());
        MemoryBenchmark.MemoryOptions o = MemoryBenchmark.MemoryOptions.parse(concat(base, "--min-pattern-length", "10"));
        assertEquals(10, o.minPatternLength());
        assertTrue(o.csvOutputFile().toString().endsWith("_minpat10_sswsi-amortized-dc3-fks.csv"),
                o.csvOutputFile().toString());
        assertThrows(IllegalArgumentException.class,
                () -> MemoryBenchmark.MemoryOptions.parse(concat(base, "--min-pattern-length", "-1")));
        MemoryBenchmark.MemoryOptions blocked = MemoryBenchmark.MemoryOptions.parse(concat(base, "--filter", "blocked"));
        assertEquals(membership.FilterKind.BLOCKED, blocked.filter());
        assertTrue(blocked.csvOutputFile().toString().endsWith("_blocked_sswsi-amortized-dc3-fks.csv"),
                blocked.csvOutputFile().toString());
    }

    private static String[] concat(String[] a, String... b) {
        String[] out = new String[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }
}
