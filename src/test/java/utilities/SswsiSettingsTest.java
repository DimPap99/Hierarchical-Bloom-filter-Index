package utilities;

import PMIndex.IPMIndexing;
import PMIndex.SlidingSuffixIndex;
import PMIndex.StreamingSlidingWindowIndex;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import tree.ssws.ChildDictionary;
import tree.ssws.Dc3SuffixTreeBuilder;
import tree.ssws.FarachColtonSuffixTreeBuilder;
import tree.ssws.SuffixTreeBuilder;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// --sswsi-mode, --sswsi-builder and --sswsi-dict: parsing, the factory, and verified harness runs.
class SswsiSettingsTest {

    private static final List<String[]> EVERY_SETTING = List.of(
            new String[]{"--sswsi-mode", "legacy"},
            new String[]{"--sswsi-mode", "amortized", "--sswsi-builder", "dc3", "--sswsi-dict", "fks"},
            new String[]{"--sswsi-mode", "amortized", "--sswsi-builder", "dc3", "--sswsi-dict", "fastutil"},
            new String[]{"--sswsi-mode", "amortized", "--sswsi-builder", "fc", "--sswsi-dict", "fks"},
            new String[]{"--sswsi-mode", "amortized", "--sswsi-builder", "fc", "--sswsi-dict", "fastutil"});

    @Test
    void defaultIsAmortizedWithDc3AndFks() {
        SswsiSettings s = SswsiSettings.parse(null, null, null);
        assertEquals(SswsiSettings.DEFAULT, s);
        assertEquals(SswsiSettings.Mode.AMORTIZED, s.mode());
        assertEquals(SuffixTreeBuilder.Kind.DC3, s.builder());
        assertEquals(ChildDictionary.Kind.FKS, s.dictionary());
        assertEquals("amortized-dc3-fks", s.label());
        // A builder or a dictionary alone keeps the other default.
        assertEquals("amortized-fc-fks", SswsiSettings.parse(null, "fc", null).label());
        assertEquals("amortized-dc3-fastutil", SswsiSettings.parse(null, null, "fastutil").label());
    }

    @Test
    void parsesEveryValueIgnoringCase() {
        assertEquals(SswsiSettings.LEGACY, SswsiSettings.parse("legacy", null, null));
        assertEquals("legacy", SswsiSettings.LEGACY.label());
        SswsiSettings s = SswsiSettings.parse("AMORTIZED", "FC", "FastUtil");
        assertEquals(new SswsiSettings(SswsiSettings.Mode.AMORTIZED, SuffixTreeBuilder.Kind.FC,
                ChildDictionary.Kind.FASTUTIL), s);
        assertEquals("amortized-fc-fastutil", s.label());
    }

    @Test
    void rejectsUnknownNamesAndChoicesLegacyDoesNotHave() {
        assertThrows(IllegalArgumentException.class, () -> SswsiSettings.parse("deamortized", null, null));
        assertThrows(IllegalArgumentException.class, () -> SswsiSettings.parse(null, "ukkonen", null));
        assertThrows(IllegalArgumentException.class, () -> SswsiSettings.parse(null, null, "hashmap"));
        assertThrows(IllegalArgumentException.class, () -> SswsiSettings.parse("legacy", "fc", null));
        assertThrows(IllegalArgumentException.class, () -> SswsiSettings.parse("legacy", null, "fks"));
    }

    @Test
    void amortizedNeedsAPowerOfTwoWindow() {
        assertDoesNotThrow(() -> SswsiSettings.DEFAULT.requireValidWindow(1 << 21));
        assertThrows(IllegalArgumentException.class, () -> SswsiSettings.DEFAULT.requireValidWindow(3 << 19));
        assertDoesNotThrow(() -> SswsiSettings.LEGACY.requireValidWindow(3 << 19));
    }

    @Test
    void harnessOptions(@TempDir Path dir) throws IOException {
        Files.createDirectories(dir.resolve("data/w10"));
        Files.createDirectories(dir.resolve("queries/w10"));
        String[] base = {"--data-root", dir.resolve("data").toString(), "--query-root", dir.resolve("queries").toString(),
                "--window", "w10"};

        assertEquals(SswsiSettings.DEFAULT, MultiBenchmarkOptions.parse(base).sswsi());
        assertEquals(Path.of(""), MultiBenchmarkOptions.parse(base).outDir());
        MultiBenchmarkOptions o = MultiBenchmarkOptions.parse(concat(base,
                "--sswsi-mode", "amortized", "--sswsi-builder", "fc", "--sswsi-dict=fastutil", "--out-dir", dir.toString()));
        assertEquals("amortized-fc-fastutil", o.sswsi().label());
        assertEquals(dir, o.outDir());
        assertEquals(SswsiSettings.LEGACY, MultiBenchmarkOptions.parse(concat(base, "--sswsi-mode", "legacy")).sswsi());
        assertEquals(0, MultiBenchmarkOptions.parse(base).minPatternLength());
        assertEquals(10, MultiBenchmarkOptions.parse(concat(base, "--min-pattern-length", "10")).minPatternLength());
        assertThrows(IllegalArgumentException.class,
                () -> MultiBenchmarkOptions.parse(concat(base, "--min-pattern-length", "-1")));

        // A window the amortized index cannot take fails at parse time, unless SSWSI does not run.
        assertThrows(IllegalArgumentException.class,
                () -> MultiBenchmarkOptions.parse(concat(base, "--window-length", "1000", "--tree-length", "1000")));
        assertDoesNotThrow(() -> MultiBenchmarkOptions.parse(concat(base, "--window-length", "1000", "--tree-length",
                "1000", "--sswsi-mode", "legacy")));
        assertDoesNotThrow(() -> MultiBenchmarkOptions.parse(concat(base, "--window-length", "1000", "--tree-length",
                "1000", "--run-suffix", "false")));
        assertThrows(IllegalArgumentException.class,
                () -> MultiBenchmarkOptions.parse(concat(base, "--sswsi-mode", "legacy", "--sswsi-dict", "fks")));
        assertThrows(IllegalArgumentException.class,
                () -> MultiBenchmarkOptions.parse(concat(base, "--out-dir", dir.resolve("missing").toString())));
    }

    @Test
    void factoryBuildsTheSelectedIndex() {
        assertInstanceOf(StreamingSlidingWindowIndex.class, IndexFactory.createSuffixIndex(1 << 10, SswsiSettings.LEGACY));
        for (SuffixTreeBuilder.Kind builder : SuffixTreeBuilder.Kind.values()) {
            for (ChildDictionary.Kind dictionary : ChildDictionary.Kind.values()) {
                SswsiSettings s = new SswsiSettings(SswsiSettings.Mode.AMORTIZED, builder, dictionary);
                IPMIndexing index = IndexFactory.createSuffixIndex(1 << 10, s);
                SlidingSuffixIndex amortized = assertInstanceOf(SlidingSuffixIndex.class, index);
                Class<? extends SuffixTreeBuilder> expected = (builder == SuffixTreeBuilder.Kind.FC)
                        ? FarachColtonSuffixTreeBuilder.class : Dc3SuffixTreeBuilder.class;
                assertEquals(expected, amortized.builder().getClass());
                assertEquals(dictionary, amortized.dictionary());
            }
        }
    }

    // The harness end to end with --verify, for every setting: HBI, SSWSI and the online suffix tree must report
    // the same occurrences for every query, with SSWSI at n = 1 against HBI at n = 2, and all three at n = 2.
    @Test
    void verifiedHarnessRunsAgreeForEverySetting(@TempDir Path dir) throws Exception {
        Random random = new Random(7);
        int w = 1 << 10;
        StringBuilder text = new StringBuilder();
        for (int i = 0; i < w; i++) text.append((char) ('a' + random.nextInt(4)));
        writeWindow(dir, "chars", text.toString(), random, false);
        List<String> tokens = new ArrayList<>();
        for (int i = 0; i < w; i++) tokens.add("10.0.0." + random.nextInt(6) + ",6");
        writeWindow(dir, "segments", String.join("\n", tokens), random, true);

        for (String mode : new String[]{"chars", "segments"}) {
            for (String[] suffixNgram : new String[][]{{"--suffix-ngram", "1"}, {"--suffix-ngram", "match"}}) {
                for (String[] setting : EVERY_SETTING) {
                    Path out = Files.createDirectories(dir.resolve("out-" + mode + "-" + suffixNgram[1] + "-"
                            + String.join("", setting).replace("--", "_")));
                    MultiBenchmarkOptions options = MultiBenchmarkOptions.parse(concat(concat(new String[]{
                            "--data-root", dir.resolve("data").toString(), "--query-root", dir.resolve("queries").toString(),
                            "--window", mode, "--mode", mode, "--window-length", "1024", "--tree-length", "1024",
                            "--ngrams", "2", "--alphabet", "8", "--fp", "0.3", "--algorithm", "bs", "--policy", "NONE",
                            "--run-suffix-tree", "true", "--verify", "true", "--warmup", "0", "--runs", "1",
                            "--type", "all", "--out-dir", out.toString()}, suffixNgram), setting));
                    BenchmarkOrchestrator.run(options);
                    String what = mode + " " + String.join(" ", suffixNgram) + " " + String.join(" ", setting);
                    assertEquals(0, ResultVerifier.mismatches(), what + ": " + ResultVerifier.summary());
                    // Two query types x two lengths x 40 queries, each compared for two more indexes.
                    assertEquals(2 * 2 * 40 * 2, ResultVerifier.checked(), what);

                    Path csv;
                    try (var files = Files.list(out)) {
                        csv = files.findFirst().orElseThrow();
                    }
                    assertTrue(csv.getFileName().toString().endsWith("_sswsi-" + options.sswsi().label() + ".csv"),
                            csv.getFileName().toString());
                    List<String> lines = Files.readAllLines(csv);
                    assertTrue(lines.get(0).contains("sswsi_mode,sswsi_builder,sswsi_dict,sswsi_ngram,suffix_tree_ngram"),
                            lines.get(0));
                    int expectedSuffixNgram = suffixNgram[1].equals("match") ? 2 : 1;
                    String settingColumns = options.sswsi().mode() == SswsiSettings.Mode.LEGACY
                            ? "legacy,,," : "amortized," + setting[3] + "," + setting[5] + ",";
                    assertTrue(lines.get(1).contains(",bs," + settingColumns + expectedSuffixNgram + ","
                            + expectedSuffixNgram + ","), lines.get(1));
                }
            }
        }
    }

    // One window of w symbols under <dir>/data/<name>/1 and 40 queries per file of lengths 10 and 20 under
    // <dir>/queries/<name>/1: uniform are substrings of the window, missing are random and mostly absent.
    private static void writeWindow(Path dir, String name, String content, Random random, boolean segments)
            throws IOException {
        Path data = Files.createDirectories(dir.resolve("data").resolve(name).resolve("1"));
        Path queries = Files.createDirectories(dir.resolve("queries").resolve(name).resolve("1"));
        Files.writeString(data.resolve("1_W10.txt"), content);
        List<String> symbols = segments ? List.of(content.split("\n")) : List.of(content.split(""));
        for (int m : new int[]{10, 20}) {
            List<String> uniform = new ArrayList<>(), missing = new ArrayList<>();
            for (int q = 0; q < 40; q++) {
                int start = random.nextInt(symbols.size() - m);
                uniform.add(String.join(segments ? " " : "", symbols.subList(start, start + m)));
                List<String> absent = new ArrayList<>();
                for (int i = 0; i < m; i++) absent.add(symbols.get(random.nextInt(symbols.size())));
                missing.add(String.join(segments ? " " : "", absent));
            }
            Files.write(queries.resolve(m + ".uniform.txt"), uniform);
            Files.write(queries.resolve(m + ".missing.txt"), missing);
        }
    }

    private static String[] concat(String[] a, String... b) {
        String[] out = new String[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }
}
