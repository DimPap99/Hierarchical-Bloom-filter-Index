package utilities;

import PMIndex.RegexIndex;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import utilities.BenchmarkEnums.QueryType;
import utilities.Workload.ReaderMode;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class WorkloadTest {

    @TempDir
    Path root;

    private Path write(String relative, String content) throws IOException {
        Path file = root.resolve(relative);
        Files.createDirectories(file.getParent());
        Files.writeString(file, content, StandardCharsets.UTF_8);
        return file;
    }

    @Test
    void charsDropLineBreaksAndKeepQuerySpaces() throws IOException {
        Path data = write("data/text.txt", "ab\ncd\r\nef");
        Path queries = write("queries/text/2.txt", " bc\n\ncd \n");
        Workload w = new Workload(null, data, queries, ReaderMode.CHARS);

        assertEquals("abcdef", new String(w.readChars()));
        assertEquals(List.of(" bc", "cd "), w.queries());
        assertEquals("text.txt", w.name());
    }

    @Test
    void segmentsSkipEmptyLinesAndSplitQueriesOnRunsOfSpaces() throws IOException {
        Path data = write("data/tok.txt", "a,b\n\nc,d\r\ne,f\n");
        Path queries = write("queries/tok/2.txt", "  a,b   c,d \n");
        Workload w = new Workload("tok", data, queries, ReaderMode.SEGMENTS);

        assertEquals(List.of("a,b", "c,d", "e,f"), w.readTokens());
        List<Workload.Query> parsed = w.parsedQueries();
        assertEquals(1, parsed.size());
        assertArrayEquals(new String[]{"a,b", "c,d"}, parsed.get(0).tokens());
    }

    @Test
    void keyMatchesHowPatternsBuildTheirNgrams() {
        RingBuffer<Object> window = new RingBuffer<>(2);
        window.append("a,b");
        window.append("c,d");
        Workload segments = new Workload("s", Path.of("x"), null, ReaderMode.SEGMENTS);
        String key = segments.key(window);
        assertEquals(new search.Pattern(new String[]{"a,b", "c,d"}, 2).nGramArr[0], key);

        RingBuffer<Object> chars = new RingBuffer<>(2);
        chars.append('x');
        chars.append('y');
        assertEquals("xy", new Workload("c", Path.of("x"), null, ReaderMode.CHARS).key(chars));
    }

    @Test
    void readerModeParsesOldSpellings() {
        assertEquals(ReaderMode.CHARS, ReaderMode.parse("chars"));
        assertEquals(ReaderMode.SEGMENTS, ReaderMode.parse("Segments"));
        assertEquals(ReaderMode.SEGMENTS, ReaderMode.parse("segments2"));
        assertThrows(IllegalArgumentException.class, () -> ReaderMode.parse("bytes"));
    }

    @Test
    void catalogResolvesFlatAndSampledLayouts() throws IOException {
        write("data/pg.txt", "hello world");
        write("queries/pg/5.txt", "hello\n");
        write("data/w5/2/2_W5.txt", "abcabc");
        write("queries/w5/2/3.rare.txt", "abc\n");

        DataCatalog catalog = DataCatalog.from(root.resolve("data"), null);
        assertEquals(root.resolve("queries").toAbsolutePath().normalize(), catalog.queryRoot());

        Workload flat = catalog.byName("pg", 1, 5, QueryType.UNIFORM, ReaderMode.CHARS);
        assertEquals(root.resolve("data/pg.txt"), flat.dataFile());
        assertEquals(root.resolve("queries/pg/5.txt"), flat.queriesFile());

        Workload sampled = catalog.byName("w5", 2, 3, QueryType.RARE, ReaderMode.CHARS);
        assertEquals(root.resolve("data/w5/2/2_W5.txt"), sampled.dataFile());
        assertEquals(root.resolve("queries/w5/2/3.rare.txt"), sampled.queriesFile());

        assertEquals(List.of(2), catalog.parts("w5"));
        assertEquals(List.of(3), catalog.patternLengths("w5", 2, QueryType.RARE));

        Workload files = catalog.files("pg.txt", "pg/5.txt", ReaderMode.CHARS);
        assertEquals(flat.dataFile(), files.dataFile());

        IllegalArgumentException missing = assertThrows(IllegalArgumentException.class,
                () -> catalog.byName("w5", 2, 40, QueryType.RARE, ReaderMode.CHARS));
        assertTrue(missing.getMessage().contains("40.rare.txt"), missing.getMessage());
    }

    @Test
    void selectionReadsFlagsAndInfersSegmentsForMawi() throws IOException {
        write("data/mawi19/1/1_W19.txt", "a\nb\n");
        write("queries/mawi19/1/10.uniform.txt", "a b\n");

        DataCatalog.Selection sel = new DataCatalog.Selection("pg2701");
        assertTrue(sel.accept("dataset", "mawi19"));
        assertTrue(sel.accept("data-root", root.resolve("data").toString()));
        assertFalse(sel.accept("fp", "0.1"));

        assertEquals(ReaderMode.SEGMENTS, sel.mode());
        assertEquals(1 << 19, sel.windowLength(1 << 21));
        Workload w = sel.resolve();
        assertEquals(root.resolve("data/mawi19/1/1_W19.txt").toAbsolutePath().normalize(), w.dataFile());
        assertEquals(ReaderMode.SEGMENTS, w.mode());
    }

    @Test
    void runFindsEveryOccurrence() throws IOException {
        Path data = write("data/t.txt", "abracadabra");
        Path queries = write("queries/t/2.txt", "abra\ncad\nzzz\n");
        Workload w = new Workload(null, data, queries, ReaderMode.CHARS);

        ExperimentRunResult r = WorkloadRunner.run(w, new RegexIndex(), 1,
                WorkloadRunner.Options.QUIET.collectMatches(true));

        assertEquals(List.of(List.of(0, 7), List.of(4), List.of()), r.matchRes());
    }

    @Test
    void streamQueriesAfterEveryWindow() throws IOException {
        Path data = write("data/t.txt", "abcdefgh");
        Path queries = write("queries/t/1.txt", "a\n");
        Workload w = new Workload(null, data, queries, ReaderMode.CHARS);

        ExperimentRunResult r = WorkloadRunner.stream(w, new RegexIndex(), 1, 4,
                WorkloadRunner.Options.QUIET.collectMatches(true));

        // 8 symbols with a window of 4: two query loads of one query each
        ArrayList<ArrayList<Integer>> loads = r.matchRes();
        assertEquals(2, loads.size());
        assertEquals(List.of(0), loads.get(0));
    }
}
