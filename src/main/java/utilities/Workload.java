package utilities;

import search.Pattern;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;
import java.util.Objects;

// One test input: a dataset file, the query file run against it, and how both are read.
// Knows how its files are read, so runners never deal with paths or reader modes:
//   - readChars() / readTokens() load the dataset as the symbols the index inserts,
//   - key(window) turns the current n-gram window into the string the index stores,
//   - queries() / patterns(ngram) load the query file and build the patterns.
// Build one with DataCatalog, or directly with the record constructor. WorkloadRunner runs it.
public record Workload(String name, Path dataFile, Path queriesFile, ReaderMode mode) {

    // How the dataset and the queries are split into symbols.
    public enum ReaderMode {
        // Every character is a symbol. Line breaks in the dataset are dropped; a query is one line.
        CHARS,
        // Every dataset line is a symbol (a token). A query is one line of space-separated tokens.
        SEGMENTS;

        // Accepts "chars" and "segments" in any case; "segments2" and other "segments*" spellings
        // count as segments.
        public static ReaderMode parse(String value) {
            String v = Objects.requireNonNull(value, "reader mode").trim().toLowerCase(Locale.ROOT);
            if (v.equals("chars") || v.equals("char")) return CHARS;
            if (v.startsWith("segment")) return SEGMENTS;
            throw new IllegalArgumentException("Unknown reader mode '" + value + "' (use chars or segments)");
        }
    }

    public Workload {
        Objects.requireNonNull(dataFile, "dataFile");
        Objects.requireNonNull(mode, "mode");
        if (name == null || name.isBlank()) {
            name = dataFile.getFileName().toString();
        }
    }

    // A workload for two explicit files, named after the dataset file.
    public static Workload of(String dataFile, String queriesFile, ReaderMode mode) {
        return new Workload(null, Path.of(dataFile), queriesFile == null ? null : Path.of(queriesFile), mode);
    }

    // The same dataset with another query file, for running several query files on one index.
    public Workload withQueries(Path otherQueriesFile) {
        return new Workload(name, dataFile, otherQueriesFile, mode);
    }

    public boolean hasQueries() {
        return queriesFile != null;
    }

    // ---------------------------------------------------------------- data side

    // The dataset as CHARS symbols: every character of the file, line breaks dropped.
    public char[] readChars() {
        requireMode(ReaderMode.CHARS);
        try {
            String raw = Files.readString(dataFile, StandardCharsets.UTF_8);
            char[] text = new char[raw.length()];
            int n = 0;
            for (int i = 0; i < raw.length(); i++) {
                char c = raw.charAt(i);
                if (c != '\n' && c != '\r') {
                    text[n++] = c;
                }
            }
            return (n == text.length) ? text : Arrays.copyOf(text, n);
        } catch (IOException e) {
            throw new UncheckedIOException("Cannot read dataset " + dataFile, e);
        }
    }

    // The dataset as SEGMENTS symbols: one token per non-empty line.
    public List<String> readTokens() {
        requireMode(ReaderMode.SEGMENTS);
        try {
            return readNonEmptyLines(dataFile);
        } catch (IOException e) {
            throw new UncheckedIOException("Cannot read dataset " + dataFile, e);
        }
    }

    // The string the index stores for the n-gram currently in window. Characters concatenate directly;
    // tokens are joined with Pattern.SEGMENT_JOINER, the same way Pattern builds query n-grams.
    public String key(RingBuffer<?> window) {
        return (mode == ReaderMode.CHARS) ? window.snapshot() : window.join(Pattern.SEGMENT_JOINER);
    }

    // ---------------------------------------------------------------- query side

    // The queries, one per non-empty line, in file order. CHARS queries are kept exactly as written,
    // since leading and trailing spaces are part of a text pattern; SEGMENTS queries are stripped.
    public List<String> queries() {
        requireQueries();
        try {
            List<String> lines = readQueries(queriesFile);
            if (mode == ReaderMode.SEGMENTS) {
                List<String> stripped = new ArrayList<>(lines.size());
                for (String line : lines) {
                    String s = line.strip();
                    if (!s.isEmpty()) stripped.add(s);
                }
                return stripped;
            }
            return lines;
        } catch (IOException e) {
            throw new UncheckedIOException("Cannot read queries " + queriesFile, e);
        }
    }

    // One query, already split into tokens in SEGMENTS mode (tokens is null in CHARS mode), so that
    // building its Pattern costs no parsing inside a timed loop.
    public record Query(String text, String[] tokens) {
        public Pattern toPattern(int ngram) {
            return (tokens == null) ? new Pattern(text, ngram) : new Pattern(tokens, ngram);
        }
    }

    // The queries, parsed for this workload's mode.
    public List<Query> parsedQueries() {
        List<String> queries = queries();
        List<Query> out = new ArrayList<>(queries.size());
        for (String q : queries) {
            out.add(new Query(q, mode == ReaderMode.SEGMENTS ? splitTokens(q) : null));
        }
        return out;
    }

    // One Pattern per query. The index fills patterns in while it searches, so this builds fresh
    // objects on every call: call it again for every index and every query load.
    public List<Pattern> patterns(int ngram) {
        List<Query> queries = parsedQueries();
        List<Pattern> out = new ArrayList<>(queries.size());
        for (Query q : queries) {
            out.add(q.toPattern(ngram));
        }
        return out;
    }

    public Pattern toPattern(String query, int ngram) {
        return toPattern(query, ngram, mode);
    }

    public static Pattern toPattern(String query, int ngram, ReaderMode mode) {
        return (mode == ReaderMode.SEGMENTS) ? new Pattern(splitTokens(query), ngram) : new Pattern(query, ngram);
    }

    // A SEGMENTS query split into its tokens; runs of spaces count as one separator.
    public static String[] splitTokens(String query) {
        return query.strip().split(" +");
    }

    // Non-empty lines of a query file, exactly as written.
    public static List<String> readQueries(Path queryFile) throws IOException {
        return readNonEmptyLines(queryFile);
    }

    // Non-empty lines; the line break ("\n", "\r\n" or "\r") is not part of the line.
    private static List<String> readNonEmptyLines(Path file) throws IOException {
        List<String> out = new ArrayList<>();
        try (BufferedReader in = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String line;
            while ((line = in.readLine()) != null) {
                if (!line.isEmpty()) out.add(line);
            }
        }
        return out;
    }

    private void requireMode(ReaderMode expected) {
        if (mode != expected) {
            throw new IllegalStateException("Workload " + name + " is read as " + mode + ", not " + expected);
        }
    }

    private void requireQueries() {
        if (queriesFile == null) {
            throw new IllegalStateException("Workload " + name + " has no query file");
        }
    }

    @Override
    public String toString() {
        return name + " [" + mode.name().toLowerCase(Locale.ROOT) + "] data=" + dataFile
                + (queriesFile == null ? "" : " queries=" + queriesFile);
    }
}
