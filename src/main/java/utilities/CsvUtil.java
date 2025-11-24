package utilities;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.Objects;

    public class CsvUtil {
    private CsvUtil() {}

    // Configuration for CSV formatting.
    public static final class Config {
        private final char delimiter;
        private final char quote;
        private final String lineSeparator;
        private final boolean alwaysQuote;
        private final Charset charset;

        private Config(char delimiter, char quote, String lineSeparator, boolean alwaysQuote, Charset charset) {
            this.delimiter = delimiter;
            this.quote = quote;
            this.lineSeparator = Objects.requireNonNull(lineSeparator, "lineSeparator");
            this.alwaysQuote = alwaysQuote;
            this.charset = Objects.requireNonNull(charset, "charset");
        }

        // Default configuration: comma, double quote, platform line separator, no forced quoting, UTF-8.
        public static Config defaults() {
            return new Config(',', '"', System.lineSeparator(), false, StandardCharsets.UTF_8);
        }

        public Config withDelimiter(char delimiter)     { return new Config(delimiter, this.quote, this.lineSeparator, this.alwaysQuote, this.charset); }
        public Config withQuote(char quote)             { return new Config(this.delimiter, quote, this.lineSeparator, this.alwaysQuote, this.charset); }
        public Config withLineSeparator(String ls)      { return new Config(this.delimiter, this.quote, ls, this.alwaysQuote, this.charset); }
        public Config withAlwaysQuote(boolean aq)       { return new Config(this.delimiter, this.quote, this.lineSeparator, aq, this.charset); }
        public Config withCharset(Charset charset)      { return new Config(this.delimiter, this.quote, this.lineSeparator, this.alwaysQuote, charset); }
    }

    // Public convenience methods for writing CSV.

    // Write a single row to a file.
    public static void writeRow(Path file, List<?> row) throws IOException {
        writeRows(file, List.of(row), Config.defaults());
    }

    // Write a single row to a file with custom configuration.
    public static void writeRow(Path file, List<?> row, Config config) throws IOException {
        writeRows(file, List.of(row), config);
    }

    // Write multiple rows (list of lists) to a file.
    public static void writeRows(Path file, List<? extends List<?>> rows) throws IOException {
        writeRows(file, rows, Config.defaults());
    }

    // Write multiple rows (list of lists) to a file with custom configuration.
    public static void writeRows(Path file, List<? extends List<?>> rows, Config config) throws IOException {
        try (BufferedWriter bw = Files.newBufferedWriter(file, config.charset)) {
            writeRows(bw, rows, config);
        }
    }

    // Write multiple rows to an OutputStream.
    public static void writeRows(OutputStream out, List<? extends List<?>> rows, Config config) throws IOException {
        try (Writer w = new BufferedWriter(new OutputStreamWriter(out, config.charset))) {
            writeRows(w, rows, config);
        }
    }

    // Write a single row to an existing Writer (caller closes).
    public static void writeRow(Writer writer, List<?> row, Config config) throws IOException {
        writeRows(writer, List.of(row), config);
    }

    // Write multiple rows to an existing Writer (caller closes).
    public static void writeRows(Writer writer, List<? extends List<?>> rows, Config config) throws IOException {
        Objects.requireNonNull(writer, "writer");
        Objects.requireNonNull(rows, "rows");
        Objects.requireNonNull(config, "config");

        for (List<?> row : rows) {
            writeOneRow(writer, row, config);
            writer.write(config.lineSeparator);
        }
        writer.flush();
    }

    // Build a CSV string from a single row.
    public static String toCsvLine(List<?> row) {
        return toCsvLine(row, Config.defaults());
    }

    // Build a CSV string from a single row with custom configuration.
    public static String toCsvLine(List<?> row, Config config) {
        StringBuilder sb = new StringBuilder();
        writeOneRowToBuilder(sb, row, config);
        return sb.toString();
    }

    // Build a CSV string for all rows.
    public static String toCsv(List<? extends List<?>> rows) {
        return toCsv(rows, Config.defaults());
    }

    // Build a CSV string for all rows with custom configuration.
    public static String toCsv(List<? extends List<?>> rows, Config config) {
        StringBuilder sb = new StringBuilder();
        try {
            for (Iterator<? extends List<?>> it = rows.iterator(); it.hasNext(); ) {
                writeOneRowToBuilder(sb, it.next(), config);
                if (it.hasNext()) sb.append(config.lineSeparator);
            }
        } catch (UncheckedIOException e) {
            // Should not happen when using StringBuilder, but kept for symmetry.
            throw e;
        }
        return sb.toString();
    }

    // Varargs convenience: write rows given as separate lists.
    @SafeVarargs
    public static void writeRows(Path file, Config config, List<?>... rows) throws IOException {
        writeRows(file, List.of(rows), config);
    }

    // Varargs convenience: write rows with default configuration.
    @SafeVarargs
    public static void writeRows(Path file, List<?>... rows) throws IOException {
        writeRows(file, List.of(rows), Config.defaults());
    }

    // Internal helpers.

    private static void writeOneRow(Writer writer, List<?> row, Config config) throws IOException {
        Objects.requireNonNull(row, "row");
        for (int i = 0; i < row.size(); i++) {
            if (i > 0) writer.write(config.delimiter);
            writer.write(quoteIfNeeded(stringify(row.get(i)), config));
        }
    }

    private static void writeOneRowToBuilder(StringBuilder sb, List<?> row, Config config) {
        Objects.requireNonNull(row, "row");
        for (int i = 0; i < row.size(); i++) {
            if (i > 0) sb.append(config.delimiter);
            sb.append(quoteIfNeeded(stringify(row.get(i)), config));
        }
    }

    // Convert any object to a field string; null becomes empty.
    private static String stringify(Object value) {
        if (value == null) return "";
        // Special case Character to avoid printing numeric code points accidentally
        if (value instanceof Character c) return String.valueOf(c);
        return String.valueOf(value);
    }

    // Apply CSV quoting rules and double inner quote characters when quoting.
    private static String quoteIfNeeded(String field, Config config) {
        boolean containsDelimiter = field.indexOf(config.delimiter) >= 0;
        boolean containsQuote = field.indexOf(config.quote) >= 0;
        boolean containsNewline = field.indexOf('\n') >= 0 || field.indexOf('\r') >= 0;

        if (config.alwaysQuote || containsDelimiter || containsQuote || containsNewline) {
            String doubled = field.replace(String.valueOf(config.quote), String.valueOf(config.quote) + config.quote);
            return config.quote + doubled + config.quote;
        }
        return field;
    }
}
