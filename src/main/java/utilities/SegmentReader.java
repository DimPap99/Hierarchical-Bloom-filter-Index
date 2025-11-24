package utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

// SegmentReader streams text split by a delimiter (e.g. '\n'), with optional query-mode switching.
public class SegmentReader implements Iterable<String>, AutoCloseable {

    private final String datasetPath;

    // ORIGINAL single-query file path (can be null if constructed with a list)
    private final String queriesFilePath;

    // NEW: optional list of query file paths (null if constructed with single file)
    private final List<String> queriesFilePaths;

    private String inputFilePath;
    private boolean opened;
    private boolean closed;

    private BufferedReader reader;

    private final char delimiter;
    private final boolean includeDelimiterInResult;
    private final boolean skipEmptySegments;

    // Constructor for one dataset file + one queries file (backwards compatible).
    public SegmentReader(String datasetPath,
                         String queriesFilePath,
                         char delimiter,
                         boolean includeDelimiterInResult,
                         boolean skipEmptySegments) {
        this.datasetPath = datasetPath;
        this.queriesFilePath = queriesFilePath;
        this.queriesFilePaths = null; // no list in this mode

        this.inputFilePath = datasetPath;
        this.delimiter = delimiter;
        this.includeDelimiterInResult = includeDelimiterInResult;
        this.skipEmptySegments = skipEmptySegments;
    }

    // Constructor for one dataset file plus multiple query files.
    public SegmentReader(String datasetPath,
                         List<String> queriesFilePaths,
                         char delimiter,
                         boolean includeDelimiterInResult,
                         boolean skipEmptySegments) {
        this.datasetPath = datasetPath;
        this.queriesFilePath = null; // not using single-path mode
        this.queriesFilePaths = List.copyOf(queriesFilePaths);

        this.inputFilePath = datasetPath;
        this.delimiter = delimiter;
        this.includeDelimiterInResult = includeDelimiterInResult;
        this.skipEmptySegments = skipEmptySegments;
    }

    // Switch to query mode; picks the first query file or the single queriesFilePath.
    public void setQueryMode() {
        if (queriesFilePaths != null && !queriesFilePaths.isEmpty()) {
            switchTo(queriesFilePaths.get(0));
        } else {
            switchTo(queriesFilePath);
        }
    }

    // Switch to a specific query file by index in the list constructor.
    public void setQueryMode(int idx) {
        if (queriesFilePaths != null) {
            if (idx < 0 || idx >= queriesFilePaths.size()) {
                throw new IndexOutOfBoundsException(
                        "Query file index " + idx + " out of range [0," + (queriesFilePaths.size() - 1) + "]");
            }
            switchTo(queriesFilePaths.get(idx));
        } else {
            // single-file mode, preserve old semantics
            if (idx != 0) {
                throw new IndexOutOfBoundsException(
                        "SegmentReader constructed with a single queries file; only index 0 is valid");
            }
            switchTo(queriesFilePath);
        }
    }

    // Go back to streaming the dataset.
    public void resetToDataset() {
        switchTo(datasetPath);
    }

    private void switchTo(String newPath) {
        try {
            close();
        } catch (Exception e) {
            throw new RuntimeException("Failed to close reader before switching sources", e);
        }
        this.inputFilePath = newPath;
        this.reader = null;
        this.opened = false;
        this.closed = false;
    }

    public void open() {
        if (opened) return;
        if (closed) throw new IllegalStateException("Already opened and closed the reader");
        try {
            this.reader = Files.newBufferedReader(Path.of(inputFilePath), StandardCharsets.UTF_8);
            this.opened = true;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() throws Exception {
        if (closed) return;
        closed = true;
        opened = false;
        if (reader != null) {
            try {
                reader.close();
            } finally {
                reader = null;
            }
        }
    }

    @Override
    public Iterator<String> iterator() {
        open();
        return new SegmentIterator(reader, delimiter, includeDelimiterInResult, skipEmptySegments, () -> {
            try {
                close();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    // =========================
    // Iterator implementation
    // =========================
    private static final class SegmentIterator implements Iterator<String> {

        private static final int EOF = -1;

        private final BufferedReader reader;
        private final char delimiter;
        private final boolean includeDelimiterInResult;
        private final boolean skipEmptySegments;
        private final Runnable onClose;

        private final char[] buf = new char[8192]; // chunk buffer for efficient reads
        private int bufLen = 0;                     // number of valid chars in buf
        private int bufPos = 0;                     // current read position in buf

        private final StringBuilder acc = new StringBuilder(8192);
        private String nextSegment = null;
        private boolean finished = false;

        SegmentIterator(BufferedReader reader,
                        char delimiter,
                        boolean includeDelimiterInResult,
                        boolean skipEmptySegments,
                        Runnable onClose) {
            this.reader = reader;
            this.delimiter = delimiter;
            this.includeDelimiterInResult = includeDelimiterInResult;
            this.skipEmptySegments = skipEmptySegments;
            this.onClose = onClose;
        }

        @Override
        public boolean hasNext() {
            if (finished) return false;
            if (nextSegment != null) return true;

            try {
                while (true) {
                    // Refill buffer if exhausted
                    if (bufPos >= bufLen) {
                        bufLen = reader.read(buf);
                        bufPos = 0;

                        if (bufLen == EOF) {
                            // End of file: flush any accumulated content as a final segment
                            if (acc.length() == 0) {
                                finish();
                                return false;
                            } else {
                                String seg = acc.toString();
                                acc.setLength(0);
                                // Handle trailing '\r' when delimiter == '\n'
                                if (delimiter == '\n' && endsWithCarriageReturn(seg)) {
                                    seg = seg.substring(0, seg.length() - 1);
                                }
                                if (skipEmptySegments && seg.isEmpty()) {
                                    // No more data anyway
                                    finish();
                                    return false;
                                }
                                nextSegment = seg;
                                return true;
                            }
                        }
                    }

                    // Consume chunk until we hit a delimiter or run out
                    while (bufPos < bufLen) {
                        char c = buf[bufPos++];
                        if (c == delimiter) {
                            // Emit segment (acc + optional delimiter)
                            String seg;
                            if (includeDelimiterInResult) {
                                acc.append(c);
                                seg = acc.toString();
                                acc.setLength(0);
                            } else {
                                seg = acc.toString();
                                acc.setLength(0);
                            }

                            // If delimiter is '\n' and the segment ends with '\r', strip it (CRLF)
                            if (delimiter == '\n' && endsWithCarriageReturn(seg)) {
                                seg = seg.substring(0, seg.length() - 1);
                            }

                            if (skipEmptySegments && seg.isEmpty()) {
                                // Continue scanning for the next non-empty segment
                                continue;
                            }
                            nextSegment = seg;
                            return true;
                        } else {
                            acc.append(c);
                        }
                    }
                    // Loop to refill buffer
                }
            } catch (IOException e) {
                finishWithIOException(e);
                throw new UncheckedIOException("Error reading segments", e);
            }
        }

        @Override
        public String next() {
            if (!hasNext()) throw new NoSuchElementException("No more segments");
            String seg = nextSegment;
            nextSegment = null;
            return seg;
        }

        private void finish() {
            if (!finished) {
                finished = true;
                try {
                    reader.close();
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                } finally {
                    onClose.run();
                }
            }
        }

        private void finishWithIOException(IOException e) {
            if (!finished) {
                finished = true;
                try {
                    reader.close();
                } catch (IOException suppressed) {
                    e.addSuppressed(suppressed);
                } finally {
                    onClose.run();
                }
            }
        }

        private static boolean endsWithCarriageReturn(String s) {
            return !s.isEmpty() && s.charAt(s.length() - 1) == '\r';
        }
    }
}
