package utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;

/**
 * SegmentReader is to String what DatasetReader is to Character:
 * it streams a file and yields segments split by a delimiter character (e.g., '\n').
 *
 * Features:
 *  - Works with large files using a BufferedReader and an internal char[] buffer.
 *  - Supports switching between dataset file and queries file (same pattern as DatasetReader).
 *  - Optionally includes the delimiter at the end of each returned segment.
 *  - Optionally skips empty segments (useful if consecutive delimiters should not produce empty strings).
 *  - Handles Windows CRLF ("\r\n") seamlessly when delimiter == '\n' (strips trailing '\r' from the segment).
 */
public class SegmentReader implements Iterable<String>, AutoCloseable {

    private final String datasetPath;
    private final String queriesFilePath;

    private String inputFilePath;
    private boolean opened;
    private boolean closed;

    private BufferedReader reader;

    private final char delimiter;
    private final boolean includeDelimiterInResult;
    private final boolean skipEmptySegments;

    // When setQueryMode() is called, we switch to queriesFilePath and assume newline-structured content,
    // but the delimiter stays whatever this instance was created with; if you want '\n', create with '\n'.
    public SegmentReader(String datasetPath,
                         String queriesFilePath,
                         char delimiter,
                         boolean includeDelimiterInResult,
                         boolean skipEmptySegments) {
        this.datasetPath = datasetPath;
        this.queriesFilePath = queriesFilePath;
        this.inputFilePath = datasetPath;
        this.delimiter = delimiter;
        this.includeDelimiterInResult = includeDelimiterInResult;
        this.skipEmptySegments = skipEmptySegments;
    }

    /**
     * Apart from datasets that are split by newline, query files are also in that format.
     * This method mirrors DatasetReader#setQueryMode() to make the switch explicit in experiment code.
     */
    public void setQueryMode() {
        switchTo(queriesFilePath);
    }

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
