package utilities;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class DatasetReader implements Iterable<Character>, AutoCloseable {

    boolean verbose;
    String inputFilePath;
    final String queriesFilePath;

    boolean lineBreakDataset;
    private BufferedReader reader;
    private boolean opened;
    private boolean closed;
    private final String datasetPath;
    private final boolean datasetLineBreak;

    public DatasetReader(String datasetPath, String queriesFilePath, boolean lineBreakDataset) {
        this.datasetPath = datasetPath;
        this.datasetLineBreak = lineBreakDataset;
        this.inputFilePath = datasetPath;
        this.queriesFilePath = queriesFilePath;
        this.lineBreakDataset = lineBreakDataset;
        this.verbose = false;
    }

    //Apart from datasets are split by \n i also keep the query files in this format.
    //This (ugly) method is used to make that explicit in the experiment code
    public void setQueryMode(){
        switchTo(queriesFilePath, true);
    }


    public void open(){
        if(opened){return;}
        if(closed){throw new IllegalStateException("Already opened the dataset buffer");}
        try{
            this.reader = Files.newBufferedReader(Path.of(inputFilePath), StandardCharsets.UTF_8);
            this.opened = true;
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public void close() throws Exception {
        if(closed){return;}
        closed = true;
        opened = false;
        if(reader != null){
            try {
                reader.close();
            } finally {
                reader = null;
            }
        }
    }

    public void resetToDataset(){
        switchTo(datasetPath, datasetLineBreak);
    }

    private void switchTo(String newPath, boolean includeLineBreaks) {
        try {
            close();
        } catch (Exception e) {
            throw new RuntimeException("Failed to close reader before switching sources", e);
        }
        this.inputFilePath = newPath;
        this.lineBreakDataset = includeLineBreaks;
        this.reader = null;
        this.opened = false;
        this.closed = false;
    }

    @Override
    public Iterator<Character> iterator() {

        this.open();
        return new CharIterator(reader, this.lineBreakDataset, () -> {
            // on iterator end, mark closed and null out the reader
            try {
                close();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
        });
    }

    private static final class CharIterator implements Iterator<Character> {

        private final BufferedReader reader;
        private final boolean includeLineBreaks;
        private final Runnable onClose; // callback to inform DatasetReader to close

        private int nextCodePoint = UNINITIALIZED;
        private static final int UNINITIALIZED = Integer.MIN_VALUE;
        private static final int EOF = -1;
        private boolean finished = false;

        CharIterator(BufferedReader reader, boolean includeLineBreaks, Runnable onClose) {
            this.reader = reader;
            this.includeLineBreaks = includeLineBreaks;
            this.onClose = onClose;
        }

        @Override
        public boolean hasNext() {
            if (finished) return false;
            if (nextCodePoint != UNINITIALIZED) return true;

            try {
                // Fetch next character
                while (true) {
                    int ch = reader.read(); // returns a UTF-16 code unit as int, or -1 at EOF
                    if (ch == EOF) {
                        finish();
                        return false;
                    }

                    char c = (char) ch;

                    if (!includeLineBreaks && (c == '\n' || c == '\r')) {
                        // skip line breaks if requested; keep reading
                        continue;
                    }

                    nextCodePoint = ch;
                    return true;
                }
            } catch (IOException e) {
                finish();
                throw new UncheckedIOException("Error reading characters", e);
            }
        }

        @Override
        public Character next() {
            if (!hasNext()) {
                throw new NoSuchElementException("No more characters");
            }
            int ch = nextCodePoint;
            nextCodePoint = UNINITIALIZED;
            return (char) ch;
        }

        private void finish() {
            if (!finished) {
                finished = true;
                try {
                    reader.close();
                } catch (IOException e) {
                    throw new UncheckedIOException(e);
                } finally {
                    // Inform the outer DatasetReader
                    onClose.run();
                }
            }
        }
    }

}
