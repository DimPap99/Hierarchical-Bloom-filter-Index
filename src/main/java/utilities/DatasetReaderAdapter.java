package utilities;

import java.util.Iterator;

/**
 * Adapter to use DatasetReader with the generic Reader<Character> interface.
 */
public final class DatasetReaderAdapter implements Reader<Character> {
    private final DatasetReader delegate;
    private Iterator<Character> it;

    public DatasetReaderAdapter(DatasetReader delegate) {
        this.delegate = delegate;
        this.it = delegate.iterator();
    }

    @Override
    public boolean hasNext() {
        return it != null && it.hasNext();
    }

    @Override
    public Character next() {
        return it.next();
    }

    @Override
    public void setQueryMode() {
        delegate.setQueryMode();
        this.it = delegate.iterator();
    }
}

