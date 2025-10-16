package utilities;

import java.util.Iterator;


 //Adapter to use SegmentReader with the generic Reader<String> interface.

 public final class SegmentReaderAdapter implements Reader<String> {
    private final SegmentReader delegate;
    private Iterator<String> it;

    public SegmentReaderAdapter(SegmentReader delegate) {
        this.delegate = delegate;
        this.it = delegate.iterator();
    }

    @Override
    public boolean hasNext() {
        return it != null && it.hasNext();
    }

    @Override
    public String next() {
        return it.next();
    }

    @Override
    public void setQueryMode() {
        delegate.setQueryMode();
        this.it = delegate.iterator();
    }
}

