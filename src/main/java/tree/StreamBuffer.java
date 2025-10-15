package tree;

/** Grow-only text store with constant-time random access backed by a fixed {@code int[]} capacity. */
public final class StreamBuffer {

    private static final int[] EMPTY = new int[0];

    private final int capacity;
    private int[] data;
    private int size;
    private long offset;
    private boolean detached;

    StreamBuffer(boolean useInts, int maxSize) {
        if (maxSize <= 0) {
            throw new IllegalArgumentException("maxSize must be positive");
        }
        this.capacity = maxSize;
        this.data = new int[maxSize];
        this.size = 0;
        this.offset = 0L;
        this.detached = false;
    }

    public void append(int value) {
        if (detached) {
            throw new IllegalStateException("StreamBuffer is detached");
        }
        if (size >= capacity) {
            throw new IllegalStateException("StreamBuffer overflow: capacity=" + capacity);
        }
        data[size++] = value;
    }

    public void append(long value) {
        if (value < Integer.MIN_VALUE || value > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Token exceeds 32-bit range: " + value);
        }
        append((int) value);
    }

    public int get(int index) {
        if (detached) {
            throw new IllegalStateException("StreamBuffer is detached");
        }
        if (index < 0 || index >= size) {
            throw new IndexOutOfBoundsException(index);
        }
        return data[index];
    }

    public int length() {
        return size;
    }

    public int lastIndex() {
        return size - 1;
    }

    public int capacity() {
        return capacity;
    }

    public Snapshot detach() {
        if (detached) {
            throw new IllegalStateException("StreamBuffer already detached");
        }
        Snapshot snapshot = new Snapshot(data, size, offset);
        data = EMPTY;
        size = 0;
        detached = true;
        return snapshot;
    }

    public void restore(Snapshot snapshot) {
        if (snapshot == null) {
            return;
        }
        if (snapshot.data.length != capacity) {
            throw new IllegalArgumentException("Snapshot does not match buffer capacity");
        }
        this.data = snapshot.data;
        this.size = snapshot.size;
        this.offset = snapshot.offset;
        this.detached = false;
    }

    public long tailPos() {
        return offset + size - 1L;
    }

    public int charAt(long globalPos) {
        int local = (int) (globalPos - offset);
        return get(local);
    }

    public static final class Snapshot {
        private final int[] data;
        private final int size;
        private final long offset;

        private Snapshot(int[] data, int size, long offset) {
            this.data = data;
            this.size = size;
            this.offset = offset;
        }
    }
}
