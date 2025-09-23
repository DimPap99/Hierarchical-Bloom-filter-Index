package utilities;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.function.Supplier;

// Generic fixed-size circular buffer preserving insertion order on iteration for ngrams.
public class RingBuffer<T> {

    private final Object[] buf;
    private int            pos    = 0;
    private int            filled = 0;

    public RingBuffer(int capacity) {
        if (capacity <= 0) {
            throw new IllegalArgumentException("capacity must be > 0");
        }
        this.buf = new Object[capacity];
    }

    public int capacity() {
        return buf.length;
    }

    public int size() {
        return Math.min(filled, buf.length);
    }

    public boolean isFilled() {
        return filled >= buf.length;
    }

    public boolean isEmpty() {
        return size() == 0;
    }

    public void append(T value) {
        Objects.requireNonNull(value, "value");
        buf[pos] = value;
        pos = (pos + 1) % buf.length;
        if (filled < buf.length) {
            filled++;
        }
    }

    public void forEach(Consumer<? super T> action) {
        Objects.requireNonNull(action, "action");
        if (isEmpty()) {
            return;
        }
        if (!isFilled()) {
            for (int i = 0; i < filled; i++) {
                action.accept(elementAt(i));
            }
            return;
        }

        int capacity = buf.length;
        int tail = capacity - pos;
        for (int i = 0; i < tail; i++) {
            action.accept(elementAt(pos + i));
        }
        for (int i = 0; i < pos; i++) {
            action.accept(elementAt(i));
        }
    }

    public String snapshot() {
        if (isEmpty()) return "";

        int expected = isFilled() ? capacity() : size();
        StringBuilder out = new StringBuilder(expected);
        forEach(ch -> out.append(ch));
        return out.toString();
    }

    public <R> R collect(Supplier<R> supplier, BiConsumer<R, ? super T> accumulator) {
        Objects.requireNonNull(supplier, "supplier");
        Objects.requireNonNull(accumulator, "accumulator");
        R container = supplier.get();
        forEach(item -> accumulator.accept(container, item));
        return container;
    }

    @SuppressWarnings("unchecked")
    private T elementAt(int index) {
        return (T) buf[index];
    }
}
