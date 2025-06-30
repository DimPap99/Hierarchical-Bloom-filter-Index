package utilities;

import java.util.Arrays;

/** Fixed–size sliding window over a character stream. */
public final class CharRingBuffer {

    public final char[] buf;          // circular buffer
    private int          pos = 0;      // next write position
    private int          filled = 0;   // how many chars we’ve seen (≤ buf.length)

    public CharRingBuffer(int k) {
        if (k <= 0) throw new IllegalArgumentException("k must be > 0");
        this.buf = new char[k];
    }

    public boolean isFilled(){
        return (filled >= buf.length);
    }
    /** Push one character from the stream. */
    public void append(char c) {
        buf[pos] = c; // a b c d      [a, b], [a, c]
        pos = (pos + 1) % buf.length;
        if (filled < buf.length) filled++;
    }

    /** Current window as a plain Java String (length ≤ k). */
    public String view() {
        if (filled == 0) return "";

        // fast path: window has not wrapped yet
        if (filled < buf.length) {
            return new String(buf, 0, filled);
        }

        // wrapped: split into two ranges and copy into tmp char[]
        int k   = buf.length;
        char[] out = new char[k];

        int tail = k - pos;                   // chars from pos..end
        System.arraycopy(buf, pos, out, 0, tail);
        System.arraycopy(buf, 0,  out, tail, pos);

        return new String(out);
    }
}
