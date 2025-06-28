package tree;

/**
 * Grow-only text store with constant-time random access.
 * Uses {@link StringBuilder} so HotSpot can scalar-replace it
 * during tight loops.
 * TODO: Instead of string buidler use int array and encode each char into an integer
 */
public final class StreamBuffer {

    public final StringBuilder data = new StringBuilder();
    /** Global position of the first char still kept in {@link #data}. */
    private long offset = 0;

    public void append(char c) { data.append(c); }

    /** 0-based global position of the last appended char. */
    public long tailPos() { return offset + data.length() - 1; }

    /** Random access by global position. */
    public char charAt(long globalPos) {
        int local = (int) (globalPos - offset);
        return data.charAt(local);
    }

    /**
     * Discard everything strictly before {@code untilExclusive}.
     * Callers must guarantee no later query needs that text.
     */
    public void dropPrefix(long untilExclusive) {
        int newStart = (int) (untilExclusive - offset);
        if (newStart <= 0) return;
        data.delete(0, newStart);
        offset = untilExclusive;
    }

    public int length() { return data.length(); }
}
