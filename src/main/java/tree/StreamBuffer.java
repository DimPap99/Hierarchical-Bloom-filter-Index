package tree;

import java.util.ArrayList;

/**
 * Grow-only text store with constant-time random access.
 * Uses {@link StringBuilder} so HotSpot can scalar-replace it
 * during tight loops.
 * TODO: Instead of string buidler use int array and encode each char into an integer
 */
public final class StreamBuffer {

    public final ArrayList<Integer> data = new ArrayList<Integer>();
    /** Global position of the first char still kept in {@link #data}. */
    private long offset = 0;

    public void append(int c) { data.add(c); }

    /** 0-based global position of the last appended char. */
    public long tailPos() { return offset + data.size() - 1; }

    /** Random access by global position. */
    public int charAt(long globalPos) {
        int local = (int) (globalPos - offset);
        return data.get(local);
    }


    public int length() { return data.size(); }
}
