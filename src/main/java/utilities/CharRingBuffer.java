package utilities;


public final class CharRingBuffer extends RingBuffer<Character> {

    public CharRingBuffer(int k) {
        super(k);
    }

    public boolean isFilled(){
        return super.isFilled();
    }
    // Push one character from the stream.
    public void append(char c) {
        super.append(c);
    }

    // Current window as a plain String (length <= k).
    public String view() {
        if (isEmpty()) return "";

        int expected = isFilled() ? capacity() : size();
        StringBuilder out = new StringBuilder(expected);
        forEach(ch -> out.append(ch));
        return out.toString();
    }
}
