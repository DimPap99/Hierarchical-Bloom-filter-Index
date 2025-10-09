package PMIndex;

import utilities.PatternResult;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * A very small “index” that just stores the original text in a StringBuilder
 * and answers pattern queries with java.util.regex.
 *
 *  • insert(String s)  – append characters/strings to the internal buffer
 *  • report(String regex) – return start positions (0-based) of every match
 *

 */
public class RegexIndex implements IPMIndexing {

    // Growable buffer that accumulates the original file contents
    private final StringBuilder text = new StringBuilder();
    // Threshold at which the internal buffer expires (clears). If set to
    // Integer.MAX_VALUE, it effectively never expires automatically.
    private final int lengthAttribute;

    /**
     * Create a RegexIndex with no practical auto-expiration threshold.
     * Existing call sites continue to work without changes.
     */
    public RegexIndex() {
        this(Integer.MAX_VALUE);
    }

    /**
     * Create a RegexIndex that automatically expires (clears) its internal
     * buffer once the accumulated text length reaches the provided threshold.
     *
     * @param lengthAttribute positive length threshold for auto-expire
     * @throws IllegalArgumentException if {@code lengthAttribute <= 0}
     */
    public RegexIndex(int lengthAttribute) {
        if (lengthAttribute <= 0) {
            throw new IllegalArgumentException("lengthAttribute must be positive");
        }
        this.lengthAttribute = lengthAttribute;
    }

    /** Append one or more characters to the indexed text. */
    @Override
    public void insert(String  s) {
        text.append(s);
        if (text.length() >= lengthAttribute) {
            expire();
        }
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(search.Pattern regex) {
        ArrayList<Integer> positions = new ArrayList<>();

        // Compile the pattern on every call – you could cache if you repeatedly
        // query the same patterns, but for the typical “many distinct queries”
        // workload the overhead is negligible compared to the scan itself.
        Pattern pattern = Pattern.compile(regex.patternTxt, Pattern.LITERAL);
        Matcher matcher = pattern.matcher(text);

        // Standard forward search – if you need overlapping matches, replace
        // start() + 1 with matcher.start() + 1 and reset matcher.region(...)
        while (matcher.find()) {
            positions.add(matcher.start());
        }
//        for (int i = 0; i < positions.size(); i++) {
//            String s = text.substring(positions.get(i), positions.get(i) + regex.patternTxt.length());
//
//            System.out.println("Pattern: " + regex.patternTxt + " Matched with: "  + s);
//
//        }
        return positions;
    }


    @Override
    public void expire() {
        text.setLength(0);
    }


    /** @return number of characters currently indexed. */
    public int length() {
        return text.length();
    }

    /** Free the internal buffer for GC when you’re done. */
    public void clear() {
        text.setLength(0);
    }

    public ArrayList<Long> getAvgTimes(search.Pattern pat){
        return null;
    }

    @Override
    public PatternResult getLatestStats() {
        return null;
    }

    @Override
    public int getTokenId(String key) {
        return -2;
    }
}
