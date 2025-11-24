package PMIndex;

import utilities.PatternResult;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

// Simple index that stores text in a StringBuilder and answers queries with java.util.regex.
public class RegexIndex implements IPMIndexing {

    // Growable buffer that accumulates the original file contents.
    private final StringBuilder text = new StringBuilder();
    // Threshold at which the internal buffer expires. Integer.MAX_VALUE disables auto-expire in practice.
    private final int lengthAttribute;

    // Create a RegexIndex with no practical auto-expiration threshold.
    public RegexIndex() {
        this(Integer.MAX_VALUE);
    }

    // Create a RegexIndex that expires once the accumulated text reaches the given threshold.
    public RegexIndex(int lengthAttribute) {
        if (lengthAttribute <= 0) {
            throw new IllegalArgumentException("lengthAttribute must be positive");
        }
        this.lengthAttribute = lengthAttribute;
    }

    // Append one or more characters to the indexed text.
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

        // Compile the pattern on every call; caching is unnecessary for many distinct queries.
        Pattern pattern = Pattern.compile(regex.patternTxt, Pattern.LITERAL);
        Matcher matcher = pattern.matcher(text);

        // Standard forward search; for overlapping matches, adjust the region and increment by 1.
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


    // Number of characters currently indexed.
    public int length() {
        return text.length();
    }

    // Free the internal buffer for GC when you are done.
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
