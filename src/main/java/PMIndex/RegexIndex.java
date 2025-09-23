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

    /** Append one or more characters to the indexed text. */
    @Override
    public void insert(String  s) {
        text.append(s);
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

        return positions;
    }

    @Override
    public void expire() {

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
}
