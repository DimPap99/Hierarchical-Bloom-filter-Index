package search;

public class MatchResult {

    public int pos;
    public boolean matched;

    public  MatchResult(boolean matched, int pos) {
        this.matched = matched;
        this.pos = pos;
    }
}
