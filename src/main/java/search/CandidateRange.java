package search;

public record CandidateRange(int startPos, int endPos) {
    public int length() { return endPos - startPos + 1; }
}
