package search;

public record Frame (
    int level,
    int intervalIdx){
    public int intervalIdx() { return intervalIdx;}
    public int level() { return level; }


}
