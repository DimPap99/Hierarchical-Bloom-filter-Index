package PMIndex;

import search.Pattern;
import utilities.PatternResult;

import java.util.ArrayList;

public interface IPMIndexing {

    void insert(String key);

    boolean exists(String key);

    ArrayList<Integer> report(Pattern key);

    public void expire();

    ArrayList<Long> getAvgTimes(Pattern pat);

    PatternResult getLatestStats();


    int getTokenId(String key);
}
