package PMIndex;

import java.util.ArrayList;

public interface IPMIndexing {

    void insert(String key);

    boolean exists(String key);

    ArrayList<Integer> report(String key);

    public void expire();
}
