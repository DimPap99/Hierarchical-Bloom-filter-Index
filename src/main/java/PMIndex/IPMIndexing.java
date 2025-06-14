package PMIndex;

import java.util.ArrayList;

public interface IPMIndexing {

    void insert(char key);

    boolean exists(String key);

    ArrayList<Integer> report(String key);

    public void expire();
}
