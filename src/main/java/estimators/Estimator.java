package estimators;

import search.Pattern;

public interface Estimator {
    public void init(int totalRecords);
    public void insert(int key);
    public double estimate(int key);

    public double[] estimateALl(Pattern p);

    public double getMin();

}
