package estimators;

import search.Pattern;

public interface Estimator {
    public void init(int totalRecords);
    public void insert(long key);
    public double estimate(long key);
    public long get(long key);

    public double[] estimateALl(Pattern p, boolean strides);

    public double getMin();

}
