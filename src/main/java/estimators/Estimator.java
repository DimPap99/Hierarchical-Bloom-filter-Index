package estimators;

public interface Estimator {
    public void init(int totalRecords);
    public void insert(char key);
    public double estimate(char key);

    public int minPruneLevel(String pattern, int streamLength, int alphabetSize, double confidence);
}
