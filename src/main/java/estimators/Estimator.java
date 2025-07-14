package estimators;

public interface Estimator {
    public void init(int totalRecords);
    public void insert(int key);
    public double estimate(int key);

    public int minPruneLevel(String pattern, int streamLength, int alphabetSize, double confidence);

    public double cost(double bloomCost, double leafCost, int Lp, int maxLevel, int patLength);
}
