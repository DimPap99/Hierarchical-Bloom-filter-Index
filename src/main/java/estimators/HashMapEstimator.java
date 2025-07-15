package estimators;

import java.util.HashMap;

public class HashMapEstimator implements Estimator {

    public HashMap<Integer, Integer> frequencies;

    public double totalRecords;

    public HashMapEstimator(int totalRecords){
        this.init(totalRecords);
    }

    @Override
    public void init(int totalRecords) {
        this.totalRecords = totalRecords;
        this.frequencies = new HashMap<>();
    }

    @Override
    public void insert(int key) {
        frequencies.merge(key, 1, Integer::sum);
    }

    @Override
    public double estimate(int key) {
        if(frequencies.containsKey(key)){
            return frequencies.get(key)/totalRecords;
        }
        return 0;
    }



}
