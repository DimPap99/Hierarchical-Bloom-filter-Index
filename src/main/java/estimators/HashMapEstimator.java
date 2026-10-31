package estimators;

import search.Pattern;

import java.util.HashMap;

public class HashMapEstimator implements Estimator {

    public HashMap<Long, Integer> frequencies;
    private long minKey = -1;
    public double totalRecords=0;

    public HashMapEstimator(int totalRecords){
        this.init(totalRecords);
    }

    @Override
    public void init(int totalRecords) {
        this.totalRecords = 0;
        this.frequencies = new HashMap<>();
    }

    @Override
    public void insert(long key) {
        frequencies.merge(key, 1, Integer::sum);
        totalRecords+=1;
    }

    @Override
    public double estimate(long key) {
        if (totalRecords <= 0) {
            return 0.0;
        }
        Integer count = frequencies.get(key);
        return (count == null) ? 0.0 : (count / totalRecords);
    }

    @Override
    public long get(long key) {
        return this.frequencies.getOrDefault(key, -1);
    }

    @Override
    public double[] estimateALl(Pattern p, boolean strides){
        long[] arr;
        if(strides){
            arr = p.effectiveNgramArr;
        }
        else{
            arr = p.nGramToLong;
        }
        double[] result = new double[arr.length];
        for(int i = 0; i < arr.length; i++){
            result[i] = this.estimate(arr[i]);
        }
        return result;
    }

    // Per-position rate of the rarest key seen so far, or 0 when nothing was inserted.
    @Override
    public double getMin() {
        if (totalRecords <= 0 || frequencies.isEmpty()) {
            return 0.0;
        }
        long rarest = -1;
        int minCount = Integer.MAX_VALUE;
        for (var e : frequencies.entrySet()) {
            if (e.getValue() < minCount) {
                minCount = e.getValue();
                rarest = e.getKey();
            }
        }
        this.minKey = rarest;
        return minCount / totalRecords;
    }


}
