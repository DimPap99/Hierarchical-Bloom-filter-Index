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

    @Override
    public double getMin() {


        long min = Long.MAX_VALUE;
        for( long k : this.frequencies.keySet()){
            if(this.frequencies.get(k)<min){ min = k;}
        }
        this.minKey = min;

        return min/this.estimate(min);
    }


}
