package estimators;

import search.Pattern;

import java.util.HashMap;

public class HashMapEstimator implements Estimator {

    public HashMap<Integer, Integer> frequencies;

    public double totalRecords;
    //TODO: Preallocate space on the hashmap to make it faster.
    public HashMapEstimator(int totalRecords){
        this.init(totalRecords);
    }

    @Override
    public void init(int totalRecords) {
        this.totalRecords = 0;
        this.frequencies = new HashMap<>();
    }

    @Override
    public void insert(int key) {
        frequencies.merge(key, 1, Integer::sum);
        totalRecords+=1;
    }

    @Override
    public double estimate(int key) {
        if(frequencies.containsKey(key)){
            return frequencies.get(key)/totalRecords;
        }
        return 0;
    }
    @Override
    public double[] estimateALl(Pattern p){
        double[] result = new double[p.nGramToInt.length];
        for(int i = 0; i < p.nGramToInt.length; i++){
            result[i] = this.estimate(p.nGramToInt[i]);
        }
        return result;
    }




}
