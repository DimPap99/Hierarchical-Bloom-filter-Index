package estimators;

import java.util.HashMap;

public class HashMapEstimator implements Estimator {

    public HashMap<Character, Integer> frequencies;

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
    public void insert(char key) {
        frequencies.merge(key, 1, Integer::sum);
    }

    @Override
    public double estimate(char key) {
        if(frequencies.containsKey(key)){
            return frequencies.get(key)/totalRecords;
        }
        return 0;
    }

    public double nodeWith_Conf(double confidence, double prob){

        return Math.log(1 - confidence) / Math.log(1 - prob);
    }

    @Override
    public int minPruneLevel(String pattern, int streamLength, int alphabetSize, double confidence) {
        double b_a;
        double p_max = 0;
        double p_temp;
        for( char c : pattern.toCharArray()){
            p_temp = this.estimate(c);
            if(p_temp > p_max){
                p_max = p_temp;
            }
        }
        b_a = nodeWith_Conf(confidence, p_max);
        return (int) (Math.floor(Math.log(streamLength/b_a)) + 1);
    }
}
