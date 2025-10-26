package estimators;


import search.Pattern;

import java.util.HashMap;

import static solvers.PatternPrunerHalley.solveB;

public class HashMapEstimatorSol implements Estimator {

    public HashMap<Long, Integer> frequencies;

    public double totalRecords;

    public HashMapEstimatorSol(int totalRecords){
        this.init(totalRecords);
    }

    @Override
    public void init(int totalRecords) {
        this.totalRecords = totalRecords;
        this.frequencies = new HashMap<>();
    }

    @Override
    public void insert(long key) {
        frequencies.merge(key, 1, Integer::sum);
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
        return 0;
    }

    @Override
    public double[] estimateALl(Pattern p, boolean strides){
        double[] result = new double[p.nGramToLong.length];
        for(int i = 0; i < p.nGramToLong.length; i++){
            result[i] = this.estimate(p.nGramToLong[i]);
        }
        return result;
    }

    @Override
    public double getMin() {
        return 0;
    }

//    public double nodeWith_Conf(double confidence, double prob){
//
//        return Math.log(1 - confidence) / Math.log(1 - prob);
//    }

//    @Override
//    public int minPruneLevel(String pattern, int streamLength, int alphabetSize, double confidence) {
//        double b_a;
//        double p_max = 0;
//        double p_temp;
//        double[] pHat = new  double[pattern.length()];
//        int i =0;
//        for( char c : pattern.toCharArray()){
//            pHat[i] = this.estimate(c);
//            i++;
//        }
//
//        b_a = solveB(pHat, confidence, 1e-10, 1e-12);
//        return (int) (Math.floor(Math.log(streamLength/b_a)) + 1);
//    }

    public double cost(double bloomCost, double leafCost, int Lp, int maxLevel, int patLength){
        return 0;
    }
}
