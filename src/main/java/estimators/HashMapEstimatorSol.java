package estimators;


import search.Pattern;

import java.util.HashMap;

import static solvers.PatternPrunerHalley.solveB;

public class HashMapEstimatorSol implements Estimator {

    public HashMap<Integer, Integer> frequencies;

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
    @Override
    public double[] estimateALl(Pattern p){
        double[] result = new double[p.nGramToInt.length];
        for(int i = 0; i < p.nGramToInt.length; i++){
            result[i] = this.estimate(p.nGramToInt[i]);
        }
        return result;
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
