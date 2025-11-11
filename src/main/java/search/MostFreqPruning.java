package search;

import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;

import static utilities.MathUtils.pruningLevel;
import static utilities.MathUtils.pruningLevelBloom;

public class MostFreqPruning implements  PruningPlan {

    public double conf;
    public double fp;
    public void PruningPlan() {
    }
    public MostFreqPruning(double conf, double fp){
       // this.conf = conf;
        this.fp = fp;
    }
    @Override
    public void prepare(Deque<Frame> stack, int lp) {
        for(int i=0; i< Math.pow(2, lp); i++){
            stack.addLast(new Frame(lp, i));
        }
    }

    @Override
        public ArrayList<Integer> pruningPlan(Pattern pattern, ImplicitTree tree, double confidence, boolean strides ) {
        ArrayList<Integer> prunedLevels = new ArrayList<>();
//
//        double b_a;
//        double p_max = 0;
//        double p_temp;
////        ArrayList<Double> probs = new ArrayList<>();
//        for( int c : pattern.nGramToInt){
//            p_temp = tree.estimator.estimate(c);
////            probs.add(p_temp);
//            if(p_temp > p_max){
//                p_max = p_temp;
//            }
//        }
//
//        int    lp   = MathUtils.pruningLevel(tree, confidence, p_max);//Math.max(0, Math.min(raw, tree.maxDepth() - 1));
        double[] pp = tree.estimator.estimateALl(pattern, strides);
        double pMin = Arrays.stream(pp).min().getAsDouble();
        int lp = Math.max(tree.effectiveRoot(), pruningLevelBloom(tree, confidence, pMin, this.fp));
        prunedLevels.add(lp);
        return prunedLevels;
    }





}
