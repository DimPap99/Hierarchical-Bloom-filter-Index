import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import algorithms.BlockSearch;
import estimators.Estimator;
import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;

import javax.swing.*;
import javax.swing.text.TableView;
import java.io.IOException;
import java.util.function.Supplier;

public class Main {
    //z00D
    //zs1000
    public static void main(String[] args) throws IOException {
        String dataFile = "/home/dimpap/IdeaProjects/Hierarchical-Bloom-filter-Index/zipf_text.txt";
        String queries = "/home/dimpap/IdeaProjects/Hierarchical-Bloom-filter-Index/substrings.txt";
        double durationHBI = 0;
        double durationIPM = 0;
        int runs = 50;
        //TODO: Adds factory for initialization of membership DS. Potentially decouple Algo code
        //do 3 warmup runs for jit
        for(int i = 0; i < 2; i++){
            Supplier<Estimator> factory = () -> new HashMapEstimator(65536);   // same instance
            Supplier<Membership> memFactory = () -> new BloomFilter();
            HBI hbi = new HBI(new BlockSearch(), 131072, 0.001, 81, 65536, factory, memFactory);
            Experiment.run(dataFile, queries, hbi, false);
            IPMIndexing index = new RegexIndex();   // switch index implementation
            Experiment.run(dataFile, queries, index, false);
        }

        for(int i = 0; i < runs; i++){

            /* build a supplier – could be a lambda or method reference */
            Supplier<Estimator> factory = () -> new HashMapEstimator(65536);   // same instance
            Supplier<Membership> memFactory = () -> new BloomFilter();

            HBI hbi = new HBI(new BlockSearch(), 131072, 0.001, 81, 65536, factory, memFactory);
            durationHBI+= Experiment.run(dataFile, queries, hbi, false);
            IPMIndexing index = new RegexIndex();   // switch index implementation
            durationIPM+= Experiment.run(dataFile, queries, index, false);
        }
        if(runs > 0){
            System.out.println("HBI avg (ms):" + (durationHBI/runs));
            System.out.println("IPM avg (ms):" + (durationIPM/runs));

        }


//        System.out.println("Holding: " + hbi);
//        System.in.read(); // Wait so you can take heap dump

    }
}
