import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import algorithms.BlockSearch;
import estimators.Estimator;
import estimators.HashMapEstimator;

import javax.swing.*;
import javax.swing.text.TableView;
import java.io.IOException;
import java.util.function.Supplier;

public class Main {
    //z00D
    //zs1000
    public static void main(String[] args) throws IOException {
        String dataFile = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/zipf_text.txt";
        String queries = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/substrings.txt";
        double durationHBI = 0;
        double durationIPM = 0;
        int runs = 40;
        //TODO: Adds factory for initialization of membership DS. Potentially decouple Algo code
        //do 3 warmup runs for jit
        for(int i = 0; i < 2; i++){
            HashMapEstimator hmEstimator = new HashMapEstimator(161072);
            Supplier<Estimator> factory = () -> hmEstimator;   // same instance
            HBI hbi = new HBI(new BlockSearch(), 262144, 0.001, 81, 65536, factory);
            Experiment.run(dataFile, queries, hbi, false);
            IPMIndexing index = new RegexIndex();   // switch index implementation
            Experiment.run(dataFile, queries, index, false);
        }
        for(int i = 0; i < runs; i++){
            HashMapEstimator hmEstimator = new HashMapEstimator(65536);

            /* build a supplier – could be a lambda or method reference */
            Supplier<Estimator> factory = () -> hmEstimator;   // same instance
            HBI hbi = new HBI(new BlockSearch(), 262144, 0.001, 81, 65536, factory);
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
