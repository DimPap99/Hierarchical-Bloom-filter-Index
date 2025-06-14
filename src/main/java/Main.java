import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import algorithms.BlockSearch;

import java.io.IOException;

public class Main {
    //16384
    public static void main(String[] args) throws IOException {
        HBI hbi = new HBI(new BlockSearch(), 131072, 0.001, 26, 131072);
        String dataFile = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/zipf_text.txt";
        String queries = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/substrings.txt";

        Experiment.run(dataFile, queries, hbi);
        IPMIndexing index = new RegexIndex();   // switch index implementation
        Experiment.run(dataFile, queries, index);

//        System.out.println("Holding: " + hbi);
//        System.in.read(); // Wait so you can take heap dump

    }
}
