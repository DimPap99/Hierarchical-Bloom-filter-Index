import PMIndex.HBI;
import algorithms.BlockSearch;

public class Main {

    public static void main(String[] args) {
        Experiment.run("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/src/main/java/simple.txt",
                "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/src/main/java/queries.txt",
                new HBI(new BlockSearch(), 32, 0.01, 26));
    }
}
