package PMIndex;

import algorithms.SearchAlgorithm;
import membership.BloomFilter;

import java.util.ArrayList;

public class HBI implements IPMIndexing {


    public int windowLength;
    public int indexedItemsCounter=-1;
    private ImplicitTree tree;
    public SearchAlgorithm searchAlgo;
    public HBI(SearchAlgorithm algo, int windowLength, double fpRate, int alphabetSize){
        this.searchAlgo = algo;
        this.windowLength = windowLength;
        this.tree = new ImplicitTree(windowLength, new BloomFilter(), fpRate, alphabetSize);
    }


    public void insert(String key){
        indexedItemsCounter++;
    }

    //Simple existence query
    public boolean exists(String key){
        return false;
    }

    public ArrayList<Integer> report(String key){
        return null;
    }
}
