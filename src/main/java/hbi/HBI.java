package hbi;

import algorithms.SearchAlgorithm;
import membership.BloomFilter;

public class HBI {


    public int windowLength;
    public int indexedItemsCounter=-1;
    private ImplicitTree tree;
    public SearchAlgorithm searchAlgo;
    public HBI(SearchAlgorithm algo, int windowLength){
        this.searchAlgo = algo;
        this.windowLength = windowLength;
        this.tree = new ImplicitTree(windowLength, new BloomFilter(24, 0.01));
    }

    public void insert(String key){
        indexedItemsCounter++;
    }

    //Simple existence query
    public void exists(String key){}

    public void report(String key){}
}
