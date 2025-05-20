package PMIndex;

import algorithms.SearchAlgorithm;
import membership.BloomFilter;
import membership.MockMembership;

import java.util.ArrayList;

public class HBI implements IPMIndexing {


    public int windowLength;
    public int indexedItemsCounter=-1;
    private ImplicitTree tree;
    public SearchAlgorithm searchAlgo;
    public HBI(SearchAlgorithm algo, int windowLength, double fpRate, int alphabetSize){
        this.searchAlgo = algo;
        this.windowLength = windowLength;
        this.tree = new ImplicitTree(windowLength, new MockMembership(), fpRate, alphabetSize);
    }


    public void insert(String key){
        indexedItemsCounter++;
        tree.insert(key);
    }

    //Simple existence query
    public boolean exists(String key){
        return false;
    }

    public ArrayList<Integer> report(String key){
        return this.searchAlgo.report(key, this.tree, false);
    }
}
