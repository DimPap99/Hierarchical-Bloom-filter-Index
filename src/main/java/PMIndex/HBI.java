package PMIndex;

import algorithms.SearchAlgorithm;
import membership.BloomFilter;
import membership.MockMembership;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Queue;

public class HBI implements IPMIndexing {


    public int windowLength;
    public int treeLength;
    public int indexedItemsCounter=-1;

    double fpRate;
    int alphabetSize;

    private ArrayList<ImplicitTree> trees;
    public SearchAlgorithm searchAlgo;
    public HBI(SearchAlgorithm algo, int windowLength, double fpRate, int alphabetSize, int tree_length){
        this.searchAlgo = algo;
        this.windowLength = windowLength;
        this.treeLength = tree_length;
        this.trees = new ArrayList<ImplicitTree>();
        this.alphabetSize = alphabetSize;
        this.fpRate = fpRate;
        this.trees.add(new ImplicitTree(treeLength, new BloomFilter(), fpRate, alphabetSize, 0));

    }


    public void expire(){
        this.trees.removeFirst();
    }
    public void insert(char key){
        indexedItemsCounter++;

        ImplicitTree lastTree = trees.getLast();
        if(lastTree.indexedItemsCounter == this.treeLength - 1){
            ImplicitTree newLastTree = new ImplicitTree(treeLength, new MockMembership(), fpRate, alphabetSize, trees.size());
            newLastTree.insert(key);
            trees.add(newLastTree);
        }
        else{
            lastTree.insert(key);
        }
    }

    //Simple existence query
    public boolean exists(String key){
        return false;
    }

    public ArrayList<Integer> report(String key){
        return this.searchAlgo.report(key.toCharArray(), this.trees, false);
    }


}
