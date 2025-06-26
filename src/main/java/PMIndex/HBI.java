package PMIndex;

import algorithms.SearchAlgorithm;
import estimators.Estimator;
import membership.BloomFilter;
import membership.Membership;
import membership.MockMembership;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Queue;
import java.util.function.Supplier;

public class HBI implements IPMIndexing {


    public int windowLength;
    public int treeLength;
    public int indexedItemsCounter=-1;

    double fpRate;
    int alphabetSize;

    private ArrayList<ImplicitTree> trees;
    public SearchAlgorithm searchAlgo;
    private Supplier<Estimator> estimatorFac;
    private Supplier<Membership> membershipFac;   // NEW

    public HBI(SearchAlgorithm algo, int windowLength, double fpRate, int alphabetSize, int tree_length, Supplier<Estimator> estimator, Supplier<Membership> membership){
        this.searchAlgo = algo;
        this.windowLength = windowLength;
        this.treeLength = tree_length;
        this.trees = new ArrayList<ImplicitTree>();
        this.alphabetSize = alphabetSize;
        this.fpRate = fpRate;
        this.estimatorFac = estimator;
        this.membershipFac = membership;
        this.trees.add(new ImplicitTree(treeLength, this.membershipFac, fpRate, alphabetSize, 0, this.estimatorFac.get()));

    }


    public void expire(){
        this.trees.removeFirst();
    }
    public void insert(char key){
        indexedItemsCounter++;

        ImplicitTree lastTree = trees.getLast();
        if(lastTree.indexedItemsCounter == this.treeLength - 1){
            ImplicitTree newLastTree = new ImplicitTree(treeLength, this.membershipFac, fpRate, alphabetSize, trees.size(), this.estimatorFac.get());
            newLastTree.estimator.insert(key);
            newLastTree.insert(key);
            trees.add(newLastTree);
        }
        else{
            lastTree.estimator.insert(key);
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
