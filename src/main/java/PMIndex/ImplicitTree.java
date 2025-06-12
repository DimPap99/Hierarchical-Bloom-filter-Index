package PMIndex;

import membership.Membership;
import utilities.HBILogger;
import utilities.Utils;

import java.util.ArrayList;
import java.util.List;

public class ImplicitTree {

    //Items will be indexed as they come. The position will be given based on total items
    //indexed for the allocated interval of the tree. Once the items indexed surpass the
    //items indexed by the tree, we will be creating a new block.
    public int intervalSize;
    public int maxDepth;
    public int indexedItemsCounter=-1;
    public int maxIdx;
    public Membership membership;
    int treeId;
    public static int ROOT_INTERVAL_IDX = -1; //The root pretty much covers everything, so we use an
    //arbitrary value to encode that. Make sure that value cannot end up at another interval

    public ImplicitTree(int intervalSize, Membership membership, double fpRate, int alphabetSize, int id){
        this.intervalSize=intervalSize;
        this.maxDepth = (int) Math.ceil(Math.log(intervalSize) / Math.log(2));   // log₂(n)
        //In my window (or block), for every position we perform d + 1 insertions (we count from 0). Therefor,
        //the number of keys to be inserted are equal to the size of the window multiplied by the depth of the tree.
        //for n grams this is min(alphabetsize^n, window-n+1)
        int distinctItems = Math.min(alphabetSize, this.intervalSize);
        int expectedInsertions = distinctItems * (this.maxDepth + 1);
        this.membership = membership;
        this.membership.init(expectedInsertions, fpRate);
        this.maxIdx = (int)Math.pow(2, this.maxDepth) - 1;
        this.treeId = id;

    }

    public int getIntervalSize(int level){
        return 1 << (this.maxDepth - level);
    }
    public int calculateDistinctItems(int alphabetSize){
        int distinctItems = 0;
        for(int i=0; i< this.maxDepth+1; i++){
            if(i == this.maxDepth){
                distinctItems+= alphabetSize * this.intervalSize; //For the levels we use ceil. Actual elements at last level
                //are <= Math.pow(2, lastDepth).
            }else{
                distinctItems+= alphabetSize * Math.pow(2, i);
            }
        }
        return distinctItems;
    }

    public String createCompositeKey(int currentDepth, int intervalIdx, String input){
        String compositeKey="";
        if(currentDepth == this.maxDepth){
            compositeKey+=currentDepth +"|"+intervalIdx%(this.maxIdx+1) +"|"+input;
        }else{
            compositeKey+=currentDepth +"|"+intervalIdx +"|"+input;
        }

        return compositeKey;
    }
    public int getLeftChild(int currentIntervalIdx){
        return currentIntervalIdx << 1;
    }

    public int getRightChild(int currentIntervalIdx){
        return (currentIntervalIdx << 1) + 1;
    }


    public void insert(String input){
        String key;
        //Increase the counter that keeps track of indexed items. This is also the position of the item we add
        this.indexedItemsCounter++;
        int intervalIdx;
        //System.out.println("Inserting item: "+input);
        for(int i=0;i<this.maxDepth+1;i++){
            //0 is root
            intervalIdx = Utils.getIntervalIndex(this.maxDepth, i, this.indexedItemsCounter);
            key = this.createCompositeKey(i, intervalIdx, input);
            //String intervalStr = Utils.intervalWithHashes(this.maxDepth, i, this.indexedItemsCounter);
            //System.out.println(intervalStr);
            //HBILogger.debug("Inserting item: "+input+" key: " + key);

            this.membership.insert(key);

        }
        //last level we index without mask
    }


}
