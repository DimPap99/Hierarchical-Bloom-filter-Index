package hbi;

import membership.Membership;
import utilities.HBILogger;
import utilities.Utils;

public class ImplicitTree {

    //Items will be indexed as they come. The position will be given based on total items
    //indexed for the allocated interval of the tree. Once the items indexed surpass the
    //items indexed by the tree, we will be creating a new block.
    int intervalSize;
    int maxDepth;
    int indexedItemsCounter=-1;

    Membership membership;

    public static int ROOT_INTERVAL_IDX = -1; //The root pretty much covers everything, so we use an
    //arbitrary value to encode that. Make sure that value cannot end up at another interval

    public void ImplictTree(int intervalSize){
           this.intervalSize=intervalSize;
           this.maxDepth= (int) Math.ceil(Math.log(intervalSize)); //We also have a final level with the actual positions but dont add it here
    }

    public String createCompositeKey(int currentDepth, int intervalIdx, String input){

        String compositeKey="";
        compositeKey+=currentDepth +"|"+intervalIdx+"|"+input;
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
        for(int i=0;i<this.maxDepth;i++){
            //0 is root
            intervalIdx = Utils.getIntervalIndex(this.maxDepth, i, this.indexedItemsCounter);
            key = this.createCompositeKey(i, intervalIdx, input);
            HBILogger.debug("Inserting item: "+input+" key: " + key);
            this.membership.insert(key);

        }

    }


}
