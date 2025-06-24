package PMIndex;

import membership.Key;
import membership.Key32;
import membership.Key64;
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
    public Key64 keyService;
    int treeId;

    public StringBuilder stream;
    private static final ThreadLocal<StringBuilder> TMP =
            ThreadLocal.withInitial(StringBuilder::new);
    public static int ROOT_INTERVAL_IDX = -1; //The root pretty much covers everything, so we use an
    //arbitrary value to encode that. Make sure that value cannot end up at another interval

    public ImplicitTree(int intervalSize, Membership membership, double fpRate, int alphabetSize, int id) {
        this.intervalSize=intervalSize;
        this.maxDepth = (int) Math.ceil(Math.log(intervalSize) / Math.log(2));   // log₂(n)
        //In my window (or block), for every position we perform d + 1 insertions (we count from 0). Therefor,
        //the number of keys to be inserted are equal to the size of the window multiplied by the depth of the tree.
        //for n grams this is min(alphabetsize^n, window-n+1)

        int expectedInsertions = this.calculateDistinctItems(alphabetSize);
        this.membership = membership;
        this.membership.init(expectedInsertions, fpRate);
        this.maxIdx = (int)Math.pow(2, this.maxDepth) - 1;
        this.treeId = id;
        this.keyService = new Key64(maxDepth, 8);
        this.stream = new StringBuilder();
    }

    public int getIntervalSize(int level){
        return 1 << (this.maxDepth - level);
    }

    public int calculateDistinctItems(int alphabetSize) {
        int total = 0;
        for (int level = 0; level < this.maxDepth; level++) {
            int nodes     = 1 << level;                 // 2^level nodes
            int interval  = this.intervalSize >> level;  // N / 2^level positions per node
            int perNode   = Math.min(alphabetSize, interval);
            total += nodes * perNode;
        }
        return total;
    }

    public String createCompositeKey(int currentDepth, int intervalIdx, char input){
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


    public void insert(char input){
        long key;
        //Increase the counter that keeps track of indexed items. This is also the position of the item we add
        this.indexedItemsCounter++;
        int intervalIdx;
        //we insert keys up to maxdepth and not at the last level (interval with all positions) as for that level
        //we store the whole stream
        for(int i=0;i<this.maxDepth;i++){
            //0 is root
            intervalIdx = Utils.getIntervalIndex(this.maxDepth, i, this.indexedItemsCounter);

            key = this.keyService.pack(i, intervalIdx, input);//this.createCompositeKey(i, intervalIdx, input);
            //String intervalStr = Utils.intervalWithHashes(this.maxDepth, i, this.indexedItemsCounter);
            //System.out.println(intervalStr);
            //HBILogger.debug("Inserting item: "+input+" key: " + key);

            this.membership.insert(key);

        }
        this.stream.append(input);
        //last level we index without mask
    }


}
