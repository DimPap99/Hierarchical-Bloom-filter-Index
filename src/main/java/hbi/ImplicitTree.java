package hbi;

public class ImplicitTree {

    //Items will be indexed as they come. The position will be given based on total items
    //indexed for the allocated interval of the tree. Once the items indexed surpass the
    //items indexed by the tree, we will be creating a new block.
    int intervalSize;
    int maxDepth;
    int indexedItems=-1;



    public void ImplictTree(int intervalSize){
           this.intervalSize=intervalSize;
           int a = 5 >>> 2
           this.maxDepth= (int) Math.ceil(Math.log(intervalSize)); //We also have a final level with the actual positions but dont add it here
    }
}
