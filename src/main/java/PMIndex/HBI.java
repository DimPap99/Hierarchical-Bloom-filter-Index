package PMIndex;

import estimators.CostFunction;
import estimators.CostFunctionMaxProb;
import membership.*;
import org.apache.commons.math3.util.Pair;
import search.*;
import estimators.Estimator;
import tree.ImplicitTree;
import tree.TreeLayout;
import utilities.PatternResult;

import java.util.*;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import static utilities.MathUtils.pruningLevel;

/**
 * Sliding-window index façade (legacy name “HBI”) adapted to the new
 * tree-layer split.  We kept the public constructor and methods unchanged
 * so the rest of the codebase keeps compiling.
 */
public final class HBI implements IPMIndexing {

    /* ---------------------------------------------------------- config */
    private final int     windowLength;
    private final int     treeLength;
    private final int     alphabetSize;
    private final double  fpRate;

    //dirty way of getting results of a query when reporting (should fix)
    private PatternResult latestQresults;

    private int BC_COST_ESTIM_ITER = 1000000;
    private int LC_COST_ESTIM_ITER = 1000000;

    //maps characters/strings (for n-grams) to an integer
    public HashMap<String, Integer> alphabetMap;

    private       long    indexedItemsCounter = -1;

    private final ArrayList<ImplicitTree<Membership>> trees
            = new ArrayList<>();

    private final SearchAlgorithm               searchAlgo;
    private final Supplier<Estimator> estimatorFac;
    private final Supplier<Membership>    membershipFac;
    private final Supplier<PruningPlan> pruningPlanFac;
    private LongKey codec;
    CostFunction cf;
    public boolean getStats = false;


    double bfCost = 96;
    double leafCost = 26;
    public ArrayList<Integer> Lp = new ArrayList<>();
    public ArrayList<Double> alphas = new ArrayList<>();
    private Verifier verifier;
    double conf;
    int maxActiveTrees;
    public HBI(SearchAlgorithm algo,
               int windowLength,
               double fpRate,
               int alphabetSize,
               int treeLength,
               Supplier<Estimator> estimator,
               Supplier<Membership> membership,
               Supplier<PruningPlan> pruningPlan,
               Verifier verifier, CostFunction cf, double conf) {

        this.searchAlgo    = algo;
        this.windowLength  = windowLength;
        this.treeLength    = treeLength;
        this.fpRate        = fpRate;
        this.alphabetSize  = alphabetSize;
        this.estimatorFac  = estimator;
        this.membershipFac = membership;
        this.pruningPlanFac = pruningPlan;
        this.verifier      = verifier;
        /* first tree */
        this.cf  = cf;
        trees.addLast(createTree());
        this.conf = conf;
        this.maxActiveTrees = (int) Math.ceil(windowLength/treeLength);
    }


    /** Evict the oldest tree unconditionally */
    public void expire() { trees.removeFirst(); }

    /** Stream a single character into the index. */
    public void insert(String c) {

        int intC = this.alphabetMap.get(c);
        indexedItemsCounter++;
        ImplicitTree lastTree = trees.getLast();
        if (lastTree.isFull()) {
            ImplicitTree<Membership> fresh = createTree();
            fresh.id = trees.size();
            fresh.estimator.insert(intC);
            fresh.append(intC, indexedItemsCounter);
            trees.add(fresh);
        } else {
            lastTree.estimator.insert(intC);
            lastTree.append(intC, indexedItemsCounter);
        }
        //think this through properly
        if(trees.getLast().indexedItemsCounter + (trees.size()-2)* this.treeLength > this.windowLength) {
            this.expire();
        }
    }
    public void  fillStackLp(int lp, Deque<Frame> fStack){
        for(int i=0; i< Math.pow(2, lp); i++){
            fStack.addLast(new Frame(lp, i));
        }
    }
    /** Simple existence query — still a TODO. */
    public boolean exists(String key) { return false; }

    /** Multi-match report delegates to the search algorithm unchanged. */
    public ArrayList<Integer> report(Pattern pat) {
    ArrayList<Integer> results = new ArrayList<>();
    //convert ngram to int representation
    for(int nIdx = 0; nIdx < pat.nGramArr.length; nIdx++) pat.nGramToInt[nIdx] = this.alphabetMap.get(pat.nGramArr[nIdx]);

     int positionOffset = -1;
     long startTime = System.currentTimeMillis();
     int lp =0;
     int lpCf = 0;
     for (int i = 0; i < this.trees.size(); i++) {
         ImplicitTree<Membership> tree = this.trees.get(i);
         tree.pruningPlan    = this.pruningPlanFac.get();
         IntervalScanner scn = new IntervalScanner(tree, pat, searchAlgo, positionOffset);
         Deque<Frame> stack = new ArrayDeque<>();
//         double[] pp = tree.estimator.estimateALl(pat);
         double pMax = Arrays.stream(tree.estimator.estimateALl(pat)).min().getAsDouble();
         lpCf = cf.minCostLp(tree, 0.001, pat, this.bfCost, this.leafCost);
         lp = pruningLevel(tree, this.conf, pMax);//cf.minCostLp(tree, 0.001, 0.5, pat, this.bfCost, this.leafCost);


         pat.charStartLp = new ArrayList<>();
         pat.charStartLp.add(lp);
         if(this.getStats) {
             this.Lp.add(lp);
             this.alphas.add(cf.getAlpha());
         }
         //fill the stack with the initial minimum Lp level
         fillStackLp(lp, stack);
         scn.seedStack(stack);

       while (scn.hasNext()) {
           CandidateRange cr = scn.next();
           if(cr == null) break;
           Pair<ArrayList<Integer>, Integer> res = this.verifier.verify(i, this.trees, cr, pat);
           for(int r : res.getFirst()){
               results.add(r);
           }
           scn.positionOffset = res.getSecond();        // advance the starting point

       }
     }
     this.latestQresults = new PatternResult(System.currentTimeMillis() - startTime, this.getAllprobes(), lp, pat, lpCf);
     return results;

    }


    public long getAvgBloomCost(Pattern pat) {
        ImplicitTree<Membership> tree = this.trees.getLast();
        int maxLvl = tree.maxDepth() -1;
//        for (int i = 0; i < maxLvl; i++) {
            long duration = 0;
            for (int z = 0; z < BC_COST_ESTIM_ITER; z++) {  // <-- BUG: tests i and increments i
                long startTime = System.nanoTime();

                long key = tree.codec.pack(maxLvl, 0, pat.nGramToInt[0]);
                tree.contains(maxLvl, key);

                duration += System.nanoTime() - startTime;
            }
            duration /= BC_COST_ESTIM_ITER;
        //}
        return duration;
    }


    public ArrayList<Long> getAvgTimes(Pattern pat) {
        long avgBloom =  getAvgBloomCost(pat);
        long avgLeaf = getAvgLeafCost(pat);
        this.bfCost = avgBloom;
        this.leafCost = avgLeaf;
        ArrayList<Long> res = new ArrayList<>();
        res.add(avgBloom);
        res.add(avgLeaf);
        return  res;
    }

    @Override
    public PatternResult getLatestStats() {
        return this.latestQresults;
    }

    public long getAvgLeafCost(Pattern pat) {
        ImplicitTree<Membership> tree = this.trees.getLast();
        long duration=0;

        for(int i = 0; i < LC_COST_ESTIM_ITER; i++){
            long startTime = System.nanoTime();
            for(int j = 0; j < pat.nGramToInt.length-1; j++){
                boolean in = tree.buffer.data.get(0) == 1;

            }
            duration += System.nanoTime() - startTime;

        }
        duration /= LC_COST_ESTIM_ITER;
        return duration;
    }

    /* -------------------------------------------------------- helpers */

    /** Builds a TreeLayout whose root spans {@code treeLength} chars. */
    private TreeLayout makeLayout() {
        // choose depth so that each deeper level halves the span until 1 char
        int levels = Integer.SIZE - Integer.numberOfLeadingZeros(treeLength);
        int leafSpan = treeLength >>> (levels - 1);   // normally 1
        return new TreeLayout(levels, leafSpan);
    }

    /** Factory for a brand-new ImplicitTree wired to our factories. */
    private ImplicitTree<Membership> createTree() {

        TreeLayout layout = makeLayout();          // same helper as before
        int maxDepth = layout.levels();
    /* --------------------------------------------------------------
       For level L:
         nodes      = 2^L                      (intervals on that level)
         interval   = treeLength / 2^L         (chars per node)
         perNode    = min(alphabetSize, interval)
         distinct   = nodes * perNode          (what you called currentLevelItems)
       -------------------------------------------------------------- */
        IntFunction<Membership> filterFactory = level -> {
            int nodes     = 1 << level;                 // 2^level
            int interval  = treeLength >> level;        // treeLength / 2^level
            int perNode   = Math.min(alphabetSize, interval);
            int distinct  = nodes * perNode;

            BloomFilter bf = new BloomFilter();
            bf.init(distinct, fpRate);
            return bf;
        };

        return new ImplicitTree<>(
                layout,
                filterFactory,          // one BloomFilter per level
                new Key64(maxDepth, alphabetSize),
                estimatorFac.get());
    }

    public int getAllprobes(){

        int sum=0;
        for(ImplicitTree<Membership> tree: this.trees){
            sum+= tree.containCounter;
            tree.containCounter = 0;
        }
        return sum;
    }

}
