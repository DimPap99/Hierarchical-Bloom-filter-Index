package PMIndex;

import membership.*;
import org.apache.commons.math3.util.Pair;
import search.*;
import estimators.Estimator;
import tree.ImplicitTree;
import tree.TreeLayout;

import java.util.ArrayList;
import java.util.ArrayDeque;
import java.util.Deque;
import java.util.HashMap;
import java.util.function.IntFunction;
import java.util.function.Supplier;

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

    //maps characters/strings (for n-grams) to an integer
    public HashMap<String, Integer> alphabetMap;

    /* ----------------------------------------------------------- state */
    private       long    indexedItemsCounter = -1;

    private final ArrayList<ImplicitTree<Membership>> trees
            = new ArrayList<>();

    /* ---------------------------------------------------------- wiring */
    private final SearchAlgorithm               searchAlgo;
    private final Supplier<Estimator> estimatorFac;
    private final Supplier<Membership>    membershipFac;
    private final Supplier<PruningPlan> pruningPlanFac;
    private LongKey codec;

    private Verifier verifier;
    public HBI(SearchAlgorithm algo,
               int windowLength,
               double fpRate,
               int alphabetSize,
               int treeLength,
               Supplier<Estimator> estimator,
               Supplier<Membership> membership,
               Supplier<PruningPlan> pruningPlan,
               Verifier verifier) {

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
        trees.addLast(createTree());
    }

    /* ----------------------------------------------------------- api */

    /** Evict the oldest tree unconditionally (same semantics as before). */
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
     for (int i = 0; i < this.trees.size(); i++) {
         ImplicitTree<Membership> tree = this.trees.get(i);
         tree.pruningPlan    = this.pruningPlanFac.get();
         IntervalScanner scn = new IntervalScanner(tree, pat, searchAlgo, positionOffset);
         Deque<Frame> stack = new ArrayDeque<>();
         int lp = tree.pruningPlan.pruningPlan(pat, tree, alphabetSize, 0.99).getFirst();
         fillStackLp(lp, stack);
         scn.seedStack(stack);

       while (scn.hasNext()) {
           CandidateRange cr = scn.next();
           if(cr == null) break;
           Pair<ArrayList<Integer>, Integer> res = this.verifier.verify(i, this.trees, cr, pat);
           for(int r : res.getFirst()){
               results.add(r);
           }
      }
     }
     return results;

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
       That is *identical* to the iniMembershipPerLevel() body you posted.
       -------------------------------------------------------------- */
        IntFunction<Membership> filterFactory = level -> {
            int nodes     = 1 << level;                 // 2^level
            int interval  = treeLength >> level;        // treeLength / 2^level
            int perNode   = Math.min(alphabetSize, interval);
            int distinct  = nodes * perNode;            // ← exact same n

            BloomFilter bf = new BloomFilter();
            bf.init(distinct, fpRate);                  // SAME CALL as before
            return bf;
        };

        return new ImplicitTree<>(
                layout,
                filterFactory,          // one ready-to-use BloomFilter per level
                new Key64(maxDepth, alphabetSize),          // your 64-bit key codec
                estimatorFac.get());
    }

}
