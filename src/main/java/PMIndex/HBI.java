package PMIndex;

import estimators.CostFunction;
import estimators.Estimator;
import membership.Key64;
import membership.Membership;
import org.apache.commons.math3.util.Pair;
import search.CandidateRange;
import search.Frame;
import search.IntervalScanner;
import search.Pattern;
import search.PruningPlan;
import search.SearchAlgorithm;
import search.Verifier;
import tree.ImplicitTree;
import tree.TreeLayout;
import utilities.AlphabetMapper;
import utilities.PatternResult;

import java.util.*;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import static utilities.MathUtils.pruningLevel;


public final class HBI implements IPMIndexing {

    private static final int DEFAULT_BC_COST_ESTIM_ITER = 1_000_000;
    private static final int DEFAULT_LC_COST_ESTIM_ITER = 1_000_000;

    private final HbiConfiguration config;
    private final HbiStats stats;

    private final int windowLength;
    private final int treeLength;
    private final int alphabetSize;
    private final double fpRate;

    private final SearchAlgorithm searchAlgo;
    private final Supplier<Estimator> estimatorFac;
    private final Supplier<Membership> membershipFac;
    private final Supplier<PruningPlan> pruningPlanFac;
    private final Verifier verifier;
    private final CostFunction cf;
    private final double conf;

    private final ArrayList<ImplicitTree<Membership>> trees = new ArrayList<>();

    private AlphabetMapper<String> alphabetMap;
    private long indexedItemsCounter = -1;

    private int bcCostEstimIter = DEFAULT_BC_COST_ESTIM_ITER;
    private int lcCostEstimIter = DEFAULT_LC_COST_ESTIM_ITER;

    private double bfCost = 97;
    private double leafCost = 26;
    private int lpOverride;
    private final int maxActiveTrees;

    public HBI(HbiConfiguration config) {
        this.config = Objects.requireNonNull(config, "config");
        this.stats = new HbiStats(config.collectStats(), config.experimentMode());
        this.windowLength = config.windowLength();
        this.treeLength = config.treeLength();
        this.alphabetSize = config.alphabetSize();
        this.fpRate = config.fpRate();
        this.searchAlgo = config.searchAlgorithm();
        this.estimatorFac = config.estimatorSupplier();
        this.membershipFac = config.membershipSupplier();
        this.pruningPlanFac = config.pruningPlanSupplier();
        this.verifier = config.verifier();
        this.cf = config.costFunction();
        this.conf = config.confidence();
        this.alphabetMap = new AlphabetMapper<>(alphabetSize);
        trees.addLast(createTree());
        this.maxActiveTrees = (int) Math.ceil((double) (windowLength / treeLength));
    }

    public HBI(SearchAlgorithm algo,
               int windowLength,
               double fpRate,
               int alphabetSize,
               int treeLength,
               Supplier<Estimator> estimator,
               Supplier<Membership> membership,
               Supplier<PruningPlan> pruningPlan,
               Verifier verifier,
               CostFunction cf,
               double conf) {
        this(HbiConfiguration.builder()
                .searchAlgorithm(algo)
                .windowLength(windowLength)
                .fpRate(fpRate)
                .alphabetSize(alphabetSize)
                .treeLength(treeLength)
                .estimatorSupplier(estimator)
                .membershipSupplier(membership)
                .pruningPlanSupplier(pruningPlan)
                .verifier(verifier)
                .costFunction(cf)
                .confidence(conf)
                .build());
    }

    public HbiStats stats() {
        return stats;
    }

    public void setAlphabetMap(AlphabetMapper<String> alphabetMap) {
        this.alphabetMap = Objects.requireNonNull(alphabetMap, "alphabetMap");
    }

    public void resetAlphabetMap(int capacity) {
        this.alphabetMap = new AlphabetMapper<>(capacity);
    }

    public void setLpOverride(int lp) {
        this.lpOverride = lp;
    }

    public int lpOverride() {
        return lpOverride;
    }

    /** Evict the oldest tree unconditionally. */
    @Override
    public void expire() {
        trees.removeFirst();
    }

    /** Stream a single character into the index. */
    @Override
    public void insert(String c) {
        int intC = alphabetMap.insert(c);
        indexedItemsCounter++;
        ImplicitTree<Membership> lastTree = trees.getLast();
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
        if (trees.getLast().indexedItemsCounter + (trees.size() - 2) * this.treeLength > this.windowLength) {
            this.expire();
        }
    }

    public void fillStackLp(int lp, Deque<Frame> fStack) {
        for (int i = 0; i < Math.pow(2, lp); i++) {
            fStack.addLast(new Frame(lp, i));
        }
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(Pattern pat) {
        ArrayList<Integer> results = new ArrayList<>();
        long queryStartNanos = System.nanoTime();
        long totalLpTimeNanos = 0L;
        for (int nIdx = 0; nIdx < pat.nGramArr.length; nIdx++) {
            int ngToInt = this.alphabetMap.getId(pat.nGramArr[nIdx]);
            pat.nGramToInt[nIdx] = ngToInt;
            if(nIdx % pat.nGram == 0) pat.effectiveNgramArr[nIdx/pat.nGram] = ngToInt;
        }
        if(pat.originalSz % pat.nGram != 0) pat.effectiveNgramArr[pat.effectiveNgramArr.length-1] = pat.nGramToInt[pat.nGramToInt.length-1];

        int positionOffset = -1;
        long startTime = System.currentTimeMillis();
        int lp = 0;
        int lpCf = 0;
        double cp_cost = 0;
        int arbitraryConfLp = 0;
        for (int i = 0; i < this.trees.size(); i++) {
            ImplicitTree<Membership> tree = this.trees.get(i);
            tree.pruningPlan = this.pruningPlanFac.get();
            IntervalScanner scn = new IntervalScanner(tree, pat, searchAlgo, positionOffset);
            Deque<Frame> stack = new ArrayDeque<>();
            double[] pp = tree.estimator.estimateALl(pat);
            double pMax = Arrays.stream(pp).min().getAsDouble();
            long lpStart = System.nanoTime();
            ArrayList<Integer> lps = null; 
//            lp = Collections.min(lps);
            totalLpTimeNanos += System.nanoTime() - lpStart;
              //cf.minCostLp(tree, 0.05, pat, 97, 26);//pruningLevel(tree, 0.99, pMax);

            if (stats.isExperimentMode()) {
                pp = tree.estimator.estimateALl(pat);
                pMax = Arrays.stream(pp).min().getAsDouble();
                lp = this.lpOverride;
                arbitraryConfLp = pruningLevel(tree, this.conf, pMax);
                int m = (int) (tree.maxDepth() - 1 - Math.ceil(Math.log(pat.nGramToInt.length) / Math.log(2)));
                cp_cost = cf.costAtLevel(tree, pp, pat.effectiveNgramArr, lp, 0.001, m);
            } else {
                lps= tree.pruningPlan.pruningPlan(pat, tree, 0.99);
//                lp = lpCf;
            }

            pat.charStartLp = lps;

            if (stats.isCollecting()) {
                stats.recordLp(lp);
                stats.recordAlpha(cf.getAlpha());
            }
            fillStackLp(lp, stack);
            scn.seedStack(stack);

            while (scn.hasNext()) {
                CandidateRange cr = scn.next();
                if (cr == null) break;
                Pair<ArrayList<Integer>, Integer> res = this.verifier.verify(i, this.trees, cr, pat);
                for (int r : res.getFirst()) {
                    results.add(r);
                }
                scn.positionOffset = res.getSecond();
            }
        }
        long queryDuration = System.nanoTime() - queryStartNanos;
        stats.recordQueryTiming(queryDuration, totalLpTimeNanos);
        if (stats.isExperimentMode()) {

            int leafProbes = this.verifier.getLeafProbes();
            int bfprobes = this.getAllprobes();
            int actualCost = bfprobes;
//            System.out.println("Pattern: " + pat.text +" Probes: " + bfprobes + " Leafs: " + bfprobes + " Actual: " + actualCost);

            this.verifier.reset();
            stats.setLatestPatternResult(new PatternResult(System.currentTimeMillis() - startTime, actualCost, lp, pat, lpCf, cp_cost, leafProbes, arbitraryConfLp));
        }

        return results;
    }

    public long getAvgBloomCost(Pattern pat) {
        ImplicitTree<Membership> tree = this.trees.getLast();
        int maxLvl = tree.maxDepth() - 1;
        long duration = 0;
        for (int z = 0; z < bcCostEstimIter; z++) {
            long startTime = System.nanoTime();
            long key = tree.codec.pack(maxLvl, 0, pat.nGramToInt[0]);
            tree.contains(maxLvl, key);
            duration += System.nanoTime() - startTime;
        }
        duration /= bcCostEstimIter;
        return duration;
    }

    @Override
    public ArrayList<Long> getAvgTimes(Pattern pat) {
        long avgBloom = getAvgBloomCost(pat);
        long avgLeaf = getAvgLeafCost(pat);
        this.bfCost = avgBloom;
        this.leafCost = avgLeaf;
        ArrayList<Long> res = new ArrayList<>();
        res.add(avgBloom);
        res.add(avgLeaf);
        return res;
    }

    @Override
    public PatternResult getLatestStats() {
        return stats.latestPatternResult();
    }

    public long getAvgLeafCost(Pattern pat) {
        ImplicitTree<Membership> tree = this.trees.getLast();
        long duration = 0;
        for (int i = 0; i < lcCostEstimIter; i++) {
            long startTime = System.nanoTime();
            for (int j = 0; j < pat.nGramToInt.length - 1; j++) {
                boolean in = tree.buffer.data.get(0) == 1;
            }
            duration += System.nanoTime() - startTime;
        }
        duration /= lcCostEstimIter;
        return duration;
    }

    private TreeLayout makeLayout() {
        int levels = Integer.SIZE - Integer.numberOfLeadingZeros(treeLength);
        int leafSpan = treeLength >>> (levels - 1);
        return new TreeLayout(levels, leafSpan);
    }

    private ImplicitTree<Membership> createTree() {
        TreeLayout layout = makeLayout();
        int maxDepth = layout.levels();
        IntFunction<Membership> filterFactory = level -> {
            int nodes = 1 << level;
            int interval = treeLength >> level;
            int perNode = Math.min(alphabetSize, interval);
            int distinct = nodes * perNode;

            Membership bf = membershipFac.get();
            bf.init(distinct, fpRate);
            return bf;
        };

        return new ImplicitTree<>(
                layout,
                filterFactory,
                new Key64(maxDepth, alphabetSize),
                estimatorFac.get());
    }

    public int getAllprobes() {
        int sum = 0;
        for (ImplicitTree<Membership> tree : this.trees) {
            sum += tree.containCounter;
            tree.containCounter = 0;
        }
        return sum;
    }
}
