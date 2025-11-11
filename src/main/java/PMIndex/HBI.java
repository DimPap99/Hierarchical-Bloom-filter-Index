package PMIndex;

import estimators.CostFunction;
import estimators.Estimator;
import estimators.HOPS;
import jdk.jshell.execution.Util;
import membership.Key64;
import membership.KeyPackingService;
import membership.Membership;
import org.apache.commons.math3.util.Pair;
import org.openjdk.jol.info.GraphLayout;
import org.openjdk.jol.vm.VM;
import search.CandidateRange;
import search.Frame;
import search.IntervalScanner;
import search.Pattern;
import search.PruningPlan;
import search.SearchAlgorithm;
import search.Verifier;
import tree.ImplicitTree;
import tree.StreamBuffer;
import tree.TreeLayout;
import utilities.*;

import java.util.*;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import static utilities.MathUtils.pruningLevel;


public final class HBI implements IPMIndexing {


    private static final int DEFAULT_BC_COST_ESTIM_ITER = 1_000_000;
    private static final int DEFAULT_LC_COST_ESTIM_ITER = 1_000_000;

    private final HbiConfiguration config;

    private int nGram = 1;
    public boolean strides = false;
    private final HbiStats stats;

    private final int windowLength;
    private final int treeLength;
    public int alphabetSize;
    public final double fpRate;

    private final SearchAlgorithm searchAlgo;
    private final Supplier<Estimator> estimatorFac;
    private final Supplier<Membership> membershipFac;
    private final Supplier<PruningPlan> pruningPlanFac;
    private final Verifier verifier;
    private final CostFunction cf;
    private boolean builtModel = false;
    private final double conf;

    public final ArrayList<ImplicitTree<Membership>> trees = new ArrayList<>();

    public StringKeyMapper keyMapper;
    private long indexedItemsCounter = -1;

    private int bcCostEstimIter = DEFAULT_BC_COST_ESTIM_ITER;
    private int lcCostEstimIter = DEFAULT_LC_COST_ESTIM_ITER;

    public HOPS pctEstimator = null;
    private double bfCost = 97;
    private double leafCost = 26;
    private int lpOverride;
    public NgramModel.Builder modelBuilder;
    public boolean memoryPolicy = false;
    public boolean isMarkov = false;

    public Utils.MemPolicy memPolicy = Utils.MemPolicy.NONE;//default policy is no policy

    private final int maxActiveTrees;
    private final HashSet<String> strhs = new HashSet<>();
    private final HashSet<Integer> assignedkeys = new HashSet<Integer>();
    // Markov config + snapshot lifecycle
    private int  desiredMaxOrder = 3;   // 2 => first-order
    private boolean markovDirty = false; // set true on insert; snapshot rebuilt lazily before query

    private void ensureMarkovBuilder() {
        if (this.modelBuilder == null) {
            int sigma = Math.max(1, this.alphabetSize);
            // maxOrder = context length + 1
            this.modelBuilder = new NgramModel.Builder(sigma, 3);
        }
    }

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
        this.keyMapper = new StringKeyMapper(this.alphabetSize, 0.0001);
        this.nGram = config.nGram();
        trees.addLast(createTree());
        this.maxActiveTrees = (int) Math.ceil((double) (windowLength / treeLength));
        if(isMarkov) {
            ensureMarkovBuilder();
        }

        this.memPolicy = config.memPolicy();
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
               double conf,
               int nGram) {
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
                .nGram(nGram)
                .build());
    }

    public HbiStats stats() {
        return stats;
    }

    public String costFunctionName() {
        return (cf == null) ? "none" : cf.getClass().getSimpleName();
    }

//    public void setAlphabetMap(AlphabetMapper<String> alphabetMap) {
//        this.alphabetMap = Objects.requireNonNull(alphabetMap, "alphabetMap");
//    }

//    public void resetAlphabetMap(int capacity) {
//        this.alphabetMap = new AlphabetMapper<>(capacity);
//    }

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
        long token = keyMapper.mapToLong(c);
//        this.strhs.add(c);
//        this.assignedkeys.add(intC);
//        if(intC == 146) System.out.println(c);
        if(isMarkov) {
            ensureMarkovBuilder();
            this.modelBuilder.observeSymbol(token); // update counts while indexing
            this.markovDirty = true;                // snapshot must be rebuilt before next query
        }

        indexedItemsCounter++;
        ImplicitTree<Membership> lastTree = trees.getLast();
        if (lastTree.isFull()) {
            ImplicitTree<Membership> fresh = createTree();
            fresh.id = trees.size();
            fresh.estimator.insert(token);
            fresh.append(token, indexedItemsCounter);
            trees.add(fresh);
//            alphabetSize = this.alphabetMap.getSize() + (int)(this.alphabetMap.getSize()*0.1);
        } else {
            lastTree.estimator.insert(token);
            lastTree.append(token, indexedItemsCounter);
        }

        if (trees.getLast().indexedItemsCounter + (trees.size() - 1) * this.treeLength > this.windowLength) {

            this.expire();
        }
    }

    public void applyMemoryPolicy() {
        if(memoryPolicy) {
            double minp = this.trees.getLast().estimator.getMin();
            int lp = pruningLevel(this.trees.getLast(), 0.95, minp);
            this.trees.getLast().dropFiltersUpToLp(lp);

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
        if (pat == null || pat.nGramArr == null || pat.nGramArr.length == 0) {
            recordDegeneratePattern(pat);
            return results;
        }

        long queryStartNanos = System.nanoTime();
        long totalLpTimeNanos = 0L;
        this.searchAlgo.setStrides(this.strides);
        for (int nIdx = 0; nIdx < pat.nGramArr.length; nIdx++) {
            long tokenVal = this.keyMapper.mapToLong(pat.nGramArr[nIdx]);
            pat.nGramToLong[nIdx] = tokenVal;
            // DO NOT register symbols or mutate the model at query time
            // if (isMarkov) { ensureMarkovBuilder(); this.modelBuilder.ensureSymbolRegistered(tokenVal); }
            if(nIdx % pat.nGram == 0) pat.effectiveNgramArr[nIdx/pat.nGram] = tokenVal;
        }
        if(pat.originalSz % pat.nGram != 0 && pat.effectiveNgramArr.length > 0 && pat.nGramToLong.length > 0) {
            pat.effectiveNgramArr[pat.effectiveNgramArr.length-1] = pat.nGramToLong[pat.nGramToLong.length-1];
        }

        // Lazily rebuild the immutable snapshot once (if new data arrived), then reuse
        if (isMarkov) {
            refreshCostModelIfNeeded();
        }

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

            long lpStart = System.nanoTime();
            ArrayList<Integer> lps = new ArrayList<>();
//            lp = Collections.min(lps);

            //cf.minCostLp(tree, 0.05, pat, 97, 26);//pruningLevel(tree, 0.99, pMax);
//
            if (stats.isExperimentMode()) {
                double[] pp = tree.estimator.estimateALl(pat, strides);
                double pMax = Arrays.stream(pp).min().getAsDouble();
                pp = tree.estimator.estimateALl(pat, this.strides);
                pMax = Arrays.stream(pp).min().getAsDouble();
                lp = this.lpOverride;
                arbitraryConfLp = pruningLevel(tree, this.conf, pMax);
                int m = (int) (tree.maxDepth() - 1 - Math.ceil(Math.log(pat.nGramToLong.length) / Math.log(2)));
                cp_cost = cf.costAtLevel(tree, pp, pat.nGramToLong, lp, 0.0, m);
                long minCostStart = System.nanoTime();
                lpCf = cf.minCostLp(tree, 0.05, pat, 97, 26, this.strides);
                stats.recordMinCostLpTime(System.nanoTime() - minCostStart);
                lps.add(lp);
            } else {
//                int s = cf.minCostLp(tree, 0.05, pat, 97, 26, false);
//                lps.add(cf.minCostLp(tree, 0.05, pat, 97, 26, false));

                if(this.cf != null) {
                    lps.add(cf.minCostLp(tree, 0.95, pat, this.bfCost, this.leafCost, this.strides));
                }
                else{
                    lps= tree.pruningPlan.pruningPlan(pat, tree, 0.99, this.strides);
                }
//                lp = lpCf;
            }
            totalLpTimeNanos += System.nanoTime() - lpStart;
            pat.charStartLp = lps;
            lp = Collections.min(lps);
            if (stats.isCollecting()) {
                stats.recordLp(lp);
                if(this.cf != null) {                stats.recordAlpha(cf.getAlpha());}
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

        if (stats.isExperimentMode()) {

            int leafProbes = this.verifier.getLeafProbes();
            int bfprobes = this.getAllprobes();
            int actualCost = bfprobes;

            this.verifier.reset();
            stats.setLatestPatternResult(new PatternResult(System.currentTimeMillis() - startTime, actualCost, lp, pat, lpCf, cp_cost, leafProbes, arbitraryConfLp));
        }
        long totalQueryTimeNanos = System.nanoTime() - queryStartNanos;
        if (stats.isCollecting()) {
            stats.recordQueryTiming(totalQueryTimeNanos, totalLpTimeNanos);
        }

        return results;
    }

    private void recordDegeneratePattern(Pattern pat) {
        if (stats != null) {
            stats.setLatestPatternResult(new PatternResult(0.0, 0, 0, pat, 0, 0.0, 0, 0));
        }
    }

    private void refreshCostModelIfNeeded() {
        if (!isMarkov) return;
        ensureMarkovBuilder();
        if (this.markovDirty) {
            this.cf.setModel(this.modelBuilder.build()); // consume prebuilt matrices
            this.markovDirty = false;
            this.builtModel = true;
        }
    }

    public long getAvgBloomCost(Pattern pat) {
        if (pat == null || pat.nGramToLong == null || pat.nGramToLong.length == 0) {
            return 0L;
        }
        ImplicitTree<Membership> tree = this.trees.getLast();
        int maxLvl = tree.maxDepth() - 1;
        long duration = 0;
        for (int z = 0; z < bcCostEstimIter; z++) {
            long startTime = System.nanoTime();
            long key = tree.codec.packWord(0, pat.nGramToLong[0]);
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

    @Override
    public int getTokenId(String key) {
        return this.keyMapper.mapToInt(key);
    }

    public long getAvgLeafCost(Pattern pat) {
        if (pat == null || pat.nGramToLong == null || pat.nGramToLong.length == 0) {
            return 0L;
        }
        ImplicitTree<Membership> tree = this.trees.getLast();
        long duration = 0;
        for (int i = 0; i < lcCostEstimIter; i++) {
            long startTime = System.nanoTime();
            StreamBuffer buffer = tree.buffer;
            for (int j = 0; j < pat.nGramToLong.length - 1; j++) {
                if (buffer.length() == 0) {
                    break;
                }
                boolean in = buffer.get(0) == 1;
            }
            duration += System.nanoTime() - startTime;
        }
        duration /= lcCostEstimIter;
        return duration;
    }

    private TreeLayout makeLayout() {
        int totalLevels = Integer.SIZE - Integer.numberOfLeadingZeros(treeLength);
        int levels = 1;
        int span = Math.max(1, treeLength);
        int target = Math.max(1, nGram);

        while (levels < totalLevels && span % 2 == 0 && span / 2 >= target) {
            span /= 2;
            levels++;
        }

        int leafSpan = Math.max(1, span);
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
                new KeyPackingService(maxDepth, alphabetSize),
                //new Key64(maxDepth, alphabetSize),
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
