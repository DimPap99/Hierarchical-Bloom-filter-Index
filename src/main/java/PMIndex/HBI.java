package PMIndex;

import estimators.CostFunction;
import estimators.CostFunctionDefaultRoot;
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

    // Monotonic global id for new trees
    private int nextTreeId = 0;

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
    public double QUANTILE = 0.05;
    private long lastPolicyQuantileKey = -1L;
    private int lastPolicyQuantileFrequency = 0;
    public boolean isRootAlg = false;
    public Utils.MemPolicy memPolicy = Utils.MemPolicy.NONE;//default policy is no policy
    public int pctEstimatorBuckets;
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
        if(this.cf instanceof CostFunctionDefaultRoot) {
            this.isRootAlg = true;
        }
        trees.addLast(createTree(this.nextTreeId++));
        this.maxActiveTrees = (int) Math.ceil((double) (windowLength / treeLength));
        if(isMarkov) {
            ensureMarkovBuilder();
        }

        Utils.MemPolicy configuredPolicy = config.memPolicy();
        this.memPolicy = configuredPolicy == null ? Utils.MemPolicy.NONE : configuredPolicy;
        this.memoryPolicy = this.memPolicy != Utils.MemPolicy.NONE;
        if(this.isRootAlg){
            this.memPolicy = Utils.MemPolicy.NONE;
            this.memoryPolicy = false;
        }
        if (this.memoryPolicy) {
            this.pctEstimatorBuckets = Math.max(1, this.alphabetSize);

            int configuredBuckets = config.buckets();
            if (configuredBuckets > 0) {
                this.pctEstimatorBuckets = configuredBuckets;
            }
            ensurePctEstimator();
        }
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
        if (!this.isRootAlg && this.memPolicy != Utils.MemPolicy.NONE) {
            this.pctEstimator.insert(token);
        }
        if(isMarkov) {
            ensureMarkovBuilder();
            this.modelBuilder.observeSymbol(token); // update counts while indexing
            this.markovDirty = true;                // snapshot must be rebuilt before next query
        }

        indexedItemsCounter++;
        ImplicitTree<Membership> lastTree = trees.getLast();
        if (lastTree.isFull()) {

            if(this.memoryPolicy && !this.isRootAlg){
                applyMemoryPolicy();
//                if(this.memPolicy == Utils.MemPolicy.PREDICTIVE){
//
//                }
            }

            ImplicitTree<Membership> fresh = createTree(nextTreeId++);
            if(!isRootAlg) {
                fresh.estimator.insert(token);
            }
            fresh.append(token, indexedItemsCounter);
            trees.add(fresh);
//            alphabetSize = this.alphabetMap.getSize() + (int)(this.alphabetMap.getSize()*0.1);
        } else {
            if(!isRootAlg) { lastTree.estimator.insert(token);}
            lastTree.append(token, indexedItemsCounter);
        }

        if (trees.getLast().indexedItemsCounter + (trees.size() - 1) * this.treeLength > this.windowLength) {

            this.expire();
        }
    }

    private HOPS ensurePctEstimator() {
        if (this.pctEstimator == null) {
            this.pctEstimator = createPctEstimator();
        }
        return this.pctEstimator;
    }

    private HOPS createPctEstimator() {
        int buckets = this.pctEstimatorBuckets;
        long seedBase = 0x9E3779B97F4A7C15L;
        seedBase = seedBase * 31 + (long) this.windowLength;
        seedBase = seedBase * 31 + (long) this.treeLength;
        seedBase = seedBase * 31 + (long) this.alphabetSize;
        SplittableRandom seedGen = new SplittableRandom(seedBase);
        long bucketSeed = seedGen.nextLong();
        long prioritySeed = seedGen.nextLong();
        return new HOPS(buckets, this.alphabetSize, bucketSeed, prioritySeed);
    }

    public void applyMemoryPolicy() {
        if (!memoryPolicy || this.trees.isEmpty()) {
            return;
        }

        ImplicitTree<Membership> latestTree = this.trees.getLast();
        HOPS hopsSampler = ensurePctEstimator();

        HOPS.QuantileEstimate estimate = hopsSampler.estimateQuantileWithKey(key -> {
            long freq = latestTree.estimator.get(key);
            if (freq <= 0) {
                return 0;
            }
            if (freq >= Integer.MAX_VALUE) {
                return Integer.MAX_VALUE;
            }
            return (int) freq;
        }, this.QUANTILE);

        this.lastPolicyQuantileKey = estimate.key;
        this.lastPolicyQuantileFrequency = estimate.frequency;
        //this always happens
        double minp = latestTree.estimator.estimate(this.lastPolicyQuantileKey);
        int lp = pruningLevel(latestTree, 0.99, minp);
        latestTree.dropFiltersUpToLp(lp);
        //clear samples
        this.pctEstimator.clear();
        // Additional policies (e.g., PREDICTIVE) can consume lastPolicyQuantileKey as needed.
    }

    public long getLastPolicyQuantileKey() {
        return lastPolicyQuantileKey;
    }

    public int getLastPolicyQuantileFrequency() {
        return lastPolicyQuantileFrequency;
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

        this.searchAlgo.setStrides(this.strides);

        // Precompute numeric tokens for the pattern once
        for (int nIdx = 0; nIdx < pat.nGramArr.length; nIdx++) {
            long tokenVal = this.keyMapper.mapToLong(pat.nGramArr[nIdx]);
            pat.nGramToLong[nIdx] = tokenVal;
            if (nIdx % pat.nGram == 0) {
                pat.effectiveNgramArr[nIdx / pat.nGram] = tokenVal;
            }
        }

        // Handle a tail n-gram for stride mode if pattern length is not divisible by nGram
        if (pat.effectiveNgramArr.length > 0 && pat.originalSz % pat.nGram != 0) {
            long tailToken = (pat.nGramToLong.length > 0)
                    ? pat.nGramToLong[pat.nGramToLong.length - 1]
                    : this.keyMapper.mapToLong(pat.patternTxt);
            pat.effectiveNgramArr[pat.effectiveNgramArr.length - 1] = tailToken;
        }

        // Rebuild Markov-aware cost model if needed
        if (isMarkov) {
            refreshCostModelIfNeeded();
        }

        // positionOffset is the global frontier in absolute coordinates
        // Any match strictly before this offset has already been verified
        int positionOffset = -1;

        long startTimeMs = System.currentTimeMillis();
        int lp = 0;
        int lpCf = 0;
        double cp_cost = 0.0;
        int arbitraryConfLp = 0;

        // We iterate trees in ascending time order (oldest to newest)
        for (int treeIdx = 0; treeIdx < this.trees.size(); treeIdx++) {
            ImplicitTree<Membership> tree = this.trees.get(treeIdx);

            // If we have already advanced our frontier past the last symbol of this tree,
            // skip it completely
            int span = tree.baseIntervalSize();
            int treeLastGlobal =
                    tree.id * span
                            + tree.buffer.lastIndex(); // last local index in this tree

            if (treeLastGlobal < positionOffset) {
                continue;
            }

            if (tree.pruningPlan == null) {
                tree.pruningPlan = this.pruningPlanFac.get();
            }

            // Build an IntervalScanner seeded with the current global frontier
            IntervalScanner scn = new IntervalScanner(tree, pat, searchAlgo, positionOffset);

            long lpStart = System.nanoTime();
            ArrayList<Integer> lps = new ArrayList<>();

            if (stats.isExperimentMode()) {
                // estimator.estimateALl gives us per-symbol estimated probabilities
                double[] pp = tree.estimator.estimateALl(pat, this.strides);
                double pMax = Arrays.stream(pp).min().getAsDouble();

                // lp from override or arbitrary confidence
                lp = this.lpOverride;
                arbitraryConfLp = pruningLevel(tree, this.conf, pMax);

                // cost of committing to a chosen level
                int m = (int) (tree.maxDepth() - 1
                        - Math.ceil(Math.log(pat.nGramToLong.length) / Math.log(2)));

                cp_cost = cf.costAtLevel(tree, pp, pat.nGramToLong, lp, 0.0, m);

                long minCostStart = System.nanoTime();
                lpCf = cf.minCostLp(tree, 0.05, pat, 97, 26, this.strides);
                stats.recordMinCostLpTime(System.nanoTime() - minCostStart);

                lps.add(lp);
            } else {
                if (this.cf != null) {
                    lps.add(cf.minCostLp(tree,
                            0.95,
                            pat,
                            this.bfCost,
                            this.leafCost,
                            this.strides));
                } else {
                    if (!this.isRootAlg) {
                        lps = tree.pruningPlan.pruningPlan(pat, tree, 0.99, this.strides);
                    } else {
                        lps.add(0);
                    }
                }
            }

            totalLpTimeNanos += System.nanoTime() - lpStart;
            pat.charStartLp = lps;
            lp = Collections.min(lps);

            if (stats.isCollecting()) {
                stats.recordLp(lp);
                if (this.cf != null) {
                    stats.recordAlpha(cf.getAlpha());
                }
            }

            scn.seedLevel(lp);

            // Walk candidate intervals from this tree
            while (scn.hasNext()) {
                CandidateRange cr = scn.next();
                if (cr == null) break;

                // Verify this candidate range, possibly crossing tree boundaries
                Pair<ArrayList<Integer>, Integer> res =
                        this.verifier.verify(treeIdx, this.trees, cr, pat);

                // Append matches, while removing immediate duplicates
                // Since verify() always advances forward, matches are nondecreasing
                for (int hit : res.getFirst()) {
                    if (results.isEmpty() || results.get(results.size() - 1) != hit) {
                        results.add(hit);
                    }
                }

                // Advance the global frontier
                positionOffset = res.getSecond();
                scn.positionOffset = positionOffset;
            }
        }

        // Experiment mode bookkeeping
        if (stats.isExperimentMode()) {
            int leafProbes = this.verifier.getLeafProbes();
            int bfprobes   = this.getAllprobes();
            int actualCost = bfprobes;

            this.verifier.reset();

            stats.setLatestPatternResult(
                    new PatternResult(
                            System.currentTimeMillis() - startTimeMs,
                            actualCost,
                            lp,
                            pat,
                            lpCf,
                            cp_cost,
                            leafProbes,
                            arbitraryConfLp
                    )
            );
        }

        if (stats.isCollecting()) {
            long totalQueryTimeNanos = System.nanoTime() - queryStartNanos;
            stats.recordQueryTiming(totalQueryTimeNanos, totalLpTimeNanos);
        }

        return results;
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

    private ImplicitTree<Membership> createTree(int treeGlobalId) {

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

        Estimator est = null;
        if (!this.isRootAlg) {
            est = estimatorFac.get();
        }

        ImplicitTree<Membership> tree = new ImplicitTree<>(
                layout,
                filterFactory,
                new KeyPackingService(maxDepth, alphabetSize),
                est
        );

        // assign the globally unique id right here
        tree.id = treeGlobalId;

        if (!isRootAlg) {
            tree.pruningPlan = pruningPlanFac.get();
        }

        return tree;
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
