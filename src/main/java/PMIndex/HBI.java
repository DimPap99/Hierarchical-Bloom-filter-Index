package PMIndex;

import estimators.CostFunction;
import estimators.CostFunctionDefaultRoot;
import estimators.Estimator;
import estimators.HOPS;
import membership.KeyPackingService;
import membership.Membership;
import org.apache.commons.math3.util.Pair;
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
import java.util.function.Predicate;
import java.util.function.Supplier;

import static utilities.MathUtils.pruningLevel;
import static utilities.MathUtils.pruningLevelBloom;


// Hierarchical Bloom filter Index (thesis Ch. 4).
//
// Coordinates: the active trees are numbered 0..trees.size()-1 (oldest first) and a position p inside
// the retained trees is tree.id * treeLength + local. When the oldest tree expires the remaining trees
// are renumbered, so these ints never grow with the stream. report() returns positions relative to the
// start of the current window, i.e. the set { i in 0..w-r : W_t[i..i+r-1] = P } of Sec. 1.2.
// windowStartPosition() converts them to absolute stream positions.
//
// Not thread-safe.
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
    // Keep levels below nGram and use token-length descent; default preserves current behaviour.
    private final boolean keepSubNgramLevels;

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
    // Absolute position of the last indexed symbol.
    private long indexedItemsCounter = -1;
    // Number of trees evicted so far; they all held exactly treeLength symbols.
    private long expiredTrees = 0;

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
    private int lastPredictiveLp = -1;
    public boolean isRootAlg = false;
    public boolean activateEstim = false;
    public Utils.MemPolicy memPolicy = Utils.MemPolicy.NONE; // default is no memory policy.
    public int pctEstimatorBuckets;
    // Markov configuration and snapshot lifecycle.
    private int  desiredMaxOrder = 3;   // 2 means first-order.
    private boolean markovDirty = false; // set true on insert; snapshot rebuilt lazily before query.

    private void ensureMarkovBuilder() {
        if (this.modelBuilder == null) {
            int sigma = Math.max(1, this.alphabetSize);
            // maxOrder = context length + 1.
            this.modelBuilder = new NgramModel.Builder(sigma, 3);
        }
    }

    private boolean estimatorsActive() {
        return this.activateEstim || !this.isRootAlg;
    }

    private boolean memoryPolicyActive() {
        return this.memoryPolicy && estimatorsActive();
    }

    public HBI(HbiConfiguration config) {
        this.config = Objects.requireNonNull(config, "config");
        this.stats = new HbiStats(config.collectStats(), config.experimentMode());
        this.windowLength = config.windowLength();
        this.treeLength = config.treeLength();
        this.alphabetSize = config.alphabetSize();
        this.fpRate = config.fpRate();
        this.keepSubNgramLevels = config.keepSubNgramLevels();
        this.searchAlgo = config.searchAlgorithm();
        if (this.searchAlgo instanceof search.BlockSearch bs) {
            bs.setUseTokenLengthGuard(this.keepSubNgramLevels);
        } else if (this.searchAlgo instanceof search.BlockSearchCharSet bscs) {
            bscs.setUseTokenLengthGuard(this.keepSubNgramLevels);
        }
        this.estimatorFac = config.estimatorSupplier();
        this.membershipFac = config.membershipSupplier();
        this.pruningPlanFac = config.pruningPlanSupplier();
        this.verifier = config.verifier();
        this.cf = config.costFunction();
        this.conf = config.confidence();
        this.keyMapper = new StringKeyMapper(this.alphabetSize, 0.0001);
        this.nGram = config.nGram();
        this.activateEstim = config.activateEstim();
        if(this.cf instanceof CostFunctionDefaultRoot) {
            this.isRootAlg = true;
        }

        Utils.MemPolicy configuredPolicy = config.memPolicy();
        this.memPolicy = configuredPolicy == null ? Utils.MemPolicy.NONE : configuredPolicy;
        this.memoryPolicy = this.memPolicy != Utils.MemPolicy.NONE;
        if (memoryPolicyActive()) {
            this.pctEstimatorBuckets = Math.max(1, this.alphabetSize);

            int configuredBuckets = config.buckets();
            if (configuredBuckets > 0) {
                this.pctEstimatorBuckets = configuredBuckets;
            }
            ensurePctEstimator();
        }
        trees.addLast(createTree(0));
    }

    public HbiStats stats() {
        return stats;
    }

    public void setLpOverride(int lp) {
        this.lpOverride = lp;
    }

    public int lpOverride() {
        return lpOverride;
    }

    public int nGram() {
        return nGram;
    }

    public int windowLength() {
        return windowLength;
    }

    public int treeLength() {
        return treeLength;
    }

    // Number of symbols indexed since construction.
    public long indexedCount() {
        return indexedItemsCounter + 1;
    }

    // Absolute stream position of the first symbol of the current window.
    public long windowStartPosition() {
        return Math.max(0L, indexedItemsCounter - windowLength + 1);
    }

    // Absolute stream position of the first symbol held by the oldest retained tree.
    private long retainedStartPosition() {
        return expiredTrees * (long) treeLength;
    }

    // Evict the oldest tree. The newest tree is never evicted, since it receives the next symbols.
    @Override
    public void expire() {
        if (trees.size() <= 1) {
            return;
        }
        trees.removeFirst();
        expiredTrees++;
        for (ImplicitTree<Membership> tree : trees) {
            tree.id--;
        }
    }

    // Drop whole trees once none of their positions lies inside the window (Sec. 4.1.3).
    private void expireOutOfWindowTrees() {
        long windowStart = windowStartPosition();
        while (trees.size() > 1 && retainedStartPosition() + treeLength <= windowStart) {
            expire();
        }
    }

    // Stream a single symbol (character or n-gram) into the index.
    @Override
    public void insert(String c) {
        insertToken(keyMapper.mapToLong(c));
    }

    // Stream a symbol that the caller has already mapped to a numeric key (Sec. 3.3).
    public void insertToken(long token) {
        if(isMarkov) {
            ensureMarkovBuilder();
            this.modelBuilder.observeSymbol(token); // update counts while indexing
            this.markovDirty = true;                // snapshot must be rebuilt before next query
        }

        indexedItemsCounter++;
        ImplicitTree<Membership> target = trees.getLast();
        if (target.isFull()) {
            // The policies decide on the tree that just filled, before the next tree starts.
            if (memoryPolicyActive()) {
                applyMemoryPolicy();
            }
            target = createTree(trees.size());
            if (this.memPolicy == Utils.MemPolicy.PREDICTIVE && this.lastPredictiveLp > 0) {
                target.dropFiltersUpToLp(this.lastPredictiveLp);
            }
            trees.add(target);
        }

        if (memoryPolicyActive()) {
            ensurePctEstimator().insert(token);
        }
        if (estimatorsActive()) {
            if (target.estimator == null) {
                target.estimator = estimatorFac.get();
            }
            target.estimator.insert(token);
        }
        target.append(token, indexedItemsCounter);

        expireOutOfWindowTrees();
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
        if (!memoryPolicyActive() || this.trees.isEmpty()) {
            return;
        }

        ImplicitTree<Membership> latestTree = this.trees.getLast();
        if (latestTree.estimator == null) {
            return;
        }
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
        // This always happens for the current policy.
        double minp = latestTree.estimator.estimate(this.lastPolicyQuantileKey);

        int lp = pruningLevelBloom(latestTree, 0.99, minp, this.fpRate);
        latestTree.dropFiltersUpToLp(lp);
        if (this.memPolicy == Utils.MemPolicy.PREDICTIVE) {
            this.lastPredictiveLp = Math.max(0, Math.min(lp, latestTree.maxDepth() - 1));
        }
        // Clear samples.
        this.pctEstimator.clear();
        // Additional policies (e.g., PREDICTIVE) can consume lastPolicyQuantileKey as needed.
    }

    // Force memory policy evaluation on the most recent tree.
    public void forceApplyMemoryPolicy() {
        if (!memoryPolicyActive() || this.trees.isEmpty()) {
            return;
        }

        ImplicitTree<Membership> latestTree = this.trees.getLast();
        if (latestTree.indexedItemsCounter < 0 || latestTree.estimator == null) {
            return; // No data or estimator available yet.
        }

        applyMemoryPolicy();
    }


    public long getLastPolicyQuantileKey() {
        return lastPolicyQuantileKey;
    }

    public int getLastPolicyQuantileFrequency() {
        return lastPolicyQuantileFrequency;
    }

    public void releaseEstimatorsMatching(Predicate<Estimator> predicate) {
        if (predicate == null) {
            return;
        }
        for (ImplicitTree<Membership> tree : this.trees) {
            Estimator estimator = tree.estimator;
            if (estimator != null && predicate.test(estimator)) {
                tree.estimator = null;
            }
        }
    }

    @Override
    public boolean exists(String key) {
        return !report(new Pattern(key, this.nGram)).isEmpty();
    }

    @Override
    public ArrayList<Integer> report(Pattern pat) {
        if (pat.nGram != this.nGram) {
            throw new IllegalArgumentException(
                    "pattern uses n-gram size " + pat.nGram + " but the index uses " + this.nGram);
        }
        // Sec. 4.3.1: a pattern needs at least one full n-gram.
        if (pat.originalSz < this.nGram) {
            throw new IllegalArgumentException(
                    "pattern length " + pat.originalSz + " is shorter than the n-gram size " + this.nGram);
        }

        ArrayList<Integer> results = new ArrayList<>();
        long queryStartNanos = System.nanoTime();
        long totalLpTimeNanos = 0L;

        this.searchAlgo.setStrides(this.strides);

        // Precompute numeric tokens for the pattern once
        if (!pat.isConverted) {
            for (int nIdx = 0; nIdx < pat.nGramArr.length; nIdx++) {
                pat.nGramToLong[nIdx] = this.keyMapper.mapToLong(pat.nGramArr[nIdx]);
            }
        }
        pat.buildStrideTokens();

        // Rebuild Markov-aware cost model if needed
        if (isMarkov) {
            refreshCostModelIfNeeded();
        }

        // positionOffset is the search frontier: every start before it is already verified or ruled out.
        // Positions of the oldest tree that slid out of the window are skipped from the start.
        final int windowStart = (int) (windowStartPosition() - retainedStartPosition());
        int positionOffset = windowStart;

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
            int treeLastPos = tree.id * span + tree.buffer.lastIndex();

            if (treeLastPos < positionOffset) {
                continue;
            }

            if (tree.pruningPlan == null) {
                tree.pruningPlan = this.pruningPlanFac.get();
            }

            // Build an IntervalScanner seeded with the current frontier
            IntervalScanner scn = new IntervalScanner(tree, pat, searchAlgo, positionOffset);

            long lpStart = System.nanoTime();
            ArrayList<Integer> lps = new ArrayList<>();

            if (stats.isExperimentMode()) {
                // estimator.estimateALl gives us per-symbol estimated probabilities
                double[] pp = tree.estimator.estimateALl(pat, this.strides);
                double pMin = Arrays.stream(pp).min().orElse(0.0);

                // lp from override or arbitrary confidence
                lp = this.lpOverride;
                arbitraryConfLp = pruningLevel(tree, this.conf, pMin);

                // cost of committing to a chosen level
                cp_cost = cf.costAtLevel(tree, pat, lp, this.strides);

                long minCostStart = System.nanoTime();
                lpCf = cf.minCostLp(tree, 0.05, pat, 97, 26, this.strides);
                long cfDur = System.nanoTime() - minCostStart;
                stats.recordMinCostLpTime(cfDur);
                if (stats.isCollecting()) {
                    stats.recordCfLp(lpCf);
                }

                lps.add(lp);
            } else {
                if (this.cf != null) {
                    long cfStart = System.nanoTime();
                    int chosen = cf.minCostLp(tree,
                            0.95,
                            pat,
                            this.bfCost,
                            this.leafCost,
                            this.strides);
                    long cfDur = System.nanoTime() - cfStart;
                    lps.add(chosen);
                    // Record CF minCost timing if stats collection is enabled
                    stats.recordMinCostLpTime(cfDur);
                    if (stats.isCollecting()) {
                        stats.recordCfLp(chosen);
                    }
                } else {
                    if (estimatorsActive() && tree.pruningPlan != null) {
                        lps = tree.pruningPlan.pruningPlan(pat, tree, this.conf, this.strides);
                    } else {
                        lps.add(0);//defaults to root/vanilla bs
                    }
                }
            }

            totalLpTimeNanos += System.nanoTime() - lpStart;
            pat.charStartLp = lps;
            lp = lps.isEmpty() ? tree.effectiveRoot() : Collections.min(lps);

            if (stats.isCollecting()) {
                stats.recordLp(lp);
                if (this.cf != null) {
                    stats.recordAlpha(cf.getAlpha());
                }
            }

            // Seed only levels that still have filters: a policy may have dropped the top of the tree.
            int seedLevel = Math.max(tree.effectiveRoot(), Math.min(lp, tree.maxDepth() - 1));
            scn.seedLevel(seedLevel);

            // Walk candidate intervals from this tree
            while (scn.hasNext()) {
                CandidateRange cr = scn.next();
                if (cr == null) break;

                // Verify this candidate range, possibly crossing tree boundaries
                Pair<ArrayList<Integer>, Integer> res =
                        this.verifier.verify(treeIdx, this.trees, cr, pat);

                // The frontier only moves forward, so every hit beyond the last one is new.
                for (int hit : res.getFirst()) {
                    if (hit < positionOffset) {
                        continue;
                    }
                    int windowPos = hit - windowStart;
                    if (results.isEmpty() || results.get(results.size() - 1) < windowPos) {
                        results.add(windowPos);
                    }
                }

                positionOffset = Math.max(positionOffset, res.getSecond());
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
        int target = this.keepSubNgramLevels ? 1 : Math.max(1, nGram);

        while (levels < totalLevels && span % 2 == 0 && span / 2 >= target) {
            span /= 2;
            levels++;
        }

        int leafSpan = Math.max(1, span);
        return new TreeLayout(levels, leafSpan);

    }

    private ImplicitTree<Membership> createTree(int treeId) {

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
        if (estimatorsActive()) {
            est = estimatorFac.get();
        }

        ImplicitTree<Membership> tree = new ImplicitTree<>(
                layout,
                filterFactory,
                new KeyPackingService(maxDepth, alphabetSize),
                est
        );

        // Slot of the tree among the retained trees; see the coordinate note at the top of the class.
        tree.id = treeId;

        if (estimatorsActive()) {
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
