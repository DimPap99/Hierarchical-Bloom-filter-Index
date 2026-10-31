package estimators;

import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import search.Pattern;
import tree.ImplicitTree;
import utilities.MathUtils;

import java.util.Arrays;

// Start-level cost model for Block Search that accounts for the pattern's own occurrences.
// Derivation and evaluation: docs/cost-function-occurrence.md.
//
// The biased estimator (CostFunctionMaxProb) multiplies per-symbol positive probabilities, as if
// the pattern's n-grams were placed in the window independently of each other. They are not: the
// n-grams of one occurrence sit next to each other, so a node that holds one of them usually holds
// all of them. On sparse alphabets (large n, packet traces) that co-location is the main reason a
// node passes, and the product predicts almost no work below the start level.
//
// Here a node holds the first m distinct pattern symbols in one of two ways:
//   signal:     an occurrence of the pattern prefix lies entirely inside the node, or
//   background: the symbols landed in the node independently (the iid model of the thesis),
// so G_m = S_m + (1 - S_m) * B_m. The prefix count behind S_m comes from the same frequency
// sketch, so the model stores nothing beyond what the index already keeps.
public class CostFunctionOccurrence extends AbstractCostFunction {

    // Below this the conditional count lambda / (1 - e^-lambda) equals 1 in double precision.
    private static final double TINY_LAMBDA = 1e-12;
    // The background is rebuilt only when the prefix count moves by more than this fraction.
    private static final double COUNT_TOLERANCE = 1e-3;

    private final double positivePrior;

    private Plan cachedPlan;
    private Pattern cachedPattern;
    private ImplicitTree<?> cachedTree;
    private int cachedRecords = -1;
    private boolean cachedStrides;

    // Realized false positive rate per level. BloomFilter.getFpRate counts every set bit, so it is
    // read once per tree state instead of once per level per query.
    private double[] betaByLevel;
    private ImplicitTree<?> betaTree;
    private int betaRecords = -1;

    public CostFunctionOccurrence() {
        this(1.0);
    }

    // positivePrior is the probability that a query whose symbols are all present in the tree
    // really occurs there. 1.0 matches workloads drawn from the window.
    public CostFunctionOccurrence(double positivePrior) {
        if (!(positivePrior >= 0.0 && positivePrior <= 1.0)) {
            throw new IllegalArgumentException("positivePrior must be in [0, 1]");
        }
        this.positivePrior = positivePrior;
    }

    @Override
    public int minCostLp(ImplicitTree tree,
                         double confInit,
                         Pattern pattern,
                         double bfCost,
                         double leafCost,
                         boolean strides) {
        Plan plan = plan(tree, pattern, strides);
        int lo = Math.min(tree.effectiveRoot(), plan.maxLevel);
        int hi = Math.max(lo, plan.deepest);
        int bestLp = lo;
        double bestCost = Double.POSITIVE_INFINITY;
        for (int lp = lo; lp <= hi; lp++) {
            double cost = plan.cost(lp);
            if (cost < bestCost) {
                bestCost = cost;
                bestLp = lp;
            }
        }
        this.predictedBloomProbeCost = bestCost;
        return bestLp;
    }

    @Override
    public double costAtLevel(ImplicitTree<?> tree, Pattern pattern, int lp, boolean strides) {
        return plan(tree, pattern, strides).cost(lp);
    }

    // Without the pattern there is no character length or stride layout, so the keys are read as
    // one overlapping n-gram per position and Block Search is assumed to descend while a child
    // spans at least as many positions as there are keys.
    @Override
    public double costAtLevel(ImplicitTree<?> tree,
                              double[] probs,
                              long[] keySeq,
                              int lp,
                              double bloomFp,
                              int stopLp) {
        if (keySeq == null || keySeq.length == 0) {
            return 0.0;
        }
        int[] offsets = new int[keySeq.length];
        for (int i = 0; i < offsets.length; i++) {
            offsets[i] = i;
        }
        return build(tree, probs, keySeq, offsets, keySeq.length, 1, false).cost(lp);
    }

    private Plan plan(ImplicitTree<?> tree, Pattern pattern, boolean strides) {
        int records = tree.buffer.length();
        if (cachedPlan != null
                && cachedPattern == pattern
                && cachedTree == tree
                && cachedRecords == records
                && cachedStrides == strides) {
            return cachedPlan;
        }
        long[] keys = strides ? pattern.effectiveNgramArr : pattern.nGramToLong;
        double[] probs = tree.estimator.estimateALl(pattern, strides);
        int[] offsets = keyOffsets(pattern, keys.length, strides);
        // BlockSearch descends while a child spans the whole pattern in characters (or tokens).
        int requiredSpan = Math.max(1, pattern.size);
        Plan plan = build(tree, probs, keys, offsets, requiredSpan, Math.max(1, pattern.nGram), strides);
        this.cachedPlan = plan;
        this.cachedPattern = pattern;
        this.cachedTree = tree;
        this.cachedRecords = records;
        this.cachedStrides = strides;
        return plan;
    }

    // Position of each tested key inside the pattern's n-gram sequence. With strides the keys sit at
    // 0, n, 2n, ... and the last one is moved to the final n-gram when the length is not a multiple
    // of n, exactly as Pattern.buildStrideTokens lays them out.
    private static int[] keyOffsets(Pattern pattern, int keyCount, boolean strides) {
        int[] offsets = new int[keyCount];
        int n = Math.max(1, pattern.nGram);
        for (int i = 0; i < keyCount; i++) {
            offsets[i] = strides ? i * n : i;
        }
        if (strides && keyCount > 0 && pattern.originalSz % n != 0) {
            offsets[keyCount - 1] = pattern.nGramToLong.length - 1;
        }
        return offsets;
    }

    private Plan build(ImplicitTree<?> tree,
                       double[] probs,
                       long[] keys,
                       int[] offsets,
                       int requiredSpan,
                       int nGram,
                       boolean strides) {
        int maxLevel = Math.max(0, tree.maxDepth() - 1);
        Model model = new Model(tree, probs, keys, offsets, nGram, strides, fpRates(tree, maxLevel),
                positivePrior);
        int deepest = Math.min(MathUtils.deepestVisitedLevel(model.width, requiredSpan), maxLevel);

        // Only the levels Block Search can visit are priced up front; a start below them (the
        // experiment harness forces those) is priced on demand.
        double[][] held = model.held(0, deepest);
        double[] start = new double[maxLevel + 1];
        Arrays.fill(start, Double.NaN);
        for (int level = 0; level <= deepest; level++) {
            start[level] = model.startCost(level, held[level]);
        }
        // Children tested at level L are the two children of every node that passes at L - 1.
        double[] below = new double[maxLevel + 2];
        for (int level = deepest; level >= 1; level--) {
            double[] parent = held[level - 1];
            double passProb = parent.length == 0 ? 1.0 : parent[parent.length - 1];
            double childTests = model.tests(held[level], parent, model.tested(level));
            below[level] = below[level + 1] + 2.0 * model.filled(level - 1) * passProb * childTests;
        }
        return new Plan(model, start, below, deepest, maxLevel);
    }

    private double[] fpRates(ImplicitTree<?> tree, int maxLevel) {
        int records = tree.buffer.length();
        if (betaByLevel == null || betaTree != tree || betaRecords != records
                || betaByLevel.length != maxLevel + 1) {
            double[] rates = new double[maxLevel + 1];
            for (int level = 0; level <= maxLevel; level++) {
                rates[level] = MathUtils.clamp01(tree.getMembershipFpRate(level));
            }
            betaByLevel = rates;
            betaTree = tree;
            betaRecords = records;
        }
        return betaByLevel;
    }

    private static int[] firstOccurrences(long[] keys) {
        LongOpenHashSet seen = new LongOpenHashSet(Math.max(4, keys.length * 2));
        int[] buffer = new int[keys.length];
        int distinct = 0;
        for (int i = 0; i < keys.length; i++) {
            if (seen.add(keys[i])) {
                buffer[distinct++] = i;
            }
        }
        return Arrays.copyOf(buffer, distinct);
    }

    // Per-pattern quantities that do not depend on the level, and G_m(L) for any level range.
    private static final class Model {
        // A background product below this is treated as zero once it can only shrink further.
        private static final double NEGLIGIBLE = 1e-15;

        final int width;
        final double records;
        private final int keyCount;
        private final int nGram;
        private final boolean strides;
        private final double[] probs;
        private final double[] beta;
        private final int[] first;
        private final double[] counts;
        private final double[] spans;
        // Count used for the residual background at m, and where that count changes.
        private final double[] countUsed;
        private final boolean[] rebuild;
        private int lastRebuild = -1;

        Model(ImplicitTree<?> tree,
              double[] probs,
              long[] keys,
              int[] offsets,
              int nGram,
              boolean strides,
              double[] beta,
              double positivePrior) {
            this.width = tree.baseIntervalSize();
            this.records = Math.max(1, tree.buffer.length());
            this.keyCount = keys.length;
            this.nGram = nGram;
            this.strides = strides;
            this.probs = probs;
            this.beta = beta;
            this.first = firstOccurrences(keys);
            int distinct = first.length;
            this.counts = new double[distinct];
            this.spans = new double[distinct];
            prefixCounts(offsets, positivePrior);
            this.countUsed = new double[distinct];
            this.rebuild = new boolean[distinct];
            double used = Double.NaN;
            for (int m = 0; m < distinct; m++) {
                double c = counts[m];
                if (Double.isNaN(used) || Math.abs(c - used) > COUNT_TOLERANCE * Math.max(1.0, used)) {
                    used = c;
                    rebuild[m] = true;
                    lastRebuild = m;
                }
                countUsed[m] = used;
            }
        }

        // c_m = E[count of the prefix | it occurs], from the iid chance count lambda of that prefix.
        // The window must hold every occurrence of the prefix, so c_m is capped by its rarest key.
        private void prefixCounts(int[] offsets, double positivePrior) {
            double logLambda = 0.0;
            double rarest = Double.POSITIVE_INFINITY;
            boolean absent = false;
            int key = 0;
            for (int m = 0; m < first.length; m++) {
                int last = first[m];
                for (; key <= last; key++) {
                    double p = MathUtils.clamp01(probs[key]);
                    if (p <= 0.0) {
                        absent = true;
                    } else {
                        logLambda += Math.log(p);
                        rarest = Math.min(rarest, p * records);
                    }
                }
                double span = offsets[last] + 1.0;
                spans[m] = span;
                if (absent) {
                    counts[m] = 0.0;
                    continue;
                }
                double starts = records - span + 1.0;
                double lambda = starts > 0.0 ? Math.exp(Math.log(starts) + logLambda) : 0.0;
                double givenPresent = lambda < TINY_LAMBDA ? 1.0 : lambda / -Math.expm1(-lambda);
                double c = positivePrior * givenPresent + (1.0 - positivePrior) * lambda;
                counts[m] = Math.min(c, rarest);
            }
        }

        // Keys Block Search tests at one node of this level.
        int tested(int level) {
            int nodeWidth = width >> level;
            int fit = (!strides || nGram <= 1) ? nodeWidth : (nodeWidth + nGram - 1) / nGram;
            return Math.min(keyCount, fit);
        }

        // Nodes of this level that hold data; the rest of the tree is still empty.
        double filled(int level) {
            double nodes = (double) (1L << level);
            return Math.min(nodes, Math.ceil(records / (double) (width >> level)));
        }

        // Expected tests over all seeded nodes when the search starts here. Empty nodes cost one
        // negative test each.
        double startCost(int level, double[] heldAtLevel) {
            double nodes = (double) (1L << level);
            double filled = filled(level);
            return filled * tests(heldAtLevel, null, tested(level)) + (nodes - filled);
        }

        // Expected tests at one node: the first always runs, test t runs while every distinct key
        // before it was positive. With parentHeld set, probabilities are conditioned on a passing parent.
        double tests(double[] held, double[] parentHeld, int tested) {
            double expected = 1.0;
            int distinct = held.length;
            for (int m = 0; m < distinct; m++) {
                int next = (m + 1 < distinct) ? first[m + 1] : tested - 1;
                int run = next - first[m];
                if (run <= 0) {
                    continue;
                }
                double g = held[m];
                if (parentHeld != null) {
                    double gp = m < parentHeld.length ? parentHeld[m] : 0.0;
                    g = gp > 0.0 ? MathUtils.clamp01(g / gp) : 0.0;
                }
                expected += run * g;
            }
            return expected;
        }

        // held[L - from][m] = G_m(L) = S_m + (1 - S_m) * B_m for levels from..to.
        double[][] held(int from, int to) {
            int levels = to - from + 1;
            double[][] out = new double[levels][];
            int[] distinctTested = new int[levels];
            double[] nodeWidth = new double[levels];
            // log of the share of positions outside one node, for placing a known count of occurrences
            double[] logOutside = new double[levels];
            double[] background = new double[levels];
            boolean[] negligible = new boolean[levels];
            int widest = 0;
            for (int i = 0; i < levels; i++) {
                int level = from + i;
                nodeWidth[i] = width >> level;
                logOutside[i] = Math.log1p(-Math.min(1.0, nodeWidth[i] / records));
                int ell = tested(level);
                int m = 0;
                while (m < first.length && first[m] < ell) {
                    m++;
                }
                distinctTested[i] = m;
                widest = Math.max(widest, m);
                out[i] = new double[m];
            }
            for (int m = 0; m < widest; m++) {
                double used = countUsed[m];
                if (rebuild[m]) {
                    // The prefix count moved, so every residual count so far changes.
                    Arrays.fill(background, 1.0);
                    Arrays.fill(negligible, false);
                    for (int j = 0; j <= m; j++) {
                        multiplyBackground(j, used, from, logOutside, distinctTested, background, negligible);
                    }
                } else {
                    multiplyBackground(m, used, from, logOutside, distinctTested, background, negligible);
                }
                double c = counts[m];
                double starts = Math.max(1.0, records - spans[m] + 1.0);
                for (int i = 0; i < levels; i++) {
                    if (m >= distinctTested[i]) {
                        continue;
                    }
                    if (!negligible[i] && m >= lastRebuild && background[i] < NEGLIGIBLE) {
                        negligible[i] = true;
                        background[i] = 0.0;
                    }
                    // c occurrences of the prefix placed among the possible starts; fit of them lie
                    // inside this node.
                    double fit = nodeWidth[i] - spans[m] + 1.0;
                    double share = Math.min(1.0, fit / starts);
                    double signal = (fit > 0.0 && c > 0.0) ? -Math.expm1(c * Math.log1p(-share)) : 0.0;
                    out[i][m] = signal + (1.0 - signal) * background[i];
                }
            }
            return out;
        }

        // One key's positive answer from occurrences outside the prefix's own, or a false positive.
        // The key's count is known, so its residual occurrences are placed uniformly and the node
        // holds one unless all of them fall outside it: 1 - (1 - b/N)^count. The iid form
        // 1 - (1 - count/N)^b treats the count as random and says 63% for a key that occurs once,
        // even for the root, which holds it for certain.
        private void multiplyBackground(int m,
                                        double prefixCount,
                                        int from,
                                        double[] logOutside,
                                        int[] distinctTested,
                                        double[] background,
                                        boolean[] negligible) {
            double residual = Math.max(0.0, MathUtils.clamp01(probs[first[m]]) * records - prefixCount);
            for (int i = 0; i < background.length; i++) {
                if (m >= distinctTested[i] || negligible[i]) {
                    continue;
                }
                double h = residual > 0.0 ? -Math.expm1(residual * logOutside[i]) : 0.0;
                double b = beta[from + i];
                background[i] *= b + (1.0 - b) * h;
            }
        }
    }

    // Costs for every start level of one pattern. Work below the start does not depend on where the
    // search started, so one suffix sum per level serves all candidates.
    private static final class Plan {
        private final Model model;
        private final double[] start;
        private final double[] below;
        final int deepest;
        final int maxLevel;

        Plan(Model model, double[] start, double[] below, int deepest, int maxLevel) {
            this.model = model;
            this.start = start;
            this.below = below;
            this.deepest = deepest;
            this.maxLevel = maxLevel;
        }

        double cost(int lp) {
            int level = Math.max(0, Math.min(lp, maxLevel));
            if (Double.isNaN(start[level])) {
                start[level] = model.startCost(level, model.held(level, level)[0]);
            }
            double deeper = level + 1 <= deepest ? below[level + 1] : 0.0;
            return start[level] + deeper;
        }
    }
}
