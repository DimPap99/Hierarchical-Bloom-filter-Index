package utilities;

import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.StreamingSlidingWindowIndex;
import PMIndex.DelayedStreamingSlidingWindowIndex;
import PMIndex.IPMIndexing;
import PMIndex.SuffixTreeIndex;
import estimators.*;
import membership.BloomFilter;
import membership.Membership;
import search.BlockSearch;
import search.Verifier;
import search.VerifierLinearLeafProbe;

import java.util.function.Supplier;

/**
 * Central place to construct indexes for benchmarks (HBI and Suffix baseline).
 */
public final class IndexFactory {

    private IndexFactory() {}

    // Configurable n-gram for the StreamingSuffix baseline (default 1).
    // Not exposed via CLI; set programmatically if you want a different value.
    private static volatile int SUFFIX_NGRAM = 1;

    public static void setSuffixNgram(int ngram) {
        SUFFIX_NGRAM = Math.max(1, ngram);
    }

    public static int getSuffixNgram() {
        return SUFFIX_NGRAM;
    }

    public static HBI createHbi(int windowLength,
                                int treeLength,
                                int alphabetSize,
                                double fpRate,
                                double confidence,
                                Utils.MemPolicy memPolicy,
                                int ngram) {
        return createHbi(windowLength, treeLength, alphabetSize, fpRate, confidence, memPolicy, ngram, "bs", 0.05, 0);
    }

    public static HBI createHbi(int windowLength,
                                int treeLength,
                                int alphabetSize,
                                double fpRate,
                                double confidence,
                                Utils.MemPolicy memPolicy,
                                int ngram,
                                String algorithm) {
        return createHbi(windowLength, treeLength, alphabetSize, fpRate, confidence, memPolicy, ngram, algorithm, 0.05, 0);
    }

    public static HBI createHbi(int windowLength,
                                int treeLength,
                                int alphabetSize,
                                double fpRate,
                                double confidence,
                                Utils.MemPolicy memPolicy,
                                int ngram,
                                String algorithm,
                                double quantile,
                                int policyBuckets) {
        //Countsketch: ε=0.05, δ=7.5e-4 → w=2048, d=8 → ~64 K
        Supplier<Estimator> estFactory = () -> new CSEstimator(2048, 8, 1);//new HashMapEstimator(treeLength);
        Supplier<Membership> memFactory = BloomFilter::new;
        Verifier verifier = new VerifierLinearLeafProbe();

        // Choose components by algorithm token
        search.SearchAlgorithm search;
        Supplier<search.PruningPlan> prFactory;
        CostFunction costFn;
        String alg = (algorithm == null) ? "bs" : algorithm.toLowerCase();

        switch (alg) {
            case "cp" -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = null;
            }
            case "cgbs" -> {
                search = new search.BlockSearchCharSet();
                prFactory = () -> new search.MultiLevelPruning(confidence, fpRate);
                costFn = null;
            }
            case "cf" -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = new CostFunctionMaxProb();
            }
            case "bs" -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = new CostFunctionDefaultRoot();
            }
            default -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = new CostFunctionDefaultRoot();
            }
        }

        boolean activateEstimators = memPolicy != Utils.MemPolicy.NONE && costFn instanceof CostFunctionDefaultRoot;

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(search)
                .windowLength(windowLength)
                .fpRate(fpRate)
                .alphabetSize(alphabetSize)
                .treeLength(treeLength)
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(verifier)
                .costFunction(costFn)
                .confidence(confidence)
                .experimentMode(false)
                .collectStats(false)
                .nGram(ngram)
                .memPolicy(memPolicy)
                .activateEstim(activateEstimators)
                .buckets(Math.max(0, policyBuckets))
                .build();

        HBI hbi = new HBI(configuration);
        hbi.QUANTILE = Math.max(0.0, Math.min(1.0, quantile));
        return hbi;
    }

    public static StreamingSlidingWindowIndex createSuffixIndex(int windowLength) {
        return new StreamingSlidingWindowIndex(windowLength);
    }

    public static StreamingSlidingWindowIndex createSuffixIndex(int windowLength, int expectedAlphabetSize) {
        return new StreamingSlidingWindowIndex(windowLength, expectedAlphabetSize);
    }

    public static IPMIndexing createSuffixTreeIndex(int expectedAlphabetSize) {
        long distinctEstimate = Math.max(1L, expectedAlphabetSize);
        double epsilon = 0.0001;
        return new SuffixTreeIndex(distinctEstimate, epsilon);
    }

    public static IPMIndexing createDelayedSuffixIndex(int windowLength,
                                                       int expectedAlphabetSize,
                                                       int delta,
                                                       boolean bufferedMode) {
        int safeDelta = Math.max(1, delta); // Delayed index requires positive delta
        return new StreamingSlidingWindowIndex(windowLength);//DelayedStreamingSlidingWindowIndex(windowLength, expectedAlphabetSize, safeDelta, bufferedMode);
    }
}
