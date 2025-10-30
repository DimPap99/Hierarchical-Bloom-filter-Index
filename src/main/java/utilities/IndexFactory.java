package utilities;

import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.StreamingSlidingWindowIndex;
import PMIndex.DelayedStreamingSlidingWindowIndex;
import PMIndex.IPMIndexing;
import estimators.CostFunction;
import estimators.CostFunctionDefaultRoot;
import estimators.CostFunctionMarkov;
import estimators.Estimator;
import estimators.HashMapEstimator;
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

    public static HBI createHbi(int windowLength,
                                int treeLength,
                                int alphabetSize,
                                double fpRate,
                                double confidence,
                                Utils.MemPolicy memPolicy,
                                int ngram) {

        Supplier<Estimator> estFactory = () -> new HashMapEstimator(treeLength);
        Supplier<Membership> memFactory = BloomFilter::new;
        Supplier<search.PruningPlan> prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
        Verifier verifier = new VerifierLinearLeafProbe();

        HbiConfiguration configuration = HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(windowLength)
                .fpRate(fpRate)
                .alphabetSize(alphabetSize)
                .treeLength(treeLength)
                .estimatorSupplier(estFactory)
                .membershipSupplier(memFactory)
                .pruningPlanSupplier(prFactory)
                .verifier(verifier)
                .costFunction(null)
                .confidence(confidence)
                .experimentMode(false)
                .collectStats(false)
                .nGram(ngram)
                .memPolicy(memPolicy)
                .build();

        return new HBI(configuration);
    }

    public static StreamingSlidingWindowIndex createSuffixIndex(int windowLength) {
        return new StreamingSlidingWindowIndex(windowLength);
    }

    public static IPMIndexing createDelayedSuffixIndex(int windowLength,
                                                       int expectedAlphabetSize,
                                                       int delta,
                                                       boolean bufferedMode) {
        int safeDelta = Math.max(1, delta); // Delayed index requires positive delta
        return new StreamingSlidingWindowIndex(windowLength);//DelayedStreamingSlidingWindowIndex(windowLength, expectedAlphabetSize, safeDelta, bufferedMode);
    }
}
