package utilities;

import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import PMIndex.StreamingSlidingWindowIndex;
import PMIndex.DelayedStreamingSlidingWindowIndex;
import PMIndex.IPMIndexing;
import PMIndex.SlidingSuffixIndex;
import PMIndex.SuffixTreeIndex;
import estimators.*;
import membership.Membership;
import search.BlockSearch;
import search.Verifier;
import search.VerifierLinearLeafProbe;

import java.util.function.Supplier;

// Central factory for constructing indexes for benchmarks (HBI and suffix baselines).
public final class IndexFactory {

    private IndexFactory() {}

    // Configurable n-gram for the StreamingSuffix baseline (default 1).
    // Not exposed via CLI; set programmatically if you want a different value.
    private static volatile int SUFFIX_NGRAM = 1;
    private static volatile int SUFFIX_TREE_NGRAM = 1;

    public static void setSuffixNgram(int ngram) {
        SUFFIX_NGRAM = Math.max(1, ngram);
    }

    public static int getSuffixNgram() {
        return SUFFIX_NGRAM;
    }

    public static void setSuffixTreeNgram(int ngram) {
        SUFFIX_TREE_NGRAM = Math.max(1, ngram);
    }

    public static int getSuffixTreeNgram() {
        return SUFFIX_TREE_NGRAM;
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
        return createHbi(windowLength, treeLength, alphabetSize, fpRate, confidence, memPolicy, ngram, algorithm,
                quantile, policyBuckets, HbiTuning.DEFAULT);
    }

    // HBI settings beyond the thesis's parameters, with their defaults: the leaf span (--min-pattern-length; 0 builds
    // every level), the Bloom filter of every level (--filter), the seed of the filters' hash functions (--seed),
    // Block Search's optional probe strategies (--query-plan, --plan-sketch-order, --plan-move-to-front,
    // --neighbour-check; all off) and the KMP verifier (--kmp-verification; off). hbi-optimizations.md has what each
    // option measured.
    public record HbiTuning(int minPatternLength, membership.FilterKind filter, long seed,
                            search.BlockSearch.Options searchOptions, boolean kmpVerification) {
        public static final HbiTuning DEFAULT =
                new HbiTuning(0, membership.FilterKind.BLOOM, 1L, search.BlockSearch.Options.DEFAULT, false);
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
                                int policyBuckets,
                                HbiTuning tuning) {
        // CountSketch: epsilon=0.05, delta=7.5e-4 -> w=2048, d=8 -> ~64K.
        Supplier<Estimator> estFactory = () -> new CSEstimator(2048, 8, 1);//new HashMapEstimator(treeLength);
        Supplier<Membership> memFactory = tuning.filter().supplier();
        Verifier verifier = tuning.kmpVerification() ? new search.VerifierKmp() : new VerifierLinearLeafProbe();

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
                //alphabetSize = 10000
            }
            case "cgbs" -> {
                search = new search.BlockSearchCharSet();
                prFactory = () -> new search.MultiLevelPruning(confidence, fpRate);
                costFn = null;
                //alphabetSize = 10000

            }
            case "cf" -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = new CostFunctionOccurrence();
//                alphabetSize = 1048576;

            }
            case "bs" -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = new CostFunctionDefaultRoot();
                //alphabetSize = 3000

            }
            default -> {
                search = new search.BlockSearch();
                prFactory = () -> new search.MostFreqPruning(confidence, fpRate);
                costFn = new CostFunctionDefaultRoot();
            }
        }

        if (search instanceof search.BlockSearch bs) {
            bs.setOptions(tuning.searchOptions());
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
                .minPatternLength(tuning.minPatternLength())
                .seed(tuning.seed())
                .build();

        HBI hbi = new HBI(configuration);
        hbi.QUANTILE = Math.max(0.0, Math.min(1.0, quantile));
        return hbi;
    }

    // SSWSI as set by --sswsi-mode, --sswsi-builder and --sswsi-dict. Neither index needs the alphabet size.
    public static IPMIndexing createSuffixIndex(int windowLength, SswsiSettings settings) {
        return switch (settings.mode()) {
            case LEGACY -> new StreamingSlidingWindowIndex(windowLength);
            case AMORTIZED -> new SlidingSuffixIndex(windowLength, settings.builder().create(), settings.dictionary());
        };
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
