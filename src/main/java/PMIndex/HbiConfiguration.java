package PMIndex;

import estimators.CostFunction;
import estimators.Estimator;
import membership.Membership;
import search.PruningPlan;
import search.SearchAlgorithm;
import search.Verifier;
import utilities.Utils;

import java.util.Objects;
import java.util.function.Supplier;

// Immutable configuration for constructing HBI instances.
public final class HbiConfiguration {

    private final SearchAlgorithm searchAlgorithm;
    private final int windowLength;
    private final double fpRate;
    private final int alphabetSize;
    private final int treeLength;
    private final Supplier<Estimator> estimatorSupplier;
    private final Supplier<Membership> membershipSupplier;
    private final Supplier<PruningPlan> pruningPlanSupplier;
    private final Verifier verifier;
    private final CostFunction costFunction;
    private final double confidence;

    private final Utils.MemPolicy memPolicy;
    private final int buckets;
    private final boolean experimentMode;
    private final boolean collectStats;
    private final boolean activateEstim;
    private final boolean keepSubNgramLevels;
    private final int minPatternLength;
    private final long seed;

    private final int nGram;

    private HbiConfiguration(Builder builder) {
        this.searchAlgorithm = Objects.requireNonNull(builder.searchAlgorithm, "searchAlgorithm");
        this.windowLength = builder.windowLength;
        this.fpRate = builder.fpRate;
        this.alphabetSize = builder.alphabetSize;
        this.treeLength = builder.treeLength;
        this.estimatorSupplier = Objects.requireNonNull(builder.estimatorSupplier, "estimatorSupplier");
        this.membershipSupplier = Objects.requireNonNull(builder.membershipSupplier, "membershipSupplier");
        this.pruningPlanSupplier = Objects.requireNonNull(builder.pruningPlanSupplier, "pruningPlanSupplier");
        this.verifier = Objects.requireNonNull(builder.verifier, "verifier");
        this.costFunction = builder.costFunction;
        this.confidence = builder.confidence;
        this.experimentMode = builder.experimentMode;
        this.collectStats = builder.collectStats;
        this.nGram = builder.nGram;
        this.memPolicy = builder.memPolicy;
        this.buckets = builder.buckets;
        this.activateEstim = builder.activateEstim;
        this.keepSubNgramLevels = builder.keepSubNgramLevels;
        this.minPatternLength = builder.minPatternLength;
        this.seed = builder.seed;
        validate();
    }

    public static Builder builder() { return new Builder(); }

    private void validate() {
        if (windowLength <= 0) {
            throw new IllegalArgumentException("windowLength must be positive");
        }
        if (treeLength <= 0) {
            throw new IllegalArgumentException("treeLength must be positive");
        }
        if (alphabetSize <= 0) {
            throw new IllegalArgumentException("alphabetSize must be positive");
        }
        if (nGram <= 0) {
            throw new IllegalArgumentException("nGram must be positive");
        }
        // Positions inside the retained trees are ints; at most windowLength + 2 * treeLength of them are kept.
        if ((long) windowLength + 2L * treeLength > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("windowLength + 2 * treeLength must fit in an int");
        }
        if (fpRate < 0.0 || fpRate >= 1.0) {
            throw new IllegalArgumentException("fpRate must be in [0,1)");
        }
        if (confidence <= 0.0 || confidence >= 1.0) {
            throw new IllegalArgumentException("confidence must be in (0,1)");
        }
        if (minPatternLength < 0) {
            throw new IllegalArgumentException("minPatternLength must be >= 0");
        }
    }
    public int buckets(){return buckets;}
    public SearchAlgorithm searchAlgorithm() { return searchAlgorithm; }
    public int windowLength() { return windowLength; }
    public double fpRate() { return fpRate; }
    public int alphabetSize() { return alphabetSize; }
    public int treeLength() { return treeLength; }
    public Supplier<Estimator> estimatorSupplier() { return estimatorSupplier; }
    public Supplier<Membership> membershipSupplier() { return membershipSupplier; }
    public Supplier<PruningPlan> pruningPlanSupplier() { return pruningPlanSupplier; }
    public Verifier verifier() { return verifier; }
    public CostFunction costFunction() { return costFunction; }
    public double confidence() { return confidence; }
    public boolean experimentMode() { return experimentMode; }
    public int nGram(){return nGram;}
    public boolean collectStats() { return collectStats; }
    public boolean activateEstim() { return activateEstim; }
    public boolean keepSubNgramLevels() { return keepSubNgramLevels; }
    public int minPatternLength() { return minPatternLength; }
    public long seed() { return seed; }

    public Utils.MemPolicy memPolicy() {
        return memPolicy;
    }

    public static final class Builder {
        private SearchAlgorithm searchAlgorithm;
        private int buckets;
        private int windowLength;
        private double fpRate;
        private int alphabetSize;
        private int treeLength;
        private Utils.MemPolicy memPolicy = Utils.MemPolicy.NONE;
        private Supplier<Estimator> estimatorSupplier;
        private Supplier<Membership> membershipSupplier;
        private Supplier<PruningPlan> pruningPlanSupplier;
        private Verifier verifier;
        private CostFunction costFunction;
        private double confidence;
        private boolean experimentMode = false;
        private boolean collectStats;
        private boolean activateEstim;
        private boolean keepSubNgramLevels;
        private int minPatternLength;
        private long seed = 1L;

        private int nGram = 1;

        private Builder() {
        }

        public Builder memPolicy(Utils.MemPolicy memPolicy) {
            this.memPolicy = (memPolicy == null) ? Utils.MemPolicy.NONE : memPolicy;
            return this;
        }
        public Builder searchAlgorithm(SearchAlgorithm searchAlgorithm) {
            this.searchAlgorithm = searchAlgorithm;
            return this;
        }

        public Builder windowLength(int windowLength) {
            this.windowLength = windowLength;
            return this;
        }

        public Builder fpRate(double fpRate) {
            this.fpRate = fpRate;
            return this;
        }

        public Builder alphabetSize(int alphabetSize) {
            this.alphabetSize = alphabetSize;
            return this;
        }

        public Builder treeLength(int treeLength) {
            this.treeLength = treeLength;
            return this;
        }

        public Builder estimatorSupplier(Supplier<Estimator> estimatorSupplier) {
            this.estimatorSupplier = estimatorSupplier;
            return this;
        }

        public Builder membershipSupplier(Supplier<Membership> membershipSupplier) {
            this.membershipSupplier = membershipSupplier;
            return this;
        }

        public Builder pruningPlanSupplier(Supplier<PruningPlan> pruningPlanSupplier) {
            this.pruningPlanSupplier = pruningPlanSupplier;
            return this;
        }

        public Builder verifier(Verifier verifier) {
            this.verifier = verifier;
            return this;
        }

        public Builder costFunction(CostFunction costFunction) {
            this.costFunction = costFunction;
            return this;
        }

        public Builder confidence(double confidence) {
            this.confidence = confidence;
            return this;
        }

        public Builder experimentMode(boolean experimentMode) {
            this.experimentMode = experimentMode;
            return this;
        }

        public Builder buckets(int buckets) {
            this.buckets = buckets;
            return this;
        }
        public Builder nGram(int nGram) {
            this.nGram = nGram;
            return this;
        }

        public Builder collectStats(boolean collectStats) {
            this.collectStats = collectStats;
            return this;
        }

        public Builder activateEstim(boolean activateEstim) {
            this.activateEstim = activateEstim;
            return this;
        }

        // Keep levels whose span is below the n-gram size (Sec. 4.3.1); off by default.
        public Builder keepSubNgramLevels(boolean keepSubNgramLevels) {
            this.keepSubNgramLevels = keepSubNgramLevels;
            return this;
        }

        // Shortest pattern length (characters, or tokens in segments mode) the queries are expected to have;
        // 0 (the default) keeps the full layout. Block Search enters a child only if the child spans at least
        // the pattern, so levels with a span below this are never probed by such patterns and are not built:
        // the leaf span becomes the smallest power of two >= minPatternLength. Shorter patterns stay exact,
        // their candidates are verified directly at the leaf level.
        public Builder minPatternLength(int minPatternLength) {
            this.minPatternLength = minPatternLength;
            return this;
        }

        // Seed of the Bloom filters' hash functions; each filter's own seed comes from it, its tree and its level,
        // so the same seed gives the same filters, false positives and probe counts in every run. Default 1.
        public Builder seed(long seed) {
            this.seed = seed;
            return this;
        }

        public HbiConfiguration build() {
            return new HbiConfiguration(this);
        }
    }
}
