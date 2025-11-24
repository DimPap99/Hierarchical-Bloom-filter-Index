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

    private int nGram;

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
        if (fpRate < 0.0 || fpRate >= 1.0) {
            throw new IllegalArgumentException("fpRate must be in [0,1)");
        }
        if (confidence <= 0.0 || confidence >= 1.0) {
            throw new IllegalArgumentException("confidence must be in (0,1)");
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
        private boolean experimentMode = true;
        private boolean collectStats;
        private boolean activateEstim;

        private int nGram;

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

        public HbiConfiguration build() {
            return new HbiConfiguration(this);
        }
    }
}
