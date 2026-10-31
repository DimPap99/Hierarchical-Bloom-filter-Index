package PMIndex;

import estimators.CSEstimator;
import estimators.CostFunction;
import estimators.CostFunctionDefaultRoot;
import estimators.CostFunctionOccurrence;
import estimators.Estimator;
import membership.BlockedBloomFilter;
import membership.BloomFilter;
import membership.Membership;
import search.BlockSearch;
import search.BlockSearchCharSet;
import search.MostFreqPruning;
import search.MultiLevelPruning;
import search.Pattern;
import search.PruningPlan;
import search.SearchAlgorithm;
import search.VerifierLinearLeafProbe;
import utilities.Utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Objects;
import java.util.function.Supplier;

/**
 * Exact pattern matching over the most recent {@code windowLength} characters of a stream.
 *
 * <pre>{@code
 * HierarchicalBloomFilterIndex index = HierarchicalBloomFilterIndex.builder()
 *         .windowLength(1 << 20)
 *         .nGram(4)
 *         .build();
 * index.append(someText);
 * long[] hits = index.search("needle");   // absolute offsets of every occurrence in the window
 * }</pre>
 *
 * <p>Each window position is summarized in a binary hierarchy of Bloom filters that prunes the search,
 * and every candidate is confirmed against the buffered window, so answers are exact: all occurrences
 * that start inside the current window and end at or before the last appended character, and nothing
 * else. Patterns must be at least {@code nGram} characters long.
 *
 * <p>Offsets are absolute: the first appended character has offset 0. Instances are not thread-safe.
 */
public final class HierarchicalBloomFilterIndex {

    /** How a query traverses the hierarchy (Sec. 4.2). All strategies return the same answers. */
    public enum SearchStrategy {
        /** Block Search: start every query at the root. */
        BLOCK_SEARCH,
        /** Confidence-Primed Block Search: start at the level where the rarest symbol stops being near-certain. */
        CONFIDENCE_PRIMED,
        /** Confidence-Gated Block Search: test each symbol only from its own gate level downward. */
        CONFIDENCE_GATED,
        /** Cost-Function Block Search: start at the level that minimizes the expected number of filter probes. */
        COST_FUNCTION
    }

    /** Per-level provisioning policy (Sec. 4.3.2) that drops top levels with little pruning value. */
    public enum ProvisioningPolicy {
        NONE, REACTIVE, PREDICTIVE
    }

    private static final long[] NO_MATCHES = new long[0];
    private static final int MAX_PACKED_NGRAM = 4; // four 16-bit chars fill one long exactly

    private final HBI hbi;
    private final int windowLength;
    private final int nGram;
    // Last nGram characters, used to emit the next n-gram symbol.
    private final char[] recent;
    // Last windowLength characters. Only kept when n-grams are hashed, to rule out hash collisions.
    private final char[] window;
    private long length;

    private HierarchicalBloomFilterIndex(Builder b) {
        this.windowLength = b.windowLength;
        this.nGram = b.nGram;
        this.recent = new char[nGram];
        this.window = (nGram > MAX_PACKED_NGRAM) ? new char[windowLength] : null;

        // The index works on n-gram symbols: n-gram i covers characters i..i+n-1, so a window of
        // windowLength characters holds windowLength - n + 1 n-grams.
        int symbolWindow = windowLength - nGram + 1;
        int treeLength = (b.treeLength > 0) ? b.treeLength : Integer.highestOneBit(symbolWindow);
        int distinct = (b.expectedDistinctSymbols > 0) ? b.expectedDistinctSymbols : treeLength;

        double confidence = b.confidence;
        double fpRate = b.fpRate;
        // Same wiring as IndexFactory: the start level comes from the cost function when there is
        // one, otherwise from the per-tree pruning plan.
        SearchAlgorithm search = (b.strategy == SearchStrategy.CONFIDENCE_GATED)
                ? new BlockSearchCharSet()
                : new BlockSearch();
        CostFunction costFunction = switch (b.strategy) {
            case BLOCK_SEARCH -> new CostFunctionDefaultRoot();
            case COST_FUNCTION -> new CostFunctionOccurrence();
            case CONFIDENCE_PRIMED, CONFIDENCE_GATED -> null;
        };
        Supplier<PruningPlan> pruningPlans = (b.strategy == SearchStrategy.CONFIDENCE_GATED)
                ? () -> new MultiLevelPruning(confidence, fpRate)
                : () -> new MostFreqPruning(confidence, fpRate);
        long sketchSeed = b.seed;
        Supplier<Estimator> estimators = () -> new CSEstimator(2048, 8, sketchSeed);
        Supplier<Membership> filters = b.blockedFilters ? BlockedBloomFilter::new : BloomFilter::new;

        Utils.MemPolicy policy = switch (b.policy) {
            case NONE -> Utils.MemPolicy.NONE;
            case REACTIVE -> Utils.MemPolicy.REACTIVE;
            case PREDICTIVE -> Utils.MemPolicy.PREDICTIVE;
        };
        int buckets = 0;
        if (policy != Utils.MemPolicy.NONE) {
            // Sec. 5.5 settings: rank error 0.025 with failure budgets of 0.05 each.
            buckets = Utils.designBucketsForRankTargetChebyshev(distinct, 0.025, 0.05, 0.05).suggestedBuckets;
        }

        HbiConfiguration config = HbiConfiguration.builder()
                .searchAlgorithm(search)
                .windowLength(symbolWindow)
                .treeLength(treeLength)
                .fpRate(fpRate)
                .alphabetSize(distinct)
                .estimatorSupplier(estimators)
                .membershipSupplier(filters)
                .pruningPlanSupplier(pruningPlans)
                .verifier(new VerifierLinearLeafProbe())
                .costFunction(costFunction)
                .confidence(confidence)
                .nGram(nGram)
                .memPolicy(policy)
                .activateEstim(policy != Utils.MemPolicy.NONE)
                .buckets(buckets)
                .minPatternLength(b.minPatternLength)
                .seed(b.seed)
                .build();
        this.hbi = new HBI(config);
        this.hbi.strides = b.strides;
    }

    public static Builder builder() {
        return new Builder();
    }

    /** Appends one character to the stream. */
    public void append(char c) {
        long pos = length++;
        recent[(int) (pos % nGram)] = c;
        if (window != null) {
            window[(int) (pos % windowLength)] = c;
        }
        if (length >= nGram) {
            hbi.insertToken(encodeRecent(pos - nGram + 1));
        }
    }

    /** Appends every character of {@code text} to the stream, in order. */
    public void append(CharSequence text) {
        for (int i = 0; i < text.length(); i++) {
            append(text.charAt(i));
        }
    }

    /**
     * Returns the absolute offsets, in increasing order, of every occurrence of {@code pattern} that
     * lies entirely inside the current window.
     *
     * @throws IllegalArgumentException if the pattern is shorter than the n-gram size
     */
    public long[] search(CharSequence pattern) {
        Pattern p = toPattern(pattern);
        if (p == null) {
            return NO_MATCHES;
        }

        ArrayList<Integer> hits = hbi.report(p);
        long windowStart = hbi.windowStartPosition();
        long[] out = new long[hits.size()];
        int count = 0;
        for (int hit : hits) {
            long offset = windowStart + hit;
            if (window == null || matchesWindow(pattern, offset)) {
                out[count++] = offset;
            }
        }
        return (count == out.length) ? out : Arrays.copyOf(out, count);
    }

    /**
     * True if {@code pattern} occurs at least once in the current window. Stops at the first occurrence.
     *
     * @throws IllegalArgumentException if the pattern is shorter than the n-gram size
     */
    public boolean contains(CharSequence pattern) {
        Pattern p = toPattern(pattern);
        if (p == null) {
            return false;
        }
        if (window == null) {
            return hbi.exists(p, null);
        }
        long windowStart = hbi.windowStartPosition();
        return hbi.exists(p, hit -> matchesWindow(pattern, windowStart + hit));
    }

    // The pattern with its tokens, or null when it cannot occur because it is longer than the window
    // or the stream.
    private Pattern toPattern(CharSequence pattern) {
        int r = pattern.length();
        if (r < nGram) {
            throw new IllegalArgumentException(
                    "pattern length " + r + " is shorter than the n-gram size " + nGram);
        }
        if (r > windowLength || r > length) {
            return null;
        }

        Pattern p = new Pattern(pattern.toString(), nGram);
        long[] tokens = new long[r - nGram + 1];
        for (int i = 0; i < tokens.length; i++) {
            tokens[i] = encode(pattern, i);
        }
        p.setTokens(tokens);
        return p;
    }

    /** Number of characters appended so far. */
    public long length() {
        return length;
    }

    /** Absolute offset of the first character in the current window. */
    public long windowStart() {
        return Math.max(0L, length - windowLength);
    }

    public int windowLength() {
        return windowLength;
    }

    public int nGram() {
        return nGram;
    }

    // n-gram starting at absolute offset start, read from the recent-character ring.
    private long encodeRecent(long start) {
        long h = seed();
        for (int i = 0; i < nGram; i++) {
            h = step(h, recent[(int) ((start + i) % nGram)]);
        }
        return h;
    }

    private long encode(CharSequence s, int from) {
        long h = seed();
        for (int i = 0; i < nGram; i++) {
            h = step(h, s.charAt(from + i));
        }
        return h;
    }

    // Up to four characters are packed losslessly; longer n-grams are hashed and every reported match
    // is re-checked against the raw characters in matchesWindow.
    private long seed() {
        return (nGram <= MAX_PACKED_NGRAM) ? 0L : 0x9E3779B97F4A7C15L ^ nGram;
    }

    private long step(long h, char c) {
        return (nGram <= MAX_PACKED_NGRAM) ? (h << 16) | c : mix64(h ^ c);
    }

    private boolean matchesWindow(CharSequence pattern, long offset) {
        for (int i = 0; i < pattern.length(); i++) {
            if (window[(int) ((offset + i) % windowLength)] != pattern.charAt(i)) {
                return false;
            }
        }
        return true;
    }

    // SplitMix64 finalizer.
    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /** Builder; only {@link #windowLength(int)} is required. Defaults follow the thesis (Ch. 5). */
    public static final class Builder {
        private int windowLength;
        private int nGram = 1;
        private int treeLength;
        private int expectedDistinctSymbols;
        private double fpRate = 0.2;
        private double confidence = 0.99;
        private boolean strides;
        private SearchStrategy strategy = SearchStrategy.BLOCK_SEARCH;
        private ProvisioningPolicy policy = ProvisioningPolicy.NONE;
        private long seed = 1L;
        private int minPatternLength;
        private boolean blockedFilters;

        private Builder() {
        }

        /** Number of most recent characters that queries search. Required. */
        public Builder windowLength(int windowLength) {
            this.windowLength = windowLength;
            return this;
        }

        /** Characters per indexed symbol (Sec. 4.3.1). Larger values prune better on low-entropy text. Default 1. */
        public Builder nGram(int nGram) {
            this.nGram = nGram;
            return this;
        }

        /**
         * Symbols per tree; must be a power of two. The window is covered by several trees and whole
         * trees expire (Sec. 4.1.3). Default: the largest power of two not above the window.
         */
        public Builder treeLength(int treeLength) {
            this.treeLength = treeLength;
            return this;
        }

        /**
         * Upper bound on distinct n-grams in one tree, used to size the Bloom filters. Underestimates
         * only cost speed; the default assumes every position is distinct.
         */
        public Builder expectedDistinctSymbols(int expectedDistinctSymbols) {
            this.expectedDistinctSymbols = expectedDistinctSymbols;
            return this;
        }

        /** Target false positive rate of each Bloom filter level. Default 0.2. */
        public Builder fpRate(double fpRate) {
            this.fpRate = fpRate;
            return this;
        }

        /** Confidence alpha for the confidence-based strategies and policies. Default 0.99. */
        public Builder confidence(double confidence) {
            this.confidence = confidence;
            return this;
        }

        /** Test only every n-th n-gram of a pattern (Sec. 4.3.1); useful for large n. Default false. */
        public Builder strides(boolean strides) {
            this.strides = strides;
            return this;
        }

        public Builder searchStrategy(SearchStrategy strategy) {
            this.strategy = Objects.requireNonNull(strategy, "strategy");
            return this;
        }

        public Builder provisioningPolicy(ProvisioningPolicy policy) {
            this.policy = Objects.requireNonNull(policy, "policy");
            return this;
        }

        /**
         * Length of the shortest pattern you expect to search. The search only enters a level whose span is at
         * least the pattern's length, so the levels finer than this are never used by such patterns; with this
         * set they are not built, which saves memory and insert time and leaves those queries unchanged.
         * Shorter patterns still get exact answers, only more slowly. Default: build every level.
         */
        public Builder minPatternLength(int minPatternLength) {
            this.minPatternLength = minPatternLength;
            return this;
        }

        /**
         * Use blocked Bloom filters: all bits of a key in one 64-bit word (one 512-bit block below a rate of about
         * 0.15), so a lookup is one memory access instead of several. Sized for the same false positive rate, with
         * up to about 1% more bits. Default false.
         */
        public Builder blockedFilters(boolean blockedFilters) {
            this.blockedFilters = blockedFilters;
            return this;
        }

        /**
         * Seed of the frequency sketch and of the Bloom filters' hash functions: the same seed gives the same
         * start-level choices and the same filter work in every run. Default 1.
         */
        public Builder seed(long seed) {
            this.seed = seed;
            return this;
        }

        public HierarchicalBloomFilterIndex build() {
            if (windowLength <= 0) {
                throw new IllegalArgumentException("windowLength must be positive");
            }
            if (nGram <= 0 || nGram > windowLength) {
                throw new IllegalArgumentException("nGram must be in [1, windowLength]");
            }
            if (treeLength < 0 || (treeLength > 0 && Integer.bitCount(treeLength) != 1)) {
                throw new IllegalArgumentException("treeLength must be a power of two");
            }
            if (expectedDistinctSymbols < 0) {
                throw new IllegalArgumentException("expectedDistinctSymbols must be positive");
            }
            if (!(fpRate > 0.0 && fpRate < 1.0)) {
                throw new IllegalArgumentException("fpRate must be in (0,1)");
            }
            if (!(confidence > 0.0 && confidence < 1.0)) {
                throw new IllegalArgumentException("confidence must be in (0,1)");
            }
            if (minPatternLength < 0) {
                throw new IllegalArgumentException("minPatternLength must be >= 0");
            }
            return new HierarchicalBloomFilterIndex(this);
        }
    }
}
