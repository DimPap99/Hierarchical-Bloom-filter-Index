package datagenerators;

import PMIndex.NgramModel;

import java.util.Arrays;
import java.util.Objects;
import java.util.Random;

/**
 * Utility builder that synthesises sequences from a variable-order Markov chain over a configurable
 * character domain. The generator keeps the {@link NgramModel.Model} snapshot in sync with the
 * emitted stream so the resulting cost function can reuse the exact chain statistics.
 */
public final class MarkovChainGenerator {

    private MarkovChainGenerator() {
    }

    public static MarkovDataset generate(int length,
                                         int minDomain,
                                         int maxDomain,
                                         int order) {
        return generate(length, minDomain, maxDomain, order, new Random());
    }

    public static MarkovDataset generate(int length,
                                         int minDomain,
                                         int maxDomain,
                                         int order,
                                         long seed) {
        return generate(length, minDomain, maxDomain, order, new Random(seed));
    }

    public static MarkovDataset generate(int length,
                                         int minDomain,
                                         int maxDomain,
                                         int order,
                                         Random random) {
        Objects.requireNonNull(random, "random");
        validateArguments(length, minDomain, maxDomain, order);

        int sigma = maxDomain - minDomain;
        int builderOrder = Math.max(1, order + 1);
        NgramModel.Builder builder = new NgramModel.Builder(sigma, builderOrder);
        StringBuilder sb = new StringBuilder(Math.max(0, length));

        if (order == 0) {
            double[] unigram = randomDistribution(sigma, random);
            for (int i = 0; i < length; i++) {
                int symbolIdx = sample(unigram, random);
                char ch = (char) (minDomain + symbolIdx);
                sb.append(ch);
                builder.observeSymbol((long) ch);
            }
            return new MarkovDataset(sb.toString(), builder.build(), order, minDomain, maxDomain, unigram, null, null);
        }

        int ctxCard = powChecked(sigma, order);
        double[] initialCtx = randomDistribution(ctxCard, random);
        double[][] transitions = new double[ctxCard][sigma];
        for (int ctx = 0; ctx < ctxCard; ctx++) {
            transitions[ctx] = randomDistribution(sigma, random);
        }

        int ctxStride = (order == 1) ? 1 : powChecked(sigma, order - 1);
        int ctx = sample(initialCtx, random);
        int[] ctxSymbols = new int[order];
        decodeContext(ctx, order, sigma, ctxSymbols);

        int produced = 0;
        for (int i = 0; i < order && produced < length; i++) {
            int symIdx = ctxSymbols[i];
            char ch = (char) (minDomain + symIdx);
            sb.append(ch);
            builder.observeSymbol((long) ch);
            produced++;
        }

        while (produced < length) {
            int nextSymbol = sample(transitions[ctx], random);
            char ch = (char) (minDomain + nextSymbol);
            sb.append(ch);
            builder.observeSymbol((long) ch);
            produced++;

            if (order == 1) {
                ctx = nextSymbol;
            } else {
                int trimmed = ctx % ctxStride;
                ctx = trimmed * sigma + nextSymbol;
            }
        }

        return new MarkovDataset(sb.toString(), builder.build(), order, minDomain, maxDomain, null, initialCtx, transitions);
    }

    private static void validateArguments(int length,
                                          int minDomain,
                                          int maxDomain,
                                          int order) {
        if (length < 0) {
            throw new IllegalArgumentException("length must be non-negative");
        }
        if (minDomain < Character.MIN_VALUE) {
            throw new IllegalArgumentException("minDomain < " + (int) Character.MIN_VALUE);
        }
        if (maxDomain > Character.MAX_VALUE + 1) {
            throw new IllegalArgumentException("maxDomain > " + ((int) Character.MAX_VALUE + 1));
        }
        if (maxDomain <= minDomain) {
            throw new IllegalArgumentException("maxDomain must be greater than minDomain");
        }
        if (order < 0) {
            throw new IllegalArgumentException("order must be >= 0");
        }
        int sigma = maxDomain - minDomain;
        if (order > 0) {
            powChecked(sigma, order); // ensure we fail early on overflow
        }
    }

    private static double[] randomDistribution(int size, Random random) {
        double[] weights = new double[size];
        double sum = 0.0;
        for (int i = 0; i < size; i++) {
            double w = random.nextDouble();
            if (w == 0.0d) {
                w = Double.MIN_NORMAL;
            }
            weights[i] = w;
            sum += w;
        }
        if (sum == 0.0) {
            Arrays.fill(weights, 1.0 / size);
            return weights;
        }
        double inv = 1.0 / sum;
        for (int i = 0; i < size; i++) {
            weights[i] *= inv;
        }
        return weights;
    }

    private static int sample(double[] distribution, Random random) {
        double r = random.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < distribution.length; i++) {
            cumulative += distribution[i];
            if (r <= cumulative || i == distribution.length - 1) {
                return i;
            }
        }
        return distribution.length - 1;
    }

    private static void decodeContext(int ctx,
                                      int order,
                                      int sigma,
                                      int[] dest) {
        int value = ctx;
        for (int i = order - 1; i >= 0; i--) {
            dest[i] = value % sigma;
            value /= sigma;
        }
    }

    private static int powChecked(int base, int exp) {
        long acc = 1L;
        for (int i = 0; i < exp; i++) {
            acc *= base;
            if (acc > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("sigma^order exceeds integer range");
            }
        }
        return (int) acc;
    }

    public static final class MarkovDataset {
        private final String data;
        private final NgramModel.Model model;
        private final int order;
        private final int minDomain;
        private final int maxDomain;
        private final double[] unigram;
        private final double[] contextInitial;
        private final double[][] transitions;

        private MarkovDataset(String data,
                              NgramModel.Model model,
                              int order,
                              int minDomain,
                              int maxDomain,
                              double[] unigram,
                              double[] contextInitial,
                              double[][] transitions) {
            this.data = data;
            this.model = model;
            this.order = order;
            this.minDomain = minDomain;
            this.maxDomain = maxDomain;
            this.unigram = copy1D(unigram);
            this.contextInitial = copy1D(contextInitial);
            this.transitions = copy2D(transitions);
        }

        public String data() {
            return data;
        }

        public NgramModel.Model model() {
            return model;
        }

        public int order() {
            return order;
        }

        public int minDomain() {
            return minDomain;
        }

        public int maxDomain() {
            return maxDomain;
        }

        public double[] unigramDistribution() {
            return copy1D(unigram);
        }

        public double[] initialContextDistribution() {
            return copy1D(contextInitial);
        }

        public double[][] transitionMatrix() {
            return copy2D(transitions);
        }

        private static double[] copy1D(double[] src) {
            return (src == null) ? null : Arrays.copyOf(src, src.length);
        }

        private static double[][] copy2D(double[][] src) {
            if (src == null) {
                return null;
            }
            double[][] copy = new double[src.length][];
            for (int i = 0; i < src.length; i++) {
                copy[i] = (src[i] == null) ? null : Arrays.copyOf(src[i], src[i].length);
            }
            return copy;
        }
    }
}
