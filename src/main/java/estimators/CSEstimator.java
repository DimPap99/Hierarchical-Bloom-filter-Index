package estimators;

import search.Pattern;


public class CSEstimator implements Estimator {

    private static final int DEFAULT_WIDTH = 2048;
    private static final int DEFAULT_DEPTH = 5;

    private final int baseWidth;
    private final int depth;
    private final long seed;

    private IntCountSketch sketch;
    private double totalRecords;

    public CSEstimator(int expectedRecords) {
        this(deriveWidth(expectedRecords), DEFAULT_DEPTH, System.nanoTime());
    }

    public CSEstimator(int width, int depth) {
        this(width, depth, System.nanoTime());
    }

    public CSEstimator(int width, int depth, long seed) {
        if (width <= 0) {
            throw new IllegalArgumentException("width must be positive");
        }
        if (depth <= 0) {
            throw new IllegalArgumentException("depth must be positive");
        }
        this.baseWidth = width;
        this.depth = depth;
        this.seed = seed;
        init(0);
    }

    @Override
    public void init(int totalRecords) {
        int width = Math.max(baseWidth, deriveWidth(totalRecords));
        this.sketch = new IntCountSketch(width, depth, seed);
        this.totalRecords = 0.0;
    }

    @Override
    public void insert(long key) {
        if (sketch == null) {
            init(0);
        }
        sketch.add(key);
        totalRecords += 1.0;
    }

    @Override
    public double estimate(long key) {
        if (sketch == null || totalRecords <= 0.0) {
            return 0.0;
        }
        long estimate = sketch.estimate(key);
        if (estimate <= 0L) {
            return 0.0;
        }
        return estimate / totalRecords;
    }

    @Override
    public long get(long key) {
        return this.sketch.estimate(key);
    }

    @Override
    public double[] estimateALl(Pattern p, boolean strides) {
        long[] tokens = strides ? p.effectiveNgramArr : p.nGramToLong;
        double[] result = new double[tokens.length];
        for (int i = 0; i < tokens.length; i++) {
            result[i] = estimate(tokens[i]);
        }
        return result;
    }

    // A Count Sketch does not track which keys it has seen, so it cannot report the rarest one.
    // The provisioning policies use the HOPS percentile sampler for that instead.
    @Override
    public double getMin() {
        throw new UnsupportedOperationException("CSEstimator cannot enumerate its keys");
    }

    private static int deriveWidth(int expectedRecords) {
        if (expectedRecords <= 0) {
            return DEFAULT_WIDTH;
        }
        int v = expectedRecords - 1;
        v |= v >>> 1;
        v |= v >>> 2;
        v |= v >>> 4;
        v |= v >>> 8;
        v |= v >>> 16;
        v = (v < 0) ? Integer.MAX_VALUE : v;
        long candidate = (long) v + 1L;
        int width = candidate > Integer.MAX_VALUE ? Integer.MAX_VALUE : (int) candidate;
        width = Math.max(2, width);
        return Math.max(DEFAULT_WIDTH, width);
    }
}
