package search;

import estimators.Estimator;

// The order in which one query tests its tokens against a node's filter (hbi-optimizations.md, Section 6; used
// with BlockSearch.Options.queryPlan). Built once per query from the pattern: each distinct token once, with the
// smallest offset (in tree positions from the start of an occurrence) at which it occurs, optionally rarest first by
// a frequency sketch. Block Search walks the plan in order and stops at the first token that tests negative; with
// moveToFront that token goes first for every later node, at every level and in every tree of the query.
public final class QueryPlan {

    private final long[] tokens;
    private final int[] offsets;
    private final int size;
    private final boolean strides;

    private QueryPlan(long[] tokens, int[] offsets, int size, boolean strides) {
        this.tokens = tokens;
        this.offsets = offsets;
        this.size = size;
        this.strides = strides;
    }

    // The tokens Block Search tests: every n-gram, or with strides every n-th one and the last
    // (Pattern.buildStrideTokens), so stride token j sits at offset min(j * n, last n-gram's offset).
    // sketch: when not null, the plan starts rarest first by its counts (ties keep the smallest offset first).
    public static QueryPlan of(Pattern p, boolean strides, Estimator sketch) {
        long[] source = strides ? p.effectiveNgramArr : p.nGramToLong;
        int lastOffset = p.nGramToLong.length - 1;
        int step = strides ? Math.max(1, p.nGram) : 1;

        long[] tokens = new long[source.length];
        int[] offsets = new int[source.length];
        // Open addressing over the tokens seen so far: 0 is empty, otherwise 1 + the token's place in the plan.
        int[] seen = new int[Integer.highestOneBit(Math.max(1, source.length)) << 2];
        int mask = seen.length - 1;
        int size = 0;
        for (int j = 0; j < source.length; j++) {
            long token = source[j];
            int slot = (int) ((token * 0x9E3779B97F4A7C15L) >>> 40) & mask;
            while (seen[slot] != 0 && tokens[seen[slot] - 1] != token) {
                slot = (slot + 1) & mask;
            }
            if (seen[slot] == 0) {   // first occurrence, so its smallest offset: offsets only grow with j
                tokens[size] = token;
                offsets[size] = Math.min(j * step, lastOffset);
                seen[slot] = ++size;
            }
        }
        if (sketch != null) {
            sortByCount(tokens, offsets, size, sketch);
        }
        return new QueryPlan(tokens, offsets, size, strides);
    }

    // Insertion sort by the sketch's count, stable, so equal counts stay in offset order.
    private static void sortByCount(long[] tokens, int[] offsets, int size, Estimator sketch) {
        long[] counts = new long[size];
        for (int i = 0; i < size; i++) counts[i] = sketch.get(tokens[i]);
        for (int i = 1; i < size; i++) {
            long c = counts[i], t = tokens[i];
            int o = offsets[i];
            int j = i - 1;
            while (j >= 0 && counts[j] > c) {
                counts[j + 1] = counts[j];
                tokens[j + 1] = tokens[j];
                offsets[j + 1] = offsets[j];
                j--;
            }
            counts[j + 1] = c;
            tokens[j + 1] = t;
            offsets[j + 1] = o;
        }
    }

    public int size() { return size; }

    public long token(int i) { return tokens[i]; }

    public int offset(int i) { return offsets[i]; }

    boolean strides() { return strides; }

    // Entry i tested negative: test it first from now on (move to front).
    public void failed(int i) {
        if (i == 0) return;
        long token = tokens[i];
        int offset = offsets[i];
        System.arraycopy(tokens, 0, tokens, 1, i);
        System.arraycopy(offsets, 0, offsets, 1, i);
        tokens[0] = token;
        offsets[0] = offset;
    }
}
