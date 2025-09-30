package PMIndex;

import java.util.*;

/** Variable-order n-gram model (order = maxOrder). No smoothing; pure backoff. */
public class  NgramModel {

    /** Immutable model used at query time. */
    public static final class Model {
        public final int sigma;
        public final int maxOrder;        // e.g., 3 => up to trigram contexts
        public final long[] uni;          // unigram counts
        public final long total;          // total symbols
        // For each context length t in [1..maxOrder-1], a map: ctxCode -> (nextSym -> count)
        public final List<Map<Long, Map<Integer, Long>>> ctxMaps;

        public Model(int sigma, int maxOrder,
                     long[] uni, long total,
                     List<Map<Long, Map<Integer, Long>>> ctxMaps) {
            this.sigma     = sigma;
            this.maxOrder  = maxOrder;
            this.uni       = uni;
            this.total     = total;
            this.ctxMaps   = ctxMaps;
        }

        /** π[v] = unigram MLE. */
        public double pi(int v) {
            long c = (v >= 0 && v < uni.length) ? uni[v] : 0L;
            return (total > 0) ? ((double)c / (double)total) : 0.0;
        }

        /** Encode a context (last t symbols) as base-sigma code. */
        private long encode(int[] seq, int startInclusive, int endExclusive) {
            long code = 0L;
            for (int i = startInclusive; i < endExclusive; i++) {
                code = code * sigma + (seq[i] & 0xffffffffL);
            }
            return code;
        }

        /**
         * Backoff conditional: try longest t in [maxOrder-1 .. 1]; if unseen, fall to unigram.
         * No smoothing; if even unigram 0, returns 0.
         */
        public double P_cond(int[] keySeq, int pos) {
            int v = keySeq[pos];
            int maxCtx = Math.min(maxOrder - 1, pos);   // how many previous symbols we have
            for (int t = maxCtx; t >= 1; t--) {
                long ctx = encode(keySeq, pos - t, pos); // last t symbols before v
                Map<Long, Map<Integer, Long>> mapT = ctxMaps.get(t - 1);
                Map<Integer, Long> row = mapT.get(ctx);
                if (row != null) {
                    long rowSum = 0L;
                    long cNext = 0L;
                    for (Map.Entry<Integer, Long> e : row.entrySet()) {
                        rowSum += e.getValue();
                        if (e.getKey() == v) cNext = e.getValue();
                    }
                    if (rowSum > 0L) {
                        return (double) cNext / (double) rowSum; // MLE for this context
                    }
                }
            }
            // backoff to unigram
            return pi(v);
        }
    }

    /** Streaming builder. Feed symbols in corpus order. */
    public static final class Builder {
        private final int sigma;
        private final int maxOrder;
        private final long[] uni;
        private long total = 0L;

        // For t=1..maxOrder-1: ctxMaps[t-1] : ctxCode -> (nextSym -> count)
        private final List<Map<Long, Map<Integer, Long>>> ctxMaps;

        // Sliding window of up to (maxOrder-1) previous symbols
        private final ArrayDeque<Integer> hist = new ArrayDeque<>();

        public Builder(int sigma, int maxOrder) {
            if (maxOrder < 1) throw new IllegalArgumentException("maxOrder must be >= 1");
            this.sigma    = sigma;
            this.maxOrder = maxOrder;
            this.uni      = new long[sigma];
            this.ctxMaps  = new ArrayList<>(Math.max(0, maxOrder - 1));
            for (int t = 1; t <= maxOrder - 1; t++) {
                ctxMaps.add(new HashMap<>());
            }
        }

        /** Observe one symbol id (same id space as your estimator). */
        public void observeSymbol(int s) {
            if (s < 0 || s >= sigma) { hist.clear(); return; }
            // update unigrams
            uni[s]++; total++;

            // encode contexts of length t = 1..maxOrder-1 (if we have at least t previous)
            if (!hist.isEmpty()) {
                // materialize history into array for simple encoding
                int[] buf = new int[hist.size()];
                int idx = 0;
                for (int x : hist) buf[idx++] = x;

                for (int t = 1; t <= Math.min(maxOrder - 1, buf.length); t++) {
                    long ctx = encodeSuffix(buf, t);  // last t of history
                    Map<Long, Map<Integer, Long>> mapT = ctxMaps.get(t - 1);
                    Map<Integer, Long> row = mapT.computeIfAbsent(ctx, k -> new HashMap<>());
                    row.put(s, row.getOrDefault(s, 0L) + 1L);
                }
            }

            // push s into history (cap length at maxOrder-1)
            hist.addLast(s);
            if (hist.size() > maxOrder - 1) hist.removeFirst();
        }

        /** Optional: break dependencies across boundaries. */
        public void resetChain() { hist.clear(); }

        private long encodeSuffix(int[] a, int t) {
            long code = 0L;
            int n = a.length;
            for (int i = n - t; i < n; i++) code = code * sigma + (a[i] & 0xffffffffL);
            return code;
        }

        public Model build() {
            return new Model(sigma, maxOrder, uni, total, ctxMaps);
        }
    }
}
