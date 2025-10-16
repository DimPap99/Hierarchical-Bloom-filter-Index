package PMIndex;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Variable-order n-gram model (order = maxOrder). No smoothing; pure backoff. */
public class  NgramModel {

    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    private static long hashContextIds(int[] seq, int startInclusive, int endExclusive) {
        long h = 0x9E3779B97F4A7C15L;
        for (int i = startInclusive; i < endExclusive; i++) {
            long v = seq[i] & 0xffffffffL;
            h = mix64(h ^ (v + 0x9E3779B97F4A7C15L));
        }
        return h;
    }

    private static long hashContextIds(long[] seq, int startInclusive, int endExclusive, Map<Long, Integer> symbolToIndex) {
        long h = 0x9E3779B97F4A7C15L;
        for (int i = startInclusive; i < endExclusive; i++) {
            Integer idx = symbolToIndex.get(seq[i]);
            if (idx == null) {
                return Long.MIN_VALUE;
            }
            long v = idx & 0xffffffffL;
            h = mix64(h ^ (v + 0x9E3779B97F4A7C15L));
        }
        return h;
    }

    /** Immutable model used at query time. */
    public static final class Model {
        public final int sigma;
        public final int maxOrder;        // e.g., 3 => up to trigram contexts
        public final long[] uni;          // unigram counts
        public final long total;          // total symbols
        // For each context length t in [1..maxOrder-1], a map: ctxCode -> (nextSym -> count)
        public final List<Map<Long, Map<Integer, Long>>> ctxMaps;
        private final Map<Long, Integer> symbolToIndex;
        private final long[] indexToSymbol;

        public Model(int sigma, int maxOrder,
                     long[] uni, long total,
                     List<Map<Long, Map<Integer, Long>>> ctxMaps,
                     Map<Long, Integer> symbolToIndex,
                     long[] indexToSymbol) {
            this.sigma     = sigma;
            this.maxOrder  = maxOrder;
            this.uni       = uni;
            this.total     = total;
            this.ctxMaps   = ctxMaps;
            this.symbolToIndex = Collections.unmodifiableMap(new HashMap<>(symbolToIndex));
            this.indexToSymbol = indexToSymbol;
        }

        /** π[v] = unigram MLE. */
        public double pi(int v) {
            long c = (v >= 0 && v < uni.length) ? uni[v] : 0L;
            return (total > 0) ? ((double)c / (double)total) : 0.0;
        }

        public Map<Long, Integer> symbolToIndex() {
            return symbolToIndex;
        }

        public int symbolIndex(long symbol) {
            Integer idx = symbolToIndex.get(symbol);
            return (idx == null) ? -1 : idx;
        }

        public long symbolForIndex(int idx) {
            if (idx < 0 || idx >= indexToSymbol.length) {
                return -1L;
            }
            return indexToSymbol[idx];
        }

        /** Encode a context (last t symbols) as base-sigma code. */
        private long encode(long[] seq, int startInclusive, int endExclusive) {
            return hashContextIds(seq, startInclusive, endExclusive, symbolToIndex);
        }

        /**
         * Backoff conditional: try longest t in [maxOrder-1 .. 1]; if unseen, fall to unigram.
         * No smoothing; if even unigram 0, returns 0.
         */
        public double P_cond(long[] keySeq, int pos) {
            int v = symbolIndex(keySeq[pos]);
            if (v < 0) {
                return 0.0;
            }
            int maxCtx = Math.min(maxOrder - 1, pos);   // how many previous symbols we have
            for (int t = maxCtx; t >= 1; t--) {
                long ctx;
                if (t == 1) {
                    int prev = symbolIndex(keySeq[pos - 1]);
                    if (prev < 0) continue;
                    ctx = prev & 0xffffffffL;
                } else {
                    ctx = encode(keySeq, pos - t, pos); // hashed context
                    if (ctx == Long.MIN_VALUE) continue;
                }
                if (ctx == Long.MIN_VALUE) {
                    continue;
                }
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
        public final int sigma;
        public final int maxOrder;
        public final long[] uni;
        public long total = 0L;

        // For t=1..maxOrder-1: ctxMaps[t-1] : ctxCode -> (nextSym -> count)
        private final List<Map<Long, Map<Integer, Long>>> ctxMaps;

        // Sliding window of up to (maxOrder-1) previous symbols
        private final ArrayDeque<Integer> hist = new ArrayDeque<>();
        private final Map<Long, Integer> symbolToIndex;
        private final long[] indexToSymbol;
        private int distinctSymbols = 0;

        public Builder(int sigma, int maxOrder) {
            if (maxOrder < 1) throw new IllegalArgumentException("maxOrder must be >= 1");
            if (sigma <= 0) throw new IllegalArgumentException("sigma must be positive");
            this.sigma    = sigma;
            this.maxOrder = maxOrder;
            this.uni      = new long[sigma];
            this.ctxMaps  = new ArrayList<>(Math.max(0, maxOrder - 1));
            for (int t = 1; t <= maxOrder - 1; t++) {
                ctxMaps.add(new HashMap<>());
            }
            this.symbolToIndex = new HashMap<>(Math.max(16, (int)Math.ceil(sigma / 0.75))); // avoid frequent resizes
            this.indexToSymbol = new long[sigma];
        }

        /** Observe one symbol id (same id space as your estimator). */
        public void observeSymbol(int s) {
            observeSymbol((long) s);
        }

        public void observeSymbol(long symbol) {
            int id = ensureSymbolIndex(symbol);

            uni[id]++;
            total++;

            if (!hist.isEmpty()) {
                int histSize = hist.size();
                int[] buf = new int[histSize];
                int idx = 0;
                for (int x : hist) buf[idx++] = x;

                for (int t = 1; t <= Math.min(maxOrder - 1, buf.length); t++) {
                    long ctx;
                    if (t == 1) {
                        ctx = buf[buf.length - 1] & 0xffffffffL;
                    } else {
                        ctx = hashContextIds(buf, buf.length - t, buf.length);
                    }
                    Map<Long, Map<Integer, Long>> mapT = ctxMaps.get(t - 1);
                    Map<Integer, Long> row = mapT.computeIfAbsent(ctx, k -> new HashMap<>());
                    row.put(id, row.getOrDefault(id, 0L) + 1L);
                }
            }

            hist.addLast(id);
            if (hist.size() > maxOrder - 1) hist.removeFirst();
        }

        /** Optional: break dependencies across boundaries. */
        public void resetChain() { hist.clear(); }

        private int ensureSymbolIndex(long symbol) {
            Integer existing = symbolToIndex.get(symbol);
            if (existing != null) {
                return existing;
            }
            int next = distinctSymbols;
            if (next >= sigma) {
                throw new IllegalStateException("Exceeded configured sigma capacity while building NgramModel: " + sigma);
            }
            symbolToIndex.put(symbol, next);
            indexToSymbol[next] = symbol;
            distinctSymbols++;
            return next;
        }

        /**
         * Ensure the provided symbol has an assigned internal index without affecting
         * unigram or context counts. Useful when queries introduce previously unseen
         * hashes that still need to be part of the Markov model.
         */
        public int ensureSymbolRegistered(long symbol) {
            return ensureSymbolIndex(symbol);
        }

        public Model build() {
            long[] uniTrim = uni;
            if (distinctSymbols < uni.length) {
                uniTrim = Arrays.copyOf(uni, distinctSymbols);
            }
            long[] idToSymbol = Arrays.copyOf(indexToSymbol, distinctSymbols);
            Map<Long, Integer> mappingSnapshot = new HashMap<>(symbolToIndex);
            return new Model(distinctSymbols, maxOrder, uniTrim, total, ctxMaps, mappingSnapshot, idToSymbol);
        }
    }
}
