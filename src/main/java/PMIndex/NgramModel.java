package PMIndex;

import java.util.*;

// Variable-order n-gram model with a prebuilt context-state operator.
// No smoothing; pure MLE from the stream.
public final class NgramModel {

    // Utilities.

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

    private static int powIntChecked(int base, int exp) {
        long acc = 1L;
        for (int i = 0; i < exp; i++) {
            acc *= base;
            if (acc > Integer.MAX_VALUE) {
                throw new IllegalArgumentException("base^exp exceeds int: " + base + "^" + exp);
            }
        }
        return (int) acc;
    }

    // Immutable snapshot.

    public static final class Model {
        public final int sigma;
        public final int maxOrder;                       // = ORDER + 1
        public final int ORDER;                          // context length t
        public final long[] uni;                         // unigram counts
        public final long total;
        public final List<Map<Long, Map<Integer, Long>>> ctxMaps;  // kept for compatibility

        private final Map<Long, Integer> symbolToIndex;
        private final long[] indexToSymbol;

        // First-order fallback (always filled).
        public final double[] PI;                        // length σ
        public final double[][] T;                       // [σ][σ]

        // Context-state operator for ORDER >= 1.
        public final int CTX_CARD;                       // = σ^ORDER, 0 if ORDER==0
        public final double[] P0_CTX;                    // length CTX_CARD, null if ORDER==0
        public final double[][] PNEXT;                   // [CTX_CARD][σ], null if ORDER==0
        public final int[][] NEXT_CTX;                   // [CTX_CARD][σ], null if ORDER==0

        public Model(int sigma,
                     int maxOrder,
                     int order,
                     long[] uni,
                     long total,
                     List<Map<Long, Map<Integer, Long>>> ctxMaps,
                     Map<Long, Integer> symbolToIndex,
                     long[] indexToSymbol,
                     double[] PI,
                     double[][] T,
                     int CTX_CARD,
                     double[] P0_CTX,
                     double[][] PNEXT,
                     int[][] NEXT_CTX) {
            this.sigma = sigma;
            this.maxOrder = maxOrder;
            this.ORDER = order;
            this.uni = uni;
            this.total = total;
            this.ctxMaps = ctxMaps;
            this.symbolToIndex = Collections.unmodifiableMap(new HashMap<>(symbolToIndex));
            this.indexToSymbol = indexToSymbol;
            this.PI = PI;
            this.T = T;
            this.CTX_CARD = CTX_CARD;
            this.P0_CTX = P0_CTX;
            this.PNEXT = PNEXT;
            this.NEXT_CTX = NEXT_CTX;
        }

        public Map<Long, Integer> symbolToIndex() { return symbolToIndex; }

        public int symbolIndex(long symbol) {
            Integer idx = symbolToIndex.get(symbol);
            return (idx == null) ? -1 : idx;
        }

        public long symbolForIndex(int idx) {
            if (idx < 0 || idx >= indexToSymbol.length) return -1L;
            return indexToSymbol[idx];
        }

        // Aggregated first-order transition probabilities derived from the observed contexts.
        public double[][] aggregatedFirstOrder() {
            return T;
        }

        public double pi(int v) {
            if (v < 0 || v >= PI.length) return 0.0;
            return PI[v];
        }

        // Hashed backoff conditional (kept for compatibility).
        public double P_cond(long[] keySeq, int pos) {
            int v = symbolIndex(keySeq[pos]);
            if (v < 0) return 0.0;
            int maxCtx = Math.min(maxOrder - 1, pos);
            for (int t = maxCtx; t >= 1; t--) {
                long ctx;
                if (t == 1) {
                    int prev = symbolIndex(keySeq[pos - 1]);
                    if (prev < 0) continue;
                    ctx = prev & 0xffffffffL;
                } else {
                    ctx = hashContextIds(keySeq, pos - t, pos, symbolToIndex);
                    if (ctx == Long.MIN_VALUE) continue;
                }
                Map<Long, Map<Integer, Long>> mapT = ctxMaps.get(t - 1);
                Map<Integer, Long> row = mapT.get(ctx);
                if (row != null) {
                    long rowSum = 0L, cNext = 0L;
                    for (Map.Entry<Integer, Long> e : row.entrySet()) {
                        rowSum += e.getValue();
                        if (e.getKey() == v) cNext = e.getValue();
                    }
                    if (rowSum > 0L) return (double) cNext / (double) rowSum;
                }
            }
            return pi(v);
        }
    }

    // ========= Streaming builder =========

    public static final class Builder {
        public final int sigma;
        public final int maxOrder;       // ORDER = maxOrder - 1
        public final int ORDER;
        public final long[] uni;
        public long total = 0L;

        // Keep the original sparse structures for compatibility
        private final List<Map<Long, Map<Integer, Long>>> ctxMaps;

        // Dense context-state structures for ORDER >= 1.
        private final int CTX_CARD;                // = sigma^ORDER, 0 if ORDER==0
        private final long[] ctxOcc;               // occurrences of each context (for P0), length CTX_CARD
        private final long[][] nextCounts;         // [CTX_CARD][sigma] counts for PNEXT rows
        private final int[][] nextCtx;             // [CTX_CARD][sigma] deterministic next context indices

        // Small history of previous symbol ids
        private final ArrayDeque<Integer> hist = new ArrayDeque<>();

        // Symbol mapping
        private final Map<Long, Integer> symbolToIndex;
        private final long[] indexToSymbol;
        private int distinctSymbols = 0;

        public Builder(int sigma, int maxOrder) {
            if (maxOrder < 1) throw new IllegalArgumentException("maxOrder must be >= 1");
            if (sigma <= 0) throw new IllegalArgumentException("sigma must be positive");
            this.sigma = sigma;
            this.maxOrder = maxOrder;
            this.ORDER = maxOrder - 1;
            this.uni = new long[sigma];

            // Sparse structures.
            this.ctxMaps = new ArrayList<>(Math.max(0, maxOrder - 1));
            for (int t = 1; t <= maxOrder - 1; t++) {
                ctxMaps.add(new HashMap<>());
            }

            // Dense context-state scaffolding.
            if (ORDER >= 1) {
                this.CTX_CARD = powIntChecked(sigma, ORDER);
                this.ctxOcc = new long[CTX_CARD];
                this.nextCounts = new long[CTX_CARD][sigma];
                this.nextCtx = new int[CTX_CARD][sigma];
                // precompute nextCtx mapping for all ctx, v
                final int powPrev = (ORDER >= 2) ? powIntChecked(sigma, ORDER - 1) : 1;
                for (int ctx = 0; ctx < CTX_CARD; ctx++) {
                    final int tail = (ORDER == 1) ? 0 : (ctx % powPrev);
                    for (int v = 0; v < sigma; v++) {
                        nextCtx[ctx][v] = tail * sigma + v;
                    }
                }
            } else {
                this.CTX_CARD = 0;
                this.ctxOcc = null;
                this.nextCounts = null;
                this.nextCtx = null;
            }

            this.symbolToIndex = new HashMap<>(Math.max(16, (int) Math.ceil(sigma / 0.75)));
            this.indexToSymbol = new long[sigma];
        }

        public void observeSymbol(int s) { observeSymbol((long) s); }

        public void observeSymbol(long symbol) {
            int id = ensureSymbolIndex(symbol);

            // unigram
            uni[id] += 1L;
            total += 1L;

            // if we have at least ORDER symbols in history, we can form a context
            if (ORDER >= 1) {
                if (hist.size() >= ORDER) {
                    int ctx = 0;
                    // oldest..newest over last ORDER ids in history
                    Iterator<Integer> it = hist.descendingIterator(); // newest..oldest
                    int[] tmp = new int[ORDER];
                    for (int i = ORDER - 1; i >= 0 && it.hasNext(); i--) {
                        tmp[i] = it.next();
                    }
                    for (int i = 0; i < ORDER; i++) ctx = ctx * sigma + tmp[i];

                    // record context occurrence
                    ctxOcc[ctx] += 1L;

                    // update context->next-symbol counts
                    nextCounts[ctx][id] += 1L;

                    // keep sparse map in sync for compatibility (optional but cheap)
                    long ctxHash = (ORDER == 1)
                            ? (tmp[ORDER - 1] & 0xffffffffL)
                            : hashContextIds(tmp, 0, ORDER);
                    Map<Long, Map<Integer, Long>> mapT = ctxMaps.get(ORDER - 1);
                    Map<Integer, Long> row = mapT.computeIfAbsent(ctxHash, k -> new HashMap<>());
                    row.put(id, row.getOrDefault(id, 0L) + 1L);
                }
            }

            // advance history
            hist.addLast(id);
            if (hist.size() > ORDER) hist.removeFirst();
        }

        // Optional: break dependencies between segments.
        public void resetChain() { hist.clear(); }

        private int ensureSymbolIndex(long symbol) {
            Integer existing = symbolToIndex.get(symbol);
            if (existing != null) return existing;
            int next = distinctSymbols;
            if (next >= sigma) {
                throw new IllegalStateException("Exceeded configured sigma capacity while building NgramModel: " + sigma);
            }
            symbolToIndex.put(symbol, next);
            indexToSymbol[next] = symbol;
            distinctSymbols++;
            return next;
        }

        // Build the immutable model; not for query time.
        public Model build() {
            // Trim symbol arrays to the actually used range, but keep sigma slots for convenience
            long[] uniTrim = uni;
            long totalSym = total;

            // π
            double[] PI = new double[sigma];
            if (totalSym > 0L) {
                double denom = (double) totalSym;
                for (int v = 0; v < sigma; v++) PI[v] = uniTrim[v] / denom;
            } else {
                double p = sigma > 0 ? 1.0 / sigma : 0.0;
                Arrays.fill(PI, p);
            }

            // First-order matrix aggregated from context counts
            double[][] T = new double[sigma][sigma];
            if (ORDER >= 1 && nextCounts != null) {
                long[][] aggregated = new long[sigma][sigma];

                if (ORDER == 1) {
                    for (int u = 0; u < sigma; u++) {
                        if (u < nextCounts.length) {
                            System.arraycopy(nextCounts[u], 0, aggregated[u], 0, sigma);
                        }
                    }
                } else {
                    final int ctxCount = Math.min(nextCounts.length, CTX_CARD);
                    for (int ctx = 0; ctx < ctxCount; ctx++) {
                        long[] row = nextCounts[ctx];
                        if (row == null) {
                            continue;
                        }
                        int prevSymbol = ctx % sigma;
                        long[] aggRow = aggregated[prevSymbol];
                        for (int v = 0; v < sigma; v++) {
                            long c = row[v];
                            if (c != 0L) {
                                aggRow[v] += c;
                            }
                        }
                    }
                }

                for (int u = 0; u < sigma; u++) {
                    long rowSum = 0L;
                    for (int v = 0; v < sigma; v++) {
                        rowSum += aggregated[u][v];
                    }
                    if (rowSum > 0L) {
                        double inv = 1.0 / (double) rowSum;
                        for (int v = 0; v < sigma; v++) {
                            T[u][v] = aggregated[u][v] * inv;
                        }
                    } else {
                        System.arraycopy(PI, 0, T[u], 0, sigma);
                    }
                }
            } else {
                for (int u = 0; u < sigma; u++) {
                    System.arraycopy(PI, 0, T[u], 0, sigma);
                }
            }

            // Context-state operator for ORDER >= 1
            int ctxCard = (ORDER >= 1) ? CTX_CARD : 0;
            double[] P0 = null;
            double[][] PNEXT = null;
            int[][] NEXT = null;

            if (ORDER >= 1) {
                P0 = new double[ctxCard];
                PNEXT = new double[ctxCard][sigma];
                NEXT = nextCtx; // already filled deterministically

                long totalCtx = 0L;
                for (int c = 0; c < ctxCard; c++) totalCtx += ctxOcc[c];
                if (totalCtx > 0L) {
                    double d = (double) totalCtx;
                    for (int c = 0; c < ctxCard; c++) P0[c] = ctxOcc[c] / d;
                } else {
                    double p = ctxCard > 0 ? 1.0 / ctxCard : 0.0;
                    Arrays.fill(P0, p);
                }

                for (int c = 0; c < ctxCard; c++) {
                    long rowSum = 0L;
                    for (int v = 0; v < sigma; v++) rowSum += nextCounts[c][v];
                    if (rowSum > 0L) {
                        double d = (double) rowSum;
                        for (int v = 0; v < sigma; v++) PNEXT[c][v] = nextCounts[c][v] / d;
                    } else {
                        // unseen context: back off to π
                        System.arraycopy(PI, 0, PNEXT[c], 0, sigma);
                    }
                }
            }

            Map<Long, Integer> mappingSnapshot = new HashMap<>(symbolToIndex);
            long[] idToSymbol = Arrays.copyOf(indexToSymbol, indexToSymbol.length);
            List<Map<Long, Map<Integer, Long>>> ctxCopy = ctxMaps;

            return new Model(
                    sigma,
                    maxOrder,
                    ORDER,
                    uniTrim,
                    totalSym,
                    ctxCopy,
                    mappingSnapshot,
                    idToSymbol,
                    PI,
                    T,
                    ctxCard,
                    P0,
                    PNEXT,
                    NEXT
            );
        }
    }
}
