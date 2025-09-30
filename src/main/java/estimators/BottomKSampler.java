package estimators;

import java.util.*;
import java.util.function.LongUnaryOperator;

/** BottomKSampler: keeps the k smallest priorities (lower is better). */
public final class BottomKSampler {

    private final int k;
    private final long prioritySeed;

    /** Max-heap by unsigned priority; root holds the largest priority among the kept k. */
    private final PriorityQueue<Entry> heap;

    /** Optional distinctness guard (to ignore duplicates); can remove if keys are known distinct. */
    private final HashSet<Integer> seen;

    public BottomKSampler(int k, long seed, boolean deduplicate) {
        if (k <= 0) throw new IllegalArgumentException("k must be > 0");
        this.k = k;
        this.prioritySeed = seed;
        this.heap = new PriorityQueue<>(k, Entry.BY_UNSIGNED_PRIORITY_DESC); // max-heap
        this.seen = deduplicate ? new HashSet<>(k * 4) : null;
    }

    /** Offer a key from the stream (distinct keys ideally). */
    public void offer(int key) {
        if (seen != null) {
            if (!seen.add(key)) return; // duplicate key → skip
        }
        long p = priorityOf(key);
        if (heap.size() < k) {
            heap.add(new Entry(key, p));
            return;
        }
        Entry worst = heap.peek(); // largest priority among the kept k
        if (unsignedLessThan(p, worst.priority)) {
            heap.poll();
            heap.add(new Entry(key, p));
        }
    }

    /** Return the current sample keys (size ≤ k). */
    public int[] sampleKeys() {
        int n = heap.size();
        int[] out = new int[n];
        int i = 0;
        for (Entry e : heap) out[i++] = e.key;
        return out;
    }

    /** Current threshold τ (unsigned). New keys must have priority < τ to enter; Long.MAX_VALUE if not full. */
    public long thresholdUnsigned() {
        if (heap.size() < k) return Long.MAX_VALUE;
        return heap.peek().priority;
    }

    public int size() { return heap.size(); }
    public int capacity() { return k; }

    // ---------- internals ----------

    private long priorityOf(int key) {
        long z = Integer.toUnsignedLong(key) ^ prioritySeed;
        return mix64(z);
    }

    private static boolean unsignedLessThan(long a, long b) {
        return (a ^ Long.MIN_VALUE) < (b ^ Long.MIN_VALUE);
    }

    /** SplitMix64 finalizer (Steele et al. 2014). */
    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    /** Heap entry */
    private static final class Entry {
        final int key;
        final long priority;
        Entry(int key, long priority) { this.key = key; this.priority = priority; }

        // Comparator: unsigned priority descending → max-heap
        static final Comparator<Entry> BY_UNSIGNED_PRIORITY_DESC =
                (a, b) -> Long.compareUnsigned(b.priority, a.priority);
    }
}
