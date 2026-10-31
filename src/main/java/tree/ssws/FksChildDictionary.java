package tree.ssws;

import java.util.Arrays;
import java.util.SplittableRandom;

import static tree.ssws.StaticSuffixTree.NONE;

/**
 * {@link ChildDictionary.Kind#FKS}: a static two-level perfect hash (Fredman, Komlós and Szemerédi, J. ACM 1984)
 * over all (node, first symbol) pairs of one tree.
 *
 * <p>The first level hashes the k keys into k buckets, and is redrawn until the buckets' squared sizes sum to at
 * most 4k (expected 2k). A bucket with b keys gets b² slots and its own second-level function, redrawn until it is
 * injective on the bucket, which succeeds with probability at least 1/2 per draw. A lookup is two hash evaluations
 * and one key comparison. Hashing is multiply-shift on the pair (first level) and a per-bucket seeded mixer of
 * that value (second level), from a fixed seed, so builds are reproducible.
 */
final class FksChildDictionary implements ChildDictionary {

    private static final long SEED = 0x2545F4914F6CDD1DL;
    private static final long GOLDEN = 0x9E3779B97F4A7C15L;
    private static final int SLOTS_PER_KEY = 4;
    /** Draws per bucket before the first level is redrawn; only needed if two keys share a first-level value. */
    private static final int DRAWS_PER_BUCKET = 64;

    private final long mulSymbol;
    private final long mulNode;
    private final int buckets;
    private final int[] bucketStart;   // buckets + 1 slot offsets; bucket b has bucketStart[b+1] - bucketStart[b] slots
    private final int[] bucketSeed;
    private final int[] slotNode;      // -1 if empty
    private final long[] slotSymbol;
    private final int[] slotChild;

    FksChildDictionary(StaticSuffixTree tree) {
        int count = 0;
        for (int p = 0; p < tree.internalCount; p++) {
            count += tree.childStart[p + 1] - firstNonTerminator(tree, p);
        }
        int[] keyNode = new int[count];
        long[] keySymbol = new long[count];
        int[] keyChild = new int[count];
        int e = 0;
        for (int p = 0; p < tree.internalCount; p++) {
            for (int q = firstNonTerminator(tree, p), to = tree.childStart[p + 1]; q < to; q++) {
                int c = tree.children[q];
                keyNode[e] = p;
                keySymbol[e] = tree.firstSymbol(p, c);
                keyChild[e] = c;
                e++;
            }
        }

        final int m = Math.max(1, count);
        SplittableRandom rnd = new SplittableRandom(SEED ^ count);
        int[] bucketOf = new int[count];
        int[] size = new int[m];
        int[] keyStart = new int[m + 1];
        int[] byBucket = new int[count];
        int[] start = new int[m + 1];
        int[] seed = new int[m];
        long mulS;
        long mulN;
        int[] sNode;
        long[] sSymbol;
        int[] sChild;
        draw:
        while (true) {
            mulS = rnd.nextLong() | 1L;
            mulN = rnd.nextLong() | 1L;
            Arrays.fill(size, 0);
            for (int k = 0; k < count; k++) {
                int b = reduce(mulS * keySymbol[k] + mulN * keyNode[k], m);
                bucketOf[k] = b;
                size[b]++;
            }
            long slots = 0;
            for (int b = 0; b < m; b++) {
                slots += (long) size[b] * size[b];
            }
            if (slots > (long) SLOTS_PER_KEY * m) {
                continue;
            }
            for (int b = 0; b < m; b++) {
                keyStart[b + 1] = keyStart[b] + size[b];
                start[b + 1] = start[b] + size[b] * size[b];
            }
            int[] cursor = Arrays.copyOf(keyStart, m);
            for (int k = 0; k < count; k++) {
                byBucket[cursor[bucketOf[k]]++] = k;
            }
            sNode = new int[(int) slots];
            Arrays.fill(sNode, -1);
            sSymbol = new long[(int) slots];
            sChild = new int[(int) slots];
            for (int b = 0; b < m; b++) {
                int keys = size[b];
                if (keys == 0) {
                    continue;
                }
                int slotCount = keys * keys;
                int draws = 0;
                while (true) {
                    if (++draws > DRAWS_PER_BUCKET) {
                        continue draw;
                    }
                    int s = (keys == 1) ? 0 : rnd.nextInt();
                    boolean injective = true;
                    for (int q = keyStart[b]; q < keyStart[b + 1]; q++) {
                        int k = byBucket[q];
                        long v = mulS * keySymbol[k] + mulN * keyNode[k];
                        int slot = start[b] + ((keys == 1) ? 0 : second(v, s, slotCount));
                        if (sNode[slot] != -1) {
                            injective = false;
                            break;
                        }
                        sNode[slot] = keyNode[k];
                        sSymbol[slot] = keySymbol[k];
                        sChild[slot] = keyChild[k];
                    }
                    if (injective) {
                        seed[b] = s;
                        break;
                    }
                    Arrays.fill(sNode, start[b], start[b] + slotCount, -1);
                }
            }
            break;
        }
        this.mulSymbol = mulS;
        this.mulNode = mulN;
        this.buckets = m;
        this.bucketStart = start;
        this.bucketSeed = seed;
        this.slotNode = sNode;
        this.slotSymbol = sSymbol;
        this.slotChild = sChild;
    }

    private static int firstNonTerminator(StaticSuffixTree tree, int p) {
        int from = tree.childStart[p];
        return (from < tree.childStart[p + 1] && tree.isTerminatorEdge(p, tree.children[from])) ? from + 1 : from;
    }

    /** Maps the high 32 bits of {@code v} to {@code [0, range)}. */
    private static int reduce(long v, int range) {
        return (int) (((v >>> 32) * range) >>> 32);
    }

    private static int second(long v, int seed, int range) {
        long z = v ^ (seed * GOLDEN);
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        z ^= z >>> 31;
        return reduce(z, range);
    }

    @Override
    public int child(int node, long symbol) {
        long v = mulSymbol * symbol + mulNode * node;
        int b = reduce(v, buckets);
        int slot = bucketStart[b];
        int slots = bucketStart[b + 1] - slot;
        if (slots == 0) {
            return NONE;
        }
        if (slots > 1) {
            slot += second(v, bucketSeed[b], slots);
        }
        return (slotNode[slot] == node && slotSymbol[slot] == symbol) ? slotChild[slot] : NONE;
    }
}
