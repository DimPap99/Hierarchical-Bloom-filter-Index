package PMIndex;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import search.Pattern;
import tree.ssws.ChildDictionary;
import tree.ssws.Dc3SuffixTreeBuilder;
import tree.ssws.RangeMax;
import tree.ssws.StaticSuffixTree;
import tree.ssws.SuffixTreeBuilder;
import utilities.PatternResult;
import utilities.TokenHasher;

import java.util.ArrayList;
import java.util.Objects;

/**
 * Sliding window string index of Bille, Fischer, Gørtz, Pedersen and Stordalen, "Sliding Window String Indexing in
 * Streams", CPM 2023: the amortized timely structure of Sec. 3.3. The mode called AMORTIZED in the benchmarks;
 * {@link StreamingSlidingWindowIndex} is LEGACY.
 *
 * <p><b>Structure (Sec. 3.1).</b> A suffix {@code S} of the stream is split into segments whose sizes are distinct
 * powers of two, at most {@code w}, decreasing to the right. Once {@code w} symbols have arrived, the leftmost
 * segment has size {@code w} and {@code w <= |S| <= 2w - 1}, so the window starts inside the leftmost segment and
 * every other segment lies entirely inside it. Each segment has a static suffix tree. Each pair of neighbours
 * {@code (a, d)} has a boundary tree over the last {@code |d|} symbols of {@code a} followed by {@code d}, with a
 * range-maximum structure over the starts of its suffixes that start in the left half (the paper's BL').
 *
 * <p><b>Arrivals (Sec. 3.3).</b> The new symbol and the rightmost {@code j} segments, of sizes
 * {@code 2^(j-1), ..., 2, 1}, become one segment {@code d} of size {@code 2^j <= w}. Only {@code d} is built, with
 * {@code B(a, d)} for its left neighbour {@code a}, never the intermediate merges of the cascade; {@code d} is the
 * rightmost segment, so it has no right boundary. When {@code d} has size {@code w}, the segment of size {@code w}
 * before it lies entirely before the window and is dropped; two segments of size {@code w} are never merged.
 *
 * <p><b>Queries (Sec. 3.2)</b> for a pattern of length {@code m}, in three disjoint parts, so no de-duplication:
 * <ol>
 *   <li>Occurrences inside a segment: every segment of size at least {@code m}.</li>
 *   <li>Occurrences that cross one boundary, between {@code s_(x+1)} and {@code s_x}, only where
 *       {@code |s_x| >= m}: the suffixes of {@code B_x} that start at or after {@code |s_x| - m + 1} in its left
 *       half. Such an occurrence ends inside {@code s_x}.</li>
 *   <li>The rest, all of which cross a boundary into a segment shorter than {@code m}: a suffix tree built at query
 *       time over {@code R}, the last {@code m - 1} symbols of the pivot (the rightmost segment of size at least
 *       {@code m}) and everything to its right. Every occurrence in {@code R} crosses the pivot's right end. With no
 *       pivot, {@code R} is everything stored.</li>
 * </ol>
 * Range-maximum structures exist only on segments of size {@code w}: the leftmost segment is the only one that can
 * start before the window, and it always has size {@code w} once the window is full.
 *
 * <p><b>Not in the paper: the window bound on boundaries.</b> Occurrences across the boundary right of the leftmost
 * segment can start before the window, because that segment is only partly inside it. The reporting threshold there
 * is {@code max(|s_x| - m + 1, window start - B.start)}; both are lower bounds on the start, so the range-maximum
 * recursion still stops correctly.
 *
 * <p>Interface behaviour as LEGACY: symbols and the n-grams of a pattern are hashed with
 * {@link TokenHasher#hashToPositiveLong}, stream positions are {@code long}, and {@link #report} returns positions
 * relative to the window start ({@link #windowStartPosition()}), in no particular order.
 *
 * <p>Defaults: FKS child dictionaries, as in the paper, and the DC3 + Kasai construction. The paper's Farach-Colton
 * construction ({@link tree.ssws.FarachColtonSuffixTreeBuilder}) is selectable and gives the same trees, but its
 * suffix and LCP arrays take 2-3 times as long, which makes inserts 1.1-1.4 times slower (docs/sswsi-results.md,
 * Phase 2). FASTUTIL is selectable too. Every segment and every boundary tree still has its own copy of its
 * symbols (docs/sswsi-baseline-handoff.md, Phase 3).
 */
public final class SlidingSuffixIndex implements IPMIndexing {

    private final int windowSize;
    private final int maxLevel;                 // log2(windowSize)
    private final SuffixTreeBuilder builder;
    private final ChildDictionary.Kind dictionary;

    private final Segment[] segments;           // left to right
    private int count;
    private long nextTokenPosition;

    private final IntArrayList hits = new IntArrayList();

    /** FKS child dictionaries (the paper's) and the DC3 construction. */
    public SlidingSuffixIndex(int windowSize) {
        this(windowSize, ChildDictionary.Kind.FKS);
    }

    public SlidingSuffixIndex(int windowSize, ChildDictionary.Kind dictionary) {
        this(windowSize, new Dc3SuffixTreeBuilder(), dictionary);
    }

    public SlidingSuffixIndex(int windowSize, SuffixTreeBuilder builder, ChildDictionary.Kind dictionary) {
        if (windowSize <= 0 || Integer.bitCount(windowSize) != 1) {
            throw new IllegalArgumentException("windowSize must be a power of two, got " + windowSize);
        }
        this.windowSize = windowSize;
        this.maxLevel = Integer.numberOfTrailingZeros(windowSize);
        this.builder = Objects.requireNonNull(builder, "builder");
        this.dictionary = Objects.requireNonNull(dictionary, "dictionary");
        this.segments = new Segment[maxLevel + 2];
    }

    @Override
    public void insert(String key) {
        if (key == null) {
            return;
        }
        long token = TokenHasher.hashToPositiveLong(key);
        long position = nextTokenPosition++;

        // The rightmost j segments have sizes 2^(j-1), ..., 2, 1; with the new symbol they become one segment of
        // size 2^j <= w, and only that segment is built.
        int j = 0;
        while (j < maxLevel && j < count && segments[count - 1 - j].size == (1 << j)) {
            j++;
        }
        int size = 1 << j;
        long[] tokens = new long[size];
        int filled = 0;
        for (int k = count - j; k < count; k++) {
            Segment s = segments[k];
            System.arraycopy(s.tokens, 0, tokens, filled, s.size);
            filled += s.size;
            segments[k] = null;
        }
        tokens[filled] = token;
        count -= j;

        long start = position - size + 1;
        if (size == windowSize && count > 0) {
            // The window is now exactly the new segment, so the old segment of size w lies entirely before it.
            if (count != 1 || segments[0].start + windowSize != start) {
                throw new IllegalStateException("a second segment of size w must follow the only other one");
            }
            segments[0] = null;
            count = 0;
        }

        StaticSuffixTree tree = builder.build(tokens, 0, size, dictionary);
        RangeMax starts = (size == windowSize) ? new RangeMax(tree) : null;
        Boundary left = (count > 0) ? buildBoundary(segments[count - 1], tokens, start) : null;
        segments[count++] = new Segment(start, size, tokens, tree, starts, left);
    }

    /** B(a, d): the last |d| symbols of a, then d, with BL' over the suffixes that start in the left half. */
    private Boundary buildBoundary(Segment a, long[] right, long rightStart) {
        int half = right.length;
        long[] text = new long[2 * half];
        System.arraycopy(a.tokens, a.size - half, text, 0, half);
        System.arraycopy(right, 0, text, half, half);
        StaticSuffixTree tree = builder.build(text, 0, text.length, dictionary);
        return new Boundary(rightStart - half, half, tree, new RangeMax(tree, half));
    }

    /** Absolute stream position of the first symbol in the current window. */
    public long windowStartPosition() {
        return Math.max(0L, nextTokenPosition - windowSize);
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(Pattern key) {
        if (key == null || key.nGramArr == null) {
            return new ArrayList<>();
        }
        long[] pattern = toPatternTokens(key);
        int m = pattern.length;
        long windowStart = windowStartPosition();
        if (m == 0 || m > nextTokenPosition - windowStart) {
            return new ArrayList<>();
        }
        IntArrayList out = hits;
        out.clear();

        int pivot = -1;
        for (int x = 0; x < count && segments[x].size >= m; x++) {
            Segment s = segments[x];
            pivot = x;

            int node = s.tree.locus(pattern);
            if (node != StaticSuffixTree.NONE) {
                int base = (int) (s.start - windowStart);
                if (s.start >= windowStart) {
                    s.tree.reportAll(node, out, base);
                } else {
                    // Only the leftmost segment, of size w, reaches before the window start.
                    s.starts.reportAtLeast(s.tree.lo(node), s.tree.hi(node), (int) (windowStart - s.start), out, base);
                }
            }

            Boundary b = s.left;
            if (b != null) {
                int bn = b.tree.locus(pattern);
                if (bn != StaticSuffixTree.NONE) {
                    long threshold = Math.max(b.half - m + 1, windowStart - b.start);
                    b.leftStarts.reportAtLeast(b.tree.lo(bn), b.tree.hi(bn), (int) threshold, out,
                            (int) (b.start - windowStart));
                }
            }
        }
        reportFromR(pattern, pivot, windowStart, out);

        ArrayList<Integer> result = new ArrayList<>(out.size());
        for (int i = 0; i < out.size(); i++) {
            result.add(out.getInt(i));
        }
        return result;
    }

    /** Part 3: a suffix tree over R, built for this query. */
    private void reportFromR(long[] pattern, int pivot, long windowStart, IntArrayList out) {
        int m = pattern.length;
        if (pivot == count - 1) {
            return;   // R would be the last m - 1 symbols of the rightmost segment
        }
        int prefix = (pivot < 0) ? 0 : m - 1;
        int length = prefix;
        for (int x = pivot + 1; x < count; x++) {
            length += segments[x].size;
        }
        if (length < m) {
            return;
        }
        long[] r = new long[length];
        long rStart;
        if (pivot < 0) {
            rStart = segments[0].start;
        } else {
            Segment p = segments[pivot];
            rStart = p.start + p.size - prefix;
            System.arraycopy(p.tokens, p.size - prefix, r, 0, prefix);
        }
        int filled = prefix;
        for (int x = pivot + 1; x < count; x++) {
            Segment s = segments[x];
            System.arraycopy(s.tokens, 0, r, filled, s.size);
            filled += s.size;
        }
        StaticSuffixTree tree = builder.build(r, 0, length, dictionary);
        int node = tree.locus(pattern);
        if (node != StaticSuffixTree.NONE) {
            // R starts before the window only when the pivot is the leftmost segment.
            tree.reportFrom(node, (int) Math.max(0L, windowStart - rStart), out, (int) (rStart - windowStart));
        }
    }

    private long[] toPatternTokens(Pattern pattern) {
        String[] grams = pattern.nGramArr;
        if (grams == null) {
            return new long[0];
        }
        long[] tokens = new long[grams.length];
        for (int i = 0; i < grams.length; i++) {
            long hashed = TokenHasher.hashToPositiveLong(grams[i]);
            tokens[i] = hashed;
            if (pattern.nGramToLong != null && i < pattern.nGramToLong.length) {
                pattern.nGramToLong[i] = hashed;
            }
        }
        return tokens;
    }

    /** Segment sizes from left to right, for tests. */
    int[] segmentSizes() {
        int[] sizes = new int[count];
        for (int x = 0; x < count; x++) {
            sizes[x] = segments[x].size;
        }
        return sizes;
    }

    /** Nothing to do: the window slides on every insert. */
    @Override
    public void expire() {
    }

    @Override
    public ArrayList<Long> getAvgTimes(Pattern pat) {
        return new ArrayList<>();
    }

    @Override
    public PatternResult getLatestStats() {
        return new PatternResult(0.0, 0, 0, null, 0, 0.0, 0, 0);
    }

    @Override
    public int getTokenId(String key) {
        Objects.requireNonNull(key, "key");
        return TokenHasher.hashToPositiveInt(key);
    }

    private static final class Segment {
        final long start;           // stream position of the first symbol
        final int size;
        final long[] tokens;
        final StaticSuffixTree tree;
        final RangeMax starts;      // size w only
        final Boundary left;        // B(left neighbour, this); null for the leftmost segment

        Segment(long start, int size, long[] tokens, StaticSuffixTree tree, RangeMax starts, Boundary left) {
            this.start = start;
            this.size = size;
            this.tokens = tokens;
            this.tree = tree;
            this.starts = starts;
            this.left = left;
        }
    }

    private static final class Boundary {
        final long start;           // stream position of the first symbol of the boundary tree's text
        final int half;             // size of the right segment; the tree covers half symbols on each side
        final StaticSuffixTree tree;
        final RangeMax leftStarts;  // BL'

        Boundary(long start, int half, StaticSuffixTree tree, RangeMax leftStarts) {
            this.start = start;
            this.half = half;
            this.tree = tree;
            this.leftStarts = leftStarts;
        }
    }
}
