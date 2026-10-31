package PMIndex;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import org.openjdk.jol.info.GraphLayout;
import search.Pattern;
import tree.ssws.ChildDictionary;
import tree.ssws.Dc3SuffixTreeBuilder;
import tree.ssws.RangeMax;
import tree.ssws.StaticSuffixTree;
import tree.ssws.SuffixTreeBuilder;
import utilities.PatternResult;
import utilities.TokenHasher;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
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
 * <p><b>Text (Sec. 2 and 3.1).</b> The paper's structure contains the suffix {@code S} once, and its suffix trees
 * represent edge labels by two indices into the text. Here {@code S} is stored once, in one array of {@code 2w}
 * symbols, always starting at index 0, and every tree (segment, boundary, and the query-time tree over R) is a view
 * {@code (offset, length)} into that array, without a copy. When a second segment of size {@code w} forms, the old
 * one is dropped and {@code S} is the new segment alone, so the new segment moves from {@code [w, 2w)} to
 * {@code [0, w)} before its tree is built: {@code w} symbols every {@code w} arrivals. No tree is alive at that
 * moment, and new symbols are written past the end of {@code S}, which no tree views.
 *
 * <p>Interface behaviour as LEGACY: symbols and the n-grams of a pattern are hashed with
 * {@link TokenHasher#hashToPositiveLong}, stream positions are {@code long}, and {@link #report} returns positions
 * relative to the window start ({@link #windowStartPosition()}), in no particular order.
 *
 * <p>Defaults: FKS child dictionaries, as in the paper, and the DC3 + Kasai construction. The paper's Farach-Colton
 * construction ({@link tree.ssws.FarachColtonSuffixTreeBuilder}) is selectable and gives the same trees, but its
 * suffix and LCP arrays take 2-3 times as long, which makes inserts 1.1-1.4 times slower (docs/sswsi-results.md,
 * Phase 2). FASTUTIL is selectable too.
 *
 * <p>With {@code -ea}, every insert ends with {@link #checkInvariants()}.
 */
public final class SlidingSuffixIndex implements IPMIndexing {

    private final int windowSize;
    private final int maxLevel;                 // log2(windowSize)
    private final SuffixTreeBuilder builder;
    private final ChildDictionary.Kind dictionary;

    private final Segment[] segments;           // left to right
    private int count;
    private long nextTokenPosition;

    private final long[] text;                  // S = text[0, nextTokenPosition - base); capacity 2w
    private long base;                          // stream position of text[0]: the leftmost segment's start

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
        this.text = new long[2 * windowSize];
    }

    @Override
    public void insert(String key) {
        if (key == null) {
            return;
        }
        long token = TokenHasher.hashToPositiveLong(key);
        long position = nextTokenPosition++;
        text[(int) (position - base)] = token;      // |S| <= 2w - 1 before this symbol

        // The rightmost j segments have sizes 2^(j-1), ..., 2, 1; with the new symbol they become one segment of
        // size 2^j <= w, and only that segment is built. Its symbols are already in place.
        int j = 0;
        while (j < maxLevel && j < count && segments[count - 1 - j].size == (1 << j)) {
            j++;
        }
        int size = 1 << j;
        for (int k = count - j; k < count; k++) {
            segments[k] = null;
        }
        count -= j;

        long start = position - size + 1;
        if (size == windowSize && count > 0) {
            // The window is now exactly the new segment, so the old segment of size w lies entirely before it.
            if (count != 1 || segments[0].start + windowSize != start) {
                throw new IllegalStateException("a second segment of size w must follow the only other one");
            }
            segments[0] = null;
            count = 0;
            // S is the new segment alone, and no tree is alive: move it to the front.
            System.arraycopy(text, windowSize, text, 0, windowSize);
            base = start;
        }

        StaticSuffixTree tree = builder.build(text, (int) (start - base), size, dictionary);
        RangeMax starts = (size == windowSize) ? new RangeMax(tree) : null;
        Boundary left = (count > 0) ? buildBoundary(start, size) : null;
        segments[count++] = new Segment(start, size, tree, starts, left);
        assert checkInvariants();
    }

    /**
     * B(a, d) for the new rightmost segment d: the last |d| symbols of its left neighbour a (which is larger), then
     * d, with BL' over the suffixes that start in the left half.
     */
    private Boundary buildBoundary(long rightStart, int half) {
        long start = rightStart - half;
        StaticSuffixTree tree = builder.build(text, (int) (start - base), 2 * half, dictionary);
        return new Boundary(start, half, tree, new RangeMax(tree, half));
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
                int shift = (int) (s.start - windowStart);
                if (s.start >= windowStart) {
                    s.tree.reportAll(node, out, shift);
                } else {
                    // Only the leftmost segment, of size w, reaches before the window start.
                    s.starts.reportAtLeast(s.tree.lo(node), s.tree.hi(node), (int) (windowStart - s.start), out, shift);
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

    /** Part 3: a suffix tree over R, built for this query. R is the end of S, so the tree is a view of it. */
    private void reportFromR(long[] pattern, int pivot, long windowStart, IntArrayList out) {
        int m = pattern.length;
        if (pivot == count - 1) {
            return;   // R would be the last m - 1 symbols of the rightmost segment
        }
        long rStart;
        if (pivot < 0) {
            rStart = segments[0].start;
        } else {
            Segment p = segments[pivot];
            rStart = p.start + p.size - (m - 1);
        }
        int length = (int) (nextTokenPosition - rStart);
        if (length < m) {
            return;
        }
        StaticSuffixTree tree = builder.build(text, (int) (rStart - base), length, dictionary);
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

    /** |S|, the number of symbols stored: all of them before {@code w} have arrived, then {@code w} to {@code 2w - 1}. */
    public long storedSymbols() {
        return nextTokenPosition - base;
    }

    /**
     * Checks the structure of Sec. 3.1 and the text layout, and throws {@link AssertionError} if anything is off;
     * returns true otherwise, so that {@code insert} can end with {@code assert checkInvariants()}, which costs
     * nothing without {@code -ea}.
     * <ul>
     *   <li>Segment sizes are distinct powers of two, at most {@code w}, decreasing to the right, and there are at
     *       most {@code log2(w) + 1} of them.</li>
     *   <li>Segments are contiguous, the first starts at {@code text[0]}, the last ends with the newest symbol, and
     *       each starts at a multiple of its size (the binary counter).</li>
     *   <li>Before {@code w} symbols have arrived, the segments hold all of them. From then on, the leftmost has size
     *       {@code w}, {@code w <= |S| <= 2w - 1}, the window starts inside the leftmost segment, and every other
     *       segment lies inside the window.</li>
     *   <li>A segment has a range-maximum structure exactly when its size is {@code w}.</li>
     *   <li>The leftmost segment has no boundary tree; every other segment {@code d} has one over the last |d|
     *       symbols of its left neighbour and d.</li>
     *   <li>Every tree is a view of its own symbols in {@code text}.</li>
     * </ul>
     */
    boolean checkInvariants() {
        long t = nextTokenPosition;
        if (count == 0) {
            check(t == 0 && base == 0, "no segments after %d symbols", t);
            return true;
        }
        check(count <= maxLevel + 1, "%d segments", count);
        check(segments[0].start == base, "the leftmost segment starts at %d, text[0] is %d", segments[0].start, base);
        long expected = base;
        for (int x = 0; x < count; x++) {
            Segment s = segments[x];
            check(s != null, "segment %d is missing", x);
            check(Integer.bitCount(s.size) == 1 && s.size <= windowSize, "segment %d has size %d", x, s.size);
            check(x == 0 || s.size < segments[x - 1].size, "segment %d is not smaller than its left neighbour", x);
            check(s.start == expected, "segment %d starts at %d, not %d", x, s.start, expected);
            check(s.start % s.size == 0, "segment %d of size %d starts at %d", x, s.size, s.start);
            check((s.starts != null) == (s.size == windowSize), "segment %d of size %d, range max %s", x, s.size,
                    s.starts);
            check(s.tree.views(text, (int) (s.start - base), s.size), "segment %d's tree is not a view of it", x);
            if (x == 0) {
                check(s.left == null, "the leftmost segment has a boundary tree");
            } else {
                Boundary b = s.left;
                check(b != null, "segment %d has no boundary tree", x);
                check(b.half == s.size && b.start == s.start - s.size, "segment %d's boundary covers [%d, %d)", x,
                        b.start, b.start + 2L * b.half);
                check(b.tree.views(text, (int) (b.start - base), 2 * b.half), "boundary %d's tree is not a view", x);
            }
            expected += s.size;
        }
        check(expected == t, "the segments end at %d, the stream at %d", expected, t);
        for (int x = count; x < segments.length; x++) {
            check(segments[x] == null, "slot %d past the last segment is not empty", x);
        }
        long stored = t - base;
        check(stored <= 2L * windowSize - 1, "%d symbols stored", stored);
        if (t < windowSize) {
            check(base == 0, "symbols before position %d were dropped before the window filled", base);
        } else {
            check(segments[0].size == windowSize, "the leftmost segment has size %d", segments[0].size);
            check(stored >= windowSize, "only %d symbols stored", stored);
            long windowStart = t - windowSize;
            check(segments[0].start <= windowStart && windowStart < segments[0].start + windowSize,
                    "the window starts at %d, outside the leftmost segment at %d", windowStart, segments[0].start);
            check(count == 1 || segments[1].start >= windowStart, "segment 1 starts before the window");
        }
        return true;
    }

    private static void check(boolean ok, String what, Object... args) {
        if (!ok) {
            throw new AssertionError("SlidingSuffixIndex: " + String.format(what, args));
        }
    }

    /**
     * Bytes reachable from this index, measured by JOL as {@code MemoryBenchmark} measures LEGACY
     * ({@link StreamingSlidingWindowIndex#estimateRetainedBytesWithoutAlphabet()}). Walks the whole object graph.
     * The builders' per-thread work arrays are not reachable from the index, so they are not included.
     */
    public long estimateRetainedBytes() {
        return GraphLayout.parseInstance(this).totalSize();
    }

    /**
     * {@link #estimateRetainedBytes()} split into the stored text, the segment trees and the boundary trees (each
     * with their child dictionaries), the range-maximum structures, and the rest. Each part is the JOL size of
     * everything reachable from it that an earlier part does not reach (a tree reaches the text, a range-maximum
     * structure its tree's suffix array).
     */
    public Map<String, Long> estimateRetainedBytesByPart() {
        List<Object> roots = new ArrayList<>();
        roots.add(text);
        long textBytes = GraphLayout.parseInstance(roots.toArray()).totalSize();
        for (int x = 0; x < count; x++) {
            roots.add(segments[x].tree);
        }
        long withSegments = GraphLayout.parseInstance(roots.toArray()).totalSize();
        for (int x = 1; x < count; x++) {
            roots.add(segments[x].left.tree);
        }
        long withBoundaries = GraphLayout.parseInstance(roots.toArray()).totalSize();
        for (int x = 0; x < count; x++) {
            Segment s = segments[x];
            if (s.starts != null) {
                roots.add(s.starts);
            }
            if (s.left != null) {
                roots.add(s.left.leftStarts);
            }
        }
        long withRangeMax = GraphLayout.parseInstance(roots.toArray()).totalSize();

        Map<String, Long> parts = new LinkedHashMap<>();
        parts.put("text", textBytes);
        parts.put("segment trees", withSegments - textBytes);
        parts.put("boundary trees", withBoundaries - withSegments);
        parts.put("range max", withRangeMax - withBoundaries);
        parts.put("other", estimateRetainedBytes() - withRangeMax);
        return parts;
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
        final StaticSuffixTree tree;
        final RangeMax starts;      // size w only
        final Boundary left;        // B(left neighbour, this); null for the leftmost segment

        Segment(long start, int size, StaticSuffixTree tree, RangeMax starts, Boundary left) {
            this.start = start;
            this.size = size;
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
