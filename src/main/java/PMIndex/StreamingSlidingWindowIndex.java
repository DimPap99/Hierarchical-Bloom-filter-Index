package PMIndex;

import org.openjdk.jol.info.GraphLayout;
import search.Pattern;
import tree.ssws.SuffixTree;
import utilities.PatternResult;
import utilities.TokenHasher;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

// Streaming sliding-window index inspired by the SSWSI hierarchy.
// Uses power-of-two segments with log-structured merges and multi-stage querying.
public class StreamingSlidingWindowIndex implements IPMIndexing {

    private final int windowSize;

    private Segment head;
    private Segment tail;
    private int segmentCount = 0;
    private int totalStoredLength = 0;
    private int nextTokenPosition = 0;

    public StreamingSlidingWindowIndex(int windowSize, int expectedAlphabetSize) {
        if (windowSize <= 0) {
            throw new IllegalArgumentException("windowSize must be positive");
        }
        this.windowSize = windowSize;
    }

    public StreamingSlidingWindowIndex(int windowSize) {
        this(windowSize, 1_024);
    }

    public long estimateRetainedBytesWithoutAlphabet() {
        long total = GraphLayout.parseInstance(this).totalSize();
        return total;
    }

    public long estimateTokenDictionaryBytes() {
        return 0L;
    }

    public long estimateSuffixTreeRemapBytes() {
        long total = 0L;
        HashSet<Boundary> visited = new HashSet<>();
        Segment current = head;
        while (current != null) {
            total += current.suffixTree.estimateTokenRemapBytes();
            Boundary boundary = current.rightBoundary;
            if (boundary != null && visited.add(boundary)) {
                total += boundary.tree.estimateTokenRemapBytes();
            }
            current = current.next;
        }
        return total;
    }

    @Override
    public void insert(String key) {
        if (key == null) {
            return;
        }
        long token = TokenHasher.hashToPositiveLong(key);
        Segment singleton = createSingletonSegment(token, nextTokenPosition);
        appendSegment(singleton);
        nextTokenPosition++;
        totalStoredLength += singleton.length;
        onSegmentInserted(singleton);
        cascadeMergesFromTail();
        trimToWindow();
    }

    private Segment createSingletonSegment(long token, int position) {
        long[] tokens = new long[]{token};
        SuffixTree tree = SuffixTree.build(tokens, windowSize);
        return new Segment(position, 0, tokens, tree);
    }

    private void appendSegment(Segment segment) {
        if (tail == null) {
            head = tail = segment;
        } else {
            tail.next = segment;
            segment.prev = tail;
            buildBoundary(tail, segment);
            tail = segment;
        }
        segmentCount++;
    }

    private void cascadeMergesFromTail() {
        Segment current = tail;
        while (current != null && current.prev != null && current.length == current.prev.length) {
            Segment left = current.prev;
            Segment merged = mergeSegments(left, current);
            replaceWithMerged(left, current, merged);
            onSegmentsMerged(merged);
            current = merged;
        }
    }

    private Segment mergeSegments(Segment left, Segment right) {
        int newLength = left.length + right.length;
        long[] tokens = new long[newLength];
        System.arraycopy(left.tokens, 0, tokens, 0, left.length);
        System.arraycopy(right.tokens, 0, tokens, left.length, right.length);
        SuffixTree tree = SuffixTree.build(tokens, windowSize);
        return new Segment(left.startPosition, left.level + 1, tokens, tree);
    }

    private void replaceWithMerged(Segment left, Segment right, Segment merged) {
        Segment prev = left.prev;
        Segment next = right.next;

        Boundary leftLeftBoundary = left.leftBoundary;
        Boundary leftRightBoundary = left.rightBoundary;
        Boundary rightRightBoundary = right.rightBoundary;

        clearBoundary(leftLeftBoundary);
        clearBoundary(leftRightBoundary); // also clears right.leftBoundary
        clearBoundary(rightRightBoundary);

        if (prev != null) {
            prev.next = merged;
            merged.prev = prev;
            buildBoundary(prev, merged);
        } else {
            merged.prev = null;
            head = merged;
        }

        if (next != null) {
            next.prev = merged;
            merged.next = next;
            buildBoundary(merged, next);
        } else {
            merged.next = null;
            tail = merged;
        }

        left.prev = left.next = null;
        right.prev = right.next = null;
        segmentCount--; // two segments became one.

        releaseBoundary(leftLeftBoundary);
        releaseBoundary(leftRightBoundary);
        releaseBoundary(rightRightBoundary);
        releaseSegment(left);
        releaseSegment(right);
    }

    private void clearBoundary(Boundary boundary) {
        if (boundary == null) {
            return;
        }
        if (boundary.left != null && boundary.left.rightBoundary == boundary) {
            boundary.left.rightBoundary = null;
        }
        if (boundary.right != null && boundary.right.leftBoundary == boundary) {
            boundary.right.leftBoundary = null;
        }
        releaseBoundary(boundary);
    }

    private void trimToWindow() {
        while (head != null && totalStoredLength - head.length >= windowSize) {
            Segment expired = head;
            totalStoredLength -= expired.length;
            removeHead();
            onSegmentExpired(expired);
            releaseSegment(expired);
        }
    }

    private void removeHead() {
        if (head == null) {
            return;
        }
        Segment newHead = head.next;
        clearBoundary(head.rightBoundary);
        head.next = null;
        if (newHead != null) {
            newHead.prev = null;
            newHead.leftBoundary = null;
        } else {
            tail = null;
        }
        head = newHead;
        segmentCount--;
    }

    private void buildBoundary(Segment left, Segment right) {
        if (left == null || right == null) {
            return;
        }
        int span = right.length;
        int leftSpan = Math.min(left.length, span);
        long[] tokens = new long[leftSpan + right.length];
        System.arraycopy(left.tokens, left.length - leftSpan, tokens, 0, leftSpan);
        System.arraycopy(right.tokens, 0, tokens, leftSpan, right.length);
        SuffixTree tree = SuffixTree.build(tokens, windowSize);
        int boundaryStart = left.startPosition + left.length - leftSpan;
        int boundaryPosition = left.startPosition + left.length;
        Boundary boundary = new Boundary(left, right, boundaryStart, boundaryPosition, tree);
        left.rightBoundary = boundary;
        right.leftBoundary = boundary;
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(Pattern key) {
        ArrayList<Integer> noResult = new ArrayList<>();
        if (key == null || key.nGramArr == null) {
            return noResult;
        }
        long[] patternTokens = toPatternTokens(key);
        if (patternTokens.length == 0) {
            return noResult;
        }
        if (tail == null) {
            return noResult;
        }
        int patternLength = patternTokens.length;
        int windowEnd = nextTokenPosition;
        int windowStart = Math.max(0, windowEnd - windowSize);

        Set<Integer> occurrences = new LinkedHashSet<>();

        collectSegmentMatches(patternTokens, patternLength, windowStart, windowEnd, occurrences);
        collectBoundaryMatches(patternTokens, patternLength, windowStart, windowEnd, occurrences);
        collectSuffixMatches(patternTokens, patternLength, windowStart, windowEnd, occurrences);

        ArrayList<Integer> result = new ArrayList<>(occurrences);
//        Collections.sort(result);
        return result;
    }

    private void collectSegmentMatches(long[] patternTokens, int patternLength, int windowStart,
                                       int windowEnd, Set<Integer> occurrences) {
        Segment current = head;
        while (current != null) {
            if (current.length >= patternLength) {
                List<Integer> matches = current.suffixTree.findOccurrences(patternTokens);
                for (int local : matches) {
                    int global = current.startPosition + local;
                    if (global >= windowStart && global + patternLength <= windowEnd) {
                        occurrences.add(global);
                    }
                }
            }
            current = current.next;
        }
    }

    private void collectBoundaryMatches(long[] patternTokens, int patternLength, int windowStart,
                                        int windowEnd, Set<Integer> occurrences) {
        Segment current = head;
        while (current != null) {
            Boundary boundary = current.rightBoundary;
            if (boundary != null) {
                List<Integer> matches = boundary.tree.findOccurrences(patternTokens);
                for (int local : matches) {
                    int global = boundary.globalStart + local;
                    if (global < boundary.boundaryPosition && global + patternLength > boundary.boundaryPosition) {
                        if (global >= windowStart && global + patternLength <= windowEnd) {
                            occurrences.add(global);
                        }
                    }
                }
            }
            current = current.next;
        }
    }

    private void collectSuffixMatches(long[] patternTokens, int patternLength, int windowStart,
                                      int windowEnd, Set<Integer> occurrences) {
        if (patternLength == 0) {
            return;
        }
        Segment pivot = findPivotSegment(patternLength);
        QuerySubstring substring;
        if (pivot == null) {
            substring = buildWindowSubstring();
        } else if (pivot.next == null) {
            return; // no smaller segments to the right; nothing crosses multiple boundaries
        } else {
            substring = buildTailSubstring(pivot, patternLength);
        }
        if (substring.tokens.length < patternLength) {
            return;
        }
        SuffixTree tree = SuffixTree.build(substring.tokens, windowSize);
        List<Integer> matches = tree.findOccurrences(patternTokens);
        for (int local : matches) {
            int global = substring.globalStart + local;
            if (global >= windowStart && global + patternLength <= windowEnd) {
                occurrences.add(global);
            }
        }
    }

    private Segment findPivotSegment(int patternLength) {
        Segment current = tail;
        while (current != null) {
            if (current.length >= patternLength) {
                return current;
            }
            current = current.prev;
        }
        return null;
    }

    private QuerySubstring buildTailSubstring(Segment pivot, int patternLength) {
        int pivotSuffix = Math.max(0, Math.min(patternLength - 1, pivot.length));
        int total = pivotSuffix;
        Segment current = pivot.next;
        while (current != null) {
            total += current.length;
            current = current.next;
        }
        if (total == 0) {
            return QuerySubstring.empty();
        }
        long[] tokens = new long[total];
        int offset = 0;
        if (pivotSuffix > 0) {
            System.arraycopy(pivot.tokens, pivot.length - pivotSuffix, tokens, 0, pivotSuffix);
            offset = pivotSuffix;
        }
        current = pivot.next;
        while (current != null) {
            System.arraycopy(current.tokens, 0, tokens, offset, current.length);
            offset += current.length;
            current = current.next;
        }
        int globalStart = pivot.startPosition + pivot.length - pivotSuffix;
        return new QuerySubstring(tokens, globalStart);
    }

    private QuerySubstring buildWindowSubstring() {
        int windowEnd = nextTokenPosition;
        int windowStart = Math.max(0, windowEnd - windowSize);
        int length = windowEnd - windowStart;
        if (length <= 0 || tail == null) {
            return QuerySubstring.empty();
        }
        long[] tokens = new long[length];
        int writeIndex = length;
        Segment current = tail;
        while (current != null && writeIndex > 0) {
            int segStart = current.startPosition;
            int segEnd = segStart + current.length;
            int copyStart = Math.max(segStart, windowStart);
            int copyEnd = Math.min(segEnd, windowEnd);
            int amount = copyEnd - copyStart;
            if (amount > 0) {
                writeIndex -= amount;
                int from = copyStart - segStart;
                System.arraycopy(current.tokens, from, tokens, writeIndex, amount);
            }
            current = current.prev;
        }
        if (writeIndex == 0) {
            return new QuerySubstring(tokens, windowStart);
        }
        int actualLength = length - writeIndex;
        if (actualLength <= 0) {
            return QuerySubstring.empty();
        }
        long[] trimmed = new long[actualLength];
        System.arraycopy(tokens, writeIndex, trimmed, 0, actualLength);
        return new QuerySubstring(trimmed, windowStart + writeIndex);
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

    @Override
    public void expire() {
        if (head == null) {
            return;
        }
        totalStoredLength -= head.length;
        Segment expired = head;
        removeHead();
        onSegmentExpired(expired);
        releaseSegment(expired);
        trimToWindow();
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

    private void onSegmentInserted(Segment segment) {
        // hook for instrumentation or logging
    }

    private void onSegmentsMerged(Segment segment) {
        // hook for instrumentation or logging
    }

    private void onSegmentExpired(Segment segment) {
        // hook for instrumentation or logging
    }

    private void releaseSegment(Segment segment) {
        if (segment == null) {
            return;
        }
        wipeArray(segment.tokens);
        if (segment.suffixTree != null) {
            segment.suffixTree.release();
        }
        releaseBoundary(segment.leftBoundary);
        releaseBoundary(segment.rightBoundary);
        segment.leftBoundary = null;
        segment.rightBoundary = null;
        segment.prev = null;
        segment.next = null;
    }

    private void releaseBoundary(Boundary boundary) {
        if (boundary == null) {
            return;
        }
        boundary.tree.release();
    }

    private void wipeArray(long[] data) {
        if (data != null) {
            Arrays.fill(data, 0L);
        }
    }

    private static final class Segment {
        private final int startPosition;
        private final int level;
        private final int length;
        private final long[] tokens;
        private final SuffixTree suffixTree;
        private Segment prev;
        private Segment next;
        private Boundary leftBoundary;
        private Boundary rightBoundary;

        private Segment(int startPosition, int level, long[] tokens, SuffixTree suffixTree) {
            this.startPosition = startPosition;
            this.level = level;
            this.tokens = tokens;
            this.length = tokens.length;
            this.suffixTree = suffixTree;
        }
    }

    private static final class Boundary {
        private final Segment left;
        private final Segment right;
        private final int globalStart;
        private final int boundaryPosition;
        private final SuffixTree tree;

        private Boundary(Segment left, Segment right, int globalStart,
                         int boundaryPosition, SuffixTree tree) {
            this.left = left;
            this.right = right;
            this.globalStart = globalStart;
            this.boundaryPosition = boundaryPosition;
            this.tree = tree;
        }
    }

    private static final class QuerySubstring {
        private static final QuerySubstring EMPTY = new QuerySubstring(new long[0], 0);
        private final long[] tokens;
        private final int globalStart;

        private QuerySubstring(long[] tokens, int globalStart) {
            this.tokens = tokens;
            this.globalStart = globalStart;
        }

        private static QuerySubstring empty() {
            return EMPTY;
        }
    }
}
