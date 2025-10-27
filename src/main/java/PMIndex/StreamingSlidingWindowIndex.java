package PMIndex;

import search.Pattern;
import tree.ssws.FarachColtonSuffixTree;
import utilities.AlphabetMapper;
import utilities.PatternResult;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;

/**
 * Streaming sliding-window index inspired by the SSWSI hierarchy. The implementation keeps
 * the amortised update policy from the paper (power-of-two segments, log-structured merges)
 * and exposes multi-stage querying: occurrences fully contained in segments, occurrences
 * crossing a single boundary, and occurrences that require the on-demand suffix rebuild for
 * the tail region.
 */
public class StreamingSlidingWindowIndex implements IPMIndexing {

    private final int windowSize;
    private final AlphabetMapper<String> alphabetMapper;

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
        this.alphabetMapper = new AlphabetMapper<>(Math.max(expectedAlphabetSize, 16));
    }

    public StreamingSlidingWindowIndex(int windowSize) {
        this(windowSize, 1_024);
    }

    @Override
    public void insert(String key) {
        if (key == null) {
            return;
        }
        int token = alphabetMapper.getId(key);
        Segment singleton = createSingletonSegment(token, nextTokenPosition);
        appendSegment(singleton);
        nextTokenPosition++;
        totalStoredLength += singleton.length;
        onSegmentInserted(singleton);
        cascadeMergesFromTail();
        trimToWindow();
    }

    private Segment createSingletonSegment(int token, int position) {
        int[] tokens = new int[]{token};
        FarachColtonSuffixTree tree = FarachColtonSuffixTree.build(tokens);
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
        int[] tokens = new int[newLength];
        System.arraycopy(left.tokens, 0, tokens, 0, left.length);
        System.arraycopy(right.tokens, 0, tokens, left.length, right.length);
        FarachColtonSuffixTree tree = FarachColtonSuffixTree.build(tokens);
        return new Segment(left.startPosition, left.level + 1, tokens, tree);
    }

    private void replaceWithMerged(Segment left, Segment right, Segment merged) {
        Segment prev = left.prev;
        Segment next = right.next;

        clearBoundary(left.leftBoundary);
        clearBoundary(left.rightBoundary); // also clears right.leftBoundary
        clearBoundary(right.rightBoundary);

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
        segmentCount--; // two segments became one
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
    }

    private void trimToWindow() {
        while (head != null && totalStoredLength - head.length >= windowSize) {
            Segment expired = head;
            totalStoredLength -= expired.length;
            removeHead();
            onSegmentExpired(expired);
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
        int[] tokens = new int[leftSpan + right.length];
        System.arraycopy(left.tokens, left.length - leftSpan, tokens, 0, leftSpan);
        System.arraycopy(right.tokens, 0, tokens, leftSpan, right.length);
        FarachColtonSuffixTree tree = FarachColtonSuffixTree.build(tokens);
        int boundaryStart = left.startPosition + left.length - leftSpan;
        int boundaryPosition = left.startPosition + left.length;
        Boundary boundary = new Boundary(left, right, tokens, boundaryStart, boundaryPosition, tree);
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
        int[] patternTokens = toPatternTokens(key);
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
        Collections.sort(result);
        return result;
    }

    private void collectSegmentMatches(int[] patternTokens, int patternLength, int windowStart,
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

    private void collectBoundaryMatches(int[] patternTokens, int patternLength, int windowStart,
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

    private void collectSuffixMatches(int[] patternTokens, int patternLength, int windowStart,
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
        FarachColtonSuffixTree tree = FarachColtonSuffixTree.build(substring.tokens);
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
        int[] tokens = new int[total];
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
        int[] tokens = new int[length];
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
        int[] trimmed = new int[actualLength];
        System.arraycopy(tokens, writeIndex, trimmed, 0, actualLength);
        return new QuerySubstring(trimmed, windowStart + writeIndex);
    }

    private int[] toPatternTokens(Pattern pattern) {
        String[] grams = pattern.nGramArr;
        if (grams == null) {
            return new int[0];
        }
        int[] tokens = new int[grams.length];
        for (int i = 0; i < grams.length; i++) {
            tokens[i] = alphabetMapper.getId(grams[i]);
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
        return alphabetMapper.getId(key);
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

    private static final class Segment {
        private final int startPosition;
        private final int level;
        private final int length;
        private final int[] tokens;
        private final FarachColtonSuffixTree suffixTree;
        private Segment prev;
        private Segment next;
        private Boundary leftBoundary;
        private Boundary rightBoundary;

        private Segment(int startPosition, int level, int[] tokens, FarachColtonSuffixTree suffixTree) {
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
        private final int[] tokens;
        private final int globalStart;
        private final int boundaryPosition;
        private final FarachColtonSuffixTree tree;

        private Boundary(Segment left, Segment right, int[] tokens, int globalStart,
                          int boundaryPosition, FarachColtonSuffixTree tree) {
            this.left = left;
            this.right = right;
            this.tokens = tokens;
            this.globalStart = globalStart;
            this.boundaryPosition = boundaryPosition;
            this.tree = tree;
        }
    }

    private static final class QuerySubstring {
        private static final QuerySubstring EMPTY = new QuerySubstring(new int[0], 0);
        private final int[] tokens;
        private final int globalStart;

        private QuerySubstring(int[] tokens, int globalStart) {
            this.tokens = tokens;
            this.globalStart = globalStart;
        }

        private static QuerySubstring empty() {
            return EMPTY;
        }
    }
}
