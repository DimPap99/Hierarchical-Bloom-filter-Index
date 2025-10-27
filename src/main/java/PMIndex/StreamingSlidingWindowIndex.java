package PMIndex;

import search.Pattern;
import tree.ssws.FarachColtonSuffixTree;
import utilities.AlphabetMapper;
import utilities.PatternResult;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Objects;

/**
 * Skeleton implementation of the streaming sliding-window string index described in the SSWSI paper.
 * At this stage we focus on the amortised segment maintenance policy: incoming tokens are grouped into
 * power-of-two segments, merged eagerly, and trimmed so the maintained suffix never drops below the
 * configured window. Query-related routines will be filled in subsequent iterations.
 */
public class StreamingSlidingWindowIndex implements IPMIndexing {

    private final int windowSize;
    private final AlphabetMapper<String> alphabetMapper;
    private final LinkedList<Segment> segments = new LinkedList<>(); // left -> oldest, right -> newest
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
        nextTokenPosition++;
        segments.addLast(singleton);
        totalStoredLength += singleton.length;
        onSegmentInserted(singleton);
        cascadeMerges();
        trimToWindow();
    }

    private Segment createSingletonSegment(int token, int position) {
        int[] tokens = new int[]{token};
        FarachColtonSuffixTree tree = FarachColtonSuffixTree.build(tokens);
        return new Segment(position, 0, tokens, tree);
    }

    private void cascadeMerges() {
        boolean merged = true;
        while (merged) {
            merged = false;
            if (segments.size() < 2) {
                return;
            }
            Segment last = segments.removeLast();
            Segment previous = segments.peekLast();
            if (previous == null || previous.length != last.length) {
                segments.addLast(last);
                return;
            }
            segments.removeLast(); // remove the previous segment we peeked
            totalStoredLength -= previous.length;
            totalStoredLength -= last.length;
            Segment mergedSegment = mergeSegments(previous, last);
            segments.addLast(mergedSegment);
            totalStoredLength += mergedSegment.length;
            merged = true;
            onSegmentsMerged(mergedSegment);
        }
    }

    private Segment mergeSegments(Segment left, Segment right) {
        int newLength = left.length + right.length;
        int[] tokens = new int[newLength];
        System.arraycopy(left.tokens, 0, tokens, 0, left.length);
        System.arraycopy(right.tokens, 0, tokens, left.length, right.length);
        FarachColtonSuffixTree tree = FarachColtonSuffixTree.build(tokens);
        Segment merged = new Segment(left.startPosition, left.level + 1, tokens, tree);
        rebuildBoundaryAround(merged);
        return merged;
    }

    private void trimToWindow() {
        while (!segments.isEmpty() && totalStoredLength - segments.peekFirst().length >= windowSize) {
            Segment expired = segments.removeFirst();
            totalStoredLength -= expired.length;
            onSegmentExpired(expired);
        }
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(Pattern key) {
        return new ArrayList<>();
    }

    @Override
    public void expire() {
        if (segments.isEmpty()) {
            return;
        }
        Segment expired = segments.removeFirst();
        totalStoredLength -= expired.length;
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
        // placeholder for instrumentation and boundary-tree maintenance
    }

    private void onSegmentsMerged(Segment segment) {
        // placeholder for instrumentation and boundary-tree maintenance
    }

    private void onSegmentExpired(Segment segment) {
        // placeholder for instrumentation and boundary-tree maintenance
    }

    private void rebuildBoundaryAround(Segment segment) {
        // placeholder for boundary tree construction shared across neighbours
    }

    private static final class Segment {
        private final int startPosition;
        private final int level;
        private final int length;
        private final int[] tokens;
        private final FarachColtonSuffixTree suffixTree;

        private Segment(int startPosition, int level, int[] tokens, FarachColtonSuffixTree suffixTree) {
            this.startPosition = startPosition;
            this.level = level;
            this.tokens = tokens;
            this.length = tokens.length;
            this.suffixTree = suffixTree;
        }
    }
}
