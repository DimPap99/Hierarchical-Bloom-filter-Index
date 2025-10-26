package PMIndex;

import search.Pattern;
import utilities.PatternResult;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;
import java.util.LinkedHashSet;

/**
 * Timely (\u03b4 = 0) streaming sliding-window string index built on top of {@link SuffixTreeIndex} segments.
 * <p>
 * This implementation follows the SSWSI hierarchy described in the project specification. Incoming n-gram “tokens”
 * are grouped into power-of-two sized segments. Each segment keeps its own suffix tree; boundaries between adjacent
 * segments are indexed by dedicated boundary trees so queries that straddle the border can be answered without
 * touching unrelated text. The hierarchy ensures that a token participates in at most \u039f(log w) merges over its
 * lifetime, yielding the desired amortised update bound while using \u039f(w) memory overall.
 * <p>
 * The index works with the existing {@link IPMIndexing} contract where {@link #insert(String)} receives the next
 * token (typically an n-gram string) and {@link #report(Pattern)} consumes search patterns built by the existing
 * query pipeline.
 */
public final class TimelySuffixWindowIndex implements IPMIndexing {

    private final int windowSize;

    /** Oldest segment first, newest segment last. */
    private final ArrayList<Segment> segments = new ArrayList<>();

    /** Parallel list of boundary trees; boundaryTrees.get(i) corresponds to segments get(i) and get(i + 1). */
    private final ArrayList<BoundaryTree> boundaryTrees = new ArrayList<>();

    private int nextGlobalIndex = 0; // monotonically increasing token counter
    private int windowStartIndex = 0; // first token index still inside the logical window

    public TimelySuffixWindowIndex(int windowSize) {
        if (windowSize <= 0) {
            throw new IllegalArgumentException("windowSize must be positive");
        }
        this.windowSize = windowSize;
    }

    @Override
    public void insert(String key) {
        if (key == null) {
            return;
        }
        appendNewToken(key);
        mergeCarriesFromTail();
        dropExpiredSegments();
    }

    private void appendNewToken(String token) {
        ArrayList<String> tokens = new ArrayList<>(1);
        tokens.add(token);
        SuffixTreeIndex suffixTree = buildSuffixTree(tokens);
        Segment segment = new Segment(tokens, nextGlobalIndex, nextGlobalIndex, suffixTree);
        segments.add(segment);
        if (segments.size() >= 2) {
            rebuildBoundaryAt(segments.size() - 2);
        }
        nextGlobalIndex++;
        windowStartIndex = Math.max(0, nextGlobalIndex - windowSize);
    }

    private void mergeCarriesFromTail() {
        boolean merged;
        do {
            merged = false;
            int size = segments.size();
            if (size < 2) {
                return;
            }

            Segment right = segments.get(size - 1);
            Segment left = segments.get(size - 2);
            if (left.length() == right.length()) {
                int mergedIndex = size - 2;
                removeBoundaryAt(mergedIndex);
                Segment mergedSegment = merge(left, right);
                segments.set(mergedIndex, mergedSegment);
                segments.remove(size - 1);
                rebuildBoundaryAt(mergedIndex);
                rebuildBoundaryAt(mergedIndex - 1);
                merged = true;
            }
        } while (merged);
        trimBoundaryList();
    }

    private Segment merge(Segment left, Segment right) {
        ArrayList<String> mergedTokens = new ArrayList<>(left.tokens.size() + right.tokens.size());
        mergedTokens.addAll(left.tokens);
        mergedTokens.addAll(right.tokens);
        SuffixTreeIndex mergedTree = buildSuffixTree(mergedTokens);
        return new Segment(mergedTokens, left.globalStart, right.globalEnd, mergedTree);
    }

    private void dropExpiredSegments() {
        while (!segments.isEmpty()) {
            Segment oldest = segments.get(0);
            if (oldest.globalEnd < windowStartIndex) {
                segments.remove(0);
                removeBoundaryAt(0);
            } else {
                break;
            }
        }
        trimBoundaryList();
        if (segments.size() >= 2) {
            rebuildBoundaryAt(0);
        }
    }

    private void rebuildBoundaryAt(int leftIndex) {
        if (leftIndex < 0 || leftIndex >= segments.size() - 1) {
            return;
        }
        Segment left = segments.get(leftIndex);
        Segment right = segments.get(leftIndex + 1);
        if (right.tokens.isEmpty()) {
            setBoundary(leftIndex, null);
            return;
        }
        int radius = right.tokens.size();
        List<String> leftSlice = sliceTail(left.tokens, radius);
        List<String> rightSlice = sliceHead(right.tokens, radius);
        if (leftSlice.isEmpty() && rightSlice.isEmpty()) {
            setBoundary(leftIndex, null);
            return;
        }
        ArrayList<String> boundaryText = new ArrayList<>(leftSlice.size() + rightSlice.size());
        boundaryText.addAll(leftSlice);
        boundaryText.addAll(rightSlice);
        if (boundaryText.isEmpty()) {
            setBoundary(leftIndex, null);
            return;
        }

        SuffixTreeIndex tree = buildSuffixTree(boundaryText);
        int leftStart = left.globalEnd - leftSlice.size() + 1;
        int splitOffset = leftSlice.size();
        setBoundary(leftIndex, new BoundaryTree(left, right, leftStart, splitOffset, tree));
    }

    private void removeBoundaryAt(int leftIndex) {
        if (leftIndex < 0 || leftIndex >= boundaryTrees.size()) {
            return;
        }
        boundaryTrees.remove(leftIndex);
    }

    private void setBoundary(int leftIndex, BoundaryTree tree) {
        while (boundaryTrees.size() < segments.size() - 1) {
            boundaryTrees.add(null);
        }
        if (leftIndex < boundaryTrees.size()) {
            boundaryTrees.set(leftIndex, tree);
        }
        trimBoundaryList();
    }

    private void trimBoundaryList() {
        while (boundaryTrees.size() > segments.size() - 1) {
            boundaryTrees.remove(boundaryTrees.size() - 1);
        }
    }

    private static List<String> sliceTail(List<String> tokens, int count) {
        if (tokens.isEmpty()) {
            return Collections.emptyList();
        }
        int start = Math.max(0, tokens.size() - count);
        return new ArrayList<>(tokens.subList(start, tokens.size()));
    }

    private static List<String> sliceHead(List<String> tokens, int count) {
        if (tokens.isEmpty()) {
            return Collections.emptyList();
        }
        int end = Math.min(tokens.size(), count);
        return new ArrayList<>(tokens.subList(0, end));
    }

    private static SuffixTreeIndex buildSuffixTree(List<String> tokens) {
        SuffixTreeIndex index = new SuffixTreeIndex(Math.max(1_024, tokens.size() * 4L), 0.0001);
        for (String token : tokens) {
            index.insert(token);
        }
        index.compactForQuerying();
        return index;
    }

    @Override
    public boolean exists(String key) {
        Objects.requireNonNull(key, "key");
        if (key.isEmpty()) {
            return false;
        }
        Pattern pattern = new Pattern(key, key.length());
        ArrayList<Integer> matches = report(pattern);
        return !matches.isEmpty();
    }

    @Override
    public ArrayList<Integer> report(Pattern pat) {
        if (pat == null || pat.nGramArr == null || pat.nGramArr.length == 0) {
            return new ArrayList<>();
        }

        List<Integer> a = reportWithinSegments(pat);
        List<Integer> b = reportCrossingBoundaries(pat);
        List<Integer> c = reportTailRegion(pat);

        LinkedHashSet<Integer> merged = new LinkedHashSet<>();
        merged.addAll(a);
        merged.addAll(b);
        merged.addAll(c);

        ArrayList<Integer> filtered = new ArrayList<>();
        int minAllowed = windowStartIndex;
        int maxAllowed = nextGlobalIndex - pat.nGramArr.length; // ensure match fits inside indexed suffix
        for (Integer pos : merged) {
            if (pos == null) {
                continue;
            }
            if (pos < minAllowed) {
                continue;
            }
            if (pos > maxAllowed) {
                continue;
            }
            filtered.add(pos);
        }
        Collections.sort(filtered);
        return filtered;
    }

    private List<Integer> reportWithinSegments(Pattern pat) {
        ArrayList<Integer> results = new ArrayList<>();
        for (Segment segment : segments) {
            ArrayList<Integer> local = segment.suffixTree.report(pat);
            if (local.isEmpty()) {
                continue;
            }
            int base = segment.globalStart;
            for (Integer idx : local) {
                int globalPos = base + idx;
                if (globalPos >= windowStartIndex) {
                    results.add(globalPos);
                }
            }
        }
        return results;
    }

    private List<Integer> reportCrossingBoundaries(Pattern pat) {
        if (boundaryTrees.isEmpty()) {
            return Collections.emptyList();
        }
        ArrayList<Integer> results = new ArrayList<>();
        int patternLength = pat.nGramArr.length;
        for (BoundaryTree boundary : boundaryTrees) {
            if (boundary == null) {
                continue;
            }
            int leftLen = boundary.leftSeg.tokens.size();
            int rightLen = boundary.rightSeg.tokens.size();
            if (Math.min(leftLen, rightLen) < patternLength) {
                continue;
            }
            ArrayList<Integer> local = boundary.suffixTree.report(pat);
            if (local.isEmpty()) {
                continue;
            }
            for (Integer offset : local) {
                int startLocal = offset;
                int endLocal = offset + patternLength;
                if (startLocal < boundary.splitOffset && endLocal > boundary.splitOffset) {
                    int globalStart = boundary.leftStartGlobal + startLocal;
                    if (globalStart >= windowStartIndex) {
                        results.add(globalStart);
                    }
                }
            }
        }
        return results;
    }

    private List<Integer> reportTailRegion(Pattern pat) {
        int patternLength = pat.nGramArr.length;
        if (segments.isEmpty()) {
            return Collections.emptyList();
        }
        int pivotIndex = -1;
        for (int i = segments.size() - 1; i >= 0; i--) {
            Segment candidate = segments.get(i);
            if (candidate.tokens.size() >= patternLength) {
                pivotIndex = i;
                break;
            }
        }
        if (pivotIndex < 0) {
            pivotIndex = 0;
        }

        Segment pivot = segments.get(pivotIndex);

        ArrayList<String> tailTokens = new ArrayList<>();
        int overlap = Math.max(0, patternLength - 1);
        int startIndexWithinPivot = Math.max(0, pivot.tokens.size() - overlap);
        tailTokens.addAll(pivot.tokens.subList(startIndexWithinPivot, pivot.tokens.size()));
        int baseGlobal = pivot.globalStart + startIndexWithinPivot;
        for (int i = pivotIndex + 1; i < segments.size(); i++) {
            tailTokens.addAll(segments.get(i).tokens);
        }
        if (tailTokens.isEmpty()) {
            return Collections.emptyList();
        }

        SuffixTreeIndex temp = buildSuffixTree(tailTokens);
        ArrayList<Integer> localMatches = temp.report(pat);
        if (localMatches.isEmpty()) {
            return Collections.emptyList();
        }
        ArrayList<Integer> results = new ArrayList<>(localMatches.size());
        for (Integer idx : localMatches) {
            int globalStart = baseGlobal + idx;
            if (globalStart >= windowStartIndex) {
                results.add(globalStart);
            }
        }
        return results;
    }

    @Override
    public void expire() {
        segments.clear();
        boundaryTrees.clear();
        nextGlobalIndex = 0;
        windowStartIndex = 0;
    }

    @Override
    public ArrayList<Long> getAvgTimes(Pattern pat) {
        return null;
    }

    @Override
    public PatternResult getLatestStats() {
        return null;
    }

    @Override
    public int getTokenId(String key) {
        return -2;
    }

    private static final class Segment {
        private final ArrayList<String> tokens;
        private final int globalStart;
        private final int globalEnd;
        private final SuffixTreeIndex suffixTree;

        private Segment(ArrayList<String> tokens, int globalStart, int globalEnd, SuffixTreeIndex suffixTree) {
            this.tokens = tokens;
            this.globalStart = globalStart;
            this.globalEnd = globalEnd;
            this.suffixTree = suffixTree;
        }

        private int length() {
            return tokens.size();
        }
    }

    private static final class BoundaryTree {
        private final Segment leftSeg;
        private final Segment rightSeg;
        private final int leftStartGlobal;
        private final int splitOffset;
        private final SuffixTreeIndex suffixTree;

        private BoundaryTree(Segment leftSeg,
                             Segment rightSeg,
                             int leftStartGlobal,
                             int splitOffset,
                             SuffixTreeIndex suffixTree) {
            this.leftSeg = leftSeg;
            this.rightSeg = rightSeg;
            this.leftStartGlobal = leftStartGlobal;
            this.splitOffset = splitOffset;
            this.suffixTree = suffixTree;
        }
    }
}
