package PMIndex;

import search.Pattern;
import tree.ssws.SuffixTree;
import utilities.AlphabetMapper;
import utilities.PatternResult;

import java.util.*;

/**
 * Delayed SSWSI variant (Section 4 in CPM 2024):
 * - Maintains only the O(log(w/delta)) largest power-of-two segments (suffix trees),
 *   so the uncovered suffix t has size < delta.
 * - Long patterns (|P| > delta/4): answer immediately by querying segment and
 *   boundary trees and by building a temporary suffix tree over R = last (m-1) of s + t.
 * - Short patterns (|P| <= delta/4): optionally buffer queries/updates up to ~delta/2
 *   ops; on flush, build suffix tree over t and answer in batch. For synchronous calls,
 *   we can also answer immediately by doing the same work eagerly.
 *
 * This class is a practical, faithful implementation on top of the segment/boundary
 * scaffolding of StreamingSlidingWindowIndex. It keeps the segments merged in powers of two;
 * tree materialization is skipped for segments shorter than minSegLen (= delta/2).
 */
public final class DelayedStreamingSlidingWindowIndex implements IPMIndexing {

    private final int windowSize;
    private final int delta; // delay parameter
    private final int minSegLen; // == delta/2 (rounded)
    private final AlphabetMapper<String> alphabetMapper;

    private Segment head;
    private Segment tail;
    private int segmentCount = 0;
    private int totalStoredLength = 0;
    private int nextTokenPosition = 0;

    private final ArrayDeque<PendingQuery> queryBuffer = new ArrayDeque<>();
    private int opsSinceLastFlush = 0; // text+pattern characters
    private final boolean bufferedMode; // true => delay short patterns; false => answer immediately

    public DelayedStreamingSlidingWindowIndex(int windowSize, int expectedAlphabetSize, int delta, boolean bufferedMode) {
        if (windowSize <= 0) throw new IllegalArgumentException("windowSize must be positive");
        if (delta <= 0) throw new IllegalArgumentException("delta must be positive");
        this.windowSize = windowSize;
        this.delta = delta;
        this.minSegLen = Math.max(1, delta / 2);
        this.alphabetMapper = new AlphabetMapper<>(Math.max(expectedAlphabetSize, 16));
        this.bufferedMode = bufferedMode;
    }

    public DelayedStreamingSlidingWindowIndex(int windowSize, int delta) {
        this(windowSize, 1024, delta, true);
    }

    @Override
    public void insert(String key) {
        if (key == null) return;
        int token = alphabetMapper.getId(key);
        Segment singleton = createSingletonSegment(token, nextTokenPosition);
        appendSegment(singleton);
        nextTokenPosition++;
        totalStoredLength += singleton.length;
        cascadeMergesFromTail();
        trimToWindow();

        // accounting for delayed mode
        opsSinceLastFlush += 1;
        if (bufferedMode && opsSinceLastFlush >= Math.max(1, delta / 2)) {
            flushBuffers();
        }
    }

    @Override
    public boolean exists(String key) { return false; }

    @Override
    public ArrayList<Integer> report(Pattern pattern) {
        int[] pat = toPatternTokens(pattern);
        ArrayList<Integer> out = new ArrayList<>();
        if (pat.length == 0) return out;

        // Short vs long
        if (pat.length <= Math.max(1, delta / 4)) {
            if (bufferedMode) {
                // Buffer query and flush opportunistically
                queryBuffer.addLast(new PendingQuery(pat.clone(), currentRightBoundary()));
                opsSinceLastFlush += pat.length;
                if (opsSinceLastFlush >= Math.max(1, delta / 2)) flushBuffers();
                // Return current best-effort answer immediately (optional); here we return the
                // eager computation to keep synchronous semantics useful.
                out.addAll(reportShortImmediate(pat));
                return out;
            } else {
                out.addAll(reportShortImmediate(pat));
                return out;
            }
        } else {
            // Long pattern: immediate answer
            out.addAll(reportLongImmediate(pat));
            return out;
        }
    }

    private List<Integer> reportLongImmediate(int[] pat) {
        LinkedHashSet<Integer> hits = new LinkedHashSet<>();
        // 1) occurrences fully in s (segments >= minSegLen) and across a single boundary
        findInSegmentsAndBoundaries(pat, hits);
        // 2) R = last (m-1) of s + t (uncovered suffix)
        int m = pat.length;
        int[] leftTail = tailFromSegments(Math.max(0, m - 1));
        int[] t = uncoveredSuffixTokens();
        int[] R = concat(leftTail, t);
        if (R.length > 0) {
            SuffixTree rTree = SuffixTree.build(R);
            List<Integer> local = rTree.findOccurrences(pat);
            int globalStart = currentRightBoundary() - R.length;
            for (int pos : local) {
                hits.add(globalStart + pos);
            }
        }
        return new ArrayList<>(hits);
    }

    private List<Integer> reportShortImmediate(int[] pat) {
        LinkedHashSet<Integer> hits = new LinkedHashSet<>();
        // 1) in s and across boundaries
        findInSegmentsAndBoundaries(pat, hits);
        // 2) across boundary (s,t) via KMP on window of size 2*(m-1)+1 centered at boundary
        int m = pat.length;
        int[] left = tailFromSegments(Math.max(0, m - 1));
        int[] right = prefix(uncoveredSuffixTokens(), Math.max(0, m - 1));
        int[] around = concat(left, right);
        int base = currentRightBoundary() - right.length - Math.max(0, m - 1);
        for (int pos : kmpAll(pat, around)) {
            hits.add(base + pos);
        }
        // 3) fully in t
        int[] t = uncoveredSuffixTokens();
        if (t.length > 0) {
            SuffixTree tree = SuffixTree.build(t);
            for (int pos : tree.findOccurrences(pat)) {
                hits.add(currentRightBoundary() - t.length + pos);
            }
        }
        return new ArrayList<>(hits);
    }

    private void flushBuffers() {
        if (queryBuffer.isEmpty()) { opsSinceLastFlush = 0; return; }
        // Build suffix tree over t once
        int[] t = uncoveredSuffixTokens();
        SuffixTree tTree = (t.length > 0) ? SuffixTree.build(t) : null;
        int rightBoundaryAtFlush = currentRightBoundary();
        while (!queryBuffer.isEmpty()) {
            PendingQuery q = queryBuffer.removeFirst();
            LinkedHashSet<Integer> hits = new LinkedHashSet<>();
            int[] pat = q.pattern;
            // in s and boundaries
            findInSegmentsAndBoundaries(pat, hits);
            // across boundary with KMP
            int m = pat.length;
            int[] left = tailFromSegments(Math.max(0, m - 1));
            int[] right = prefix(uncoveredSuffixTokens(), Math.max(0, m - 1));
            int[] around = concat(left, right);
            int base = rightBoundaryAtFlush - right.length - Math.max(0, m - 1);
            for (int pos : kmpAll(pat, around)) hits.add(base + pos);
            // fully in t at the time of flush
            if (tTree != null) {
                for (int pos : tTree.findOccurrences(pat)) {
                    hits.add(rightBoundaryAtFlush - t.length + pos);
                }
            }
            // Note: results are not returned here (buffered). This method is for maintenance;
            // the caller's synchronous report already returned best-effort answers.
        }
        opsSinceLastFlush = 0;
    }

    private void findInSegmentsAndBoundaries(int[] pat, Set<Integer> out) {
        // segments
        for (Segment s = head; s != null; s = s.next) {
            if (s.suffixTree != null) {
                for (int pos : s.suffixTree.findOccurrences(pat)) {
                    out.add(s.startPosition + pos);
                }
            }
        }
        // boundaries
        for (Segment s = head; s != null; s = s.next) {
            if (s.rightBoundary != null && s.rightBoundary.tree != null) {
                Boundary b = s.rightBoundary;
                for (int pos : b.tree.findOccurrences(pat)) {
                    out.add(b.globalStart + pos);
                }
            }
        }
    }

    // Removed explicit sorting to match Streaming index behavior and reduce overhead.

    private int currentRightBoundary() {
        return nextTokenPosition; // last appended position + 1
    }

    // ===== Segment maintenance (power-of-two merges), only materialize trees for len >= minSegLen =====

    private Segment createSingletonSegment(int token, int position) {
        int[] tokens = new int[]{token};
        SuffixTree tree = (tokens.length >= minSegLen) ? SuffixTree.build(tokens) : null;
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
            current = merged;
        }
    }

    private Segment mergeSegments(Segment left, Segment right) {
        int newLength = left.length + right.length;
        int[] tokens = new int[newLength];
        System.arraycopy(left.tokens, 0, tokens, 0, left.length);
        System.arraycopy(right.tokens, 0, tokens, left.length, right.length);
        SuffixTree tree = (newLength >= minSegLen) ? SuffixTree.build(tokens) : null;
        return new Segment(left.startPosition, left.level + 1, tokens, tree);
    }

    private void replaceWithMerged(Segment left, Segment right, Segment merged) {
        Segment prev = left.prev;
        Segment next = right.next;

        clearBoundary(left.leftBoundary);
        clearBoundary(left.rightBoundary);
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
        segmentCount--;
    }

    private void buildBoundary(Segment left, Segment right) {
        // Only materialize boundary trees when both sides are covered by segment trees
        if (left == null || right == null) return;
        left.rightBoundary = null;
        right.leftBoundary = null;
        if (left.suffixTree == null || right.suffixTree == null) return;
        int span = right.length;
        int leftSpan = Math.min(left.length, span);
        int[] tokens = new int[leftSpan + right.length];
        System.arraycopy(left.tokens, left.length - leftSpan, tokens, 0, leftSpan);
        System.arraycopy(right.tokens, 0, tokens, leftSpan, right.length);
        SuffixTree tree = SuffixTree.build(tokens);
        int boundaryStart = left.startPosition + left.length - leftSpan;
        int boundaryPosition = left.startPosition + left.length;
        Boundary boundary = new Boundary(left, right, tokens, boundaryStart, boundaryPosition, tree);
        left.rightBoundary = boundary;
        right.leftBoundary = boundary;
    }

    private void clearBoundary(Boundary b) {
        if (b == null) return;
        if (b.left != null && b.left.rightBoundary == b) b.left.rightBoundary = null;
        if (b.right != null && b.right.leftBoundary == b) b.right.leftBoundary = null;
    }

    private void trimToWindow() {
        while (head != null && totalStoredLength - head.length >= windowSize) {
            totalStoredLength -= head.length;
            removeHead();
        }
    }

    private void removeHead() {
        if (head == null) return;
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

    private int[] toPatternTokens(Pattern pattern) {
        String[] grams = pattern.nGramArr;
        if (grams == null) return new int[0];
        int[] tokens = new int[grams.length];
        for (int i = 0; i < grams.length; i++) tokens[i] = alphabetMapper.getId(grams[i]);
        return tokens;
    }

    private int[] uncoveredSuffixTokens() {
        // Concatenate tokens from the tail made of segments < minSegLen
        ArrayDeque<int[]> parts = new ArrayDeque<>();
        int total = 0;
        for (Segment s = tail; s != null; s = s.prev) {
            if (s.length >= minSegLen) break;
            parts.addFirst(s.tokens);
            total += s.length;
            if (total >= delta) break; // cap at delta for safety
        }
        int[] out = new int[total];
        int w = 0;
        for (int[] arr : parts) { System.arraycopy(arr, 0, out, w, arr.length); w += arr.length; }
        return out;
    }

    private int[] tailFromSegments(int need) {
        if (need <= 0) return new int[0];
        int have = 0;
        ArrayDeque<int[]> parts = new ArrayDeque<>();
        for (Segment s = tail; s != null && have < need; s = s.prev) {
            if (s.length < minSegLen || s.suffixTree == null) continue;
            int take = Math.min(s.length, need - have);
            if (take == s.length) {
                parts.addFirst(s.tokens);
                have += take;
            } else {
                int[] slice = new int[take];
                System.arraycopy(s.tokens, s.length - take, slice, 0, take);
                parts.addFirst(slice);
                have += take;
            }
        }
        int[] out = new int[have];
        int w = 0;
        for (int[] p : parts) { System.arraycopy(p, 0, out, w, p.length); w += p.length; }
        return out;
    }

    private static int[] concat(int[] a, int[] b) {
        if (a.length == 0) return b.clone();
        if (b.length == 0) return a.clone();
        int[] out = new int[a.length + b.length];
        System.arraycopy(a, 0, out, 0, a.length);
        System.arraycopy(b, 0, out, a.length, b.length);
        return out;
    }

    private static int[] prefix(int[] a, int len) {
        int k = Math.min(len, a.length);
        int[] out = new int[k];
        System.arraycopy(a, 0, out, 0, k);
        return out;
    }

    private static List<Integer> kmpAll(int[] pat, int[] text) {
        ArrayList<Integer> out = new ArrayList<>();
        if (pat.length == 0 || text.length == 0) return out;
        int[] pi = new int[pat.length];
        for (int i = 1, j = 0; i < pat.length; i++) {
            while (j > 0 && pat[i] != pat[j]) j = pi[j - 1];
            if (pat[i] == pat[j]) j++;
            pi[i] = j;
        }
        for (int i = 0, j = 0; i < text.length; i++) {
            while (j > 0 && text[i] != pat[j]) j = pi[j - 1];
            if (text[i] == pat[j]) {
                j++;
                if (j == pat.length) {
                    out.add(i - pat.length + 1);
                    j = pi[j - 1];
                }
            }
        }
        return out;
    }

    @Override
    public void expire() { /* not used in this variant */ }

    @Override
    public ArrayList<Long> getAvgTimes(Pattern pat) { return new ArrayList<>(); }

    @Override
    public PatternResult getLatestStats() { return new PatternResult(0.0, 0, 0, null, 0, 0.0, 0, 0); }

    @Override
    public int getTokenId(String key) { return alphabetMapper.getId(key); }

    // ===== data holders =====
    private static final class Segment {
        private final int startPosition;
        private final int level;
        private final int length;
        private final int[] tokens;
        private final SuffixTree suffixTree; // null if not materialized (< minSegLen)
        private Segment prev;
        private Segment next;
        private Boundary leftBoundary;
        private Boundary rightBoundary;

        private Segment(int startPosition, int level, int[] tokens, SuffixTree tree) {
            this.startPosition = startPosition;
            this.level = level;
            this.tokens = tokens;
            this.length = tokens.length;
            this.suffixTree = tree;
        }
    }

    private static final class Boundary {
        private final Segment left;
        private final Segment right;
        private final int[] tokens;
        private final int globalStart;
        private final int boundaryPosition;
        private final SuffixTree tree;

        private Boundary(Segment left, Segment right, int[] tokens, int globalStart,
                         int boundaryPosition, SuffixTree tree) {
            this.left = left;
            this.right = right;
            this.tokens = tokens;
            this.globalStart = globalStart;
            this.boundaryPosition = boundaryPosition;
            this.tree = tree;
        }
    }

    private static final class PendingQuery {
        final int[] pattern;
        final int rightBoundaryAtArrival;
        PendingQuery(int[] pattern, int rightBoundaryAtArrival) {
            this.pattern = pattern;
            this.rightBoundaryAtArrival = rightBoundaryAtArrival;
        }
    }
}
