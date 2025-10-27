package tree.ssws;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;

/**
 * SuffixTreeDC3
 *
 * This class builds a compressed suffix tree for an integer token sequence.
 *
 * Public contract
 *
 *   SuffixTreeDC3.build(int[] tokens)
 *     Build a suffix tree for the given array of integer tokens. We append a unique
 *     sentinel token that is strictly smaller than any real token so that all suffixes
 *     are distinct. The caller must guarantee that no real token collides with that
 *     sentinel after shifting, which we enforce by construction.
 *
 *   findOccurrences(int[] pattern)
 *     Return all start offsets (0-based) where the given pattern occurs in the
 *     original text (without counting the sentinel). Offsets are filtered so that
 *     matches that extend past the end of the real text are discarded.
 *
 * Internals
 *
 *   The build pipeline is:
 *
 *     1. Append a sentinel and simultaneously generate a dense non-negative
 *        "normalized" alphabet using one radix sort pass and rank assignment.
 *        This removes the need for normalizeWithSentinel(…) plus extra copies.
 *
 *     2. Build the suffix array in linear time using DC3 (also called Skew).
 *
 *     3. Build the Longest Common Prefix (LCP) array in linear time with Kasai.
 *
 *     4. Reconstruct the explicit compressed suffix tree in linear time from
 *        (suffix array, LCP array), using an array-based stack rather than
 *        allocating StackEntry objects. Node children are stored in a primitive
 *        Int2ObjectOpenHashMap<Edge> to avoid boxing overhead from HashMap<Integer,Edge>.
 *
 *  Complexity
 *
 *    Build time is overall O(n) where n is the text length including sentinel,
 *    under the standard integer alphabet assumptions we already enforce.
 *
 *    Memory is O(n) nodes and edges plus the text array.
 *
 *    Query time for findOccurrences is O(m + occ) where m is pattern length
 *    and occ is the number of reported matches. We descend once over the
 *    pattern and then collect leaf suffix indices.
 */
public final class SuffixTreeDC3 {

    private final Node root;
    private final int[] text; // includes the sentinel at the end
    private final int originalLength; // length before appending sentinel

    private SuffixTreeDC3(Node root, int[] text, int originalLength) {
        this.root = root;
        this.text = text;
        this.originalLength = originalLength;
    }

    /**
     * Build a suffix tree for the provided integer token array. The builder appends a
     * unique sentinel value that is strictly smaller than any existing token so that
     * each suffix becomes unique. We then:
     *
     *   - Produce a normalized "rank" array of the same length with values in
     *     {0, 1, ..., distinct-1} using a single radix sort on 64-bit keys
     *     (RadixSorter.sortParallel). This both enforces a dense alphabet and
     *     guarantees that the sentinel is globally smallest.
     *
     *   - Run DC3 (Skew) on that normalized array to get a suffix array in O(n).
     *
     *   - Run Kasai to get the Longest Common Prefix (LCP) array in O(n).
     *
     *   - Reconstruct the explicit suffix tree using a stack-based linear pass.
     *
     * Returned edges refer to the *original* token array with the sentinel,
     * not the normalized ranks. That means findOccurrences compares real tokens.
     */
    public static SuffixTreeDC3 build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }

        // ------------------------------------------------------------
        // Step 0. Append sentinel and in the same pass prepare data
        // for normalization and suffix array construction.
        // We also avoid a separate normalizeWithSentinel(…).
        //
        // We choose sentinel = (minToken - 1) unless minToken is Integer.MIN_VALUE,
        // in which case we just reuse Integer.MIN_VALUE. The sentinel is strictly
        // smaller than every normal token. This guarantees that the last suffix,
        // which starts at that sentinel, is lexicographically smallest.
        // ------------------------------------------------------------
        final int n0 = alphabetMappedText.length;
        final int[] terminated = new int[n0 + 1];

        int minVal = (n0 == 0) ? 0 : alphabetMappedText[0];
        for (int v : alphabetMappedText) {
            if (v < minVal) {
                minVal = v;
            }
        }
        final int sentinel = (n0 == 0)
                ? Integer.MIN_VALUE
                : (minVal == Integer.MIN_VALUE ? Integer.MIN_VALUE : (minVal - 1));

        System.arraycopy(alphabetMappedText, 0, terminated, 0, n0);
        terminated[n0] = sentinel;

        final int n = n0 + 1;

        // We now build a dense non-negative normalized alphabet for DC3.
        //
        // We do this in-place here instead of using a separate normalizeWithSentinel.
        // We assign each position i a 64-bit key of the form:
        //
        //   key = ((terminated[i] - sentinel) << 1) | flag
        //
        // where flag = 1 for all real symbols and flag = 0 for the sentinel
        // at the end. Because sentinel == minimum value, terminated[i] - sentinel
        // is always >= 0. Shifting by one bit leaves room to force the sentinel's
        // "flag" to 0 so it becomes globally smallest.
        //
        // We then radix-sort these keys, stable, and assign consecutive integer
        // ranks 0,1,2,... to equal keys, which is exactly what DC3 expects:
        //
        //   - integer alphabet
        //   - non-negative
        //   - alphabet size O(n)
        //
        long[] keys = new long[n];
        int[] order = new int[n];
        final long base = (long) sentinel; // subtract sentinel so every value is >= 0

        for (int i = 0; i < n; i++) {
            long adjusted = ((long) terminated[i]) - base; // non-negative
            long k = (adjusted << 1) | 1L; // default flag = 1
            if (i == n - 1) {
                // make the very last position (the sentinel suffix) lexicographically smallest
                k = (adjusted << 1);
            }
            keys[i] = k;
            order[i] = i;
        }

        // radix sort keys[ ] and reorder order[ ] in parallel
        RadixSorter.sortParallel(keys, order);

        // assign dense ranks: equal keys get the same rank
        int[] normalized = new int[n];
        int rank = 0;
        long lastKey = keys[0];
        normalized[order[0]] = 0;
        for (int i = 1; i < n; i++) {
            long k = keys[i];
            if (k != lastKey) {
                rank++;
                lastKey = k;
            }
            normalized[order[i]] = rank;
        }

        // ------------------------------------------------------------
        // Step 1. Build suffix array on normalized ranks via DC3 / Skew.
        // This is O(n) time. We keep normalized non-negative so DC3 can do
        // radix/counting sorts without comparison sorting.
        // ------------------------------------------------------------
        int[] sa = SuffixArrayBuilder.buildSuffixArray(normalized);

        // ------------------------------------------------------------
        // Step 2. Build LCP array with Kasai's algorithm in O(n) time.
        // ------------------------------------------------------------
        int[] lcp = (n > 1)
                ? SuffixArrayBuilder.buildLcpArray(normalized, sa)
                : new int[0];

        // ------------------------------------------------------------
        // Step 3. Reconstruct the explicit compressed suffix tree in O(n)
        // using a stack-based algorithm. We switched to an array-backed
        // stack and primitive maps to drastically reduce allocation.
        // ------------------------------------------------------------
        Node root = LinearBuilder.buildFromSuffixArray(terminated, sa, lcp);

        return new SuffixTreeDC3(root, terminated, n0);
    }

    /**
     * Return the root node of the suffix tree.
     */
    public Node getRoot() {
        return root;
    }

    /**
     * Return a defensive copy of the underlying text, including the sentinel.
     */
    public int[] getText() {
        return Arrays.copyOf(text, text.length);
    }

    /**
     * Return the original text length, without counting the sentinel that was appended.
     */
    public int getOriginalLength() {
        return originalLength;
    }

    /**
     * Find all starting offsets where pattern occurs in the original text.
     * The pattern is given as an array of integer tokens in the same alphabet
     * as the text (without the sentinel).
     *
     * We walk down the suffix tree by following edges and comparing along
     * each edge's substring. When we consume the entire pattern we collect
     * all leaves beneath that point. This takes time proportional to
     * (pattern length + number of matches).
     */
    public List<Integer> findOccurrences(int[] pattern) {
        if (pattern == null || pattern.length == 0) {
            return Collections.emptyList();
        }

        Node current = root;
        int patternIndex = 0;

        while (patternIndex < pattern.length) {

            int symbol = pattern[patternIndex];
            Edge edge = current.edges.get(symbol);
            if (edge == null) {
                return Collections.emptyList();
            }

            int edgeStart = edge.start;
            int edgeEndExclusive = edge.end;
            int edgeLen = edgeEndExclusive - edgeStart;

            int consumed = 0;
            while (consumed < edgeLen && patternIndex < pattern.length) {
                if (text[edgeStart + consumed] != pattern[patternIndex]) {
                    return Collections.emptyList();
                }
                consumed++;
                patternIndex++;
            }

            // If we finished matching the pattern in the middle (or end) of this edge,
            // then every suffix hanging below edge.child is a full match.
            if (patternIndex == pattern.length) {
                return collectOccurrences(edge.child, pattern.length);
            }

            // Otherwise we consumed the full edge label but still have pattern left,
            // so descend to the child node and continue.
            current = edge.child;
        }

        // If we exited naturally, then pattern ended exactly at a node.
        // Collect leaves under that node.
        return collectOccurrences(current, pattern.length);
    }

    private List<Integer> collectOccurrences(Node node, int patternLength) {
        if (node == null) {
            return Collections.emptyList();
        }
        List<Integer> matches = new ArrayList<>();
        collectDFS(node, patternLength, matches);
        return matches;
    }

    private void collectDFS(Node node, int patternLength, List<Integer> out) {
        if (node.isLeaf()) {
            int idx = node.suffixIndex;
            // Filter out matches that would run past the original text length
            if (idx >= 0 && idx + patternLength <= originalLength) {
                out.add(idx);
            }
            return;
        }
        for (Edge e : node.edges.values()) {
            collectDFS(e.child, patternLength, out);
        }
    }

    /**
     * Node in the suffix tree.
     *
     * We store outgoing edges in an Int2ObjectOpenHashMap<Edge> (from the FastUtil
     * library) instead of a java.util.HashMap<Integer, Edge>. This avoids boxing
     * integers into java.lang.Integer objects and avoids HashMap's node overhead.
     *
     * suffixIndex is only meaningful for leaves. Internal nodes have suffixIndex = -1.
     */
    public static final class Node {
        private final Int2ObjectOpenHashMap<Edge> edges;
        private int suffixIndex = -1;

        private Node() {
            this.edges = new Int2ObjectOpenHashMap<>();
        }

        private static Node leaf(int suffixIndex) {
            Node node = new Node();
            node.suffixIndex = suffixIndex;
            return node;
        }

        /**
         * Expose the edges as an unmodifiable Map<Integer, Edge> for callers that expect
         * java.util.Map. We pay boxing only in this accessor, not in the hot loops.
         */
        public Map<Integer, Edge> getEdges() {
            return Collections.unmodifiableMap(edges);
        }

        public boolean isLeaf() {
            return suffixIndex >= 0;
        }

        public int getSuffixIndex() {
            return suffixIndex;
        }
    }

    /**
     * Edge in the suffix tree.
     *
     * Each edge is labeled implicitly by text[start .. end),
     * where end is exclusive. We never copy substrings.
     */
    public static final class Edge {
        private final int start;
        private int end;
        private final Node child;

        private Edge(int start, int end, Node child) {
            this.start = start;
            this.end = end;
            this.child = child;
        }

        public int getStart() {
            return start;
        }

        public int getEnd() {
            return end;
        }

        public Node getChild() {
            return child;
        }
    }

    /**
     * LinearBuilder
     *
     * Convert (suffix array, Longest Common Prefix array) into an explicit compressed
     * suffix tree in O(n) time using a stack. We removed ArrayDeque<StackEntry> and
     * replaced it with parallel primitive stacks:
     *
     *   Node[] nodeStack
     *   int[] depthStack
     *   Edge[] edgeStack
     *   int top
     *
     * which avoids allocating StackEntry objects per suffix and is much friendlier
     * for the Java Virtual Machine.
     *
     * depthStack[k] is the total string depth (number of symbols from root) of nodeStack[k].
     * edgeStack[k] is the incoming edge that led to nodeStack[k] from its parent, or null for root.
     */
    private static final class LinearBuilder {

        private LinearBuilder() {
        }

        static Node buildFromSuffixArray(int[] text, int[] sa, int[] lcp) {
            final int n = sa.length;
            final Node root = new Node();

            // Maximum possible depth in symbols is at most text.length.
            // We allocate arrays of size text.length + 1 to be safe.
            final int maxDepth = text.length + 1;
            final Node[] nodeStack = new Node[maxDepth];
            final int[] depthStack = new int[maxDepth];
            final Edge[] edgeStack = new Edge[maxDepth];

            int top = 0;
            nodeStack[0] = root;
            depthStack[0] = 0;
            edgeStack[0] = null;

            for (int i = 0; i < n; i++) {

                int suffixStart = sa[i];
                int lcpValue = (i == 0) ? 0 : lcp[i - 1];

                // We will need to know the incoming edge that got popped last,
                // in case we must split it to create a new internal node.
                Edge lastIncomingEdge = null;

                // Pop stack until we are at depth <= lcpValue.
                while (depthStack[top] > lcpValue) {
                    lastIncomingEdge = edgeStack[top];
                    top--;
                }

                // If the current top depth is still less than lcpValue,
                // then the new Longest Common Prefix boundary falls in
                // the middle of the edge that we just popped. We must
                // split that edge and create an internal node.
                if (depthStack[top] < lcpValue) {
                    if (lastIncomingEdge == null) {
                        throw new IllegalStateException("Inconsistent suffix array / LCP combination");
                    }

                    Node parentNode = nodeStack[top];
                    Node oldChild = lastIncomingEdge.child;

                    // Split point inside lastIncomingEdge's label.
                    int splitOffset = lcpValue - depthStack[top];
                    int splitPoint = lastIncomingEdge.start + splitOffset;

                    Node internal = new Node();

                    // Tail edge from the split point down to the original child.
                    Edge tail = new Edge(splitPoint, lastIncomingEdge.end, oldChild);
                    internal.edges.put(text[tail.start], tail);

                    // Shorten the incoming edge so that it now ends at the split
                    // point and points to our new internal node.
                    lastIncomingEdge.end = splitPoint;
                    // We must replace the child in the parentNode's map as well.
                    // Because we do not track the parent edge key separately,
                    // we can safely re-put using the original start symbol.
                    parentNode.edges.put(text[lastIncomingEdge.start],
                            new Edge(lastIncomingEdge.start, lastIncomingEdge.end, internal));

                    // Push the internal node on the stack.
                    top++;
                    nodeStack[top] = internal;
                    depthStack[top] = lcpValue;
                    edgeStack[top] = new Edge(lastIncomingEdge.start, lastIncomingEdge.end, internal);
                }

                // Now we are at the correct parent node for this suffix.
                Node parent = nodeStack[top];

                // Create a fresh leaf node that represents suffix starting at suffixStart.
                Node leaf = Node.leaf(suffixStart);
                int edgeStart = suffixStart + depthStack[top];
                Edge leafEdge = new Edge(edgeStart, text.length, leaf);

                parent.edges.put(text[edgeStart], leafEdge);

                // Push leaf on stack so the next suffix can reuse its path depth.
                top++;
                nodeStack[top] = leaf;
                depthStack[top] = text.length - suffixStart; // full suffix depth
                edgeStack[top] = leafEdge;
            }

            return root;
        }
    }
}
