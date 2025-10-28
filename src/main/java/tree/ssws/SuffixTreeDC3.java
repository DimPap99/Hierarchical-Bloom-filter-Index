package tree.ssws;

import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * SuffixTreeDC3
 *
 * Compressed suffix tree built from an integer token sequence.
 * Build pipeline:
 *   1. Append sentinel smaller than any real token.
 *   2. Normalize tokens to a dense non-negative alphabet using radix/ranking.
 *   3. Build suffix array in O(n) time using DC3 (Skew).
 *   4. Build Longest Common Prefix (LCP) with Kasai in O(n).
 *   5. Build explicit compressed suffix tree in O(n) from suffix array + LCP.
 *
 * Query:
 *   findOccurrences(int[] pattern) returns all start offsets in the original text
 *   (no sentinel) where 'pattern' occurs.
 */
public final class SuffixTreeDC3 {

    private final Node root;
    private final int[] text;           // includes sentinel at end
    private final int originalLength;   // length before sentinel

    private SuffixTreeDC3(Node root, int[] text, int originalLength) {
        this.root = root;
        this.text = text;
        this.originalLength = originalLength;
    }

    /**
     * Build a suffix tree from an integer token array. We:
     *  - append a unique sentinel strictly smaller than any input token,
     *  - create a dense non-negative alphabet via radix ranking,
     *  - run DC3 / Skew to get the suffix array,
     *  - run Kasai to get the LCP,
     *  - build a compressed suffix tree in linear time.
     */
    public static SuffixTreeDC3 build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }

        // Step 0. Append sentinel that is strictly smaller than all real symbols.
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

        // Step 0.5. Build dense non-negative ranks so DC3 can work efficiently.
        //
        // We encode each position i (including the sentinel) into a sortable 64-bit key.
        // keys[i] is monotonic with respect to (token value, sentinel-flag).
        // Then we radix-sort those keys and assign consecutive integer ranks.
        long[] keys = new long[n];
        int[] order = new int[n];

        final long base = (long) sentinel;
        for (int i = 0; i < n; i++) {
            long adjusted = ((long) terminated[i]) - base; // >= 0
            long k = (adjusted << 1) | 1L; // default: "real" symbol
            if (i == n - 1) {
                // Force the last position (the sentinel suffix) to be globally smallest.
                k = (adjusted << 1);
            }
            keys[i] = k;
            order[i] = i;
        }

        // Radix sort the (key,index) pairs. This is now optimized to skip unused high bytes.
        RadixSorter.sortParallel(keys, order);

        // Assign dense ranks: equal keys get same rank, monotone increasing otherwise.
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

        // Step 1. Build suffix array using DC3 / Skew on normalized array.
        int[] sa = SuffixArrayBuilder.buildSuffixArray(normalized);

        // Step 2. Build LCP via Kasai.
        int[] lcp = (n > 1)
                ? SuffixArrayBuilder.buildLcpArray(normalized, sa)
                : new int[0];

        // Step 3. Build explicit compressed suffix tree.
        Node root = LinearBuilder.buildFromSuffixArray(terminated, sa, lcp);

        return new SuffixTreeDC3(root, terminated, n0);
    }

    /**
     * Return the root node.
     */
    public Node getRoot() {
        return root;
    }

    /**
     * Return a copy of the underlying text (with sentinel).
     */
    public int[] getText() {
        return Arrays.copyOf(text, text.length);
    }

    /**
     * Return the logical text length, excluding the sentinel.
     */
    public int getOriginalLength() {
        return originalLength;
    }

    /**
     * Find all starting offsets where the pattern occurs in the original text.
     * Pattern is an array of tokens from the same alphabet (without the sentinel).
     * We descend the tree. If we finish matching the pattern in the middle of an edge,
     * we gather all leaf suffix indices below that point. We filter out matches that
     * would run past the original-length boundary.
     */
    public List<Integer> findOccurrences(int[] pattern) {
        if (pattern == null || pattern.length == 0) {
            return Collections.emptyList();
        }

        Node current = root;
        int patternIndex = 0;

        while (patternIndex < pattern.length) {

            int symbol = pattern[patternIndex];
            Edge edge = current.getEdge(symbol);
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

            // Pattern ended in the middle or at the end of this edge.
            if (patternIndex == pattern.length) {
                return collectOccurrences(edge.child, pattern.length);
            }

            // Otherwise we fully consumed this edge and still have pattern left.
            current = edge.child;
        }

        // If we exit the loop naturally, pattern ended exactly at a node.
        return collectOccurrences(current, pattern.length);
    }

    /**
     * Iteratively gather all leaves under 'node'. Each leaf's suffixIndex is a match start.
     * We filter out matches that would extend beyond originalLength.
     *
     * We do this iteratively rather than recursively to avoid deep recursion overhead
     * and Java Virtual Machine stack pressure on large repetitive inputs.
     */
    private List<Integer> collectOccurrences(Node node, int patternLength) {
        if (node == null) {
            return Collections.emptyList();
        }

        List<Integer> matches = new ArrayList<>();
        ArrayList<Node> stack = new ArrayList<>();
        stack.add(node);

        while (!stack.isEmpty()) {
            Node cur = stack.remove(stack.size() - 1);
            if (cur.isLeaf()) {
                int idx = cur.getSuffixIndex();
                if (idx >= 0 && idx + patternLength <= originalLength) {
                    matches.add(idx);
                }
                continue;
            }
            for (Edge e : cur.edgesIterable()) {
                stack.add(e.child);
            }
        }

        return matches;
    }

    /**
     * A node in the suffix tree.
     *
     * We use an adaptive child storage strategy for performance.
     *
     *   - For up to two children we store (key, edge) pairs directly in fields
     *     without allocating any map.
     *
     *   - If a node ever needs more than two distinct children, we "promote"
     *     it to an Int2ObjectOpenHashMap<Edge> from FastUtil. We then move
     *     those inline children into the map once, and from that point on we
     *     just use the map.
     *
     * This saves a huge amount of allocation and rehash churn compared to
     * creating a brand new map for every node up front.
     */
    public static final class Node {

        // Inline child slot A
        private int keyA = Integer.MIN_VALUE;
        private Edge edgeA = null;

        // Inline child slot B
        private int keyB = Integer.MIN_VALUE;
        private Edge edgeB = null;

        // Promoted map for >2 children
        private Int2ObjectOpenHashMap<Edge> edgeMap = null;

        // Leaf payload
        private int suffixIndex = -1; // only meaningful for leaves

        private Node() {
            // no allocation up front
        }

        private static Node leaf(int suffixIndex) {
            Node node = new Node();
            node.suffixIndex = suffixIndex;
            return node;
        }

        public boolean isLeaf() {
            return suffixIndex >= 0;
        }

        public int getSuffixIndex() {
            return suffixIndex;
        }

        /**
         * Return the edge starting with 'symbol', or null if none.
         */
        public Edge getEdge(int symbol) {
            if (edgeMap != null) {
                return edgeMap.get(symbol);
            }
            if (edgeA != null && keyA == symbol) {
                return edgeA;
            }
            if (edgeB != null && keyB == symbol) {
                return edgeB;
            }
            return null;
        }

        /**
         * Insert or replace an outgoing edge under the given first symbol.
         */
        public void putEdge(int symbol, Edge e) {
            if (edgeMap != null) {
                edgeMap.put(symbol, e);
                return;
            }

            // Try slot A
            if (edgeA == null || keyA == symbol) {
                keyA = symbol;
                edgeA = e;
                return;
            }

            // Try slot B
            if (edgeB == null || keyB == symbol) {
                keyB = symbol;
                edgeB = e;
                return;
            }

            // Need to promote to a map
            edgeMap = new Int2ObjectOpenHashMap<>(4);
            edgeMap.put(keyA, edgeA);
            edgeMap.put(keyB, edgeB);
            edgeMap.put(symbol, e);

            // Optional: clear inline slots to save a tiny bit of memory.
            keyA = Integer.MIN_VALUE;
            keyB = Integer.MIN_VALUE;
            edgeA = null;
            edgeB = null;
        }

        /**
         * An iterable view of all outgoing edges. Callers do not need to know if we are
         * using inline slots or a map.
         */
        public Iterable<Edge> edgesIterable() {
            if (edgeMap != null) {
                return edgeMap.values();
            }

            if (edgeA != null && edgeB != null) {
                return Arrays.asList(edgeA, edgeB);
            } else if (edgeA != null) {
                return Collections.singletonList(edgeA);
            } else if (edgeB != null) {
                return Collections.singletonList(edgeB);
            } else {
                return Collections.emptyList();
            }
        }
    }

    /**
     * A compressed edge in the suffix tree. The label is text[start .. end),
     * end exclusive. We do not copy substrings.
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
     * Convert a suffix array + Longest Common Prefix array into an explicit compressed suffix tree
     * in O(n). We keep an explicit stack of:
     *
     *   nodeStack[k]  = the node at stack depth k
     *   depthStack[k] = the string depth (number of symbols from root) for that node
     *   edgeStack[k]  = the incoming edge that leads to nodeStack[k] from its parent
     *
     * No user code needs to call this directly.
     */
    private static final class LinearBuilder {

        private LinearBuilder() {
        }

        static Node buildFromSuffixArray(int[] text, int[] sa, int[] lcp) {
            final int n = sa.length;
            final Node root = new Node();

            // Worst-case depth in symbols is <= text.length.
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

                // lastIncomingEdge remembers the edge we popped most recently
                // so that if we need to split it to create an internal node,
                // we know which edge to rewrite.
                Edge lastIncomingEdge = null;

                // Pop until the top of stack has depth <= current LCP.
                while (depthStack[top] > lcpValue) {
                    lastIncomingEdge = edgeStack[top];
                    top--;
                }

                // If the current top depth is still < lcpValue,
                // we need to split lastIncomingEdge to create an internal node.
                if (depthStack[top] < lcpValue) {
                    if (lastIncomingEdge == null) {
                        throw new IllegalStateException("Inconsistent suffix array / LCP combination");
                    }

                    Node parentNode = nodeStack[top];
                    Node oldChild = lastIncomingEdge.child;

                    // Split inside lastIncomingEdge's label.
                    int splitOffset = lcpValue - depthStack[top];
                    int splitPoint = lastIncomingEdge.start + splitOffset;

                    // New internal node.
                    Node internal = new Node();

                    // Tail edge from splitPoint down to the old child.
                    Edge tail = new Edge(splitPoint, lastIncomingEdge.end, oldChild);
                    internal.putEdge(text[tail.start], tail);

                    // Shorten the original incoming edge to end at the split.
                    lastIncomingEdge.end = splitPoint;

                    // Replace in parentNode: parentNode now points to 'internal' via that shortened edge.
                    parentNode.putEdge(text[lastIncomingEdge.start],
                            new Edge(lastIncomingEdge.start, lastIncomingEdge.end, internal));

                    // Push the internal node.
                    top++;
                    nodeStack[top] = internal;
                    depthStack[top] = lcpValue;
                    edgeStack[top] = new Edge(lastIncomingEdge.start, lastIncomingEdge.end, internal);
                }

                // Attach the leaf for this suffix under nodeStack[top].
                Node parent = nodeStack[top];

                Node leaf = Node.leaf(suffixStart);
                int edgeStart = suffixStart + depthStack[top];
                Edge leafEdge = new Edge(edgeStart, text.length, leaf);

                parent.putEdge(text[edgeStart], leafEdge);

                // Push the leaf so that the next suffix can build on it.
                top++;
                nodeStack[top] = leaf;
                depthStack[top] = text.length - suffixStart; // full suffix depth
                edgeStack[top] = leafEdge;
            }

            return root;
        }
    }
}
