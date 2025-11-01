package tree.ssws;

import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * SuffixTree
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
public final class SuffixTree {

    private final Node root;
    private static final  int SENTINEL =0; //reserve 0 as sentinel
    private final int[] text;           // includes sentinel at end (remapped per-tree)
    private final int originalLength;   // length before sentinel
    // Map original token -> local remapped token in [1..distinct]. Missing => not present.
    private final Int2IntOpenHashMap tokenRemap;
    private boolean released = false;

    private SuffixTree(Node root, int[] text, int originalLength, Int2IntOpenHashMap tokenRemap) {
        this.root = root;
        this.text = text;
        this.originalLength = originalLength;
        this.tokenRemap = tokenRemap;
    }

    /**
     * Build a suffix tree from an integer token array. We:
     *  - append a unique sentinel strictly smaller than any input token,
     *  - create a dense non-negative alphabet via radix ranking,
     *  - run DC3 / Skew to get the suffix array,
     *  - run Kasai to get the LCP,
     *  - build a compressed suffix tree in linear time.
     */
    /**
     * Build a suffix tree from an integer token array that already obeys:
     *   - All real tokens are >= 1.
     *   - The integer 0 is globally reserved for the sentinel only.
     *
     * Construction steps in this optimized version:
     *   1. Append SENTINEL (which must be 0) to the end of the array.
     *   2. Run DC3 / Skew directly on the terminated array to get the suffix array.
     *   3. Run Kasai on the same terminated array to get the LCP array.
     *   4. Build the compressed suffix tree from suffix array + LCP.
     *
     * This skips the expensive rank-normalization and 64-bit radix prepass,
     * because the caller guarantees that 0 is already the unique global minimum.
     */
    public static SuffixTree build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }

        // Precondition on incoming symbols from the global AlphabetMapper:
        //  - All real symbols are >= 1 (0 is reserved for the sentinel globally).
        // However, their numeric IDs can be very large (global alphabet growth),
        // which would blow up counting-sort buckets inside DC3 if used directly.
        // To keep construction and queries independent of global alphabet size,
        // we remap the per-text symbols to a dense local range 1..distinct.

        final int n0 = alphabetMappedText.length;

        // Build a dense remap: collect and sort a copy, assign ranks in [1..distinct].
        // We avoid allocating an enormous direct-address table by using a compact map.
        int[] copy = Arrays.copyOf(alphabetMappedText, n0);
        // Sort in-place (LSD radix) to get non-decreasing order.
        RadixSorter.sortInPlace(copy);
        Int2IntOpenHashMap remap = new Int2IntOpenHashMap(Math.max(4, n0 * 2));
        remap.defaultReturnValue(-1);
        int rank = 1; // start at 1; 0 is reserved for the sentinel
        int prev = Integer.MIN_VALUE;
        for (int v : copy) {
            if (v != prev) {
                remap.put(v, rank++);
                prev = v;
            }
        }

        // Apply remap to the text and append sentinel 0 at the end.
        final int[] terminated = new int[n0 + 1];
        for (int i = 0; i < n0; i++) {
            int mapped = remap.getOrDefault(alphabetMappedText[i], -1);
            if (mapped <= 0) {
                // Should never happen as all symbols in alphabetMappedText were inserted above.
                // Guard defensively to keep the structure valid.
                throw new IllegalStateException("Missing remap for symbol " + alphabetMappedText[i]);
            }
            terminated[i] = mapped;
        }
        terminated[n0] = SENTINEL; // 0 is unique sentinel

        final int n = n0 + 1;

        // Step 1. Build suffix array directly on the terminated integer array.
        // SuffixArrayBuilder.buildSuffixArray requires:
        //   - All values are non-negative.
        //   - The smallest value is the unique sentinel at the end.
        // With SENTINEL = 0 and all real tokens >= 1, these hold.
        int[] sa = SuffixArrayBuilder.buildSuffixArray(terminated);

        // Step 2. Build Longest Common Prefix (LCP) via Kasai on the SAME terminated array.
        int[] lcp = (n > 1)
                ? SuffixArrayBuilder.buildLcpArray(terminated, sa)
                : new int[0];

        // Step 3. Build explicit compressed suffix tree from suffix array + LCP.
        Node root = LinearBuilder.buildFromSuffixArray(terminated, sa, lcp);

        // originalLength is the logical text length before the sentinel was appended.
        return new SuffixTree(root, terminated, n0, remap);
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

        // Translate pattern symbols via per-tree remap. If any symbol
        // does not exist in this tree, there can be no occurrences.
        int[] q = new int[pattern.length];
        for (int i = 0; i < pattern.length; i++) {
            int mapped = tokenRemap.getOrDefault(pattern[i], -1);
            if (mapped <= 0) {
                return Collections.emptyList();
            }
            q[i] = mapped;
        }

        Node current = root;
        int patternIndex = 0;

        while (patternIndex < pattern.length) {

            int symbol = q[patternIndex];
            Edge edge = current.getEdge(symbol);
            if (edge == null) {
                return Collections.emptyList();
            }

            int edgeStart = edge.start;
            int edgeEndExclusive = edge.end;
            int edgeLen = edgeEndExclusive - edgeStart;

            int consumed = 0;
            while (consumed < edgeLen && patternIndex < pattern.length) {
                if (text[edgeStart + consumed] != q[patternIndex]) {
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
     * Release internal structures eagerly so that detached trees relinquish memory sooner.
     * Safe to call multiple times.
     */
    public void release() {
        if (released) {
            return;
        }
        released = true;
        Arrays.fill(text, 0);
        tokenRemap.clear();
        clearNode(root);
    }

    private void clearNode(Node node) {
        if (node == null) {
            return;
        }
        if (node.edgeA != null) {
            clearNode(node.edgeA.child);
            node.edgeA = null;
            node.keyA = Integer.MIN_VALUE;
        }
        if (node.edgeB != null) {
            clearNode(node.edgeB.child);
            node.edgeB = null;
            node.keyB = Integer.MIN_VALUE;
        }
        if (node.edgeMap != null) {
            for (Edge edge : node.edgeMap.values()) {
                clearNode(edge.child);
            }
            node.edgeMap.clear();
            node.edgeMap = null;
        }
        node.suffixIndex = -1;
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
