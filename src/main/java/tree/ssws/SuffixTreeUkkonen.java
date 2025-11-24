package tree.ssws;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

// Suffix tree over an integer alphabet with a unique sentinel appended.
public final class SuffixTreeUkkonen {

    // Sentinel value appended to every text. Must be smaller than any real token.
    private static final int SENTINEL = -1;

    // Root of the suffix tree.
    private final Node root;

    // Full text including the sentinel.
    private final int[] text;

    // Length of the original text before appending the sentinel.
    private final int originalLength;

    private SuffixTreeUkkonen(Node root, int[] text, int originalLength) {
        this.root = root;
        this.text = text;
        this.originalLength = originalLength;
    }

    // Build a suffix tree for the given integer token array using Ukkonen's algorithm.
    public static SuffixTreeUkkonen build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }
        // Append sentinel
        int[] terminated = Arrays.copyOf(alphabetMappedText, alphabetMappedText.length + 1);
        terminated[terminated.length - 1] = SENTINEL;

        // Build using Ukkonen
        UkkonenBuilder builder = new UkkonenBuilder(terminated);
        Node builtRoot = builder.build();

        // Wrap in final tree object
        return new SuffixTreeUkkonen(builtRoot, terminated, alphabetMappedText.length);
    }

    // Return the root node.
    public Node getRoot() {
        return root;
    }

    // Return a defensive copy of the underlying text, including the sentinel.
    public int[] getText() {
        return Arrays.copyOf(text, text.length);
    }

    // Return the number of original symbols (without the sentinel).
    public int getOriginalLength() {
        return originalLength;
    }

    // Find all starting indices where pattern occurs in the original text.
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
                // No outgoing edge that starts with this symbol.
                return Collections.emptyList();
            }

            int edgeStart = edge.getStart();
            int edgeEndExclusive = edge.getEnd();
            int edgeLen = edgeEndExclusive - edgeStart;

            int consumed = 0;
            while (consumed < edgeLen && patternIndex < pattern.length) {
                if (text[edgeStart + consumed] != pattern[patternIndex]) {
                    return Collections.emptyList();
                }
                consumed++;
                patternIndex++;
            }

            // If we consumed the entire pattern, every leaf below this child is an occurrence.
            if (patternIndex == pattern.length) {
                return collectOccurrences(edge.getChild(), pattern.length);
            }

            // Otherwise we matched the whole edge and still have pattern left; descend.
            current = edge.getChild();
        }

        // If we exit the loop naturally, we matched the entire pattern at a node.
        return collectOccurrences(current, pattern.length);
    }

    // Depth-first collection of all suffix positions under a node, filtered by originalLength.
    private List<Integer> collectOccurrences(Node node, int patternLength) {
        if (node == null) {
            return Collections.emptyList();
        }
        List<Integer> matches = new ArrayList<>();
        collect(node, patternLength, matches);
        return matches;
    }

    private void collect(Node node, int patternLength, List<Integer> out) {
        if (node.isLeaf()) {
            int index = node.suffixIndex;
            // Filter out suffixes that start beyond the original range or hit the sentinel.
            if (index >= 0 && index + patternLength <= originalLength) {
                out.add(index);
            }
            return;
        }
        for (Edge edge : node.edges.values()) {
            collect(edge.child, patternLength, out);
        }
    }

    // Node in the suffix tree.
    public static final class Node {
        private final Map<Integer, Edge> edges = new HashMap<>();

        // suffixIndex is meaningful for leaves.
        private int suffixIndex = -1;

        // suffixLink is used during construction, not needed for queries.
        private Node suffixLink = null;

        private Node() {
        }

        private static Node leaf(int suffixIndex) {
            Node node = new Node();
            node.suffixIndex = suffixIndex;
            return node;
        }

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

    // Edge in the suffix tree stored in implicit [start, end) form.
    public static final class Edge {
        private final int start;
        private final EndRef endRef;
        private final Node child;

        private Edge(int start, EndRef endRef, Node child) {
            this.start = start;
            this.endRef = endRef;
            this.child = child;
        }

        public int getStart() {
            return start;
        }

        /**
         * Return the exclusive end index of this edge's label in the text.
         */
        public int getEnd() {
            return endRef.value;
        }

        public Node getChild() {
            return child;
        }
    }

    /**
     * Tiny mutable integer holder so that all current leaf edges can share
     * the same end index in O(1) time. Ukkonen's algorithm updates exactly
     * one EndRef (the global "leaf end") per phase.
     */
    private static final class EndRef {
        int value;
        EndRef(int v) { this.value = v; }
    }

    /**
     * Internal builder that runs Ukkonen's algorithm on an integer array.
     * This constructs the compressed suffix tree in overall O(n) time.
     */
    private static final class UkkonenBuilder {
        private final int[] text;        // full text with sentinel
        private final Node root;         // root of the suffix tree
        private final EndRef leafEndRef; // global end reference for all current leaves

        // Ukkonen state
        private Node activeNode;
        private int activeEdgeIdx;       // index in 'text' that names the active edge
        private int activeLength;        // how many symbols we have matched on that edge

        private int remainingSuffixCount;
        private Node lastNewInternalNode;

        UkkonenBuilder(int[] text) {
            this.text = text;
            this.root = new Node();
            this.leafEndRef = new EndRef(0);

            this.activeNode = root;
            this.activeEdgeIdx = -1;
            this.activeLength = 0;
            this.remainingSuffixCount = 0;
            this.lastNewInternalNode = null;
        }

        /**
         * Run Ukkonen's algorithm over all characters in 'text'.
         * Then assign suffix indices to leaves.
         */
        Node build() {
            for (int pos = 0; pos < text.length; pos++) {
                extend(pos);
            }
            // After the full scan, assign suffixIndex values to leaves
            setSuffixIndicesDFS(root, 0);
            return root;
        }

        /**
         * Process text[pos]. This is one Ukkonen phase.
         *
         * We
         *   1. Increase the shared leaf end
         *   2. Add all pending suffixes (remainingSuffixCount)
         *   3. Maintain suffix links between newly created internal nodes
         */
        private void extend(int pos) {
            // Step 1. The new character at position pos is now part of all active leaves.
            leafEndRef.value = pos + 1;

            // We will need to insert all suffixes that end at 'pos'
            remainingSuffixCount++;
            lastNewInternalNode = null;

            while (remainingSuffixCount > 0) {

                // If we are not currently in the middle of an edge, the next suffix
                // we need to insert starts at position pos itself.
                if (activeLength == 0) {
                    activeEdgeIdx = pos;
                }

                int currentSymbol = text[activeEdgeIdx];
                Edge edge = activeNode.edges.get(currentSymbol);

                if (edge == null) {
                    // Rule 2 (no edge starting with currentSymbol)
                    // We create a new leaf edge from activeNode to a new leaf node.
                    Node leafNode = Node.leaf(-1);
                    Edge leafEdge = new Edge(pos, leafEndRef, leafNode);
                    activeNode.edges.put(currentSymbol, leafEdge);

                    // If we created an internal node in the previous iteration,
                    // its suffix link should point to activeNode.
                    if (lastNewInternalNode != null) {
                        lastNewInternalNode.suffixLink = activeNode;
                        lastNewInternalNode = null;
                    }
                } else {
                    // There is an edge starting with currentSymbol.
                    // Try to walk down this edge if our activeLength is large enough.
                    if (walkDown(edge)) {
                        continue;
                    }

                    // Check whether the next symbol on this edge matches text[pos].
                    int edgeNextSymbol = text[edge.start + activeLength];
                    if (edgeNextSymbol == text[pos]) {
                        // Rule 3 (current character is already on the edge)
                        // We just move forward one character in the active edge and stop.
                        if (lastNewInternalNode != null && activeNode != root) {
                            lastNewInternalNode.suffixLink = activeNode;
                            lastNewInternalNode = null;
                        }

                        activeLength++;
                        // Crucial: Rule 3 ends the phase early.
                        return;
                    }

                    // Otherwise we must split the edge.
                    // We create a new internal node 'splitNode' in the middle of 'edge'.
                    Node splitNode = new Node();

                    // Edge prefix from activeNode to splitNode.
                    // Its label is text[edge.start .. edge.start+activeLength)
                    EndRef splitEndRef = new EndRef(edge.start + activeLength);
                    Edge splitEdge = new Edge(edge.start, splitEndRef, splitNode);

                    // Remainder of the old edge becomes child of splitNode.
                    // Its label is text[edge.start+activeLength .. edge.getEnd())
                    Edge remainderEdge = new Edge(edge.start + activeLength, edge.endRef, edge.child);
                    splitNode.edges.put(text[remainderEdge.start], remainderEdge);

                    // New leaf edge for the new suffix starting at 'pos'
                    Node leafNode = Node.leaf(-1);
                    Edge leafEdge = new Edge(pos, leafEndRef, leafNode);
                    splitNode.edges.put(text[pos], leafEdge);

                    // Attach splitEdge to activeNode, replacing the old 'edge'
                    activeNode.edges.put(text[splitEdge.start], splitEdge);

                    // Suffix link maintenance
                    if (lastNewInternalNode != null) {
                        lastNewInternalNode.suffixLink = splitNode;
                    }
                    lastNewInternalNode = splitNode;
                }

                // One suffix has been added
                remainingSuffixCount--;

                // Update the active point
                if (activeNode == root && activeLength > 0) {
                    // If we are at the root and we have matched some characters
                    // on an edge, we shrink the activeLength by one and advance
                    // activeEdgeIdx to the next suffix to insert.
                    activeLength--;
                    activeEdgeIdx = pos - remainingSuffixCount + 1;
                } else if (activeNode != root) {
                    // Follow the suffix link if possible
                    activeNode = (activeNode.suffixLink != null) ? activeNode.suffixLink : root;
                }
            }
        }

        /**
         * Try to skip down the given edge if activeLength is long enough.
         * If we consume the entire edge, we move the active point down to the child.
         *
         * Return true if we moved down and should continue without finishing this extension step.
         * Return false if we are now somewhere in the middle of this edge.
         */
        private boolean walkDown(Edge edge) {
            int edgeLen = edge.getEnd() - edge.start;
            if (activeLength >= edgeLen) {
                activeEdgeIdx += edgeLen;
                activeLength -= edgeLen;
                activeNode = edge.child;
                return true;
            }
            return false;
        }

        /**
         * After the tree is fully built we assign each leaf node a suffixIndex,
         * which is the starting position in 'text' of the suffix represented by
         * that leaf.
         *
         * We do a depth first search. We track 'depth', the number of symbols
         * consumed from the root to reach the current node. If the total text
         * length is N, then a leaf reached at depth d represents the suffix
         * that starts at N - d.
         */
        private void setSuffixIndicesDFS(Node node, int depth) {
            if (node.edges.isEmpty()) {
                // Leaf
                node.suffixIndex = text.length - depth;
                return;
            }
            for (Edge e : node.edges.values()) {
                int edgeLen = e.getEnd() - e.start;
                setSuffixIndicesDFS(e.child, depth + edgeLen);
            }
        }
    }
}
