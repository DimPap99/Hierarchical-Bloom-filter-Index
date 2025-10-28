package tree.ssws;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * Minimal suffix tree representation that keeps the API surface compatible with the
 * Farach-Colton construction while relying on a simple quadratic builder for now. Tokens
 * are assumed to be integers produced by the {@link utilities.AlphabetMapper}; a unique
 * sentinel value is appended during construction to preserve suffix uniqueness.
 */
public final class SuffixTreeQuadratic {

    private static final int SENTINEL = -1;

    private final Node root;
    private final int[] text; // includes the sentinel at the end
    private final int originalLength;

    private SuffixTreeQuadratic(Node root, int[] text, int originalLength) {
        this.root = root;
        this.text = text;
        this.originalLength = originalLength;
    }

    public static SuffixTreeQuadratic build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }
        int[] terminated = Arrays.copyOf(alphabetMappedText, alphabetMappedText.length + 1);
        terminated[terminated.length - 1] = SENTINEL;
        Node root = new Node();
        for (int i = 0; i < terminated.length; i++) {
            insertSuffix(root, terminated, i);
        }
        return new SuffixTreeQuadratic(root, terminated, alphabetMappedText.length);
    }

    public Node getRoot() {
        return root;
    }

    public int[] getText() {
        return Arrays.copyOf(text, text.length);
    }

    public int getOriginalLength() {
        return originalLength;
    }

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
            int edgeEnd = edge.end;
            int remaining = Math.min(edgeEnd, text.length) - edgeStart;
            int consumed = 0;
            while (consumed < remaining && patternIndex < pattern.length) {
                if (text[edgeStart + consumed] != pattern[patternIndex]) {
                    return Collections.emptyList();
                }
                consumed++;
                patternIndex++;
            }
            if (patternIndex == pattern.length) {
                return collectOccurrences(edge.child, pattern.length);
            }
            current = edge.child;
        }
        return collectOccurrences(current, pattern.length);
    }

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
            if (index >= 0 && index + patternLength <= originalLength) {
                out.add(index);
            }
            return;
        }
        for (Edge edge : node.edges.values()) {
            collect(edge.child, patternLength, out);
        }
    }

    private static void insertSuffix(Node root, int[] terminatedText, int start) {
        Node current = root;
        int index = start;
        while (index < terminatedText.length) {
            int symbol = terminatedText[index];
            Edge edge = current.edges.get(symbol);
            if (edge == null) {
                Node leaf = Node.leaf(start);
                current.edges.put(symbol, new Edge(index, terminatedText.length, leaf));
                return;
            }
            int edgeIndex = edge.start;
            int edgeEnd = edge.end;
            int matchLength = 0;
            while (edgeIndex + matchLength < edgeEnd && index + matchLength < terminatedText.length) {
                if (terminatedText[edgeIndex + matchLength] != terminatedText[index + matchLength]) {
                    break;
                }
                matchLength++;
            }
            if (edgeIndex + matchLength == edgeEnd) {
                current = edge.child;
                index += matchLength;
                continue;
            }
            if (index + matchLength == terminatedText.length) {
                // incoming suffix is already represented
                return;
            }
            Node split = new Node();
            // remainder of existing edge becomes child of the split node
            Edge tail = new Edge(edge.start + matchLength, edge.end, edge.child);
            split.edges.put(terminatedText[tail.start], tail);
            // new leaf for the remainder of the incoming suffix
            Node leaf = Node.leaf(start);
            Edge leafEdge = new Edge(index + matchLength, terminatedText.length, leaf);
            split.edges.put(terminatedText[leafEdge.start], leafEdge);
            // shorten original edge to the common prefix and point to the split node
            edge.end = edge.start + matchLength;
            edge.child = split;
            return;
        }
        // empty suffix; ensure a leaf exists off the root for the sentinel alone
        if (!current.edges.containsKey(SENTINEL)) {
            Node leaf = Node.leaf(start);
            current.edges.put(SENTINEL, new Edge(terminatedText.length - 1, terminatedText.length, leaf));
        }
    }

    public static final class Node {
        private final Map<Integer, Edge> edges = new HashMap<>();
        private int suffixIndex = -1;

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

    public static final class Edge {
        private int start;
        private int end;
        private Node child;

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
}