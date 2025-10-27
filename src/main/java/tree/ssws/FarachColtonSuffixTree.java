package tree.ssws;

import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * Minimal suffix tree representation built via a radix-sorted rank transformation. The class exposes the
 * skeleton required for the Farach-Colton construction while keeping the implementation simple enough to
 * iterate on. Future revisions will replace the naive insertion routine with the full randomized linear-time
 * construction described in the paper; the surrounding API already aligns with that effort.
 */
public final class FarachColtonSuffixTree {

    private static final int SENTINEL = -1;

    private final Node root;
    private final int[] text; // includes the sentinel at the end

    private FarachColtonSuffixTree(Node root, int[] text) {
        this.root = root;
        this.text = text;
    }

    public static FarachColtonSuffixTree build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }
        int[] ranked = RadixSorter.rankTransform(alphabetMappedText);
        int[] terminated = Arrays.copyOf(ranked, ranked.length + 1);
        terminated[terminated.length - 1] = SENTINEL;
        Node root = new Node();
        for (int i = 0; i < terminated.length; i++) {
            insertSuffix(root, terminated, i);
        }
        return new FarachColtonSuffixTree(root, terminated);
    }

    public Node getRoot() {
        return root;
    }

    public int[] getText() {
        return Arrays.copyOf(text, text.length);
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
