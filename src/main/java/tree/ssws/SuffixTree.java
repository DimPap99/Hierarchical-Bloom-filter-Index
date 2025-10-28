package tree.ssws;

import java.util.*;
import java.util.function.Consumer;

/**
 * SuffixTree
 *
 * This class builds a suffix tree over an integer token array using Ukkonen's algorithm
 * in linear time in the length of the input in the sense of big O notation O(n).
 *
 * Public contract
 *
 *   SuffixTree.build(int[] tokens)
 *     Build a compressed suffix tree for the given integer token array.
 *     A single sentinel value is appended automatically so that all suffixes
 *     become explicit leaves. The caller must guarantee that no real token is -1.
 *
 *   findOccurrences(int[] pattern)
 *     Return all start offsets where "pattern" occurs in the original text
 *     that is, before the sentinel was appended. Offsets are zero based.
 *
 * Why this is faster and more memory friendly than the naive version
 *
 *   1. We do not create Edge objects. Each Node stores the edge label from its
 *      parent directly as (start, end). The "end" position is inclusive, so the
 *      edge covers text[start .. end], both indices included.
 *
 *   2. Leaf edges all share a single global integer end called leafEndValue,
 *      instead of each leaf edge having its own EndRef object. A leaf simply
 *      stores end = LEAF_SENTINEL. During traversal we interpret LEAF_SENTINEL
 *      as leafEndValue. This matches what your SuffixTreeIndex does for the
 *      streaming case.
 *
 *   3. Children of a node are stored in a custom IntChildMap, which uses four
 *      inline slots before upgrading to a tiny open addressed hash table called
 *      IntNodeMap. This avoids allocating java.util.HashMap for almost every
 *      node, which is a very big win for suffix trees, because most internal
 *      nodes only have one or two outgoing children.
 *
 *   4. We assign suffixIndex at leaf creation time. That means we never need
 *      a cleanup traversal just to figure out which suffix each leaf represents.
 *
 *   5. After we finish building the tree, we run a light compaction pass that
 *      trims maps to minimal capacity and drops suffix links since suffix links
 *      are only needed while running Ukkonen. This matches the compactForQuerying
 *      idea in your SuffixTreeIndex.
 *
 * Complexity summary
 *
 *   Build time is linear in the number of tokens, plus the sentinel,
 *   in the sense of big O notation O(n).
 *   Query time for findOccurrences is linear in the pattern length
 *   plus the number of matches that are reported, which is standard
 *   for suffix tree search.
 */
public final class SuffixTree {

    /**
     * Sentinel value appended to the text. The caller must ensure this does not
     * appear in the original token stream. It must compare strictly less than
     * any real token if you rely on lexicographic uniqueness of suffixes.
     * For most tokenizers that map natural tokens to non negative integers,
     * -1 is safe.
     */
    private static final int SENTINEL = -1;

    /**
     * Root of the suffix tree.
     */
    private final Node root;

    /**
     * The full text including the sentinel at the end.
     */
    private final int[] text;

    /**
     * Length of the original text before appending the sentinel.
     */
    private final int originalLength;

    /**
     * Final inclusive end index for all open leaf edges.
     * After build completes this is text.length - 1.
     * We keep it so that effectiveEnd works for leaves.
     */
    private final int finalLeafEndValue;

    private SuffixTree(Node root,
                               int[] text,
                               int originalLength,
                               int finalLeafEndValue) {
        this.root = root;
        this.text = text;
        this.originalLength = originalLength;
        this.finalLeafEndValue = finalLeafEndValue;
    }

    /**
     * Build a suffix tree for the given integer token array using Ukkonen's algorithm.
     * A unique sentinel is appended so that all suffixes become explicit and terminate
     * at leaves. This guarantees that every suffix corresponds to exactly one leaf.
     *
     * The resulting tree uses the same allocation strategy as SuffixTreeIndex
     * in other words, Node holds (start, end) for the incoming edge and keeps its
     * children in a compact IntChildMap rather than a java.util.HashMap.
     */
    public static SuffixTree build(int[] alphabetMappedText) {
        if (alphabetMappedText == null) {
            throw new IllegalArgumentException("text cannot be null");
        }

        // Append the sentinel
        int n0 = alphabetMappedText.length;
        int[] terminated = Arrays.copyOf(alphabetMappedText, n0 + 1);
        terminated[n0] = SENTINEL;

        // Run Ukkonen's algorithm in a builder that uses the compact node layout
        UkkonenBuilder builder = new UkkonenBuilder(terminated);
        Node builtRoot = builder.build();

        // After building we optionally compact suffix links and trim maps
        compactForQuerying(builtRoot);

        // Wrap in final tree object
        return new SuffixTree(
                builtRoot,
                terminated,
                n0,
                builder.getLeafEndValue()
        );
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
     * Return the number of original symbols without the sentinel.
     */
    public int getOriginalLength() {
        return originalLength;
    }

    /**
     * Pattern search
     *
     * We walk the pattern down the tree starting at the root.
     * At each step we choose the child edge whose first token matches
     * the next pattern token. That child node stores the edge label
     * boundaries in the fields start and end. We then compare the pattern
     * tokens to text[start .. end] directly. We never allocate new strings
     * and we never rescan the same pattern position twice.
     *
     * If we consume the entire pattern in the middle of an edge, that is
     * still a match, because compressed suffix tree edges are path compressed
     * and there is no branching in the middle of an edge label. All suffixes
     * hanging under the child node correspond to that match.
     *
     * Once the pattern is fully matched, we collect all leaves reachable
     * under the current node. For each leaf we report its suffixIndex,
     * as long as that match stays within the originalLength and does not
     * run into the sentinel.
     *
     * Time cost is linear in the length of the pattern plus the number
     * of matches reported. That is the classic suffix tree query bound.
     */
    public List<Integer> findOccurrences(int[] pattern) {
        if (pattern == null || pattern.length == 0) {
            return Collections.emptyList();
        }

        Node current = root;
        int patternIndex = 0;

        while (patternIndex < pattern.length) {
            int symbol = pattern[patternIndex];
            Node next = getChild(current, symbol);
            if (next == null) {
                return Collections.emptyList();
            }

            int edgeStart = next.start;
            int edgeEndInclusive = next.effectiveEnd(finalLeafEndValue);
            int cursor = edgeStart;

            // Walk along this edge while tokens match
            while (patternIndex < pattern.length && cursor <= edgeEndInclusive) {
                if (text[cursor] != pattern[patternIndex]) {
                    return Collections.emptyList();
                }
                cursor++;
                patternIndex++;
            }

            // If we consumed the full pattern, collect leaves below 'next'
            if (patternIndex == pattern.length) {
                ArrayList<Integer> matches = new ArrayList<>();
                collectLeafStarts(next, pattern.length, matches);
                return matches;
            }

            // Otherwise we must have consumed the entire edge
            if (cursor > edgeEndInclusive) {
                current = next;
                continue;
            }

            // We got stuck in the middle of the edge without consuming
            // the full pattern. That means mismatch
            return Collections.emptyList();
        }

        // If we exit naturally, we matched the pattern exactly at an explicit node
        ArrayList<Integer> matches = new ArrayList<>();
        collectLeafStarts(current, pattern.length, matches);
        return matches;
    }

    /**
     * Depth first traversal of the subtree under node, collecting starting
     * positions for all suffixes at leaves. We filter out suffixes that
     * would extend past the originalLength, so matches that include the
     * sentinel are not reported.
     */
    private void collectLeafStarts(Node node, int patternLength, List<Integer> out) {
        if (node == null) {
            return;
        }

        if (!hasChildren(node)) {
            // Leaf case
            // suffixIndex is set during construction for leaves
            int startIndex = (node.suffixIndex >= 0) ? node.suffixIndex : node.start;
            if (startIndex >= 0 && startIndex + patternLength <= originalLength) {
                out.add(startIndex);
            }
            return;
        }

        forEachChild(node, child -> collectLeafStarts(child, patternLength, out));
    }

    /**
     * After building we want to drop construction only pointers like suffixLink
     * and trim the child maps so they do not waste capacity.
     *
     * This is similar to compactForQuerying in SuffixTreeIndex.
     * It makes the final tree cheaper to keep resident in memory.
     * Under the landmark model and to be fair since we did not implement Farach Colton as in the original paper
     * before any memory bench marking this is used
     */
    public static void compactForQuerying(Node node) {
        if (node == null) {
            return;
        }
        if (hasChildren(node)) {
            // Trim the map itself first
            node.children.trimToSize();
            // Then recurse
            forEachChild(node, SuffixTree::compactForQuerying);
        } else {
            // No children at all, drop the structure to null to save memory
            node.children = null;
        }
        // Suffix links are only needed during construction
        node.suffixLink = null;
    }

    /* =========================================================================
       Internal node, child map, and Ukkonen builder
       ========================================================================= */

    /**
     * Node represents both
     *   1. an edge label from its parent, given by the inclusive [start, end]
     *      indices into the global text array, and
     *   2. a branching point with zero or more outgoing children.
     *
     * This is the exact same trick you are using in SuffixTreeIndex.
     *
     * Fields
     *   children       Either null, or an IntChildMap for outgoing edges
     *   suffixLink     Used only during Ukkonen construction, cleared later
     *   start          Inclusive start index in text[] of the edge label
     *   end            Inclusive end index in text[] of the edge label
     *                  or LEAF_SENTINEL for "use the global leaf end"
     *   suffixIndex    For leaves, the starting index in the text for that suffix
     *
     * LEAF_SENTINEL means "this node is currently a leaf whose edge keeps growing"
     * and we should read its effective end from the shared leafEndValue
     */
    public static final class Node {
        private IntChildMap children;
        private Node suffixLink;
        private int start;
        private int end; // inclusive, or LEAF_SENTINEL
        private int suffixIndex = -1;

        private static final int LEAF_SENTINEL = Integer.MIN_VALUE;

        private Node(int start, int end) {
            this.start = start;
            this.end = end;
        }

        private int effectiveEnd(int leafEndValue) {
            return (end == LEAF_SENTINEL) ? leafEndValue : end;
        }

        private int edgeLength(int leafEndValue) {
            if (start == -1) {
                return 0;
            }
            return effectiveEnd(leafEndValue) - start + 1;
        }
    }

    /**
     * IntChildMap is a compact associative container from "first token on outgoing edge"
     * to the child Node that represents that edge.
     *
     * It starts with four inline slots, then upgrades to an IntNodeMap which is a
     * small open addressed hash table on primitive int keys. This avoids creating
     * java.util.HashMap and it avoids boxing keys into Integer objects.
     *
     * This matches your ChildMap but specialized to int instead of long, so that
     * it can be used with int token arrays.
     */
    private static final class IntChildMap {
        // mode 0 empty, 1..4 number of inline entries, 5 means using hash table
        private byte mode = 0;

        private int k1, k2, k3, k4;
        private Node v1, v2, v3, v4;

        private IntNodeMap map;

        boolean isEmpty() {
            return mode == 0;
        }

        Node get(int key) {
            switch (mode) {
                case 0:  return null;
                case 1:  return (k1 == key ? v1 : null);
                case 2:  return (k1 == key ? v1 : (k2 == key ? v2 : null));
                case 3:  return (k1 == key ? v1 :
                        (k2 == key ? v2 :
                                (k3 == key ? v3 : null)));
                case 4:  return (k1 == key ? v1 :
                        (k2 == key ? v2 :
                                (k3 == key ? v3 :
                                        (k4 == key ? v4 : null))));
                default:
                    return map.get(key);
            }
        }

        void put(int key, Node value) {
            switch (mode) {
                case 0:
                    k1 = key; v1 = value; mode = 1; return;
                case 1:
                    if (k1 == key) { v1 = value; return; }
                    k2 = key; v2 = value; mode = 2; return;
                case 2:
                    if (k1 == key) { v1 = value; return; }
                    if (k2 == key) { v2 = value; return; }
                    k3 = key; v3 = value; mode = 3; return;
                case 3:
                    if (k1 == key) { v1 = value; return; }
                    if (k2 == key) { v2 = value; return; }
                    if (k3 == key) { v3 = value; return; }
                    k4 = key; v4 = value; mode = 4; return;
                case 4:
                    if (k1 == key) { v1 = value; return; }
                    if (k2 == key) { v2 = value; return; }
                    if (k3 == key) { v3 = value; return; }
                    if (k4 == key) { v4 = value; return; }
                    // upgrade to hash table
                    map = new IntNodeMap(8);
                    map.put(k1, v1);
                    map.put(k2, v2);
                    map.put(k3, v3);
                    map.put(k4, v4);
                    // free inline slots
                    v1 = v2 = v3 = v4 = null;
                    mode = 5;
                    map.put(key, value);
                    return;
                default:
                    map.put(key, value);
            }
        }

        void trimToSize() {
            if (mode == 5 && map != null) {
                map.trimToSize();
            }
        }

        void forEach(Consumer<Node> consumer) {
            switch (mode) {
                case 0:  return;
                case 1:  consumer.accept(v1); return;
                case 2:  consumer.accept(v1); consumer.accept(v2); return;
                case 3:  consumer.accept(v1); consumer.accept(v2); consumer.accept(v3); return;
                case 4:  consumer.accept(v1); consumer.accept(v2); consumer.accept(v3); consumer.accept(v4); return;
                default: map.forEach(consumer);
            }
        }
    }

    /**
     * IntNodeMap is the open addressed fallback map for IntChildMap once we exceed
     * four children. It uses linear probing, keeps primitive int keys and Node values,
     * and supports trimming after the build is complete.
     */
    private static final class IntNodeMap {
        private static final float LOAD_FACTOR = 0.75f;

        private int[] keys;
        private Node[] values;
        private byte[] states; // 0 empty, 1 occupied, 2 deleted (2 is unused here)
        private int mask;
        private int size;
        private int threshold;

        IntNodeMap(int initialCapacity) {
            initialise(Math.max(4, initialCapacity));
        }

        Node get(int key) {
            int idx = findIndex(key);
            return (idx >= 0) ? values[idx] : null;
        }

        void put(int key, Node value) {
            int slot = findSlot(key);
            if (states[slot] == 1) {
                values[slot] = value;
                return;
            }
            keys[slot] = key;
            values[slot] = value;
            states[slot] = 1;
            size++;
            if (size >= threshold) {
                rehash(keys.length << 1);
            }
        }

        void trimToSize() {
            if (size == 0) {
                initialise(4);
                return;
            }
            int minNeeded = (int) Math.ceil(size / (double) LOAD_FACTOR);
            int target = 1;
            while (target < minNeeded) {
                target <<= 1;
            }
            if (keys.length != target) {
                rehash(target);
            }
        }

        void forEach(Consumer<Node> consumer) {
            for (int i = 0; i < states.length; i++) {
                if (states[i] == 1) {
                    consumer.accept(values[i]);
                }
            }
        }

        private int findIndex(int key) {
            int idx = mix(key) & mask;
            while (states[idx] != 0) {
                if (states[idx] == 1 && keys[idx] == key) {
                    return idx;
                }
                idx = (idx + 1) & mask;
            }
            return -1;
        }

        private int findSlot(int key) {
            int idx = mix(key) & mask;
            while (states[idx] == 1 && keys[idx] != key) {
                idx = (idx + 1) & mask;
            }
            return idx;
        }

        private void rehash(int newCapacity) {
            int[] oldKeys = keys;
            Node[] oldValues = values;
            byte[] oldStates = states;

            initialise(newCapacity);
            for (int i = 0; i < oldStates.length; i++) {
                if (oldStates[i] == 1) {
                    put(oldKeys[i], oldValues[i]);
                }
            }
        }

        private void initialise(int capacity) {
            int pow2 = 1;
            while (pow2 < capacity) {
                pow2 <<= 1;
            }
            keys = new int[pow2];
            values = new Node[pow2];
            states = new byte[pow2];
            mask = pow2 - 1;
            threshold = (int) (pow2 * LOAD_FACTOR);
            size = 0;
        }

        private int mix(int key) {
            // 32 bit mix similar to Murmur style avalanche
            int z = key;
            z ^= (z >>> 16);
            z *= 0x85ebca6b;
            z ^= (z >>> 13);
            z *= 0xc2b2ae35;
            z ^= (z >>> 16);
            return z;
        }
    }

    /**
     * UkkonenBuilder implements Ukkonen's algorithm (which is an online linear time suffix tree
     * construction algorithm) for an integer array. It constructs the final Node graph directly
     * using the compact Node and IntChildMap structures above. It never creates intermediate
     * Edge objects and never creates a java.util.HashMap per node.
     *
     * The code structure is deliberately very close to SuffixTreeIndex.extendSuffixTree.
     * The difference is that here we are iterating across a fixed int[] called text,
     * not across a growing LongSequence, but the mechanics are the same.
     */
    private static final class UkkonenBuilder {
        private final int[] text;       // full text including sentinel
        private final Node root;        // suffix tree root

        // Active point and global state from Ukkonen's algorithm
        private Node activeNode;
        private int activeEdge = -1;    // index in text naming the active edge
        private int activeLength = 0;

        private int remainingSuffixCount = 0;
        private int leafEndValue = -1;  // inclusive end shared by all current leaves
        private Node lastCreatedInternalNode = null;

        UkkonenBuilder(int[] text) {
            this.text = text;
            this.root = new Node(-1, -1);
            this.root.suffixLink = this.root;
            this.activeNode = this.root;
        }

        Node build() {
            for (int pos = 0; pos < text.length; pos++) {
                extend(pos);
            }
            return root;
        }

        int getLeafEndValue() {
            return leafEndValue;
        }

        private void extend(int pos) {
            // The new symbol at position pos extends all current leaves
            leafEndValue = pos;
            remainingSuffixCount++;
            lastCreatedInternalNode = null;

            while (remainingSuffixCount > 0) {
                if (activeLength == 0) {
                    activeEdge = pos;
                }

                int currentEdgeToken = text[activeEdge];
                Node next = getChild(activeNode, currentEdgeToken);

                if (next == null) {
                    // Create a new leaf
                    Node leaf = new Node(pos, Node.LEAF_SENTINEL);
                    // This suffix starts at (pos - remainingSuffixCount + 1)
                    leaf.suffixIndex = pos - remainingSuffixCount + 1;
                    putChild(activeNode, currentEdgeToken, leaf);

                    // Suffix link hookup for the last created internal node
                    if (lastCreatedInternalNode != null) {
                        lastCreatedInternalNode.suffixLink = activeNode;
                        lastCreatedInternalNode = null;
                    }
                } else {
                    if (walkDown(next)) {
                        continue;
                    }

                    int nextToken = text[next.start + activeLength];
                    int currentToken = text[pos];
                    if (nextToken == currentToken) {
                        // Rule 3 in Ukkonen
                        if (lastCreatedInternalNode != null && activeNode != root) {
                            lastCreatedInternalNode.suffixLink = activeNode;
                            lastCreatedInternalNode = null;
                        }
                        activeLength++;
                        break;
                    }

                    // We have a mismatch in the middle of an edge.
                    // We split that edge and create a new internal node.
                    int splitEndVal = next.start + activeLength - 1;
                    Node split = new Node(next.start, splitEndVal);
                    putChild(activeNode, currentEdgeToken, split);

                    // New leaf from current position
                    Node leaf = new Node(pos, Node.LEAF_SENTINEL);
                    leaf.suffixIndex = pos - remainingSuffixCount + 1;
                    putChild(split, currentToken, leaf);

                    // Advance 'next' to start after the split and hook it under split
                    next.start = next.start + activeLength;
                    int nextEdgeToken = text[next.start];
                    putChild(split, nextEdgeToken, next);

                    // Suffix link stitching
                    if (lastCreatedInternalNode != null) {
                        lastCreatedInternalNode.suffixLink = split;
                    }
                    lastCreatedInternalNode = split;
                    split.suffixLink = root;
                }

                remainingSuffixCount--;

                if (activeNode == root && activeLength > 0) {
                    activeLength--;
                    activeEdge = pos - remainingSuffixCount + 1;
                } else if (activeNode != root) {
                    activeNode = (activeNode.suffixLink != null) ? activeNode.suffixLink : root;
                }
            }
        }

        private boolean walkDown(Node next) {
            int edgeLength = next.edgeLength(leafEndValue);
            if (activeLength >= edgeLength) {
                activeEdge += edgeLength;
                activeLength -= edgeLength;
                activeNode = next;
                return true;
            }
            return false;
        }
    }

    /* =========================================================================
       Tiny helpers mirroring SuffixTreeIndex helpers but for int keys
       ========================================================================= */

    private static Node getChild(Node n, int token) {
        return (n.children == null) ? null : n.children.get(token);
    }

    private static void putChild(Node n, int token, Node child) {
        if (n.children == null) {
            n.children = new IntChildMap();
        }
        n.children.put(token, child);
    }

    private static boolean hasChildren(Node n) {
        return n.children != null && !n.children.isEmpty();
    }

    private static void forEachChild(Node n, Consumer<Node> c) {
        if (n.children != null) {
            n.children.forEach(c);
        }
    }
}
