package PMIndex;

import search.Pattern;
import utilities.PatternResult;
import utilities.StringKeyMapper;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Objects;
import java.util.function.Consumer;

/**
 * Suffix tree that stores hashed n-gram tokens instead of raw characters. Incoming
 * strings are mapped to longs via {@link StringKeyMapper} so the representation
 * matches the numeric token stream used by HBI. Ukkonen's algorithm keeps the tree
 * in sync as tokens arrive, allowing pattern queries to execute in linear time with
 * respect to the number of tokens in the query.
 */
public class SuffixTreeIndex implements IPMIndexing {

    private final StringKeyMapper keyMapper;
    private final LongSequence tokenStream = new LongSequence();
    private final Node root;

    private Node activeNode;
    private int activeEdge = -1;
    private int activeLength = 0;
    private int remainingSuffixCount = 0;
    private End leafEnd = new End(-1);
    private Node lastCreatedInternalNode = null;

    public SuffixTreeIndex() {
        this(new StringKeyMapper(1_000_000L, 0.0001));
    }

    public SuffixTreeIndex(long expectedDistinctTokens, double epsilon) {
        this(new StringKeyMapper(expectedDistinctTokens, epsilon));
    }

    public SuffixTreeIndex(StringKeyMapper mapper) {
        this.keyMapper = Objects.requireNonNull(mapper, "mapper");
        this.root = new Node(-1, new End(-1), null);
        this.root.suffixLink = this.root;
        this.activeNode = this.root;
    }

    @Override
    public void insert(String key) {
        if (key == null || key.isEmpty()) {
            return;
        }
        long token = keyMapper.mapToLong(key);
        int pos = tokenStream.add(token);
        extendSuffixTree(pos);
    }

    private void extendSuffixTree(int pos) {
        leafEnd.value = pos;
        remainingSuffixCount++;
        lastCreatedInternalNode = null;

        while (remainingSuffixCount > 0) {
            if (activeLength == 0) {
                activeEdge = pos;
            }
            long currentEdgeToken = tokenStream.get(activeEdge);
            Node next = activeNode.children.get(currentEdgeToken);
            if (next == null) {
                Node leaf = new Node(pos, leafEnd, activeNode);
                leaf.suffixIndex = pos - remainingSuffixCount + 1;
                activeNode.children.put(currentEdgeToken, leaf);

                if (lastCreatedInternalNode != null) {
                    lastCreatedInternalNode.suffixLink = activeNode;
                    lastCreatedInternalNode = null;
                }
            } else {
                if (walkDown(next)) {
                    continue;
                }
                long nextToken = tokenStream.get(next.start + activeLength);
                long currentToken = tokenStream.get(pos);
                if (nextToken == currentToken) {
                    if (lastCreatedInternalNode != null && activeNode != root) {
                        lastCreatedInternalNode.suffixLink = activeNode;
                        lastCreatedInternalNode = null;
                    }
                    activeLength++;
                    break;
                }

                End splitEnd = new End(next.start + activeLength - 1);
                Node split = new Node(next.start, splitEnd, activeNode);
                activeNode.children.put(currentEdgeToken, split);

                Node leaf = new Node(pos, leafEnd, split);
                leaf.suffixIndex = pos - remainingSuffixCount + 1;
                split.children.put(currentToken, leaf);

                next.start = next.start + activeLength;
                next.parent = split;
                long nextEdgeToken = tokenStream.get(next.start);
                split.children.put(nextEdgeToken, next);

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
                activeNode = activeNode.suffixLink != null ? activeNode.suffixLink : root;
            }
        }
    }

    private boolean walkDown(Node next) {
        int edgeLength = next.edgeLength();
        if (activeLength >= edgeLength) {
            activeEdge += edgeLength;
            activeLength -= edgeLength;
            activeNode = next;
            return true;
        }
        return false;
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(Pattern key) {
        if (key == null) {
            return new ArrayList<>();
        }
        long[] patternTokens = toTokenSequence(key);
        return reportInternal(patternTokens);
    }

    private long[] toTokenSequence(Pattern pattern) {
        if (pattern.nGramArr == null || pattern.nGramArr.length == 0) {
            return new long[0];
        }
        int len = pattern.nGramArr.length;
        long[] tokens = new long[len];
        for (int i = 0; i < len; i++) {
            long token = keyMapper.mapToLong(pattern.nGramArr[i]);
            tokens[i] = token;
            if (pattern.nGramToLong != null && i < pattern.nGramToLong.length) {
                pattern.nGramToLong[i] = token;
            }
            if (pattern.effectiveNgramArr != null && pattern.nGram > 0 && i % pattern.nGram == 0) {
                int idx = i / pattern.nGram;
                if (idx < pattern.effectiveNgramArr.length) {
                    pattern.effectiveNgramArr[idx] = token;
                }
            }
        }
        if (pattern.effectiveNgramArr != null && pattern.effectiveNgramArr.length > 0
                && pattern.nGramToLong != null && pattern.nGramToLong.length > 0
                && pattern.nGram > 0 && pattern.originalSz % pattern.nGram != 0) {
            pattern.effectiveNgramArr[pattern.effectiveNgramArr.length - 1] =
                    pattern.nGramToLong[pattern.nGramToLong.length - 1];
        }
        return tokens;
    }

    private ArrayList<Integer> reportInternal(long[] patternTokens) {
        ArrayList<Integer> result = new ArrayList<>();
        if (patternTokens.length == 0) {
            for (int i = 0; i < tokenStream.size(); i++) {
                result.add(i);
            }
            return result;
        }
        Node current = root;
        int patternIndex = 0;

        while (patternIndex < patternTokens.length) {
            long searchToken = patternTokens[patternIndex];
            Node next = current.children.get(searchToken);
            if (next == null) {
                return result;
            }

            int edgeStart = next.start;
            int edgeEnd = next.end.value;
            int edgeIndex = edgeStart;
            while (patternIndex < patternTokens.length && edgeIndex <= edgeEnd) {
                if (patternTokens[patternIndex] != tokenStream.get(edgeIndex)) {
                    return new ArrayList<>();
                }
                patternIndex++;
                edgeIndex++;
            }

            if (patternIndex == patternTokens.length) {
                collectSuffixIndices(next, result);
                Collections.sort(result);
                return result;
            }

            if (edgeIndex > edgeEnd) {
                current = next;
            } else {
                return new ArrayList<>();
            }
        }
        Collections.sort(result);
        return result;
    }

    private void collectSuffixIndices(Node node, ArrayList<Integer> result) {
        if (node == null) {
            return;
        }
        if (node.suffixIndex >= 0) {
            result.add(node.suffixIndex);
            return;
        }
        node.children.forEach(child -> collectSuffixIndices(child, result));
    }

    @Override
    public void expire() {
        tokenStream.clear();
        root.children.clear();
        root.start = -1;
        root.end = new End(-1);
        root.suffixLink = root;
        root.parent = null;
        activeNode = root;
        activeEdge = -1;
        activeLength = 0;
        remainingSuffixCount = 0;
        leafEnd = new End(-1);
        lastCreatedInternalNode = null;
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

    private static final class Node {
        private final LongNodeMap children = new LongNodeMap();
        private Node suffixLink;
        private int start;
        private End end;
        private Node parent;
        private int suffixIndex = -1;

        private Node(int start, End end, Node parent) {
            this.start = start;
            this.end = end;
            this.parent = parent;
        }

        private int edgeLength() {
            if (start == -1) {
                return 0;
            }
            return end.value - start + 1;
        }
    }

    private static final class End {
        private int value;

        private End(int value) {
            this.value = value;
        }
    }

    private static final class LongSequence {
        private long[] data = new long[16];
        private int size = 0;

        int add(long value) {
            ensureCapacity(size + 1);
            data[size] = value;
            return size++;
        }

        long get(int index) {
            if (index < 0 || index >= size) {
                throw new IndexOutOfBoundsException("Index: " + index + ", size: " + size);
            }
            return data[index];
        }

        int size() {
            return size;
        }

        void clear() {
            size = 0;
        }

        private void ensureCapacity(int capacity) {
            if (capacity <= data.length) {
                return;
            }
            int newCapacity = data.length << 1;
            while (newCapacity < capacity) {
                newCapacity <<= 1;
            }
            data = Arrays.copyOf(data, newCapacity);
        }
    }

    private static final class LongNodeMap {
        private static final float LOAD_FACTOR = 0.75f;
        private long[] keys;
        private Node[] values;
        private byte[] states; // 0 = empty, 1 = occupied, 2 = deleted
        private int mask;
        private int size;
        private int threshold;

        LongNodeMap() {
            initialise(16);
        }

        Node get(long key) {
            int idx = findIndex(key);
            return idx >= 0 ? values[idx] : null;
        }

        void put(long key, Node value) {
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

        void clear() {
            Arrays.fill(states, (byte) 0);
            Arrays.fill(values, null);
            size = 0;
        }

        void forEach(Consumer<Node> consumer) {
            for (int i = 0; i < states.length; i++) {
                if (states[i] == 1) {
                    consumer.accept(values[i]);
                }
            }
        }

        private int findIndex(long key) {
            int idx = mix(key) & mask;
            while (states[idx] != 0) {
                if (states[idx] == 1 && keys[idx] == key) {
                    return idx;
                }
                idx = (idx + 1) & mask;
            }
            return -1;
        }

        private int findSlot(long key) {
            int idx = mix(key) & mask;
            while (states[idx] == 1 && keys[idx] != key) {
                idx = (idx + 1) & mask;
            }
            return idx;
        }

        private void rehash(int newCapacity) {
            long[] oldKeys = keys;
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
            keys = new long[pow2];
            values = new Node[pow2];
            states = new byte[pow2];
            mask = pow2 - 1;
            threshold = (int) (pow2 * LOAD_FACTOR);
            size = 0;
        }

        private int mix(long key) {
            long z = key;
            z = (z ^ (z >>> 33)) * 0xff51afd7ed558ccdL;
            z = (z ^ (z >>> 33)) * 0xc4ceb9fe1a85ec53L;
            z = z ^ (z >>> 33);
            return (int) z;
        }
    }
}
