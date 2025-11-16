package PMIndex;

import org.openjdk.jol.info.GraphLayout;
import org.openjdk.jol.vm.VM;
import search.Pattern;
import utilities.PatternResult;
import utilities.StringKeyMapper;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Objects;
import java.util.function.Consumer;

/**
 * Suffix tree that stores hashed n-gram tokens as longs. Incoming strings are mapped to longs via
 * {@link StringKeyMapper} so the representation matches the numeric token stream used by HBI.
 * Ukkonen's algorithm keeps the tree in sync as tokens arrive, allowing pattern queries to execute
 * in linear time w.r.t. number of tokens in the query.
 *
 * Memory optimizations included:
 *  - Lazy children: Node.children is null until the first child is inserted (big leaf savings)
 *  - Four-inline-slot ChildMap before upgrading to a hash table
 *  - Integer edge ends with a single global leafEndValue (no per-node End objects)
 *  - Optional compaction pass after build
 *  - JOL memory reports (full and partitioned)
 */
public class SuffixTreeIndex implements IPMIndexing {

    private final StringKeyMapper keyMapper;
    private final LongSequence tokenStream;
    private final Node root;

    // Ukkonen active point
    private Node activeNode;
    private int activeEdge = -1;
    private int activeLength = 0;

    // Ukkonen bookkeeping
    private int remainingSuffixCount = 0;
    private int leafEndValue = -1;                // global, shared end for all leaves
    private Node lastCreatedInternalNode = null;

    public SuffixTreeIndex() {
        this(new StringKeyMapper(1 << 21, 0.0001));
    }

    public SuffixTreeIndex(long expectedDistinctTokens, double epsilon) {
        this(new StringKeyMapper(expectedDistinctTokens, epsilon));
    }

    /** Construct with a provided mapper and a dynamically-resizable token stream. */
    public SuffixTreeIndex(StringKeyMapper mapper) {
        this.keyMapper = Objects.requireNonNull(mapper, "mapper");
        this.tokenStream = new LongSequence();
        this.root = new Node(-1, -1);
        this.root.suffixLink = this.root;
        this.activeNode = this.root;
    }

    /** Fixed-capacity token stream, preallocating space for the full text to index. */
    public SuffixTreeIndex(long expectedDistinctTokens, double epsilon, int totalTokens) {
        this(new StringKeyMapper(expectedDistinctTokens, epsilon), totalTokens);
    }

    /** Construct with a provided mapper and a fixed-capacity token stream. */
    public SuffixTreeIndex(StringKeyMapper mapper, int totalTokens) {
        this.keyMapper = Objects.requireNonNull(mapper, "mapper");
        if (totalTokens <= 0) {
            throw new IllegalArgumentException("totalTokens must be positive");
        }
        this.tokenStream = new LongSequence(totalTokens, true);
        this.root = new Node(-1, -1);
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
        // advance the global leaf end
        leafEndValue = pos;
        remainingSuffixCount++;
        lastCreatedInternalNode = null;

        while (remainingSuffixCount > 0) {
            if (activeLength == 0) {
                activeEdge = pos;
            }
            long currentEdgeToken = tokenStream.get(activeEdge);
            Node next = getChild(activeNode, currentEdgeToken);

            if (next == null) {
                // create a new leaf with LEAF_SENTINEL end
                Node leaf = new Node(pos, Node.LEAF_SENTINEL);
                leaf.suffixIndex = pos - remainingSuffixCount + 1;
                putChild(activeNode, currentEdgeToken, leaf);

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

                // split the edge
                int splitEndVal = next.start + activeLength - 1;
                Node split = new Node(next.start, splitEndVal);
                putChild(activeNode, currentEdgeToken, split);

                // new leaf from current position
                Node leaf = new Node(pos, Node.LEAF_SENTINEL);
                leaf.suffixIndex = pos - remainingSuffixCount + 1;
                putChild(split, currentToken, leaf);

                // advance 'next' to start after the split
                next.start = next.start + activeLength;
                long nextEdgeToken = tokenStream.get(next.start);
                putChild(split, nextEdgeToken, next);

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

    @Override
    public boolean exists(String key) {
        // Not implemented in this index; report() is the primary query op.
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

    /**
     * Aligns with Java regex semantics for empty patterns:
     * reports matches at all N+1 boundaries when the pattern has zero tokens.
     * Keeps exact path-length accounting when we finish mid-edge.
     */
    private ArrayList<Integer> reportInternal(long[] patternTokens) {
        ArrayList<Integer> result = new ArrayList<>();

        // empty pattern: report every boundary including the one after the last token
        if (patternTokens.length == 0) {
            int n = tokenStream.size();
            for (int i = 0; i <= n; i++) {
                result.add(i);
            }
            return result;
        }

        Node current = root;
        int patternIndex = 0;
        int pathLength = 0;

        while (patternIndex < patternTokens.length) {
            long searchToken = patternTokens[patternIndex];
            Node next = getChild(current, searchToken);
            if (next == null) {
                return result;
            }

            int edgeStart = next.start;
            int edgeEnd   = next.effectiveEnd(leafEndValue);
            int edgeIndex = edgeStart;

            // walk along this edge while tokens match
            while (patternIndex < patternTokens.length && edgeIndex <= edgeEnd) {
                if (patternTokens[patternIndex] != tokenStream.get(edgeIndex)) {
                    return new ArrayList<>();
                }
                patternIndex++;
                edgeIndex++;
            }

            // if we consumed the entire pattern, collect answers under 'next'
            if (patternIndex == patternTokens.length) {
                int consumedOnThisEdge = edgeIndex - edgeStart;
                int nextPathLength = pathLength + consumedOnThisEdge;

                collectSuffixIndices(next, nextPathLength, result);
                return result;
            }

            // otherwise, if we consumed the whole edge, descend to child
            if (edgeIndex > edgeEnd) {
                int consumedOnThisEdge = edgeEnd - edgeStart + 1;
                pathLength += consumedOnThisEdge;
                current = next;
            } else {
                return new ArrayList<>();
            }
        }
        return result;
    }

    private void collectSuffixIndices(Node node, int pathLength, ArrayList<Integer> result) {
        if (node == null) return;

        ArrayDeque<Node> stack = new ArrayDeque<>();
        ArrayDeque<Integer> pathLengths = new ArrayDeque<>();
        stack.push(node);
        pathLengths.push(pathLength);

        while (!stack.isEmpty()) {
            Node current = stack.pop();
            int currentPathLength = pathLengths.pop();

            if (!hasChildren(current)) {
                // For leaves, the suffix start is the leaf's start position. Prefer recorded suffixIndex if present.
                int start = (current.suffixIndex >= 0) ? current.suffixIndex : current.start;
                if (start >= 0) {
                    current.suffixIndex = start;
                    result.add(start);
                }
                continue;
            }

            final int nextPathLengthBase = currentPathLength;
            forEachChild(current, child -> {
                stack.push(child);
                pathLengths.push(nextPathLengthBase + child.edgeLength(leafEndValue));
            });
        }
    }

    @Override
    public void expire() {
        tokenStream.clear();
        if (root.children != null) {
            root.children.clear();
            root.children = null;
        }
        root.start = -1;
        root.end = -1;
        root.suffixLink = root;
        activeNode = root;
        activeEdge = -1;
        activeLength = 0;
        remainingSuffixCount = 0;
        leafEndValue = -1;
        lastCreatedInternalNode = null;
    }

    /**
     * Optional compaction step to reduce memory after all insertions are complete.
     * Trims all child maps to minimal capacity (or drops if leaf), clears construction-only pointers,
     * and shrinks the token stream backing array to exact size.
     * Invoke only when you will not insert more tokens.
     */
    public void compactForQuerying() {
        tokenStream.shrinkToFit();
        trimDfs(root);
        lastCreatedInternalNode = null;
    }

    private void trimDfs(Node node) {
        if (node == null) return;
        if (!hasChildren(node)) {
            node.children = null; // ensure leaves hold no ChildMap at all
        } else {
            node.children.trimToSize();
            forEachChild(node, this::trimDfs);
        }
        node.dropConstructionOnlyPointers();
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

    // =========================
    // === JOL memory reports ==
    // =========================

    /**
     * Produce a memory footprint report for the suffix tree index using JOL.
     *
     * @param includeFootprintTable when true, append the full class histogram for the index root.
     * @return human-readable report string in mebibytes
     */
    public String jolMemoryReport(boolean includeFootprintTable) {
        StringBuilder sb = new StringBuilder(4_096);

        sb.append("=== JOL / VM details ===\n");
        sb.append(VM.current().details()).append('\n');

        GraphLayout totalLayout = GraphLayout.parseInstance(this);
        sb.append("\n=== SuffixTree total (index as root) ===\n");
        appendLayout(sb, "Total", totalLayout.totalSize());

        GraphLayout nodesLayout = GraphLayout.parseInstance(root);
        appendLayout(sb, "Node graph (root)", nodesLayout.totalSize());

        GraphLayout tokenLayout = GraphLayout.parseInstance(tokenStream);
        appendLayout(sb, "Token stream", tokenLayout.totalSize());

        GraphLayout mapperLayout = GraphLayout.parseInstance(keyMapper);
        appendLayout(sb, "StringKeyMapper", mapperLayout.totalSize());

        long remainder = totalLayout.totalSize()
                - nodesLayout.totalSize()
                - tokenLayout.totalSize()
                - mapperLayout.totalSize();
        if (remainder < 0) {
            remainder = 0;
        }
        appendLayout(sb, "Other fields", remainder);

        if (includeFootprintTable) {
            sb.append("\n--- Class footprint (SuffixTree root) ---\n");
            sb.append(totalLayout.toFootprint()).append('\n');
        }

        return sb.toString();
    }

    /** Same report as {@link #jolMemoryReport(boolean)} plus the numeric total MiB. */
    public utilities.MemoryUsageReport jolMemoryReportWithTotal(boolean includeFootprintTable) {
        String txt = jolMemoryReport(includeFootprintTable);
        long totalBytes = GraphLayout.parseInstance(this).totalSize();
        double totalMiB = totalBytes / (1024.0 * 1024.0);
        return new utilities.MemoryUsageReport(txt, totalMiB);
    }

    private static void appendLayout(StringBuilder sb, String label, long bytes) {
        sb.append(String.format(Locale.ROOT,
                "%s: %d B (%.3f MiB)%n",
                label,
                bytes,
                bytes / (1024.0 * 1024.0)));
    }

    /**
     * Compact, partitioned report to compare directly with HBI's partitioned JOL report.
     * Returns exactly four lines: Total, Nodes, TokenStream, Mapper.
     */
    public String jolMemoryReportPartitioned() {
        long totalBytes  = GraphLayout.parseInstance(this).totalSize();
        long nodesBytes  = GraphLayout.parseInstance(root).totalSize();
        long streamBytes = GraphLayout.parseInstance(tokenStream).totalSize();
        long mapperBytes = GraphLayout.parseInstance(keyMapper).totalSize();

        Locale L = Locale.ROOT;
        String totalLine  = String.format(L, "Total: %d B (%.3f MiB)",  totalBytes,  totalBytes  / (1024.0 * 1024.0));
        String nodesLine  = String.format(L, "Nodes: %d B (%.3f MiB)",  nodesBytes,  nodesBytes  / (1024.0 * 1024.0));
        String streamLine = String.format(L, "TokenStream: %d B (%.3f MiB)", streamBytes, streamBytes / (1024.0 * 1024.0));
        String mapLine    = String.format(L, "Mapper: %d B (%.3f MiB)", mapperBytes, mapperBytes / (1024.0 * 1024.0));

        return totalLine + "\n" + nodesLine + "\n" + streamLine + "\n" + mapLine;
    }

    /** Same as {@link #jolMemoryReportPartitioned()} but also returns the numeric total MiB. */
    public utilities.MemoryUsageReport jolMemoryReportPartitionedWithTotal() {
        String txt = jolMemoryReportPartitioned();
        long totalBytes = GraphLayout.parseInstance(this).totalSize();
        double totalMiB = totalBytes / (1024.0 * 1024.0);
        return new utilities.MemoryUsageReport(txt, totalMiB);
    }

    // =====================
    // === Debug helpers ===
    // =====================

    /**
     * Debug: Walk the tree for the given query, printing traversed edges and the leaf starts
     * under the matched node. Limits leaf printing to the tailWindow of the stream for focus.
     */
    public void debugDumpForQuery(String query, int nGram, int tailWindow) {
        long[] patternTokens = mapQueryToTokens(query, nGram);
        System.out.printf("DEBUG[SuffixTree]: query=\"%s\" tokens=%d tailWindow=%d%n",
                query, patternTokens.length, tailWindow);

        Node current = root;
        int i = 0;
        while (i < patternTokens.length) {
            long token = patternTokens[i];
            Node next = getChild(current, token);
            if (next == null) {
                System.out.printf("DEBUG: no edge for token[%d], abort.%n", i);
                return;
            }

            int s = next.start;
            int e = next.effectiveEnd(leafEndValue);
            int elen = next.edgeLength(leafEndValue);
            System.out.printf("DEBUG: edge start=%d end=%d len=%d%n", s, e, elen);

            int idx = s;
            while (i < patternTokens.length && idx <= e) {
                long t = tokenStream.get(idx);
                if (t != patternTokens[i]) {
                    System.out.printf("DEBUG: mismatch at token[%d] streamIdx=%d%n", i, idx);
                    return;
                }
                i++; idx++;
            }

            if (i == patternTokens.length) {
                System.out.println("DEBUG: pattern fully consumed; collecting leaves under matched node...");
                debugPrintLeaves(next, tailWindow);
                return;
            }

            if (idx > e) {
                current = next; // descend
            } else {
                System.out.println("DEBUG: stopped mid-edge unexpectedly.");
                return;
            }
        }
    }

    private long[] mapQueryToTokens(String query, int nGram) {
        search.Pattern p = new search.Pattern(query, nGram);
        long[] tokens = new long[p.nGramArr.length];
        for (int i = 0; i < p.nGramArr.length; i++) {
            tokens[i] = keyMapper.mapToLong(p.nGramArr[i]);
        }
        return tokens;
    }

    private void debugPrintLeaves(Node node, int tailWindow) {
        int n = tokenStream.size();
        int tailCutoff = Math.max(0, n - Math.max(0, tailWindow));
        int[] counter = new int[1];
        System.out.printf("DEBUG: streamSize=%d tailCutoff=%d%n", n, tailCutoff);
        debugPrintLeavesDfs(node, tailCutoff, counter);
        System.out.printf("DEBUG: totalLeavesVisited=%d%n", counter[0]);
    }

    private void debugPrintLeavesDfs(Node node, int tailCutoff, int[] counter) {
        if (!hasChildren(node)) {
            int start = (node.suffixIndex >= 0) ? node.suffixIndex : node.start;
            int end = node.effectiveEnd(leafEndValue);
            if (start >= tailCutoff) {
                System.out.printf("DEBUG: leaf start=%d end=%d%n", start, end);
            }
            counter[0]++;
            return;
        }
        forEachChild(node, child -> debugPrintLeavesDfs(child, tailCutoff, counter));
    }

    /** Quick check: does the token stream match the query tokens at a given start? */
    public boolean debugCheckAt(String query, int nGram, int start) {
        long[] tokens = mapQueryToTokens(query, nGram);
        int n = tokenStream.size();
        if (start < 0 || start + tokens.length > n) {
            System.out.printf("DEBUG: start %d out of range [0,%d) for %d tokens%n", start, n, tokens.length);
            return false;
        }
        for (int j = 0; j < tokens.length; j++) {
            long t = tokenStream.get(start + j);
            if (t != tokens[j]) {
                System.out.printf("DEBUG: mismatch at start=%d, j=%d (stream=%d, queryToken=%d)\n",
                        start, j, t, tokens[j]);
                return false;
            }
        }
        System.out.printf("DEBUG: exact token match at start=%d for query=\"%s\" (nGram=%d)\n", start, query, nGram);
        return true;
    }

    /** Scan a range [from, to] and print all token-exact matches for the query. */
    public void debugScanRange(String query, int nGram, int from, int to) {
        long[] tokens = mapQueryToTokens(query, nGram);
        int n = tokenStream.size();
        int m = tokens.length;
        from = Math.max(0, from);
        to = Math.min(n - m, to);
        System.out.printf("DEBUG: scanRange query=\"%s\" nGram=%d from=%d to=%d (n=%d, m=%d)\n",
                query, nGram, from, to, n, m);
        for (int s = from; s <= to; s++) {
            boolean ok = true;
            for (int j = 0; j < m; j++) {
                if (tokenStream.get(s + j) != tokens[j]) { ok = false; break; }
            }
            if (ok) {
                System.out.printf("DEBUG: scan match at start=%d\n", s);
            }
        }
    }

    // =========================
    // ==== Compact node type ==
    // =========================
    private static final class Node {
        // LAZY: null until first child is inserted
        private ChildMap children;
        private Node suffixLink;
        private int start;
        private int end;                 // explicit end or LEAF_SENTINEL
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

        /** Release construction-only pointer to reduce memory after build. */
        private void dropConstructionOnlyPointers() {
            this.suffixLink = null;
        }
    }

    // =========================
    // === Child helpers =======
    // =========================
    private static Node getChild(Node n, long key) {
        return (n.children == null) ? null : n.children.get(key);
    }

    private static void putChild(Node n, long key, Node child) {
        if (n.children == null) n.children = new ChildMap();
        n.children.put(key, child);
    }

    private static boolean hasChildren(Node n) {
        return n.children != null && !n.children.isEmpty();
    }

    private static void forEachChild(Node n, Consumer<Node> c) {
        if (n.children != null) n.children.forEach(c);
    }

    // =========================
    // === Token storage =======
    // =========================
    private static final class LongSequence {
        private long[] data;
        private int size = 0;
        private final boolean fixedCapacity;

        LongSequence() {
            this(16, false);
        }

        LongSequence(int capacity) {
            this(capacity, false);
        }

        LongSequence(int capacity, boolean fixed) {
            if (capacity <= 0) {
                throw new IllegalArgumentException("capacity must be positive");
            }
            this.data = new long[capacity];
            this.fixedCapacity = fixed;
        }

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

        /** Shrink the backing array to the exact logical size. Call after all insertions. */
        void shrinkToFit() {
            if (fixedCapacity) return; // preserve single-allocation guarantee
            if (data.length != size) {
                data = Arrays.copyOf(data, size);
            }
        }

        private void ensureCapacity(int capacity) {
            if (capacity <= data.length) return;
            if (fixedCapacity) {
                throw new IllegalStateException("LongSequence capacity exceeded: requested " + capacity +
                        " > allocated " + data.length);
            }
            int newCapacity = data.length << 1;
            while (newCapacity < capacity) newCapacity <<= 1;
            data = Arrays.copyOf(data, newCapacity);
        }
    }

    // ===========================================================
    // === Compact child container with 4-inline fast path     ===
    // ===========================================================
    private static final class ChildMap {
        // mode 0 empty, 1..4 number of inline entries, 5 upgraded hash table
        private byte mode = 0;

        private long k1, k2, k3, k4;
        private Node v1, v2, v3, v4;

        private LongNodeMap map;

        boolean isEmpty() {
            return mode == 0;
        }

        Node get(long key) {
            switch (mode) {
                case 0:  return null;
                case 1:  return k1 == key ? v1 : null;
                case 2:  return (k1 == key ? v1 : (k2 == key ? v2 : null));
                case 3:  return (k1 == key ? v1 : (k2 == key ? v2 : (k3 == key ? v3 : null)));
                case 4:  return (k1 == key ? v1 : (k2 == key ? v2 : (k3 == key ? v3 : (k4 == key ? v4 : null))));
                default: return map.get(key);
            }
        }

        void put(long key, Node value) {
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
                    // upgrade
                    map = new LongNodeMap(8);
                    map.put(k1, v1); map.put(k2, v2); map.put(k3, v3); map.put(k4, v4);
                    // free inline slots
                    v1 = v2 = v3 = v4 = null;
                    mode = 5;
                    map.put(key, value);
                    return;
                default:
                    map.put(key, value);
            }
        }

        void clear() {
            mode = 0;
            v1 = v2 = v3 = v4 = null;
            map = null;
        }

        void trimToSize() {
            if (mode == 5 && map != null) {
                map.trimToSize();
            }
        }

        void forEach(Consumer<Node> consumer) {
            switch (mode) {
                case 0: return;
                case 1: consumer.accept(v1); return;
                case 2: consumer.accept(v1); consumer.accept(v2); return;
                case 3: consumer.accept(v1); consumer.accept(v2); consumer.accept(v3); return;
                case 4: consumer.accept(v1); consumer.accept(v2); consumer.accept(v3); consumer.accept(v4); return;
                default:
                    map.forEach(consumer);
            }
        }
    }

    // ============================================
    // === Hash map for the upgraded child table ===
    // ============================================
    private static final class LongNodeMap {
        private static final float LOAD_FACTOR = 0.75f;

        private long[] keys;
        private Node[] values;
        private byte[] states; // 0 empty, 1 occupied, 2 deleted
        private int mask;
        private int size;
        private int threshold;

        LongNodeMap() {
            initialise(16);
        }

        LongNodeMap(int initialCapacity) {
            initialise(Math.max(4, initialCapacity));
        }

        boolean isEmpty() {
            return size == 0;
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

        /** Reduce arrays to the smallest power-of-two capacity that satisfies the load factor. */
        void trimToSize() {
            if (size == 0) {
                initialise(4);
                return;
            }
            int minNeeded = (int) Math.ceil(size / (double) LOAD_FACTOR);
            int target = 1;
            while (target < minNeeded) target <<= 1;
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
