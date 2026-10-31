package tree.ssws;

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.Arrays;

/**
 * Suffix tree of {@code text[offset, offset + n)} followed by a unique terminator, stored in arrays and immutable
 * after construction. Built by a {@link SuffixTreeBuilder}.
 *
 * <p>The tree is a view: it keeps a reference to {@code text} and reads {@code text[offset, offset + n)} whenever it
 * is searched, so those entries must not change while the tree is in use. The rest of {@code text} may change;
 * {@code SlidingSuffixIndex} keeps every tree in one shared array.
 *
 * <p><b>Nodes.</b> Leaf {@code i} is suffix {@code sa[i]}, the i-th suffix in suffix-array order, and has node id
 * {@code ~i}, which is negative. Internal nodes have ids {@code 0 .. internalNodeCount() - 1}, and the root is 0.
 * An internal node stores its string depth and its leaf range {@code [lo, hi]} in suffix-array order. Children are
 * stored in suffix-array order, which is the order of the ranks of their first symbols (names follow first
 * occurrence, see {@link RankReducer}), so a depth-first traversal visits the leaves in suffix-array order.
 *
 * <p><b>No stored edge labels.</b> The edge into child {@code c} of {@code p} is
 * {@code text[sa[lo(c)] + depth(p), sa[lo(c)] + depth(c))}, where position {@code n} is the terminator. A leaf's
 * depth counts the terminator: the leaf of suffix {@code s} has depth {@code n + 1 - s}. The terminator is not a
 * symbol: it never matches a pattern symbol, and the child dictionary leaves out edges that start with it.
 *
 * <p>Nothing here recurses on tree depth: a run of one symbol gives a tree as deep as it is long.
 */
public final class StaticSuffixTree {

    /** No such node. Leaves are {@code ~i} for {@code 0 <= i <= n}, so this is never a node id. */
    public static final int NONE = Integer.MIN_VALUE;

    final long[] text;
    final int offset;
    final int n;
    final int[] sa;              // n + 1 suffixes in sorted order; sa[0] = n, the terminator alone
    final int internalCount;
    final int[] depth;           // per internal node: string depth
    final int[] lo;              // per internal node: first leaf, in suffix-array order
    final int[] hi;              // per internal node: last leaf
    final int[] childStart;      // internalCount + 1 offsets into children
    final int[] children;        // node ids; each parent's children in suffix-array order
    final ChildDictionary dictionary;

    private StaticSuffixTree(long[] text, int offset, int n, int[] sa, int internalCount, int[] depth, int[] lo,
                             int[] hi, int[] childStart, int[] children, ChildDictionary.Kind dictionary) {
        this.text = text;
        this.offset = offset;
        this.n = n;
        this.sa = sa;
        this.internalCount = internalCount;
        this.depth = depth;
        this.lo = lo;
        this.hi = hi;
        this.childStart = childStart;
        this.children = children;
        this.dictionary = dictionary.build(this);   // last: it reads the arrays above
    }

    /** Number of symbols, without the terminator. */
    public int length() {
        return n;
    }

    /** {@code length() + 1}: one leaf per suffix, including the terminator on its own. */
    public int leafCount() {
        return n + 1;
    }

    public int internalNodeCount() {
        return internalCount;
    }

    public static boolean isLeaf(int node) {
        return node < 0;
    }

    /** String depth; a leaf's depth counts the terminator. */
    public int depth(int node) {
        return (node >= 0) ? depth[node] : n + 1 - sa[~node];
    }

    /** First leaf of the node's subtree, in suffix-array order. */
    public int lo(int node) {
        return (node >= 0) ? lo[node] : ~node;
    }

    /** Last leaf of the node's subtree, in suffix-array order. */
    public int hi(int node) {
        return (node >= 0) ? hi[node] : ~node;
    }

    /** Start of the suffix at position {@code rank} of the suffix array. */
    public int suffix(int rank) {
        return sa[rank];
    }

    public int childCount(int node) {
        return (node >= 0) ? childStart[node + 1] - childStart[node] : 0;
    }

    /** The {@code k}-th child of an internal node, in suffix-array order. */
    public int child(int node, int k) {
        return children[childStart[node] + k];
    }

    /** The child whose edge starts with {@code symbol}, through the child dictionary, or {@link #NONE}. */
    public int childBySymbol(int node, long symbol) {
        return dictionary.child(node, symbol);
    }

    /** True if this is the tree of {@code text[offset, offset + length)}: the same array, not a copy of it. */
    public boolean views(long[] text, int offset, int length) {
        return this.text == text && this.offset == offset && this.n == length;
    }

    /** The symbol at {@code position < length()}. */
    public long symbolAt(int position) {
        if (position < 0 || position >= n) {
            throw new IndexOutOfBoundsException(position);
        }
        return text[offset + position];
    }

    /** Position of the first symbol of the edge from internal node {@code parent} into its child. */
    int firstSymbolPosition(int parent, int child) {
        return sa[lo(child)] + depth[parent];
    }

    /** True if the edge from {@code parent} into {@code child} starts with the terminator. */
    boolean isTerminatorEdge(int parent, int child) {
        return firstSymbolPosition(parent, child) == n;
    }

    /** First symbol of the edge from {@code parent} into {@code child}; not for a terminator edge. */
    long firstSymbol(int parent, int child) {
        return text[offset + firstSymbolPosition(parent, child)];
    }

    /**
     * The locus of {@code pattern}: the highest node whose path label starts with the pattern, so that its leaf
     * range holds exactly the pattern's occurrences. {@link #NONE} if the pattern does not occur. The pattern must
     * not be empty.
     */
    public int locus(long[] pattern) {
        final int m = pattern.length;
        final long[] text = this.text;
        int node = 0;
        int matched = 0;          // = depth(node)
        while (true) {
            int c = dictionary.child(node, pattern[matched]);
            if (c == NONE) {
                return NONE;
            }
            int start;
            int end;
            if (c < 0) {
                start = sa[~c];
                if (m > n - start) {
                    return NONE;  // the pattern runs into the terminator
                }
                end = m;
            } else {
                start = sa[lo[c]];
                end = Math.min(depth[c], m);
            }
            // The dictionary matched the edge's first symbol; compare the rest of it.
            int base = offset + start;
            for (int k = matched + 1; k < end; k++) {
                if (text[base + k] != pattern[k]) {
                    return NONE;
                }
            }
            if (end == m) {
                return c;
            }
            node = c;
            matched = end;
        }
    }

    /** Appends {@code base + s} for the start {@code s} of every suffix below {@code node}. */
    public void reportAll(int node, IntArrayList out, int base) {
        for (int i = lo(node), last = hi(node); i <= last; i++) {
            out.add(base + sa[i]);
        }
    }

    /** Appends {@code base + s} for every suffix start {@code s >= minStart} below {@code node}. */
    public void reportFrom(int node, int minStart, IntArrayList out, int base) {
        for (int i = lo(node), last = hi(node); i <= last; i++) {
            int s = sa[i];
            if (s >= minStart) {
                out.add(base + s);
            }
        }
    }

    /**
     * Builds the tree from the suffix array of {@code text[offset, offset + n)} plus terminator ({@code n + 1}
     * entries) and its LCP array ({@code lcp[i]} = LCP of the suffixes {@code sa[i]} and {@code sa[i + 1]}).
     *
     * <p>Stack algorithm: the stack holds the rightmost path. For each suffix in sorted order, pop the nodes deeper
     * than its LCP with the previous suffix, attaching each to the node below it; if no node has exactly that
     * depth, create one between them. Each node is attached to its parent when it is popped, so every parent
     * receives its children in suffix-array order.
     */
    static StaticSuffixTree fromSuffixArray(long[] text, int offset, int n, int[] sa, int[] lcp,
                                            ChildDictionary.Kind dictionary) {
        final int leaves = n + 1;
        BuildWorkspace ws = WORKSPACE.get();
        ws.ensure(leaves);
        final int[] depth = ws.depth;
        final int[] lo = ws.lo;
        final int[] hi = ws.hi;
        final int[] stack = ws.stack;
        final int[] stackDepth = ws.stackDepth;
        final int[] edgeParent = ws.edgeParent;
        final int[] edgeChild = ws.edgeChild;

        int internal = 1;                    // node 0, the root: depth 0, leaves [0, n]
        depth[0] = 0;
        lo[0] = 0;
        int top = 0;
        stack[0] = 0;
        stackDepth[0] = 0;
        int edges = 0;
        for (int i = 0; i < leaves; i++) {
            int l = (i == 0) ? 0 : lcp[i - 1];
            int last = NONE;
            while (stackDepth[top] > l) {
                last = stack[top--];
                if (last >= 0) {
                    hi[last] = i - 1;
                }
                if (stackDepth[top] >= l) {
                    edgeParent[edges] = stack[top];
                    edgeChild[edges++] = last;
                }
            }
            if (stackDepth[top] < l) {
                // The previous leaf is always popped (the terminator is unique), so last is set here.
                int v = internal++;
                depth[v] = l;
                lo[v] = (last >= 0) ? lo[last] : ~last;
                edgeParent[edges] = v;
                edgeChild[edges++] = last;
                stack[++top] = v;
                stackDepth[top] = l;
            }
            stack[++top] = ~i;
            stackDepth[top] = leaves - sa[i];
        }
        while (top > 0) {
            int x = stack[top--];
            if (x >= 0) {
                hi[x] = leaves - 1;
            }
            edgeParent[edges] = stack[top];
            edgeChild[edges++] = x;
        }
        hi[0] = leaves - 1;

        // Group the edges by parent; a stable counting sort keeps each parent's children in order.
        int[] childStart = new int[internal + 1];
        for (int e = 0; e < edges; e++) {
            childStart[edgeParent[e] + 1]++;
        }
        for (int v = 0; v < internal; v++) {
            childStart[v + 1] += childStart[v];
        }
        int[] cursor = ws.cursor;
        System.arraycopy(childStart, 0, cursor, 0, internal);
        int[] children = new int[edges];
        for (int e = 0; e < edges; e++) {
            children[cursor[edgeParent[e]]++] = edgeChild[e];
        }
        return new StaticSuffixTree(text, offset, n, sa, internal, Arrays.copyOf(depth, internal),
                Arrays.copyOf(lo, internal), Arrays.copyOf(hi, internal), childStart, children, dictionary);
    }

    private static final ThreadLocal<BuildWorkspace> WORKSPACE = ThreadLocal.withInitial(BuildWorkspace::new);

    /** Scratch arrays reused across builds on one thread; the tree keeps exact-size copies. */
    private static final class BuildWorkspace {
        int[] depth = new int[0];
        int[] lo = new int[0];
        int[] hi = new int[0];
        int[] cursor = new int[0];
        int[] stack = new int[1];
        int[] stackDepth = new int[1];
        int[] edgeParent = new int[0];
        int[] edgeChild = new int[0];

        void ensure(int leaves) {
            if (depth.length >= leaves) {
                return;
            }
            int cap = Math.max(leaves, depth.length + (depth.length >> 1));
            depth = new int[cap];                // internal nodes: at most max(1, n) < leaves + 1
            lo = new int[cap];
            hi = new int[cap];
            cursor = new int[cap];
            stack = new int[cap + 1];            // root, internal nodes and one leaf on the rightmost path
            stackDepth = new int[cap + 1];
            edgeParent = new int[2 * cap];       // internal - 1 + leaves edges
            edgeChild = new int[2 * cap];
        }
    }
}
