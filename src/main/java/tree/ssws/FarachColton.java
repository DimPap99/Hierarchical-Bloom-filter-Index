package tree.ssws;

import java.util.Arrays;

/**
 * Suffix tree construction of Farach (FOCS 1997) and Farach-Colton, Ferragina and Muthukrishnan (J. ACM 47(6),
 * 2000), linear for an integer alphabet {@code [1, n]}; Bille et al. build their suffix trees with it ([8] there).
 *
 * <p>Each level returns its tree as a suffix array plus LCP array: the depth-first leaf order of the tree and the
 * string depth of the lowest common ancestor of each pair of neighbouring leaves, which determine a suffix tree
 * with sorted children. The level above rebuilds the tries it needs from them, and the top level hands them to
 * {@link StaticSuffixTree#fromSuffixArray}.
 *
 * <p>The papers count from 1: the "odd" suffixes start at the 0-based even positions 0, 2, 4, ... and the
 * "even" suffixes at the 0-based odd positions. One level, for {@code s} of length {@code n}:
 * <ol>
 *   <li><b>Odd tree.</b> Rank the pairs {@code (s[2i], s[2i+1])} of {@code s$} (padded with one 0) into
 *       {@code s'} of length {@code n/2 + 1} and recurse. Suffix {@code i} of {@code s'} is the odd suffix
 *       {@code 2i}. Two odd suffixes that agree on {@code l'} pairs share {@code 2l'} symbols, or {@code 2l' + 1}
 *       if their next pairs start with the same symbol.</li>
 *   <li><b>Even tree.</b> The even suffix {@code j} is {@code s[j]} followed by the odd suffix {@code j + 1}, so
 *       one bucket pass over the odd order sorts the even suffixes. Neighbours with the same first symbol share
 *       1 + the LCP of their odd suffixes, a range minimum over the odd LCP array.</li>
 *   <li><b>Over-merge.</b> Walk both tries at once and merge edges that start with the same symbol over the
 *       length of the shorter one, without looking at the rest: some paths get merged too far.</li>
 *   <li><b>True depths.</b> For a merged node {@code v}, take an odd leaf {@code x} and an even leaf {@code y}
 *       whose lowest common ancestor is {@code v}. Every such pair has the same LCP, {@code L(v)}, and it is
 *       1 + the LCP of {@code x + 1} and {@code y + 1}, which is {@code L} of their lowest common ancestor
 *       {@code d(v)}: {@code L(v) = 1 + L(d(v))}, {@code L(root) = 0}.</li>
 *   <li><b>Unmerge.</b> Top down, a merged node with {@code L(v)} equal to its merged depth is a true node. If
 *       not, the two paths part at depth {@code L(v)}, strictly inside both merged edges, and the node there gets
 *       the original odd and even subtrees below that point, in the order of their next symbol. The result is
 *       written out as suffix array plus LCP array.</li>
 * </ol>
 * Only the levels recurse (log n deep); every walk over a tree uses an explicit stack, because a run of one
 * symbol gives trees as deep as they are long.
 */
final class FarachColton {

    /** Inputs up to this length are sorted directly. */
    static final int BASE = 32;

    private static final byte MERGED = 0;
    private static final byte REF_ODD = 1;
    private static final byte REF_EVEN = 2;

    private FarachColton() {
    }

    static void build(int[] s, int n, int sigma, int[] sa, int[] lcp) {
        build(s, n, sigma, sa, lcp, BASE);
    }

    /**
     * Suffix array and LCP array of {@code s[0, n)} followed by the terminator.
     *
     * @param s     symbols {@code s[0, n)} in {@code [1, sigma]}, then {@code s[n] = s[n + 1] = 0}
     * @param sa    receives the {@code n + 1} suffixes in sorted order; {@code sa[0] = n}
     * @param lcp   receives {@code n} values: {@code lcp[i]} is the LCP of suffixes {@code sa[i]} and
     *              {@code sa[i + 1]}
     * @param base  inputs of length at most {@code max(base, 2)} are sorted directly (tests use 2, so that every
     *              step runs on small inputs)
     */
    static void build(int[] s, int n, int sigma, int[] sa, int[] lcp, int base) {
        if (n <= Math.max(2, base)) {
            sortDirectly(s, n, sa, lcp);
        } else {
            new Level(s, n, sigma, base).run(sa, lcp);
        }
    }

    /** Insertion sort of the suffixes; the unique terminator ends every comparison. */
    static void sortDirectly(int[] s, int n, int[] sa, int[] lcp) {
        for (int i = 0; i <= n; i++) {
            int j = i - 1;
            while (j >= 0 && less(s, i, sa[j])) {
                sa[j + 1] = sa[j];
                j--;
            }
            sa[j + 1] = i;
        }
        for (int i = 0; i < n; i++) {
            int a = sa[i];
            int b = sa[i + 1];
            int h = 0;
            while (s[a + h] == s[b + h]) {
                h++;
            }
            lcp[i] = h;
        }
    }

    private static boolean less(int[] s, int a, int b) {
        while (s[a] == s[b]) {
            a++;
            b++;
        }
        return s[a] < s[b];
    }

    /**
     * A compacted trie of some of the suffixes, rebuilt from their sorted order and LCPs with the stack algorithm
     * of {@link StaticSuffixTree#fromSuffixArray}. Leaf {@code k} is suffix {@code sa[k]} and has id {@code ~k};
     * internal nodes have ids from 0 (the root, which may have a single child here).
     */
    private static final class Trie {
        final int[] sa;
        final int[] lcp;
        final int leaves;
        final int n;
        final int internal;
        final int[] depth;
        final int[] lo;
        final int[] hi;
        final int[] childStart;
        final int[] children;

        Trie(int[] sa, int[] lcp, int leaves, int n) {
            this.sa = sa;
            this.lcp = lcp;
            this.leaves = leaves;
            this.n = n;
            int cap = leaves + 1;
            int[] depth = new int[cap];
            int[] lo = new int[cap];
            int[] hi = new int[cap];
            int[] stack = new int[cap + 1];
            int[] stackDepth = new int[cap + 1];
            int[] edgeParent = new int[2 * cap];
            int[] edgeChild = new int[2 * cap];
            int internal = 1;
            int top = 0;
            int edges = 0;
            for (int i = 0; i < leaves; i++) {
                int l = (i == 0) ? 0 : lcp[i - 1];
                int last = StaticSuffixTree.NONE;
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
                    int v = internal++;
                    depth[v] = l;
                    lo[v] = (last >= 0) ? lo[last] : ~last;
                    edgeParent[edges] = v;
                    edgeChild[edges++] = last;
                    stack[++top] = v;
                    stackDepth[top] = l;
                }
                stack[++top] = ~i;
                stackDepth[top] = n + 1 - sa[i];
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
            int[] childStart = new int[internal + 1];
            for (int e = 0; e < edges; e++) {
                childStart[edgeParent[e] + 1]++;
            }
            for (int v = 0; v < internal; v++) {
                childStart[v + 1] += childStart[v];
            }
            int[] cursor = Arrays.copyOf(childStart, internal);
            int[] children = new int[edges];
            for (int e = 0; e < edges; e++) {
                children[cursor[edgeParent[e]]++] = edgeChild[e];
            }
            this.internal = internal;
            this.depth = depth;
            this.lo = lo;
            this.hi = hi;
            this.childStart = childStart;
            this.children = children;
        }

        int depthOf(int c) {
            return (c >= 0) ? depth[c] : n + 1 - sa[~c];
        }

        int loOf(int c) {
            return (c >= 0) ? lo[c] : ~c;
        }

        int hiOf(int c) {
            return (c >= 0) ? hi[c] : ~c;
        }

        /** A suffix in the subtree of {@code c}: its path label starts with the path to {@code c}. */
        int someSuffix(int c) {
            return sa[loOf(c)];
        }

        /** True if {@code c} is a leaf whose whole suffix, terminator included, ends at depth {@code d}. */
        boolean exhaustedLeaf(int c, int d) {
            return c < 0 && depthOf(c) == d;
        }
    }

    /** One level of the recursion. */
    private static final class Level {
        private final int[] s;
        private final int n;
        private final int sigma;
        private final int base;

        // Merged nodes (T_m). Node 0 is the root. For each: the string depth after over-merging, the odd and the
        // even node at or just below the merged point, the parent, and a range of entries: its children, in the
        // order of their first symbol. Children have larger ids than their parent.
        private int merged;
        private int[] mDepth;
        private int[] mOdd;
        private int[] mEven;
        private int[] mParent;
        private int[] mFirst;
        private int[] mEnd;
        // An entry is a merged child, or an original odd or even subtree hanging below the merged node as it is.
        private byte[] entKind;
        private int[] entId;
        private int[] trueDepth;   // L(v)
        private int out;           // suffixes emitted so far

        Level(int[] s, int n, int sigma, int base) {
            this.s = s;
            this.n = n;
            this.sigma = sigma;
            this.base = base;
        }

        void run(int[] sa, int[] lcp) {
            // Step 1: the odd suffixes, from the suffix array of the pairs.
            int np = n / 2 + 1;
            int[] pairs = new int[np + 2];
            int sigmaPairs = rankPairs(pairs, np);
            int[] saPairs = new int[np + 1];
            int[] lcpPairs = new int[np];
            build(pairs, np, sigmaPairs, saPairs, lcpPairs, base);
            int[] saOdd = new int[np];
            int[] lcpOdd = new int[np - 1];
            for (int k = 0; k < np; k++) {
                saOdd[k] = 2 * saPairs[k + 1];          // saPairs[0] is the terminator of s'
            }
            for (int k = 0; k + 1 < np; k++) {
                int l = lcpPairs[k + 1];
                lcpOdd[k] = 2 * l + ((s[saOdd[k] + 2 * l] == s[saOdd[k + 1] + 2 * l]) ? 1 : 0);
            }
            int[] oddRank = new int[n + 1];
            for (int k = 0; k < np; k++) {
                oddRank[saOdd[k]] = k;
            }

            // Step 2: the even suffixes, by first symbol over the odd order of j + 1.
            int ne = (n + 1) / 2;
            int[] saEven = new int[ne];
            int[] next = new int[sigma + 1];
            for (int j = 1; j <= n; j += 2) {
                next[s[j]]++;
            }
            for (int c = 0, sum = 0; c <= sigma; c++) {
                int t = next[c];
                next[c] = sum;
                sum += t;
            }
            if ((n & 1) == 1) {
                saEven[next[0]++] = n;                   // the terminator alone, first
            }
            for (int k = 0; k < np; k++) {
                int p = saOdd[k];
                if (p > 0) {
                    saEven[next[s[p - 1]]++] = p - 1;
                }
            }
            int[] lcpEven = new int[ne - 1];
            RangeMin oddLcp = new RangeMin(lcpOdd, np - 1);
            for (int k = 0; k + 1 < ne; k++) {
                int j1 = saEven[k];
                int j2 = saEven[k + 1];
                lcpEven[k] = (s[j1] != s[j2]) ? 0 : 1 + oddLcp.min(oddRank[j1 + 1], oddRank[j2 + 1] - 1);
            }

            // Steps 3 to 5.
            Trie odd = new Trie(saOdd, lcpOdd, np, n);
            Trie even = new Trie(saEven, lcpEven, ne, n);
            overMerge(odd, even);
            trueDepths(odd, even);
            emit(odd, even, sa, lcp);
            if (out != n + 1) {
                throw new IllegalStateException("emitted " + out + " of " + (n + 1) + " suffixes");
            }
        }

        /** Dense ranks, from 1, of the pairs (s[2i], s[2i+1]), i < np, by two counting-sort passes. */
        private int rankPairs(int[] pairs, int np) {
            int[] count = new int[sigma + 2];
            int[] bySecond = new int[np];
            int[] order = new int[np];
            for (int i = 0; i < np; i++) {
                count[s[2 * i + 1] + 1]++;
            }
            for (int c = 1; c <= sigma + 1; c++) {
                count[c] += count[c - 1];
            }
            for (int i = 0; i < np; i++) {
                bySecond[count[s[2 * i + 1]]++] = i;
            }
            Arrays.fill(count, 0);
            for (int i = 0; i < np; i++) {
                count[s[2 * i] + 1]++;
            }
            for (int c = 1; c <= sigma + 1; c++) {
                count[c] += count[c - 1];
            }
            for (int k = 0; k < np; k++) {
                int i = bySecond[k];
                order[count[s[2 * i]]++] = i;
            }
            int rank = 0;
            int prevA = -1;
            int prevB = -1;
            for (int k = 0; k < np; k++) {
                int i = order[k];
                int a = s[2 * i];
                int b = s[2 * i + 1];
                if (a != prevA || b != prevB) {
                    rank++;
                    prevA = a;
                    prevB = b;
                }
                pairs[i] = rank;
            }
            return rank;   // pairs[np] = pairs[np + 1] = 0 from the allocation
        }

        /** Step 3: merges the two tries, comparing only the first symbol of each edge. */
        private void overMerge(Trie odd, Trie even) {
            int nodes = odd.internal + odd.leaves + even.internal + even.leaves;
            int cap = nodes + 1;
            mDepth = new int[cap];
            mOdd = new int[cap];
            mEven = new int[cap];
            mParent = new int[cap];
            mFirst = new int[cap];
            mEnd = new int[cap];
            entKind = new byte[cap + nodes];
            entId = new int[cap + nodes];
            int[] todo = new int[cap];
            int count = 1;                // node 0: the roots of both tries, depth 0
            mParent[0] = -1;
            int sp = 0;
            todo[sp++] = 0;
            int e = 0;
            while (sp > 0) {
                int w = todo[--sp];
                int d = mDepth[w];
                int co = mOdd[w];
                int ce = mEven[w];
                // The branches below the merged point: the node's children, or the one edge the point lies on.
                boolean oddAtNode = odd.depthOf(co) == d;
                boolean evenAtNode = even.depthOf(ce) == d;
                int oBase = (oddAtNode && co >= 0) ? odd.childStart[co] : 0;
                int oCount = !oddAtNode ? 1 : (co >= 0) ? odd.childStart[co + 1] - oBase : 0;
                int eBase = (evenAtNode && ce >= 0) ? even.childStart[ce] : 0;
                int eCount = !evenAtNode ? 1 : (ce >= 0) ? even.childStart[ce + 1] - eBase : 0;
                int bo = (oCount == 0) ? 0 : oddAtNode ? odd.children[oBase] : co;
                int be = (eCount == 0) ? 0 : evenAtNode ? even.children[eBase] : ce;
                int fo = (oCount == 0) ? Integer.MAX_VALUE : s[odd.someSuffix(bo) + d];
                int fe = (eCount == 0) ? Integer.MAX_VALUE : s[even.someSuffix(be) + d];
                mFirst[w] = e;
                int i = 0;
                int j = 0;
                while (i < oCount || j < eCount) {
                    boolean takeOdd = fo <= fe;
                    boolean takeEven = fe <= fo;
                    if (takeOdd && takeEven) {
                        int v = count++;
                        mDepth[v] = Math.min(odd.depthOf(bo), even.depthOf(be));
                        mOdd[v] = bo;
                        mEven[v] = be;
                        mParent[v] = w;
                        entKind[e] = MERGED;
                        entId[e++] = v;
                        todo[sp++] = v;
                    } else if (takeOdd) {
                        entKind[e] = REF_ODD;
                        entId[e++] = bo;
                    } else {
                        entKind[e] = REF_EVEN;
                        entId[e++] = be;
                    }
                    if (takeOdd) {
                        if (++i < oCount) {
                            bo = odd.children[oBase + i];
                            fo = s[odd.someSuffix(bo) + d];
                        } else {
                            fo = Integer.MAX_VALUE;
                        }
                    }
                    if (takeEven) {
                        if (++j < eCount) {
                            be = even.children[eBase + j];
                            fe = s[even.someSuffix(be) + d];
                        } else {
                            fe = Integer.MAX_VALUE;
                        }
                    }
                }
                mEnd[w] = e;
            }
            merged = count;
        }

        /** Step 4: L(v) for every merged node, through the d-links. */
        private void trueDepths(Trie odd, Trie even) {
            int m = merged;
            // attOdd[k], attEven[k]: the merged node that the k-th odd or even leaf hangs from in T_m (or that it
            // is the exhausted leaf of). A referenced subtree is a range of leaves, so these are range fills.
            int[] att = new int[n + 1];
            for (int w = 0; w < m; w++) {
                for (int q = mFirst[w]; q < mEnd[w]; q++) {
                    if (entKind[q] != MERGED) {
                        Trie t = (entKind[q] == REF_ODD) ? odd : even;
                        int c = entId[q];
                        for (int k = t.loOf(c), hi = t.hiOf(c); k <= hi; k++) {
                            att[t.sa[k]] = w;
                        }
                    }
                }
                if (odd.exhaustedLeaf(mOdd[w], mDepth[w])) {
                    att[odd.sa[~mOdd[w]]] = w;
                }
                if (even.exhaustedLeaf(mEven[w], mDepth[w])) {
                    att[even.sa[~mEven[w]]] = w;
                }
            }

            // An odd leaf x and an even leaf y with LCA(x, y) = w, bottom up. repOdd/repEven: any leaf below w.
            int[] repOdd = new int[m];
            int[] repEven = new int[m];
            int[] pairX = new int[m];
            int[] pairY = new int[m];
            for (int w = m - 1; w >= 0; w--) {
                int selfOdd = odd.exhaustedLeaf(mOdd[w], mDepth[w]) ? odd.sa[~mOdd[w]] : -1;
                int selfEven = even.exhaustedLeaf(mEven[w], mDepth[w]) ? even.sa[~mEven[w]] : -1;
                int ro = selfOdd;
                int re = selfEven;
                int firstOddEntry = -1;
                for (int q = mFirst[w]; q < mEnd[w]; q++) {
                    int o = entryOdd(q, odd, repOdd);
                    int ev = entryEven(q, even, repEven);
                    if (ro < 0) {
                        ro = o;
                    }
                    if (re < 0) {
                        re = ev;
                    }
                    if (firstOddEntry < 0 && o >= 0) {
                        firstOddEntry = q;
                    }
                }
                repOdd[w] = ro;
                repEven[w] = re;
                if (w == 0) {
                    continue;
                }
                int x;
                int y;
                if (selfOdd >= 0) {
                    x = selfOdd;
                    y = re;
                } else if (selfEven >= 0) {
                    x = ro;
                    y = selfEven;
                } else {
                    x = entryOdd(firstOddEntry, odd, repOdd);
                    y = -1;
                    for (int q = mFirst[w]; q < mEnd[w] && y < 0; q++) {
                        if (q != firstOddEntry) {
                            y = entryEven(q, even, repEven);
                        }
                    }
                    if (y < 0) {
                        // Every even leaf is below the first odd entry: take the odd leaf from another entry.
                        y = entryEven(firstOddEntry, even, repEven);
                        x = -1;
                        for (int q = mFirst[w]; q < mEnd[w] && x < 0; q++) {
                            if (q != firstOddEntry) {
                                x = entryOdd(q, odd, repOdd);
                            }
                        }
                    }
                }
                if (x < 0 || y < 0) {
                    throw new IllegalStateException("merged node " + w + " has no odd/even leaf pair");
                }
                pairX[w] = x;
                pairY[w] = y;
            }

            // LCA in T_m by preorder and range minimum: for pre(u) < pre(v), the node with the fewest edges to the
            // root in (pre(u), pre(v)] is a child of LCA(u, v). (Not the smallest string depth: the range also
            // holds the later descendants of u, which can be shallower in symbols than the child towards v.)
            int[] pre = new int[m];
            int[] byPre = new int[m];
            int[] key = new int[m];
            int[] edgesToRoot = new int[m];
            int[] stack = new int[m];
            int sp = 0;
            stack[sp++] = 0;
            int c = 0;
            while (sp > 0) {
                int w = stack[--sp];
                pre[w] = c;
                byPre[c] = w;
                key[c] = edgesToRoot[w];
                c++;
                for (int q = mEnd[w] - 1; q >= mFirst[w]; q--) {
                    if (entKind[q] == MERGED) {
                        edgesToRoot[entId[q]] = edgesToRoot[w] + 1;
                        stack[sp++] = entId[q];
                    }
                }
            }
            RangeMin shallowest = new RangeMin(key, m);
            int[] link = new int[m];
            for (int w = 1; w < m; w++) {
                int u = att[pairX[w] + 1];
                int v = att[pairY[w] + 1];
                if (u == v) {
                    link[w] = u;
                } else {
                    int a = pre[u];
                    int b = pre[v];
                    link[w] = mParent[byPre[shallowest.argMin(Math.min(a, b) + 1, Math.max(a, b))]];
                }
            }

            // L(v) = 1 + L(d(v)): follow the links to a known value, then fill the chain back in.
            trueDepth = new int[m];
            Arrays.fill(trueDepth, -1);
            trueDepth[0] = 0;
            int[] chain = stack;
            for (int w = 1; w < m; w++) {
                int len = 0;
                int v = w;
                while (trueDepth[v] < 0) {
                    if (len == m) {
                        throw new IllegalStateException("d-links do not reach the root");
                    }
                    chain[len++] = v;
                    v = link[v];
                }
                int depth = trueDepth[v];
                while (len > 0) {
                    trueDepth[chain[--len]] = ++depth;
                }
            }
        }

        private int entryOdd(int q, Trie odd, int[] repOdd) {
            return switch (entKind[q]) {
                case REF_ODD -> odd.someSuffix(entId[q]);
                case REF_EVEN -> -1;
                default -> repOdd[entId[q]];
            };
        }

        private int entryEven(int q, Trie even, int[] repEven) {
            return switch (entKind[q]) {
                case REF_EVEN -> even.someSuffix(entId[q]);
                case REF_ODD -> -1;
                default -> repEven[entId[q]];
            };
        }

        /** Step 5: walks T_m top down, unmerges where L(v) is short of the merged depth, writes SA and LCP. */
        private void emit(Trie odd, Trie even, int[] sa, int[] lcp) {
            int m = merged;
            int[] frameNode = new int[m];
            int[] frameNext = new int[m];
            int[] frameLcp = new int[m];   // LCP between the node's first leaf and the leaf before it
            int top = 0;
            frameNode[0] = 0;
            frameNext[0] = mFirst[0];
            top = 1;
            out = 0;
            while (top > 0) {
                int f = top - 1;
                int w = frameNode[f];
                int q = frameNext[f];
                if (q == mEnd[w]) {
                    top--;
                    continue;
                }
                frameNext[f] = q + 1;
                int before = (q == mFirst[w]) ? frameLcp[f] : mDepth[w];
                int id = entId[q];
                if (entKind[q] == REF_ODD) {
                    copy(odd, id, before, sa, lcp);
                } else if (entKind[q] == REF_EVEN) {
                    copy(even, id, before, sa, lcp);
                } else if (trueDepth[id] == mDepth[id]) {
                    frameNode[top] = id;
                    frameNext[top] = mFirst[id];
                    frameLcp[top] = before;
                    top++;
                } else {
                    // Merged too far: the odd and even paths part at depth l, inside both merged edges.
                    int l = trueDepth[id];
                    int co = mOdd[id];
                    int ce = mEven[id];
                    if (s[odd.someSuffix(co) + l] < s[even.someSuffix(ce) + l]) {
                        copy(odd, co, before, sa, lcp);
                        copy(even, ce, l, sa, lcp);
                    } else {
                        copy(even, ce, before, sa, lcp);
                        copy(odd, co, l, sa, lcp);
                    }
                }
            }
        }

        /** Appends the leaves below {@code c}, with their LCPs; the first one shares {@code before} with its predecessor. */
        private void copy(Trie t, int c, int before, int[] sa, int[] lcp) {
            if (out > 0) {
                lcp[out - 1] = before;
            }
            for (int k = t.loOf(c), hi = t.hiOf(c); k <= hi; k++) {
                sa[out] = t.sa[k];
                if (k < hi) {
                    lcp[out] = t.lcp[k];
                }
                out++;
            }
        }
    }
}
