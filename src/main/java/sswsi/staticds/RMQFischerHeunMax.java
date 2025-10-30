package sswsi.staticds;

import it.unimi.dsi.fastutil.ints.IntArrayList;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

/**
 * Fischer–Heun RMQ for maximum on a static integer array.
 * Preprocessing O(n), queries O(1), space O(n).
 *
 * Implementation via max-Cartesian tree + Euler tour on depths (±1 property),
 * and block decomposition with microstructure tables and a sparse table on block minima.
 */
final class RMQFischerHeunMax {
    private final int n;
    private final int[] arr;

    // Cartesian tree (max)
    private final int[] left;
    private final int[] right;
    private final int[] parent;
    private final int root;

    // Euler tour arrays
    private final int m;              // 2n-1
    private final int[] E;            // nodes in tour
    private final int[] L;            // depth at each tour step
    private final int[] first;        // first occurrence of node in E

    // Block decomposition on L
    private final int b;              // block size
    private final int B;              // number of blocks
    private final int[] blockType;    // type code per block
    private final int[] blockMinPos;  // absolute pos in L of min depth in block
    private final MinSparseTable blockMinSt; // over depths at blockMinPos
    private final Map<Long, int[]> microTables = new HashMap<>(); // key=(type<<6)|len

    RMQFischerHeunMax(int[] arr) {
        this.n = arr.length;
        this.arr = Arrays.copyOf(arr, n);
        if (n == 0) {
            this.left = this.right = this.parent = new int[0];
            this.root = -1;
            this.m = 0; this.E = this.L = this.first = new int[0];
            this.b = 1; this.B = 0;
            this.blockType = this.blockMinPos = new int[0];
            this.blockMinSt = new MinSparseTable(new int[0]);
            return;
        }

        // 1) Build max-cartesian tree in O(n)
        this.left = new int[n];
        this.right = new int[n];
        this.parent = new int[n];
        Arrays.fill(left, -1);
        Arrays.fill(right, -1);
        Arrays.fill(parent, -1);

        IntArrayList stack = new IntArrayList();
        for (int i = 0; i < n; i++) {
            int last = -1;
            while (!stack.isEmpty() && arr[stack.getInt(stack.size() - 1)] <= arr[i]) {
                last = stack.removeInt(stack.size() - 1);
            }
            if (!stack.isEmpty()) {
                int p = stack.getInt(stack.size() - 1);
                right[p] = i;
                parent[i] = p;
            }
            if (last != -1) {
                left[i] = last;
                parent[last] = i;
            }
            stack.add(i);
        }
        int r = stack.getInt(stack.size() - 1);
        while (parent[r] != -1) r = parent[r];
        this.root = r;

        // 2) Euler tour and depths
        this.m = 2 * n - 1;
        this.E = new int[m];
        this.L = new int[m];
        this.first = new int[n];
        Arrays.fill(first, -1);

        int idx = 0;
        // iterative DFS to avoid recursion limits
        class Frame { int u, state, depth; Frame(int u, int s, int d){this.u=u;this.state=s;this.depth=d;} }
        java.util.ArrayDeque<Frame> st = new java.util.ArrayDeque<>();
        st.push(new Frame(root, 0, 0));
        while (!st.isEmpty()) {
            Frame f = st.peek();
            if (f.state == 0) {
                // pre
                E[idx] = f.u; L[idx] = f.depth; if (first[f.u] == -1) first[f.u] = idx; idx++;
                f.state = 1;
                if (left[f.u] != -1) {
                    st.push(new Frame(left[f.u], 0, f.depth + 1));
                }
            } else if (f.state == 1) {
                // between children
                E[idx] = f.u; L[idx] = f.depth; idx++;
                f.state = 2;
                if (right[f.u] != -1) {
                    st.push(new Frame(right[f.u], 0, f.depth + 1));
                }
            } else {
                // exit
                E[idx] = f.u; L[idx] = f.depth; idx++;
                st.pop();
            }
        }
        // idx should be m; guard if tree had missing children
        if (idx != m) {
            // trim
            // but leave arrays as-is for safety
        }

        // 3) Fischer–Heun over L (min RMQ on depths)
        int mm = m;
        int log2 = (mm <= 1) ? 1 : (31 - Integer.numberOfLeadingZeros(mm));
        int bs = Math.max(1, log2 / 2);
        this.b = bs;
        this.B = (mm + b - 1) / b;
        this.blockType = new int[B];
        this.blockMinPos = new int[B];

        int[] blockMinDepth = new int[Math.max(B,1)];
        for (int bi = 0; bi < B; bi++) {
            int start = bi * b;
            int end = Math.min(mm - 1, start + b - 1);
            int len = end - start + 1;
            int code = buildBlockType(L, start, len);
            blockType[bi] = codeFor(code, len);
            // precompute micro-table for this (type,len) if absent
            long key = key(code, len);
            if (!microTables.containsKey(key)) {
                int[] tbl = buildMicroTable(L, start, len);
                microTables.put(key, tbl);
            }
            // block min
            int minPos = start;
            for (int i2 = start + 1; i2 <= end; i2++) {
                if (L[i2] < L[minPos]) minPos = i2;
            }
            blockMinPos[bi] = minPos;
            blockMinDepth[bi] = L[minPos];
        }
        this.blockMinSt = new MinSparseTable(blockMinDepth);
    }

    // Return index of maximum in arr[l..r]
    int argmax(int l, int r) {
        if (l > r) { int t = l; l = r; r = t; }
        int i = first[l];
        int j = first[r];
        if (i > j) { int t = i; i = j; j = t; }
        int pos = rmqDepthPos(i, j);
        return E[pos]; // node id equals position in arr
    }

    // RMQ on L depths, return absolute position in L of minimal depth
    private int rmqDepthPos(int i, int j) {
        int bi = i / b, bj = j / b;
        if (bi == bj) {
            return bi * b + inBlockArgMin(bi, i % b, j % b);
        }
        int leftPos = bi * b + inBlockArgMin(bi, i % b, b - 1);
        int rightPos = bj * b + inBlockArgMin(bj, 0, j % b);
        if (bi + 1 > bj - 1) {
            return (L[leftPos] <= L[rightPos]) ? leftPos : rightPos;
        }
        int midBlock = blockMinSt.argminBlock(bi + 1, bj - 1);
        int midPos = blockMinPos[midBlock];
        int best = leftPos;
        if (L[rightPos] < L[best]) best = rightPos;
        if (L[midPos] < L[best]) best = midPos;
        return best;
    }

    private int inBlockArgMin(int blockIndex, int offL, int offR) {
        int start = blockIndex * b;
        int len = Math.min(b, m - start);
        int code = decodeType(blockType[blockIndex]);
        long key = key(code, len);
        int[] tbl = microTables.get(key);
        int idx = offL * len + offR;
        return tbl[idx]; // returns offset within block
    }

    private static long key(int code, int len) {
        return (((long) code) << 6) ^ len;
    }
    private static int codeFor(int code, int len) { return (code << 6) ^ len; }
    private static int decodeType(int packed) { return (packed >>> 6); }

    private int buildBlockType(int[] L, int start, int len) {
        int code = 0;
        for (int k = 0; k < len - 1; k++) {
            int d = L[start + k + 1] - L[start + k];
            int bit = (d > 0) ? 1 : 0; // +1 => 1, -1 => 0
            code |= (bit << k);
        }
        return code;
    }

    // Build in-block RMQ table for positions [0..len-1] using brute force over L
    private int[] buildMicroTable(int[] L, int start, int len) {
        int[] tbl = new int[len * len];
        for (int i = 0; i < len; i++) {
            int minPos = i;
            for (int j = i; j < len; j++) {
                if (L[start + j] < L[start + minPos]) minPos = j;
                tbl[i * len + j] = minPos -  i + i; // store absolute within block
            }
        }
        return tbl;
    }

    // Min sparse table over block minima depths
    private static final class MinSparseTable {
        private final int n;
        private final int K;
        private final int[][] st; // store argmin block index
        private final int[] depth;
        private final int[] lg;

        MinSparseTable(int[] depth) {
            this.depth = depth;
            this.n = depth.length;
            if (n == 0) { K = 0; st = new int[1][0]; lg = new int[1]; return; }
            this.lg = new int[n + 1];
            for (int i = 2; i <= n; i++) lg[i] = lg[i >> 1] + 1;
            this.K = lg[n];
            this.st = new int[K + 1][n];
            for (int i = 0; i < n; i++) st[0][i] = i;
            for (int k = 1; k <= K; k++) {
                int len = 1 << k;
                int half = len >> 1;
                for (int i = 0; i + len <= n; i++) {
                    int a = st[k - 1][i];
                    int b = st[k - 1][i + half];
                    st[k][i] = (depth[a] <= depth[b]) ? a : b;
                }
            }
        }

        int argminBlock(int l, int r) {
            if (l > r) return -1;
            int len = r - l + 1;
            int k = lg[len];
            int a = st[k][l];
            int b = st[k][r - (1 << k) + 1];
            return (depth[a] <= depth[b]) ? a : b;
        }
    }
}

