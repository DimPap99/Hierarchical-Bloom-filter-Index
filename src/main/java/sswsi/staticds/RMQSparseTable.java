package sswsi.staticds;

/**
 * Sparse table for range-maximum queries returning the index of the maximum.
 * Preprocessing O(n log n), queries O(1), space O(n log n).
 */
final class RMQSparseTable {
    private final int n;
    private final int[] log2;
    private final int[][] st; // indices into 'arr' of the max in the block
    private final int[] arr;

    RMQSparseTable(int[] arr) {
        this.arr = arr;
        this.n = arr.length;
        this.log2 = new int[n + 1];
        for (int i = 2; i <= n; i++) {
            log2[i] = log2[i >> 1] + 1;
        }
        int K = (n == 0) ? 0 : (31 - Integer.numberOfLeadingZeros(n));
        this.st = new int[K + 1][n];
        for (int i = 0; i < n; i++) {
            st[0][i] = i;
        }
        for (int k = 1; k <= K; k++) {
            int len = 1 << k;
            int half = len >> 1;
            for (int i = 0; i + len <= n; i++) {
                int leftIdx = st[k - 1][i];
                int rightIdx = st[k - 1][i + half];
                st[k][i] = (arr[leftIdx] >= arr[rightIdx]) ? leftIdx : rightIdx;
            }
        }
    }

    // Return index of max in [l, r]
    int argmax(int l, int r) {
        if (l > r) return -1;
        int len = r - l + 1;
        int k = log2[len];
        int leftIdx = st[k][l];
        int rightIdx = st[k][r - (1 << k) + 1];
        return (arr[leftIdx] >= arr[rightIdx]) ? leftIdx : rightIdx;
    }
}

