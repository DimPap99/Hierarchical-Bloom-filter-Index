package tree.ssws;

import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;

import static tree.ssws.StaticSuffixTree.NONE;

/**
 * {@link ChildDictionary.Kind#FASTUTIL}: up to two children per node are compared inline, as in the legacy tree;
 * a node with more children gets its own {@code Long2IntOpenHashMap}.
 */
final class FastutilChildDictionary implements ChildDictionary {

    private final long[] keys;                 // two inline slots per node: first symbols
    private final int[] kids;                  // two inline slots per node: children, NONE if empty
    private final Long2IntOpenHashMap[] maps;  // per node with more than two children; null if there is none

    FastutilChildDictionary(StaticSuffixTree tree) {
        int nodes = tree.internalCount;
        keys = new long[2 * nodes];
        kids = new int[2 * nodes];
        Long2IntOpenHashMap[] big = null;
        for (int p = 0; p < nodes; p++) {
            int from = tree.childStart[p];
            int to = tree.childStart[p + 1];
            // The terminator sorts first, so only the first child can start with it.
            if (from < to && tree.isTerminatorEdge(p, tree.children[from])) {
                from++;
            }
            int count = to - from;
            if (count <= 2) {
                for (int k = 0; k < 2; k++) {
                    int c = (k < count) ? tree.children[from + k] : NONE;
                    kids[2 * p + k] = c;
                    keys[2 * p + k] = (c != NONE) ? tree.firstSymbol(p, c) : 0L;
                }
            } else {
                if (big == null) {
                    big = new Long2IntOpenHashMap[nodes];
                }
                Long2IntOpenHashMap map = new Long2IntOpenHashMap(count);
                map.defaultReturnValue(NONE);
                for (int q = from; q < to; q++) {
                    map.put(tree.firstSymbol(p, tree.children[q]), tree.children[q]);
                }
                big[p] = map;
                kids[2 * p] = NONE;
                kids[2 * p + 1] = NONE;
            }
        }
        maps = big;
    }

    @Override
    public int child(int node, long symbol) {
        if (maps != null) {
            Long2IntOpenHashMap map = maps[node];
            if (map != null) {
                return map.get(symbol);
            }
        }
        int s = node << 1;
        // An empty slot holds NONE, so returning it on a key match is still correct.
        if (keys[s] == symbol) {
            return kids[s];
        }
        if (keys[s + 1] == symbol) {
            return kids[s + 1];
        }
        return NONE;
    }
}
