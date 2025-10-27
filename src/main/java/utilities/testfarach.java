package utilities;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;
import java.util.function.Consumer;

/**
 * Farach–Colton suffix-tree construction in O(n) time for integer alphabets.
 * Call FarachColton.build(String) or FarachColton.build(int[]) to obtain the root node.
 */
public final class testfarach {
    public static Map<Character, Integer> map = new HashMap<>();

    private testfarach() {
    }

    public static Node build(String input) {
        int[] remapped = remapToIntegers(input);
        return build(remapped);
    }

    public static Node build(int[] alphabet) {
        int[] copy = Arrays.copyOf(alphabet, alphabet.length);
        return constructSuffixTree(copy);
    }

    static Node constructSuffixTree(int[] input) {
        int[] S = appendUniqueChar(input);
        int originalLength = S.length - 1;

        if (originalLength == 0) {
            Node root = new Node(0, 0);
            root.addChild(new Node(1, 1));
            root.updateLeafList();
            return root;
        }
        if (originalLength == 1) {
            Node root = new Node(0, 0);
            root.addChild(new Node(2, 1));
            root.addChild(new Node(1, 2));
            root.updateLeafList();
            return root;
        }

        Node odd = buildOddTree(S);
        Node even = buildEvenTree(odd, S);
        Node merged = overmerge(even, odd, S);
        computeLcpTree(merged);
        adjustOvermerge(merged, S);
        cleanup(merged);
        return merged;
    }

    

    private static Node buildOddTree(int[] S) {
        int n = S.length;
        int pairCount = n / 2;
        int[][] pairs = new int[pairCount][2];
        for (int i = 0; i < pairCount; i++) {
            pairs[i][0] = S[2 * i];
            pairs[i][1] = S[2 * i + 1];
        }
        int[][] sortedPairs =  countingSortPairs(pairs);

        Map<Long, Integer> rankMap = new HashMap<>();
        int rank = 1;
        long lastKey = Long.MIN_VALUE;
        for (int[] p : sortedPairs) {
            long key = packPair(p[0], p[1]);
            if (key != lastKey) {
                rankMap.put(key, rank++);
                lastKey = key;
            }
        }

        int[] reduced = new int[pairCount];
        for (int i = 0; i < pairCount; i++) {
            reduced[i] = rankMap.get(packPair(S[2 * i], S[2 * i + 1]));
        }

        Node tree = constructSuffixTree(reduced);
        tree.updateLeafList();
        extendOddTree(tree, n);
        tree.updateLeafList();
        repairOddTree(tree, S);
        tree.updateLeafList();
        return tree;
    }

    private static void extendOddTree(Node root, int totalLength) {
        root.dfs(node -> {
            if (node.isLeaf()) {
                node.suffixIndex = node.suffixIndex * 2 - 1;
                node.depth = totalLength - node.suffixIndex + 1;
            } else {
                node.depth *= 2;
            }
        });
    }

    private static void repairOddTree(Node root, int[] S) {
        Deque<Node> stack = new ArrayDeque<>();
        stack.push(root);
        while (!stack.isEmpty()) {
            Node node = stack.pop();
            if (node.isLeaf()) {
                continue;
            }
            List<Node> children = node.children;
            List<Node> newChildren = new ArrayList<>();
            List<Node> bucket = new ArrayList<>();
            int currentChar = Integer.MIN_VALUE;

            for (Node child : children) {
                if (child.depth == node.depth) {
                    continue;
                }
                Node sample = child.isLeaf() ? child : child.leafList.get(0);
                int pos = sample.suffixIndex + node.depth - 1;
                if (pos < 0 || pos >= S.length) {
                    continue;
                }
                int firstChar = S[pos];
                if (currentChar == Integer.MIN_VALUE || firstChar == currentChar) {
                    bucket.add(child);
                    currentChar = firstChar;
                } else {
                    flushBucket(node, newChildren, bucket);
                    bucket.clear();
                    bucket.add(child);
                    currentChar = firstChar;
                }
            }
            flushBucket(node, newChildren, bucket);
            node.children.clear();
            newChildren.forEach(node::addChild);

            if (node.children.size() == 1 && node.suffixIndex != 0) {
                Node only = node.children.get(0);
                node.depth = only.depth;
                node.suffixIndex = only.suffixIndex;
                node.children = only.children;
                node.children.forEach(c -> c.parent = node);
            }
            for (Node child : node.children) {
                stack.push(child);
            }
        }
    }

    private static void flushBucket(Node parent, List<Node> output, List<Node> bucket) {
        if (bucket.isEmpty()) {
            return;
        }
        if (bucket.size() == 1) {
            output.add(bucket.get(0));
        } else {
            Node inner = new Node(parent.depth + 1, -1);
            bucket.forEach(inner::addChild);
            output.add(inner);
        }
    }

    private static Node buildEvenTree(Node oddRoot, int[] S) {
        List<Node> leaves = new ArrayList<>();
        oddRoot.dfs(node -> {
            if (node.isLeaf()) {
                leaves.add(node);
            }
        });

        List<Integer> sortedOdd = new ArrayList<>();
        for (Node leaf : leaves) {
            sortedOdd.add(leaf.suffixIndex);
        }

        List<int[]> evenWithPrefix = new ArrayList<>();
        int n = S.length - 1;
        for (int id : sortedOdd) {
            if (id <= 1) {
                continue;
            }
            int pos = id - 2;
            if (pos < 0 || pos >= S.length) {
                continue;
            }
            evenWithPrefix.add(new int[]{S[pos], id});
        }

        evenWithPrefix.sort(Comparator.comparingInt(a -> a[0]));
        List<Integer> evenSuffixes = new ArrayList<>();
        for (int[] tuple : evenWithPrefix) {
            evenSuffixes.add(tuple[1] - 1);
        }
        if ((n & 1) == 0) {
            evenSuffixes.add(n);
        }

        Map<Integer, Node> idToNode = new HashMap<>();
        oddRoot.dfs(node -> {
            if (node.isLeaf()) {
                idToNode.put(node.suffixIndex, node);
            }
        });

        Map<Pair, Integer> lcp = new HashMap<>();
        if (!evenSuffixes.isEmpty()) {
            List<Pair> queries = new ArrayList<>();
            for (int i = 0; i < evenSuffixes.size() - 1; i++) {
                int a = evenSuffixes.get(i);
                int b = evenSuffixes.get(i + 1);
                int aheadA = a - 1;
                int aheadB = b - 1;
                boolean charsEqual = aheadA >= 0 && aheadA < S.length
                        && aheadB >= 0 && aheadB < S.length
                        && S[aheadA] == S[aheadB];
                if (charsEqual && a < n && b < n
                        && idToNode.containsKey(a + 1) && idToNode.containsKey(b + 1)) {
                    queries.add(new Pair(a + 1, b + 1));
                } else {
                    int lcpVal = charsEqual ? 1 : 0;
                    lcp.put(new Pair(a, b), lcpVal);
                }
            }
            if (!queries.isEmpty()) {
                Map<Pair, Node> lcaMap = tarjanLca(oddRoot, idToNode, queries);
                for (Pair q : queries) {
                    int evenA = q.first - 1;
                    int evenB = q.second - 1;
                    int value = lcaMap.containsKey(q)
                            ? lcaMap.get(q).depth + 1
                            : 0;
                    lcp.put(new Pair(evenA, evenB), value);
                }
            }
        }

        Node root = new Node(0, 0);
        Deque<Node> stack = new ArrayDeque<>();
        stack.push(root);
        Integer prev = null;
        for (int idx = 0; idx < evenSuffixes.size(); idx++) {
            int suffix = evenSuffixes.get(idx);
            int currLcp = (prev == null) ? 0 : lcp.getOrDefault(new Pair(prev, suffix), 0);

            while (!stack.isEmpty() && stack.peek().depth > currLcp) {
                stack.pop();
            }
            Node parent = stack.peek();
            if (parent.depth < currLcp) {
                Node inner = new Node(currLcp, -1);
                if (!parent.children.isEmpty()) {
                    Node last = parent.children.remove(parent.children.size() - 1);
                    parent.addChild(inner);
                    inner.addChild(last);
                } else {
                    parent.addChild(inner);
                }
                stack.push(inner);
                parent = inner;
            }

            Node leaf = new Node((S.length - suffix), suffix);
            parent.addChild(leaf);
            stack.push(leaf);
            prev = suffix;
        }
        root.updateLeafList();
        return root;
    }

    private static Node overmerge(Node even, Node odd, int[] S) {
        Node merged = new Node(0, 0);
        mergeDfs(merged, even, odd, S);
        merged.updateLeafList();
        return merged;
    }

    private static void mergeDfs(Node target, Node even, Node odd, int[] S) {
        List<Node> evenChildren = even.children;
        List<Node> oddChildren = odd.children;
        int e = 0;
        int o = 0;
        while (e < evenChildren.size() || o < oddChildren.size()) {
            Node eChild = e < evenChildren.size() ? evenChildren.get(e) : null;
            Node oChild = o < oddChildren.size() ? oddChildren.get(o) : null;

            Integer eChar = (eChild == null) ? null : edgeFirstChar(eChild, S);
            Integer oChar = (oChild == null) ? null : edgeFirstChar(oChild, S);

            if (eChild == null) {
                transferChild(target, oChild, false);
                o++;
                continue;
            }
            if (oChild == null) {
                transferChild(target, eChild, true);
                e++;
                continue;
            }

            if (!eChar.equals(oChar)) {
                if (oChar < eChar) {
                    transferChild(target, oChild, false);
                    o++;
                } else {
                    transferChild(target, eChild, true);
                    e++;
                }
            } else {
                int eEdgeLen = eChild.depth - target.depth;
                int oEdgeLen = oChild.depth - target.depth;
                if (eEdgeLen != oEdgeLen) {
                    Node shortChild = (eEdgeLen < oEdgeLen) ? eChild : oChild;
                    Node longChild = (shortChild == eChild) ? oChild : eChild;
                    boolean shortIsEven = (shortChild == eChild);

                    Node inner = new Node(shortChild.depth, shortChild.suffixIndex);
                    target.addChild(inner);
                    inner.evenSubtree = eChild;
                    inner.oddSubtree = oChild;

                    Node substitute = new Node(shortChild.depth, shortChild.suffixIndex);
                    substitute.addChild(longChild);

                    if (shortChild.children.isEmpty()) {
                        if (shortIsEven) {
                            inner.lcaEven = shortChild;
                        } else {
                            inner.lcaOdd = shortChild;
                        }
                    }

                    if (shortIsEven) {
                        mergeDfs(inner, shortChild, substitute, S);
                    } else {
                        mergeDfs(inner, substitute, shortChild, S);
                    }

                    bubbleLca(inner, target);
                } else {
                    Node inner = new Node(eChild.depth, eChild.suffixIndex);
                    inner.evenSubtree = eChild;
                    inner.oddSubtree = oChild;
                    if (eChild.children.isEmpty()) {
                        inner.lcaEven = eChild;
                    }
                    if (oChild.children.isEmpty()) {
                        inner.lcaOdd = oChild;
                    }
                    target.addChild(inner);
                    mergeDfs(inner, eChild, oChild, S);
                    bubbleLca(inner, target);
                }
                e++;
                o++;
            }
        }
    }

    private static void bubbleLca(Node child, Node parent) {
        boolean overwrite = parent.overwriteLca;
        if (child.lcaEven != null && (!parent.hasLcaEven || overwrite)) {
            parent.lcaEven = child.lcaEven;
            parent.hasLcaEven = true;
            overwrite = parent.overwriteLca;
        }
        if (child.lcaOdd != null && (!parent.hasLcaOdd || overwrite)) {
            parent.lcaOdd = child.lcaOdd;
            parent.hasLcaOdd = true;
            parent.overwriteLca = parent.hasLcaEven && parent.hasLcaOdd;
        }
    }

    private static void transferChild(Node target, Node child, boolean isEven) {
        if (isEven) {
            target.hasLcaEven = true;
            target.lcaEven = child.isLeaf() ? child : child.leafList.get(0);
        } else {
            target.hasLcaOdd = true;
            target.lcaOdd = child.isLeaf() ? child : child.leafList.get(0);
        }
        target.addChild(child);
    }

    private static Integer edgeFirstChar(Node node, int[] S) {
        Node sample = node.isLeaf() ? node : node.leafList.get(0);
        int pos = sample.suffixIndex + node.parent.depth - 1;
        if (pos < 0) {
            pos = 0;
        } else if (pos >= S.length) {
            pos = S.length - 1;
        }
        return S[pos];
    }

    private static void computeLcpTree(Node root) {
        List<Pair> basePairs = new ArrayList<>();
        root.dfs(node -> {
            if (node.lcaEven != null && node.lcaOdd != null) {
                basePairs.add(new Pair(node.lcaEven.suffixIndex, node.lcaOdd.suffixIndex));
            }
        });
        if (basePairs.isEmpty()) {
            root.lcpDepth = 0;
            return;
        }

        Map<Integer, Node> idToNode = new HashMap<>();
        root.dfs(node -> {
            if (node.isLeaf()) {
                idToNode.put(node.suffixIndex, node);
            }
        });

        Set<Pair> queries = new HashSet<>(basePairs.size() * 2);
        for (Pair pair : basePairs) {
            queries.add(pair);
            if (idToNode.containsKey(pair.first + 1) && idToNode.containsKey(pair.second + 1)) {
                queries.add(new Pair(pair.first + 1, pair.second + 1));
            }
        }

        Map<Pair, Node> lcaMap = tarjanLca(root, idToNode, queries);
        for (Pair pair : basePairs) {
            Node lca = lcaMap.get(pair);
            if (lca == null || lca == root) {
                continue;
            }
            if (!idToNode.containsKey(pair.first + 1) || !idToNode.containsKey(pair.second + 1)) {
                continue;
            }
            Pair shifted = new Pair(pair.first + 1, pair.second + 1);
            Node nextLca = lcaMap.get(shifted);
            if (nextLca == null || nextLca == lca) {
                continue;
            }
            lca.suffixLink = nextLca;
        }

        root.lcpDepth = 0;
        ensureLcpDepth(root);
        root.dfs(testfarach::ensureLcpDepth);
    }

    private static void ensureLcpDepth(Node node) {
        if (node.lcpDepth != null) {
            return;
        }
        if (node.parent == null) {
            node.lcpDepth = 0;
            return;
        }
        if (node.suffixLink == null) {
            return;
        }
        ensureLcpDepth(node.suffixLink);
        if (node.suffixLink.lcpDepth != null) {
            node.lcpDepth = node.suffixLink.lcpDepth + 1;
        }
    }

    private static void adjustOvermerge(Node root, int[] S) {
        root.bfs(node -> {
            if (node.lcpDepth != null && node.depth != node.lcpDepth) {
                node.children.clear();
                node.depth = node.lcpDepth;
                node.suffixIndex = -1;

                Node evenSub = node.evenSubtree;
                Node oddSub = node.oddSubtree;

                reparent(evenSub);
                reparent(oddSub);

                int evenChar = S[node.lcaEven.suffixIndex + node.lcpDepth - 1];
                int oddChar = S[node.lcaOdd.suffixIndex + node.lcpDepth - 1];

                if (evenChar < oddChar) {
                    node.addChild(evenSub);
                    node.addChild(oddSub);
                } else {
                    node.addChild(oddSub);
                    node.addChild(evenSub);
                }
            }
        });
        root.updateLeafList();
    }

    private static void reparent(Node subtree) {
        subtree.bfs(node -> node.children.forEach(child -> child.parent = node));
    }

    private static void cleanup(Node root) {
        root.dfs(node -> {
            node.suffixLink = null;
            node.evenSubtree = null;
            node.oddSubtree = null;
            node.lcaEven = null;
            node.lcaOdd = null;
            node.hasLcaEven = false;
            node.hasLcaOdd = false;
            node.overwriteLca = false;
            node.lcpDepth = null;
        });
    }

    private static Map<Pair, Node> tarjanLca(Node root,
                                             Map<Integer, Node> idToNode,
                                             Collection<Pair> queries) {
        Map<Node, List<Node>> q = new IdentityHashMap<>();
        Map<Pair, Node> answers = new HashMap<>();
        for (Pair pair : queries) {
            Node u = idToNode.get(pair.first);
            Node v = idToNode.get(pair.second);
            if (u == null || v == null) {
                continue;
            }
            q.computeIfAbsent(u, k -> new ArrayList<>()).add(v);
            q.computeIfAbsent(v, k -> new ArrayList<>()).add(u);
        }

        DisjointSet ds = new DisjointSet();
        root.dfs(node -> node.visited = false);
        tarjanDfs(root, q, answers, ds, new IdentityHashMap<>());
        return answers;
    }

    private static void tarjanDfs(Node u,
                                  Map<Node, List<Node>> queries,
                                  Map<Pair, Node> answers,
                                  DisjointSet ds,
                                  IdentityHashMap<Node, Node> ancestor) {
        ds.makeSet(u);
        ancestor.put(u, u);
        for (Node child : u.children) {
            tarjanDfs(child, queries, answers, ds, ancestor);
            ds.union(u, child);
            ancestor.put(ds.find(u), u);
        }
        u.visited = true;
        for (Node v : queries.getOrDefault(u, Collections.emptyList())) {
            if (v.visited) {
                Node lca = ancestor.get(ds.find(v));
                answers.put(new Pair(u.suffixIndex, v.suffixIndex), lca);
            }
        }
    }

    private static int[] remapToIntegers(String s) {

        int next = 1;
        int[] out = new int[s.length()];
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);
            map.putIfAbsent(c, next++);
            out[i] = map.get(c);
        }
        return out;
    }

    static int[] appendUniqueChar(int[] input) {
        int min = Integer.MAX_VALUE;
        for (int v : input) {
            if (v < min) {
                min = v;
            }
        }
        long candidate = (long) min - 1L;
        int sentinel;
        if (candidate < Integer.MIN_VALUE) {
            sentinel = Integer.MIN_VALUE;
        } else {
            sentinel = (int) candidate;
        }
        int[] result = Arrays.copyOf(input, input.length + 1);
        result[input.length] = sentinel;
        return result;
    }

    private static int[][] countingSortPairs(int[][] pairs) {
        if (pairs.length == 0) {
            return pairs;
        }
        int[][] sorted = countingSort(pairs, 1);
        sorted = countingSort(sorted, 0);
        return sorted;
    }

    private static int[][] countingSort(int[][] pairs, int index) {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (int[] pair : pairs) {
            min = Math.min(min, pair[index]);
            max = Math.max(max, pair[index]);
        }
        int offset = (min < 0) ? -min : 0;
        int size = max + offset + 1;
        int[] count = new int[size];
        for (int[] pair : pairs) {
            count[pair[index] + offset]++;
        }
        for (int i = 1; i < size; i++) {
            count[i] += count[i - 1];
        }
        int[][] out = new int[pairs.length][2];
        for (int i = pairs.length - 1; i >= 0; i--) {
            int[] pair = pairs[i];
            int key = pair[index] + offset;
            out[--count[key]] = pair;
        }
        return out;
    }

    private static long packPair(int a, int b) {
        return ((long) a << 32) ^ (b & 0xffffffffL);
    }

    private static final class DisjointSet {
        private final IdentityHashMap<Node, Node> parent = new IdentityHashMap<>();

        void makeSet(Node node) {
            parent.put(node, node);
        }

        Node find(Node node) {
            Node p = parent.get(node);
            if (p != node) {
                parent.put(node, p = find(p));
            }
            return p;
        }

        void union(Node a, Node b) {
            Node pa = find(a);
            Node pb = find(b);
            if (pa != pb) {
                parent.put(pb, pa);
            }
        }
    }

    public static final class Node {
        private List<Node> children = new ArrayList<>();
        private List<Node> leafList = new ArrayList<>();
        private Node parent;
        private Node evenSubtree;
        private Node oddSubtree;
        private Node suffixLink;
        private Node lcaEven;
        private Node lcaOdd;
        private boolean overwriteLca;
        private boolean hasLcaEven;
        private boolean hasLcaOdd;
        private boolean visited;
        private Integer lcpDepth;
        private int depth;
        private int suffixIndex;

        Node(int depth, int suffixIndex) {
            this.depth = depth;
            this.suffixIndex = suffixIndex;
        }

        void addChild(Node child) {
            children.add(child);
            child.parent = this;
        }

        boolean isLeaf() {
            return children.isEmpty();
        }

        void bfs(java.util.function.Consumer<Node> consumer) {
            java.util.Deque<Node> deque = new java.util.ArrayDeque<>();
            deque.add(this);
            while (!deque.isEmpty()) {
                Node node = deque.poll();
                consumer.accept(node);
                for (Node child : node.children) {
                    deque.add(child);
                }
            }
        }

        void dfs(java.util.function.Consumer<Node> consumer) {
            java.util.Deque<Node> stack = new java.util.ArrayDeque<>();
            stack.push(this);
            while (!stack.isEmpty()) {
                Node node = stack.pop();
                consumer.accept(node);
                for (int i = node.children.size() - 1; i >= 0; i--) {
                    stack.push(node.children.get(i));
                }
            }
        }

        List<Node> updateLeafList() {
            if (isLeaf()) {
                leafList = new ArrayList<>();
                leafList.add(this);
                return leafList;
            }
            leafList = new ArrayList<>();
            for (Node child : children) {
                leafList.addAll(child.updateLeafList());
            }
            return leafList;
        }

        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            bfs(node -> {
                sb.append("depth=").append(node.depth)
                        .append(" suffix=").append(node.suffixIndex)
                        .append(" children=").append(node.children.size())
                        .append('\n');
            });
            return sb.toString();
        }

        // ===== Added public getters for the adapter =====
        public List<Node> getChildren() { return java.util.Collections.unmodifiableList(children); }
        public List<Node> getLeafList() { return java.util.Collections.unmodifiableList(leafList); }
        public Node getParent() { return parent; }
        public int getDepth() { return depth; }
        public int getSuffixIndex() { return suffixIndex; }
    }


    private static final class Pair {
        final int first;
        final int second;

        Pair(int f, int s) {
            this.first = f;
            this.second = s;
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof Pair)) {
                return false;
            }
            Pair other = (Pair) o;
            return first == other.first && second == other.second;
        }

        @Override
        public int hashCode() {
            return (first * 31) ^ second;
        }
    }



    public static String readWholeFile(String absolutePath) throws IOException {
        Path p;
        p = Path.of(absolutePath);
        return Files.readString(p, StandardCharsets.UTF_8);
    }

    public static void main(String[] args) throws IOException {
        String t = readWholeFile("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/uniform_text_21_experiment.txt") +"$";
        long s = System.currentTimeMillis();
        testfarach.Node root = testfarach.build(t);
        System.out.println("Build time: " + (System.currentTimeMillis() - s) + "ms");


        String ss = "00S00T";
        int[] sss = remapToIntegers(ss);

        List<Node> children =  root.getChildren();

        for(Node c : children){
            List<Node> leafs = c.getLeafList();
            for(Node leaf : leafs){
                System.out.println(leaf.toString());
            }
            int a = 2;
        }
    }
}
