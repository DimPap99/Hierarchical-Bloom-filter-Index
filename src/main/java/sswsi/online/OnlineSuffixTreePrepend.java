package sswsi.online;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

// Online suffix tree with prepend and rollback, wrapping a reversed Ukkonen core.
public interface OnlineSuffixTreePrepend {

    final class InsertResult {
        public final int parentV;
        public final int newLeafW;
        public final int firstSymbolOnNewEdge;
        public final List<Integer> matches;

        InsertResult(int parentV, int newLeafW, int firstSymbolOnNewEdge, List<Integer> matches) {
            this.parentV = parentV;
            this.newLeafW = newLeafW;
            this.firstSymbolOnNewEdge = firstSymbolOnNewEdge;
            this.matches = matches;
        }
    }

    void resetWithText(IntList revS3);

    void beginTxn(int queryId);

    void prepend(int symbol);

    InsertResult finishTxn();

    void rollbackTxn();

    List<Integer> enumerateMatches();

    int logicalLength();

    default void updateBaseRevS3(IntList revS3) {
        resetWithText(revS3);
    }

    final class RandomizedKopelowitzStyle implements OnlineSuffixTreePrepend {

        // Physical storage; logical text is reversed (charAt maps logical to physical).
        private final IntArrayList text = new IntArrayList();

        // Ukkonen active point
        private Node root;
        private Node activeNode;
        private int activeEdgeIndex;
        private int activeLength;
        private int remainingSuffixCount;
        private int leafEnd; // logical last index
        private Node lastCreatedInternal;

        // Query-time
        private final ArrayList<LogEntry> writeLog = new ArrayList<>();
        private final List<Integer> lastMatches = new ArrayList<>();
        private final IntArrayList txnStream = new IntArrayList();
        private boolean inTxn;
        private int currentTxnId;

        @Override
        public void resetWithText(IntList revS3) {
            clearAll();
            if (revS3 == null) return;
            // Build by logically prepending symbols (physically append in reverse order)
            for (int i = revS3.size() - 1; i >= 0; i--) {
                append(revS3.getInt(i));
            }
        }

        @Override
        public void updateBaseRevS3(IntList revS3) {
            if (revS3 == null) return;
            int n = revS3.size();
            if (n == logicalLength() + 1) {
                // True left extension by one symbol
                int newSym = revS3.getInt(0);
                append(newSym);
            } else {
                // Fallback rebuild on block boundary shift
                resetWithText(revS3);
            }
        }

        @Override
        public void beginTxn(int queryId) {
            inTxn = true;
            currentTxnId = queryId;
            txnStream.clear();
            lastMatches.clear();
            writeLog.clear();
        }

        @Override
        public void prepend(int symbol) {
            txnStream.add(symbol);
        }

        @Override
        public InsertResult finishTxn() {
            if (!inTxn) {
                return new InsertResult(0, 0, 0, Collections.emptyList());
            }
            if (txnStream.isEmpty()) {
                inTxn = false;
                return new InsertResult(0, 0, 0, Collections.emptyList());
            }
            int sentinel = txnStream.getInt(0);
            int m = Math.max(0, txnStream.size() - 1);
            int[] revP = new int[m];
            for (int i = 0; i < m; i++) revP[i] = txnStream.getInt(m - i);

            lastMatches.clear();
            LocusResult loc = findLocusOrMismatch(revP);
            int firstSymbol;
            if (loc.mismatch) {
                firstSymbol = loc.mismatchSymbol;
            } else {
                Node v = ensureExplicitNodeAtLocus(loc);
                Node dollarLeaf = new Node(leafEnd + 1, LEAF_SENTINEL);
                dollarLeaf.suffixIndex = -1;
                addChildWithLog(v, 0, dollarLeaf);
                enumerateFromNodeExcludingKeyBounded(v, 0, m, lastMatches);
                firstSymbol = 0;
            }

            inTxn = false;
            txnStream.clear();
            return new InsertResult(0, leafEnd + 1, firstSymbol, lastMatches);
        }

        @Override
        public void rollbackTxn() {
            inTxn = false;
            txnStream.clear();
            lastMatches.clear();
            rollbackWriteLog();
        }

        @Override
        public List<Integer> enumerateMatches() {
            return Collections.unmodifiableList(lastMatches);
        }

        @Override
        public int logicalLength() {
            return text.size();
        }

        // ===== Ukkonen over logical reversal =====
        private void clearAll() {
            text.clear();
            root = new Node(-1, -1);
            root.suffixLink = root;
            activeNode = root;
            activeEdgeIndex = -1;
            activeLength = 0;
            remainingSuffixCount = 0;
            leafEnd = -1;
            lastCreatedInternal = null;
            inTxn = false;
            currentTxnId = -1;
            lastMatches.clear();
            txnStream.clear();
            writeLog.clear();
        }

        private void append(int symbol) {
            text.add(symbol);
            leafEnd = text.size() - 1;
            remainingSuffixCount++;
            lastCreatedInternal = null;

            while (remainingSuffixCount > 0) {
                if (activeLength == 0) {
                    activeEdgeIndex = leafEnd;
                }
                int edgeSymbol = charAt(activeEdgeIndex);
                Node next = activeNode.child(edgeSymbol);

                if (next == null) {
                    Node leaf = new Node(leafEnd, LEAF_SENTINEL);
                    leaf.suffixIndex = leafEnd - remainingSuffixCount + 1;
                    activeNode.putChild(edgeSymbol, leaf);
                    if (activeNode.leafPtr == null) activeNode.leafPtr = leaf;
                    if (lastCreatedInternal != null) {
                        lastCreatedInternal.suffixLink = activeNode;
                        lastCreatedInternal = null;
                    }
                } else {
                    if (walkDown(next)) continue;
                    int nextSymbol = charAt(next.start + activeLength);
                    if (nextSymbol == symbol) {
                        if (lastCreatedInternal != null && activeNode != root) {
                            lastCreatedInternal.suffixLink = activeNode;
                            lastCreatedInternal = null;
                        }
                        activeLength++;
                        break;
                    }
                    int splitEnd = next.start + activeLength - 1;
                    Node split = new Node(next.start, splitEnd);
                    activeNode.putChild(edgeSymbol, split);

                    Node leaf = new Node(leafEnd, LEAF_SENTINEL);
                    leaf.suffixIndex = leafEnd - remainingSuffixCount + 1;
                    split.putChild(symbol, leaf);
                    split.leafPtr = leaf;

                    next.start = next.start + activeLength;
                    int nextEdgeSymbol = charAt(next.start);
                    split.putChild(nextEdgeSymbol, next);

                    if (lastCreatedInternal != null) {
                        lastCreatedInternal.suffixLink = split;
                    }
                    lastCreatedInternal = split;
                    split.suffixLink = root;
                }

                remainingSuffixCount--;
                if (activeNode == root && activeLength > 0) {
                    activeLength--;
                    activeEdgeIndex = leafEnd - remainingSuffixCount + 1;
                } else if (activeNode != root) {
                    activeNode = activeNode.suffixLink != null ? activeNode.suffixLink : root;
                }
            }
        }

        private boolean walkDown(Node next) {
            int edgeLen = next.edgeLength(leafEnd);
            if (activeLength >= edgeLen) {
                activeEdgeIndex += edgeLen;
                activeLength -= edgeLen;
                activeNode = next;
                return true;
            }
            return false;
        }

        private int charAt(int logicalIndex) {
            int size = text.size();
            int physical = size - 1 - logicalIndex;
            return text.getInt(physical);
        }

        // ===== Query helpers =====
        private LocusResult findLocusOrMismatch(int[] pattern) {
            Node current = root;
            Node parent = null;
            int edgeFirstSymbol = -1;
            int idx = 0;
            while (idx < pattern.length) {
                int sym = pattern[idx];
                Node next = current.child(sym);
                if (next == null) return LocusResult.mismatch(sym);
                int i = next.start;
                int end = next.effectiveEnd(leafEnd);
                int consumed = 0;
                while (i + consumed <= end && idx + consumed < pattern.length) {
                    if (charAt(i + consumed) != pattern[idx + consumed]) {
                        return LocusResult.mismatch(pattern[idx]);
                    }
                    consumed++;
                }
                if (idx + consumed == pattern.length) {
                    int edgeLen = end - i + 1;
                    boolean atNode = (consumed == edgeLen);
                    return LocusResult.at(current, next, charAt(next.start), consumed, atNode);
                }
                parent = current;
                edgeFirstSymbol = charAt(next.start);
                current = next;
                idx += consumed;
            }
            return LocusResult.at(parent, current, edgeFirstSymbol, 0, true);
        }

        private Node ensureExplicitNodeAtLocus(LocusResult loc) {
            if (loc.atNode) return loc.child;
            Node parent = loc.parent;
            Node child = loc.child;
            int splitOffset = loc.consumedWithinEdge;
            int oldStart = child.start;
            int splitPoint = oldStart + splitOffset;
            Node split = new Node(oldStart, splitPoint - 1);
            int edgeKey = loc.edgeFirstSymbol;
            replaceChildWithLog(parent, edgeKey, child, split);
            logModifyStart(child, oldStart);
            child.start = splitPoint;
            int childKey = charAt(child.start);
            addChildWithLog(split, childKey, child);
            // Set leaf pointer conservatively
            if (child.isLeaf()) split.leafPtr = child;
            return split;
        }

        private void enumerateFromNodeExcludingKeyBounded(Node v, int excludedKey, int patLen, List<Integer> out) {
            ArrayDeque<Node> stack = new ArrayDeque<>();
            if (v.children != null) {
                for (Int2ObjectMap.Entry<Node> e : v.children.int2ObjectEntrySet()) {
                    if (e.getIntKey() == excludedKey) continue;
                    stack.push(e.getValue());
                }
            }
            while (!stack.isEmpty()) {
                int budget = 4;
                while (budget-- > 0 && !stack.isEmpty()) {
                    Node cur = stack.pop();
                    if (cur.isLeaf()) {
                        int start = cur.suffixIndex;
                        if (start >= 0 && start + patLen <= logicalLength()) out.add(start);
                        continue;
                    }
                    if (cur.leafPtr != null) {
                        int start = cur.leafPtr.suffixIndex;
                        if (start >= 0 && start + patLen <= logicalLength()) out.add(start);
                    }
                    if (cur.children != null) {
                        for (Node ch : cur.children.values()) {
                            stack.push(ch);
                        }
                    }
                }
            }
        }

        private void addChildWithLog(Node parent, int key, Node child) {
            parent.children.put(key, child);
            writeLog.add(LogEntry.addChild(parent, key));
        }

        private void replaceChildWithLog(Node parent, int key, Node oldChild, Node newChild) {
            parent.children.put(key, newChild);
            writeLog.add(LogEntry.replaceChild(parent, key, oldChild));
        }

        private void logModifyStart(Node node, int oldStart) {
            writeLog.add(LogEntry.modifyStart(node, oldStart));
        }

        private void rollbackWriteLog() {
            for (int i = writeLog.size() - 1; i >= 0; i--) {
                LogEntry e = writeLog.get(i);
                switch (e.kind) {
                    case ADD_CHILD -> e.parent.children.remove(e.key);
                    case REPLACE_CHILD -> e.parent.children.put(e.key, e.oldChild);
                    case MODIFY_START -> e.node.start = e.oldStart;
                }
            }
            writeLog.clear();
        }

        private static final int LEAF_SENTINEL = -1;

        private static final class LocusResult {
            final boolean mismatch;
            final int mismatchSymbol;
            final Node parent;
            final Node child;
            final int edgeFirstSymbol;
            final int consumedWithinEdge;
            final boolean atNode;
            private LocusResult(boolean mismatch, int mismatchSymbol, Node parent, Node child,
                                int edgeFirstSymbol, int consumedWithinEdge, boolean atNode) {
                this.mismatch = mismatch;
                this.mismatchSymbol = mismatchSymbol;
                this.parent = parent;
                this.child = child;
                this.edgeFirstSymbol = edgeFirstSymbol;
                this.consumedWithinEdge = consumedWithinEdge;
                this.atNode = atNode;
            }
            static LocusResult mismatch(int sym) { return new LocusResult(true, sym, null, null, -1, 0, false); }
            static LocusResult at(Node parent, Node child, int edgeFirstSymbol, int consumedWithinEdge, boolean atNode) {
                return new LocusResult(false, -1, parent, child, edgeFirstSymbol, consumedWithinEdge, atNode);
            }
        }

        private static final class LogEntry {
            enum Kind { ADD_CHILD, REPLACE_CHILD, MODIFY_START }
            final Kind kind; final Node parent; final int key; final Node oldChild; final Node node; final int oldStart;
            private LogEntry(Kind kind, Node parent, int key, Node oldChild, Node node, int oldStart) {
                this.kind = kind; this.parent = parent; this.key = key; this.oldChild = oldChild; this.node = node; this.oldStart = oldStart;
            }
            static LogEntry addChild(Node parent, int key) { return new LogEntry(Kind.ADD_CHILD, parent, key, null, null, 0); }
            static LogEntry replaceChild(Node parent, int key, Node oldChild) { return new LogEntry(Kind.REPLACE_CHILD, parent, key, oldChild, null, 0); }
            static LogEntry modifyStart(Node node, int oldStart) { return new LogEntry(Kind.MODIFY_START, null, 0, null, node, oldStart); }
        }

        private final class Node {
            int start; int end; Int2ObjectOpenHashMap<Node> children = new Int2ObjectOpenHashMap<>();
            Node suffixLink; int suffixIndex = -1; Node leafPtr;
            Node(int start, int end) { this.start = start; this.end = end; this.children.defaultReturnValue(null); }
            boolean isLeaf() { return children.isEmpty(); }
            int effectiveEnd(int logicalEnd) { return (end == LEAF_SENTINEL) ? logicalEnd : end; }
            int edgeLength(int logicalEnd) { return (start < 0) ? 0 : effectiveEnd(logicalEnd) - start + 1; }
            Node child(int symbol) { return children.get(symbol); }
            void putChild(int symbol, Node child) { children.put(symbol, child); if (leafPtr == null && child.isLeaf()) leafPtr = child; }
        }
    }
}
