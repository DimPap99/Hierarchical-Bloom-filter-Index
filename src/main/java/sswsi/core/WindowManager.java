package sswsi.core;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

// Maintains the sliding window and block layout for the stream.
public final class WindowManager {

    public static final class BlockView {
        public final int id;
        public final int globalStart;
        public final int globalEndExclusive;
        public final IntList tokens;

        BlockView(int id, int globalStart, int globalEndExclusive, IntList tokens) {
            this.id = id;
            this.globalStart = globalStart;
            this.globalEndExclusive = globalEndExclusive;
            this.tokens = tokens;
        }

        public int length() {
            return tokens.size();
        }
    }

    private static final class Block {
        final int id;
        int globalStart;
        int globalEndExclusive;
        final IntArrayList tokens = new IntArrayList();

        Block(int id, int globalStart) {
            this.id = id;
            this.globalStart = globalStart;
            this.globalEndExclusive = globalStart;
        }

        BlockView view() {
            return new BlockView(id, globalStart, globalEndExclusive, tokens);
        }
    }

    private final int windowSize;
    private final int blockLength;
    private final ArrayDeque<Block> blocks = new ArrayDeque<>();
    private int totalInserted = 0;
    private int windowStartGlobal = 0;

    public WindowManager(int windowSize, int blockLength) {
        if (windowSize <= 0) {
            throw new IllegalArgumentException("windowSize must be positive");
        }
        if (blockLength <= 0) {
            throw new IllegalArgumentException("blockLength must be positive");
        }
        this.windowSize = windowSize;
        this.blockLength = blockLength;
    }

    public int windowSize() {
        return windowSize;
    }

    public int blockLength() {
        return blockLength;
    }

    public int totalInserted() {
        return totalInserted;
    }

    public int windowStartGlobal() {
        return windowStartGlobal;
    }

    public int windowEndExclusive() {
        return totalInserted;
    }

    public int windowLength() {
        return Math.min(windowSize, totalInserted);
    }

    public int toWindowOffset(int globalPosition) {
        return globalPosition - windowStartGlobal;
    }

    // Append one token (by its mapped integer id) and trim to the window.
    public void append(int token) {
        int globalIndex = totalInserted;
        ensureBlockFor(globalIndex);
        Block block = blocks.peekLast();
        block.tokens.add(token);
        block.globalEndExclusive = globalIndex + 1;
        totalInserted++;
        trimToWindow();
    }

    private void ensureBlockFor(int globalIndex) {
        int blockId = Math.floorDiv(globalIndex, blockLength);
        if (blocks.isEmpty() || blocks.peekLast().id != blockId) {
            Block fresh = new Block(blockId, globalIndex);
            blocks.addLast(fresh);
        }
    }

    private void trimToWindow() {
        int cutoff = Math.max(0, totalInserted - windowSize);
        windowStartGlobal = cutoff;
        while (!blocks.isEmpty() && blocks.peekFirst().globalEndExclusive <= cutoff) {
            blocks.removeFirst();
        }
        if (blocks.isEmpty()) {
            return;
        }
        Block first = blocks.peekFirst();
        if (first.globalStart < cutoff) {
            int remove = cutoff - first.globalStart;
            if (remove >= first.tokens.size()) {
                first.tokens.clear();
                first.globalStart = cutoff;
            } else {
                first.tokens.removeElements(0, remove);
                first.globalStart = cutoff;
            }
        }
    }

    public List<BlockView> snapshotBlocks() {
        if (blocks.isEmpty()) {
            return Collections.emptyList();
        }
        List<BlockView> out = new ArrayList<>(blocks.size());
        for (Block block : blocks) {
            out.add(block.view());
        }
        return Collections.unmodifiableList(out);
    }

    public IntArrayList gatherFromTail(int startFromTailInclusive,
                                       int endFromTailInclusive,
                                       IntArrayList reuse) {
        if (reuse == null) {
            reuse = new IntArrayList();
        } else {
            reuse.clear();
        }
        if (blocks.isEmpty()) {
            return reuse;
        }
        int size = blocks.size();
        int lower = Math.max(0, startFromTailInclusive);
        int upper = Math.max(lower, endFromTailInclusive);
        for (int offset = upper; offset >= lower; offset--) {
            if (offset >= size) {
                continue;
            }
            Block block = getFromTail(offset);
            reuse.addAll(block.tokens);
        }
        return reuse;
    }

    public BlockView getBlockFromTail(int offset) {
        if (blocks.isEmpty()) {
            return null;
        }
        if (offset < 0 || offset >= blocks.size()) {
            return null;
        }
        return getFromTail(offset).view();
    }

    private Block getFromTail(int offset) {
        int idx = 0;
        for (var it = blocks.descendingIterator(); it.hasNext(); ) {
            Block candidate = it.next();
            if (idx == offset) {
                return candidate;
            }
            idx++;
        }
        throw new IllegalArgumentException("Offset " + offset + " outside current block range");
    }

    public IntArrayList windowAsList(IntArrayList reuse) {
        if (reuse == null) {
            reuse = new IntArrayList(windowLength());
        } else {
            reuse.clear();
        }
        for (Block block : blocks) {
            reuse.addAll(block.tokens);
        }
        return reuse;
    }
}
