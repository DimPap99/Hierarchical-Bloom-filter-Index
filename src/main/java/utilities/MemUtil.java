package utilities;

import PMIndex.HBI;
import org.openjdk.jol.info.GraphLayout;
import org.openjdk.jol.vm.VM;
import tree.StreamBuffer;

public class MemUtil {

    // Detailed JOL report for an HBI instance, optionally including per-tree and class footprint.
    public String jolMemoryReport(boolean includePerTree, boolean includeFootprintTable, HBI hbi) {
        StringBuilder sb = new StringBuilder(16_384);

        // VM details (useful to interpret alignment, header sizes, and compressed oops status)
        sb.append("=== JOL / VM details ===\n");
        sb.append(VM.current().details()).append('\n');

        // Total retained size from the HBI instance as the single root
        GraphLayout total = GraphLayout.parseInstance(hbi);
        sb.append("\n=== HBI total (hbi as root) ===\n");
        sb.append("Total bytes       : ").append(total.totalSize()).append(" B\n");
        sb.append("Total bytes (MiB) : ")
                .append(String.format(java.util.Locale.ROOT, "%.3f", total.totalSize() / (1024.0 * 1024.0)))
                .append(" MiB\n");

        if (includeFootprintTable) {
            // Class-by-class histogram (very helpful to spot large contributors, e.g., BitSet, arrays)
            sb.append("\n--- Class footprint (HBI root) ---\n");
            sb.append(total.toFootprint()).append('\n');
        }

        if (includePerTree) {
            sb.append("\n=== Per-tree breakdown (each tree as its own root; do NOT sum) ===\n");
            for (int i = 0; i < hbi.trees.size(); i++) {
                tree.ImplicitTree<membership.Membership> t = hbi.trees.get(i);
                GraphLayout gl = GraphLayout.parseInstance(t);

                sb.append("Tree #").append(i)
                        .append("  bytes=").append(gl.totalSize()).append(" B")
                        .append("  (").append(String.format(java.util.Locale.ROOT, "%.3f",
                                gl.totalSize() / (1024.0 * 1024.0))).append(" MiB)\n");

                // You can comment hbi out if the output becomes too verbose.
                sb.append(gl.toFootprint()).append('\n');
            }
        }

        return sb.toString();
    }

    // Like jolMemoryReport but also returns the total MiB value.
    public MemoryUsageReport jolMemoryReportWithTotal(boolean includePerTree, boolean includeFootprintTable, HBI hbi) {
        String txt = jolMemoryReport(includePerTree, includeFootprintTable, hbi);
        long totalBytes = GraphLayout.parseInstance(hbi).totalSize();
        double totalMiB = totalBytes / (1024.0 * 1024.0);
        return new MemoryUsageReport(txt, totalMiB);
    }

    // Partitioned report for trees, AlphabetMapper, estimators, and total HBI memory.
    public String jolMemoryReportPartitioned(HBI hbi) {
        // Build a short report with four lines for easier parsing.
        StringBuilder sb = new StringBuilder(1024);

        // 0) Collect and detach estimators to avoid double counting in the "trees" bucket.
        java.util.List<estimators.Estimator> savedEstimators = new java.util.ArrayList<>(hbi.trees.size());
        for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
            savedEstimators.add(t.estimator); // may be null for empty trees; that is fine
            t.estimator = null;               // detach to exclude from "trees" bucket
        }

        try {
            // Trees (excluding estimators)
            long treesBytes = org.openjdk.jol.info.GraphLayout.parseInstance(hbi.trees).totalSize();

            // AlphabetMapper
            long alphabetBytes = 0;

            // Estimators (the ones that perform insertions)
            long estimatorsBytes = 0L;
            if (!savedEstimators.isEmpty()) {
                // Measure all estimator objects as independent roots.
                Object[] roots = savedEstimators.toArray();
                estimatorsBytes = org.openjdk.jol.info.GraphLayout.parseInstance(roots).totalSize();
            }

            // Reattach estimators before measuring the grand total.
            for (int i = 0; i < hbi.trees.size(); i++) {
                hbi.trees.get(i).estimator = savedEstimators.get(i);
            }

            // Grand total (entire HBI object graph)
            long totalBytes = org.openjdk.jol.info.GraphLayout.parseInstance(hbi).totalSize();

            // Format as bytes and MiB.
            java.util.Locale L = java.util.Locale.ROOT;
            String totalLine      = String.format(L, "Total: %d B (%.3f MiB)", totalBytes,      totalBytes      / (1024.0 * 1024.0));
            String treesLine      = String.format(L, "Trees (excl. estimators): %d B (%.3f MiB)", treesBytes,   treesBytes     / (1024.0 * 1024.0));
            String alphabetLine   = String.format(L, "AlphabetMapper: %d B (%.3f MiB)",           alphabetBytes, alphabetBytes / (1024.0 * 1024.0));
            String estimatorsLine = String.format(L, "Estimators: %d B (%.3f MiB)",               estimatorsBytes, estimatorsBytes / (1024.0 * 1024.0));

            // Exactly four lines, in the requested order.
            sb.append(totalLine).append('\n');
            sb.append(treesLine).append('\n');
            sb.append(alphabetLine).append('\n');
            sb.append(estimatorsLine);

            return sb.toString();
        } finally {
            // Safety: ensure estimators are reattached even if an exception occurs.
            for (int i = 0; i < hbi.trees.size(); i++) {
                if (hbi.trees.get(i).estimator == null) {
                    hbi.trees.get(i).estimator = (i < savedEstimators.size()) ? savedEstimators.get(i) : null;
                }
            }
        }
    }

    // Like jolMemoryReportPartitioned but also returns the total MiB value.
    public MemoryUsageReport jolMemoryReportPartitionedWithTotal(HBI hbi) {
        String txt = jolMemoryReportPartitioned(hbi);
        long totalBytes = org.openjdk.jol.info.GraphLayout.parseInstance(hbi).totalSize();
        double totalMiB = totalBytes / (1024.0 * 1024.0);
        return new MemoryUsageReport(txt, totalMiB);
    }




    public String fastMemoryEstimate(boolean includePayloadApprox, int approxIntegerBytes, HBI hbi) {
        // approxIntegerBytes parameter retained for backwards compatibility; long-backed buffers ignore it.
        final double LN2 = Math.log(2.0);
        final int ALIGN  = VM.current().objectAlignment();    // often 8 or 16
        final int REF    = inferReferenceSizeFromVmDetails();

        // Tunable overhead constants (ballpark, kept conservative)
        final int BITSET_OBJ_OVERHEAD      = 24;  // shell object fields (BitSet itself)
        final int LONG_ARRAY_HDR           = 16;  // long[] header
        final int BLOOM_FILTER_OBJ_OVERHEAD= 64;  // BloomFilter fields (n,m,k,p,seeds refs, temp arrays, etc.)
        final int AL_LIST_OBJ_OVERHEAD     = 24;  // ArrayList object header/fields (for levelFilters)
        final int OBJ_ARRAY_HDR            = 16;  // Object[] header
        final int STREAM_BUFFER_OBJ_OVERHEAD = 24; // StreamBuffer object shell

        long treesBytes = 0L;

        // Sum across trees
        for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
            //  Bloom filters per level (constructed exactly like your factory)
            int levels = t.effectiveDepth();                     // total levels in hbi tree
            int root = t.effectiveRoot();
            for (int level = root; level < levels; level++) {
                int nodes    = 1 << level;
                int interval = t.intervalSize(level);      // equals (treeLength >> level)
                int perNode  = Math.min(hbi.alphabetSize, interval);
                long nDistinct = (long) nodes * (long) perNode;

                if (nDistinct <= 0) continue;

                double mBitsD = - (nDistinct * Math.log(hbi.fpRate)) / (LN2 * LN2);
                long   mBits  = (long) Math.ceil(mBitsD);
                long   words  = (mBits + 63L) >>> 6;       // ceil(mBits / 64)
                long   bitsetBytes = alignUp(LONG_ARRAY_HDR + words * 8L, ALIGN);
                long   bfBytes     = alignUp(BITSET_OBJ_OVERHEAD + bitsetBytes + BLOOM_FILTER_OBJ_OVERHEAD, ALIGN);

                treesBytes += bfBytes;
            }

            //  LevelDirectory (list of per-level filters)
            long levelDirRefs = alignUp(AL_LIST_OBJ_OVERHEAD, ALIGN)
                    + alignUp(OBJ_ARRAY_HDR + (long) levels * REF, ALIGN);
            treesBytes += levelDirRefs;

            //  Codec + small objects (very small; lumped into a constant per tree)
            final int SMALLS = 128;  // Key64, TreeLayout, a few small fields
            treesBytes += SMALLS;

            //  StreamBuffer (backed by int[])
            int size = t.buffer.length();  // number of symbols in hbi tree buffer
            long bufferShell = alignUp(STREAM_BUFFER_OBJ_OVERHEAD, ALIGN);
            long payloadBytes = includePayloadApprox ? (long) size * Long.BYTES : 0L;
            long arrayBytes = alignUp(LONG_ARRAY_HDR + payloadBytes, ALIGN);
            treesBytes += bufferShell + arrayBytes;
        }

        //  AlphabetMapper measured with a small JOL traversal
        long alphabetBytes = 0L;

        //  Estimators measured with a small JOL traversal
        long estimatorsBytes = 0L;
        if (!hbi.trees.isEmpty()) {
            java.util.ArrayList<Object> roots = new java.util.ArrayList<>(hbi.trees.size());
            for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
                if (t.estimator != null) roots.add(t.estimator);
            }
            if (!roots.isEmpty()) {
                estimatorsBytes = GraphLayout.parseInstance(roots.toArray()).totalSize();
            }
        }

        long totalApprox = treesBytes + alphabetBytes + estimatorsBytes;

        java.util.Locale L = java.util.Locale.ROOT;
        String totalLine    = String.format(L, "Total (approx): %d B (%.3f MiB)", totalApprox, totalApprox / (1024.0 * 1024.0));
        String treesLabel   = includePayloadApprox
                ? "Trees (approx, incl. payload)"
                : "Trees (approx, excl. payload)";
        String treesLine    = String.format(L, "%s: %d B (%.3f MiB)", treesLabel, treesBytes, treesBytes / (1024.0 * 1024.0));
        String alphabetLine = String.format(L, "AlphabetMapper (JOL): %d B (%.3f MiB)", alphabetBytes, alphabetBytes / (1024.0 * 1024.0));
        String estimLine    = String.format(L, "Estimators (JOL): %d B (%.3f MiB)", estimatorsBytes, estimatorsBytes / (1024.0 * 1024.0));

        return totalLine + "\n" + treesLine + "\n" + alphabetLine + "\n" + estimLine;
    }

    // Align size up to the nearest multiple of alignment (power of two).
    private static long alignUp(long size, int alignment) {
        long a = alignment;
        return (size + (a - 1)) & ~(a - 1);
    }

    // Align size upwards with overflow safety.
    private static long alignUpSafe(long size, int alignment) {
        if (size <= 0) return 0L;
        long a = alignment;
        long capped = (size > Long.MAX_VALUE - (a - 1)) ? Long.MAX_VALUE - (a - 1) : size;
        long aligned = (capped + (a - 1)) & ~(a - 1);
        return aligned < 0 ? Long.MAX_VALUE : aligned;
    }

    // Safe addition that saturates at Long.MAX_VALUE.
    private static long safeAdd(long x, long y) {
        long r = x + y;
        if (((x ^ r) & (y ^ r)) < 0) return Long.MAX_VALUE;
        return r;
    }

    // Infer reference size (4 with compressed oops, else 8) from JOL VM details.
    private static int inferReferenceSizeFromVmDetails() {
        try {
            String d = VM.current().details().toLowerCase(java.util.Locale.ROOT);
            if (d.contains("compressed oops") || d.contains("compressed ordinary object pointers")) return 4;
        } catch (Throwable ignored) { }
        return 8;
    }

    // Exact JOL size of all membership filters across trees.
    public long jolExactMembershipBytes(HBI hbi) {
        java.util.ArrayList<Object> roots = new java.util.ArrayList<>();
        for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
            int levels = t.totalFilters();               // pass-through to LevelDirectory.depth()
            for (int L = 0; L < levels; L++) {
                membership.Membership m = t.filterAtLevel(L);
                if (m != null) roots.add(m);
            }
        }
        if (roots.isEmpty()) return 0L;
        return GraphLayout.parseInstance(roots.toArray()).totalSize();
    }
    // Exact JOL size of tree structures, excluding payload and estimators.
    public long jolExactTreesNoPayload(HBI hbi) {
        // Detach estimators to avoid counting them here
        java.util.ArrayList<estimators.Estimator> saved = new java.util.ArrayList<>(hbi.trees.size());
        for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
            saved.add(t.estimator);
            t.estimator = null;
        }

        // Temporarily detach each buffer's payload (O(1))
        java.util.ArrayList<StreamBuffer.Snapshot> snaps = new java.util.ArrayList<>(hbi.trees.size());
        for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
            snaps.add(t.buffer.detach());
        }

        try {
            return GraphLayout.parseInstance(hbi.trees).totalSize();
        } finally {
            // Restore buffers and estimators
            for (int i = 0; i < hbi.trees.size(); i++) {
                hbi.trees.get(i).buffer.restore(snaps.get(i));
            }
            for (int i = 0; i < hbi.trees.size(); i++) {
                if (hbi.trees.get(i).estimator == null) hbi.trees.get(i).estimator = saved.get(i);
            }
        }
    }

    // Hybrid report: exact membership and estimators, approximate shells and buffers.
    public String jolHybridMemoryReport(int approxIntegerBytes, HBI hbi) {
        // approxIntegerBytes retained for API compatibility; stream buffers are long-backed.
        // Exact membership bytes
        long membershipBytes = jolExactMembershipBytes(hbi);

        // Approximate: the small tree structures around the filters
        final int ALIGN = VM.current().objectAlignment();
        final int REF   = inferReferenceSizeFromVmDetails();
        final int AL_LIST_OBJ_OVERHEAD = 24;
        final int OBJ_ARRAY_HDR        = 16;
        final int SMALLS_PER_TREE      = 128;
        final int STREAM_BUFFER_OBJ_OVERHEAD = 24;
        final int LONG_ARRAY_HDR       = 16;

        long treeStructuresMinusMembership = 0L;
        long streamBufferBytes = 0L;

        for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
            int levels =  t.totalFilters();

            // LevelDirectory shells (ArrayList + Object[] of references)
            long levelList = alignUpSafe(AL_LIST_OBJ_OVERHEAD, ALIGN);
            long levelArr  = alignUpSafe(OBJ_ARRAY_HDR + (long) levels * REF, ALIGN);
            treeStructuresMinusMembership = safeAdd(treeStructuresMinusMembership, safeAdd(levelList, levelArr));

            // Small per-tree structures (codec/layout/etc.)
            treeStructuresMinusMembership = safeAdd(treeStructuresMinusMembership, SMALLS_PER_TREE);

            // StreamBuffer: full (structure + payload) so you can see data cost
            int size = t.buffer.length();
            long bufShell = alignUpSafe(STREAM_BUFFER_OBJ_OVERHEAD, ALIGN);
            long bufArray = alignUpSafe(LONG_ARRAY_HDR + (long) size * Long.BYTES, ALIGN);
            streamBufferBytes = safeAdd(streamBufferBytes, safeAdd(bufShell, bufArray));
        }

        // Exact with JOL: AlphabetMapper and Estimators (both small graphs)
        long alphabetBytes = 0L;

        long estimatorsBytes = 0L;
        if (!hbi.trees.isEmpty()) {
            java.util.ArrayList<Object> roots = new java.util.ArrayList<>(hbi.trees.size());
            for (tree.ImplicitTree<membership.Membership> t : hbi.trees) {
                if (t.estimator != null) roots.add(t.estimator);
            }
            if (!roots.isEmpty()) estimatorsBytes = GraphLayout.parseInstance(roots.toArray()).totalSize();
        }

        long treeStructuresBytes = safeAdd(membershipBytes, treeStructuresMinusMembership);
        long total = safeAdd(treeStructuresBytes, safeAdd(streamBufferBytes, safeAdd(alphabetBytes, estimatorsBytes)));

        java.util.Locale L = java.util.Locale.ROOT;
        return String.format(L, "Total (hybrid): %d B (%.3f MiB)", total, total/(1024.0*1024.0)) + "\n"
                + String.format(L, "TreeStructures (exact membership + approx shells): %d B (%.3f MiB)", treeStructuresBytes, treeStructuresBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "  |_ Membership (JOL exact): %d B (%.3f MiB)", membershipBytes, membershipBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "  |_ Non-membership shells (approx): %d B (%.3f MiB)", treeStructuresMinusMembership, treeStructuresMinusMembership/(1024.0*1024.0)) + "\n"
                + String.format(L, "StreamBuffer (approx): %d B (%.3f MiB)", streamBufferBytes, streamBufferBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "AlphabetMapper (JOL): %d B (%.3f MiB)", alphabetBytes, alphabetBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "Estimators (JOL): %d B (%.3f MiB)", estimatorsBytes, estimatorsBytes/(1024.0*1024.0));
    }


}
