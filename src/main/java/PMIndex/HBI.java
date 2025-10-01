package PMIndex;

import estimators.CostFunction;
import estimators.Estimator;
import jdk.jshell.execution.Util;
import membership.Key64;
import membership.Membership;
import org.apache.commons.math3.util.Pair;
import org.openjdk.jol.info.GraphLayout;
import org.openjdk.jol.vm.VM;
import search.CandidateRange;
import search.Frame;
import search.IntervalScanner;
import search.Pattern;
import search.PruningPlan;
import search.SearchAlgorithm;
import search.Verifier;
import tree.ImplicitTree;
import tree.TreeLayout;
import utilities.AlphabetMapper;
import utilities.PatternResult;

import java.util.*;
import java.util.function.IntFunction;
import java.util.function.Supplier;

import static utilities.MathUtils.pruningLevel;


public final class HBI implements IPMIndexing {

    private static final int DEFAULT_BC_COST_ESTIM_ITER = 1_000_000;
    private static final int DEFAULT_LC_COST_ESTIM_ITER = 1_000_000;

    private final HbiConfiguration config;

    private int nGram = 1;
    private final HbiStats stats;

    private final int windowLength;
    private final int treeLength;
    private final int alphabetSize;
    private final double fpRate;

    private final SearchAlgorithm searchAlgo;
    private final Supplier<Estimator> estimatorFac;
    private final Supplier<Membership> membershipFac;
    private final Supplier<PruningPlan> pruningPlanFac;
    private final Verifier verifier;
    private final CostFunction cf;
    private boolean builtModel = false;
    private final double conf;

    private final ArrayList<ImplicitTree<Membership>> trees = new ArrayList<>();

    private AlphabetMapper<String> alphabetMap;
    private long indexedItemsCounter = -1;

    private int bcCostEstimIter = DEFAULT_BC_COST_ESTIM_ITER;
    private int lcCostEstimIter = DEFAULT_LC_COST_ESTIM_ITER;


    private double bfCost = 97;
    private double leafCost = 26;
    private int lpOverride;
    public NgramModel.Builder modelBuilder;
    public boolean memoryPolicy = false;
    private final int maxActiveTrees;

    public HBI(HbiConfiguration config) {
        this.config = Objects.requireNonNull(config, "config");
        this.stats = new HbiStats(config.collectStats(), config.experimentMode());
        this.windowLength = config.windowLength();
        this.treeLength = config.treeLength();
        this.alphabetSize = config.alphabetSize();
        this.fpRate = config.fpRate();
        this.searchAlgo = config.searchAlgorithm();
        this.estimatorFac = config.estimatorSupplier();
        this.membershipFac = config.membershipSupplier();
        this.pruningPlanFac = config.pruningPlanSupplier();
        this.verifier = config.verifier();
        this.cf = config.costFunction();
        this.conf = config.confidence();
        this.alphabetMap = new AlphabetMapper<>(alphabetSize);
        this.nGram = config.nGram();

        trees.addLast(createTree());
        this.maxActiveTrees = (int) Math.ceil((double) (windowLength / treeLength));
        //this.modelBuilder = new NgramModel.Builder(89*89, 2);
    }

    public HBI(SearchAlgorithm algo,
               int windowLength,
               double fpRate,
               int alphabetSize,
               int treeLength,
               Supplier<Estimator> estimator,
               Supplier<Membership> membership,
               Supplier<PruningPlan> pruningPlan,
               Verifier verifier,
               CostFunction cf,
               double conf,
               int nGram) {
        this(HbiConfiguration.builder()
                .searchAlgorithm(algo)
                .windowLength(windowLength)
                .fpRate(fpRate)
                .alphabetSize(alphabetSize)
                .treeLength(treeLength)
                .estimatorSupplier(estimator)
                .membershipSupplier(membership)
                .pruningPlanSupplier(pruningPlan)
                .verifier(verifier)
                .costFunction(cf)
                .confidence(conf)
                .nGram(nGram)
                .build());
    }

    public HbiStats stats() {
        return stats;
    }

    public void setAlphabetMap(AlphabetMapper<String> alphabetMap) {
        this.alphabetMap = Objects.requireNonNull(alphabetMap, "alphabetMap");
    }

    public void resetAlphabetMap(int capacity) {
        this.alphabetMap = new AlphabetMapper<>(capacity);
    }

    public void setLpOverride(int lp) {
        this.lpOverride = lp;
    }

    public int lpOverride() {
        return lpOverride;
    }

    /** Evict the oldest tree unconditionally. */
    @Override
    public void expire() {
        trees.removeFirst();
    }

    /** Stream a single character into the index. */
    @Override
    public void insert(String c) {
        int intC = alphabetMap.insert(c);
//        this.modelBuilder.observeSymbol(intC);
        indexedItemsCounter++;
        ImplicitTree<Membership> lastTree = trees.getLast();
        if (lastTree.isFull()) {
            ImplicitTree<Membership> fresh = createTree();
            fresh.id = trees.size();
            fresh.estimator.insert(intC);
            fresh.append(intC, indexedItemsCounter);
            trees.add(fresh);
        } else {
            lastTree.estimator.insert(intC);
            lastTree.append(intC, indexedItemsCounter);
        }
        if (trees.getLast().indexedItemsCounter + (trees.size() - 2) * this.treeLength > this.windowLength) {
            this.expire();
        }
    }

    public void applyMemoryPolicy() {
        if(memoryPolicy) {
            double minp = this.trees.getLast().estimator.getMin();
            int lp = pruningLevel(this.trees.getLast(), 0.95, minp);
            this.trees.getLast().dropFiltersUpToLp(lp);

        }
    }

    public void fillStackLp(int lp, Deque<Frame> fStack) {
        for (int i = 0; i < Math.pow(2, lp); i++) {
            fStack.addLast(new Frame(lp, i));
        }
    }

    @Override
    public boolean exists(String key) {
        return false;
    }

    @Override
    public ArrayList<Integer> report(Pattern pat) {

        if(!this.builtModel) {
//            this.modelBuilder.Sigma
//            cf.setModel(this.modelBuilder.build());

        };
        ArrayList<Integer> results = new ArrayList<>();
        long queryStartNanos = System.nanoTime();
        long totalLpTimeNanos = 0L;
        for (int nIdx = 0; nIdx < pat.nGramArr.length; nIdx++) {
            int ngToInt = this.alphabetMap.getId(pat.nGramArr[nIdx]);
            pat.nGramToInt[nIdx] = ngToInt;
            if(nIdx % pat.nGram == 0) pat.effectiveNgramArr[nIdx/pat.nGram] = ngToInt;
        }
        if(pat.originalSz % pat.nGram != 0) pat.effectiveNgramArr[pat.effectiveNgramArr.length-1] = pat.nGramToInt[pat.nGramToInt.length-1];

        int positionOffset = -1;
        long startTime = System.currentTimeMillis();
        int lp = 0;
        int lpCf = 0;
        double cp_cost = 0;
        int arbitraryConfLp = 0;
        for (int i = 0; i < this.trees.size(); i++) {
            ImplicitTree<Membership> tree = this.trees.get(i);
            tree.pruningPlan = this.pruningPlanFac.get();
            IntervalScanner scn = new IntervalScanner(tree, pat, searchAlgo, positionOffset);
            Deque<Frame> stack = new ArrayDeque<>();
            double[] pp = tree.estimator.estimateALl(pat);
            double pMax = Arrays.stream(pp).min().getAsDouble();
            long lpStart = System.nanoTime();
            ArrayList<Integer> lps = null;
//            lp = Collections.min(lps);
            totalLpTimeNanos += System.nanoTime() - lpStart;
              //cf.minCostLp(tree, 0.05, pat, 97, 26);//pruningLevel(tree, 0.99, pMax);

            if (stats.isExperimentMode()) {
                pp = tree.estimator.estimateALl(pat);
                pMax = Arrays.stream(pp).min().getAsDouble();
                lp = this.lpOverride;
                arbitraryConfLp = pruningLevel(tree, this.conf, pMax);
                int m = (int) (tree.maxDepth() - 1 - Math.ceil(Math.log(pat.nGramToInt.length) / Math.log(2)));
                cp_cost = cf.costAtLevel(tree, pp, pat.nGramToInt, lp, 0.0, m);
                lpCf = cf.minCostLp(tree, 0.05, pat, 97, 26);
            } else {
                lps= tree.pruningPlan.pruningPlan(pat, tree, 0.99);
//                lp = lpCf;
            }

            pat.charStartLp = lps;

            if (stats.isCollecting()) {
                stats.recordLp(lp);
                stats.recordAlpha(cf.getAlpha());
            }
            fillStackLp(lp, stack);
            scn.seedStack(stack);

            while (scn.hasNext()) {
                CandidateRange cr = scn.next();
                if (cr == null) break;
                Pair<ArrayList<Integer>, Integer> res = this.verifier.verify(i, this.trees, cr, pat);
                for (int r : res.getFirst()) {
                    results.add(r);
                }
                scn.positionOffset = res.getSecond();
            }
        }
        long queryDuration = System.nanoTime() - queryStartNanos;
        stats.recordQueryTiming(queryDuration, totalLpTimeNanos);
        if (stats.isExperimentMode()) {

            int leafProbes = this.verifier.getLeafProbes();
            int bfprobes = this.getAllprobes();
            int actualCost = bfprobes;
//            System.out.println("Pattern: " + pat.text +" Probes: " + bfprobes + " Leafs: " + bfprobes + " Actual: " + actualCost);

            this.verifier.reset();
            stats.setLatestPatternResult(new PatternResult(System.currentTimeMillis() - startTime, actualCost, lp, pat, lpCf, cp_cost, leafProbes, arbitraryConfLp));
        }

        return results;
    }

    public long getAvgBloomCost(Pattern pat) {
        ImplicitTree<Membership> tree = this.trees.getLast();
        int maxLvl = tree.maxDepth() - 1;
        long duration = 0;
        for (int z = 0; z < bcCostEstimIter; z++) {
            long startTime = System.nanoTime();
            long key = tree.codec.pack(maxLvl, 0, pat.nGramToInt[0]);
            tree.contains(maxLvl, key);
            duration += System.nanoTime() - startTime;
        }
        duration /= bcCostEstimIter;
        return duration;
    }

    @Override
    public ArrayList<Long> getAvgTimes(Pattern pat) {
        long avgBloom = getAvgBloomCost(pat);
        long avgLeaf = getAvgLeafCost(pat);
        this.bfCost = avgBloom;
        this.leafCost = avgLeaf;
        ArrayList<Long> res = new ArrayList<>();
        res.add(avgBloom);
        res.add(avgLeaf);
        return res;
    }

    @Override
    public PatternResult getLatestStats() {
        return stats.latestPatternResult();
    }

    public long getAvgLeafCost(Pattern pat) {
        ImplicitTree<Membership> tree = this.trees.getLast();
        long duration = 0;
        for (int i = 0; i < lcCostEstimIter; i++) {
            long startTime = System.nanoTime();
            for (int j = 0; j < pat.nGramToInt.length - 1; j++) {
                boolean in = tree.buffer.data.get(0) == 1;
            }
            duration += System.nanoTime() - startTime;
        }
        duration /= lcCostEstimIter;
        return duration;
    }

    private TreeLayout makeLayout() {
        int levels = Integer.SIZE - Integer.numberOfLeadingZeros(treeLength);
        levels = levels - (nGram/2);

        int leafSpan = treeLength >>> (levels - 1);
        return new TreeLayout(levels, leafSpan);

    }

    private ImplicitTree<Membership> createTree() {
        TreeLayout layout = makeLayout();
        int maxDepth = layout.levels();
        IntFunction<Membership> filterFactory = level -> {
            int nodes = 1 << level;
            int interval = treeLength >> level;
            int perNode = Math.min(alphabetSize, interval);
            int distinct = nodes * perNode;

            Membership bf = membershipFac.get();
            bf.init(distinct, fpRate);
            return bf;
        };

        return new ImplicitTree<>(
                layout,
                filterFactory,
                new Key64(maxDepth, alphabetSize),
                estimatorFac.get());
    }

    public int getAllprobes() {
        int sum = 0;
        for (ImplicitTree<Membership> tree : this.trees) {
            sum += tree.containCounter;
            tree.containCounter = 0;
        }
        return sum;
    }

    public String jolMemoryReport(boolean includePerTree, boolean includeFootprintTable) {
        StringBuilder sb = new StringBuilder(16_384);

        // VM details (useful to interpret alignment, header sizes, and compressed oops status)
        sb.append("=== JOL / VM details ===\n");
        sb.append(VM.current().details()).append('\n');

        // Total retained size from the HBI instance as the single root
        GraphLayout total = GraphLayout.parseInstance(this);
        sb.append("\n=== HBI total (this as root) ===\n");
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
            for (int i = 0; i < this.trees.size(); i++) {
                tree.ImplicitTree<membership.Membership> t = this.trees.get(i);
                GraphLayout gl = GraphLayout.parseInstance(t);

                sb.append("Tree #").append(i)
                        .append("  bytes=").append(gl.totalSize()).append(" B")
                        .append("  (").append(String.format(java.util.Locale.ROOT, "%.3f",
                                gl.totalSize() / (1024.0 * 1024.0))).append(" MiB)\n");

                // You can comment this out if the output becomes too verbose.
                sb.append(gl.toFootprint()).append('\n');
            }
        }

        return sb.toString();
    }

    public String jolMemoryReportPartitioned() {
        // Build a report; keep exactly four lines of results for easy parsing/printing.
        StringBuilder sb = new StringBuilder(1024);

        // 0) Collect and detach estimators to avoid double counting in the "trees" bucket.
        java.util.List<estimators.Estimator> savedEstimators = new java.util.ArrayList<>(this.trees.size());
        for (tree.ImplicitTree<membership.Membership> t : this.trees) {
            savedEstimators.add(t.estimator); // may be null for empty trees; that is fine
            t.estimator = null;               // detach to exclude from "trees" bucket
        }

        try {
            // 1) Trees (excluding estimators)
            long treesBytes = org.openjdk.jol.info.GraphLayout.parseInstance(this.trees).totalSize();

            // 2) AlphabetMapper
            long alphabetBytes = (this.alphabetMap != null)
                    ? org.openjdk.jol.info.GraphLayout.parseInstance(this.alphabetMap).totalSize()
                    : 0L;

            // 3) Estimators (the ones that perform insertions)
            long estimatorsBytes = 0L;
            if (!savedEstimators.isEmpty()) {
                // Measure all estimator objects as independent roots.
                Object[] roots = savedEstimators.toArray();
                estimatorsBytes = org.openjdk.jol.info.GraphLayout.parseInstance(roots).totalSize();
            }

            // Reattach estimators before measuring the grand total.
            for (int i = 0; i < this.trees.size(); i++) {
                this.trees.get(i).estimator = savedEstimators.get(i);
            }

            // 4) Grand total (entire HBI object graph)
            long totalBytes = org.openjdk.jol.info.GraphLayout.parseInstance(this).totalSize();

            // Format as bytes and mebibytes (MiB = 1,048,576 bytes).
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
            for (int i = 0; i < this.trees.size(); i++) {
                if (this.trees.get(i).estimator == null) {
                    this.trees.get(i).estimator = (i < savedEstimators.size()) ? savedEstimators.get(i) : null;
                }
            }
        }
    }



    /**
     * FAST (approximate) memory report in four lines:
     *   1) Total (approx) = Trees + AlphabetMapper + Estimators
     *   2) Trees (approx, excludes estimators; payload optional and approximated)
     *   3) AlphabetMapper (measured with a small JOL traversal)
     *   4) Estimators (measured with a small JOL traversal)
     *
     * Definitions (explained):
     * - Retained size: total bytes of heap objects reachable from a root. Here we avoid
     *   a full traversal for trees and instead compute an approximation from the known layout.
     * - Bloom filter size model:
     *      m_bits  = ceil( -(n * ln(p)) / (ln(2)^2) ), where n = nodes * min(Σ, interval)
     *      words   = ceil(m_bits / 64)
     *      bitset_bytes ≈ 16 (long[] header) + 8 * words
     *      bloom_obj_overhead ≈ 64 bytes (BloomFilter fields) + 24 bytes (BitSet shell) + 16 bytes (int[2] seeds) + 16 bytes (byte[8])
     *   These constants are ballpark values on a 64-bit Java Virtual Machine (JVM) with compressed ordinary object pointers.
     * - StreamBuffer payload (optional and approximate): we do NOT walk each Integer.
     *   We model it as:
     *      refs_array  ≈ 16 (Object[] header) + ref_size * size
     *      integers    ≈ size * approxIntegerBytes   (each Integer object + alignment)
     *   By default, we assume ref_size = VM.current().oopSize() (usually 4) and approxIntegerBytes = 16.
     *
     * IMPORTANT:
     * - This is designed to be QUICK, not exact. It will typically be within a reasonable factor of a full JOL walk,
     *   and it completely avoids the "hang" you saw on large windows.
     * - You do NOT need Serviceability Agent (SA). If you see JOL warnings, run with:
     *     -Djdk.attach.allowAttachSelf=true  -XX:+EnableDynamicAgentLoading  -Djol.skipSA=true
     *
     * @return four lines: Total (approx) / Trees (approx) / AlphabetMapper (JOL) / Estimators (JOL)
     * @return four lines: Total (approx) / Trees (approx) / AlphabetMapper (JOL) / Estimators (JOL)
     */

    public String fastMemoryEstimate(boolean includePayloadApprox, int approxIntegerBytes) {
        final double LN2 = Math.log(2.0);
        final int ALIGN  = VM.current().objectAlignment();    // often 8 or 16
        final int REF    = inferReferenceSizeFromVmDetails();

        // Tunable overhead constants (ballpark, kept conservative)
        final int BITSET_OBJ_OVERHEAD      = 24;  // shell object fields (BitSet itself)
        final int LONG_ARRAY_HDR           = 16;  // long[] header
        final int BLOOM_FILTER_OBJ_OVERHEAD= 64;  // BloomFilter fields (n,m,k,p,seeds refs, temp arrays, etc.)
        final int AL_LIST_OBJ_OVERHEAD     = 24;  // ArrayList object header/fields (for levelFilters and StreamBuffer)
        final int OBJ_ARRAY_HDR            = 16;  // Object[] header

        long treesBytes = 0L;

        // Sum across trees
        for (tree.ImplicitTree<membership.Membership> t : this.trees) {
            // ---- Bloom filters per level (constructed exactly like your factory) ----
            int levels = t.effectiveDepth();                     // total levels in this tree
            int root = t.effectiveRoot();
            for (int level = root; level < levels; level++) {
                int nodes    = 1 << level;
                int interval = t.intervalSize(level);      // equals (treeLength >> level)
                int perNode  = Math.min(this.alphabetSize, interval);
                long nDistinct = (long) nodes * (long) perNode;

                if (nDistinct <= 0) continue;

                double mBitsD = - (nDistinct * Math.log(this.fpRate)) / (LN2 * LN2);
                long   mBits  = (long) Math.ceil(mBitsD);
                long   words  = (mBits + 63L) >>> 6;       // ceil(mBits / 64)
                long   bitsetBytes = alignUp(LONG_ARRAY_HDR + words * 8L, ALIGN);
                long   bfBytes     = alignUp(BITSET_OBJ_OVERHEAD + bitsetBytes + BLOOM_FILTER_OBJ_OVERHEAD, ALIGN);

                treesBytes += bfBytes;
            }

            // ---- LevelDirectory (list of per-level filters) ----
            long levelDirRefs = alignUp(AL_LIST_OBJ_OVERHEAD, ALIGN)
                    + alignUp(OBJ_ARRAY_HDR + (long) levels * REF, ALIGN);
            treesBytes += levelDirRefs;

            // ---- Codec + small objects (very small; lumped into a constant per tree) ----
            final int SMALLS = 128;  // Key64, TreeLayout, a few small fields
            treesBytes += SMALLS;

            // ---- StreamBuffer (exclude heavy payload unless requested) ----
            int size = t.buffer.length();  // number of symbols in this tree buffer
            // Reference array (Object[]) to hold Integer refs
            long refsArray = alignUp(OBJ_ARRAY_HDR + (long) size * REF, ALIGN)
                    + alignUp(AL_LIST_OBJ_OVERHEAD, ALIGN);   // ArrayList shell
            treesBytes += refsArray;

            if (includePayloadApprox && size > 0) {
                // Each Integer object (very rough) — default ~16 bytes with compressed OOPs
                long integers = (long) size * approxIntegerBytes;
                treesBytes += integers;
            }
        }

        // ---- AlphabetMapper measured with a small JOL traversal ----
        long alphabetBytes = (this.alphabetMap != null)
                ? GraphLayout.parseInstance(this.alphabetMap).totalSize()
                : 0L;

        // ---- Estimators measured with a small JOL traversal ----
        long estimatorsBytes = 0L;
        if (!this.trees.isEmpty()) {
            java.util.ArrayList<Object> roots = new java.util.ArrayList<>(this.trees.size());
            for (tree.ImplicitTree<membership.Membership> t : this.trees) {
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

    /** Aligns size up to the nearest multiple of 'alignment' (which must be a power of two). */
    private static long alignUp(long size, int alignment) {
        long a = alignment;
        return (size + (a - 1)) & ~(a - 1);
    }

    /** Align size upwards to the nearest multiple of 'alignment' (power of two), with overflow safety. */
    private static long alignUpSafe(long size, int alignment) {
        if (size <= 0) return 0L;
        long a = alignment;
        long capped = (size > Long.MAX_VALUE - (a - 1)) ? Long.MAX_VALUE - (a - 1) : size;
        long aligned = (capped + (a - 1)) & ~(a - 1);
        return aligned < 0 ? Long.MAX_VALUE : aligned;
    }

    /** Safe addition that saturates at Long.MAX_VALUE (avoids wrap-around). */
    private static long safeAdd(long x, long y) {
        long r = x + y;
        if (((x ^ r) & (y ^ r)) < 0) return Long.MAX_VALUE;
        return r;
    }

    /** Infer reference size (4 with compressed ordinary object pointers, else 8) from JOL VM details. */
    private static int inferReferenceSizeFromVmDetails() {
        try {
            String d = VM.current().details().toLowerCase(java.util.Locale.ROOT);
            if (d.contains("compressed oops") || d.contains("compressed ordinary object pointers")) return 4;
        } catch (Throwable ignored) { }
        return 8;
    }

    /* ===== Reflection helpers to hide StreamBuffer payloads during JOL ===== */

    private static final class ArrayListSnapshot {
        final Object[] elementData;
        final int size;
        ArrayListSnapshot(Object[] elementData, int size) { this.elementData = elementData; this.size = size; }
    }

    private static java.lang.reflect.Field AL_ELEMENT_DATA;
    private static java.lang.reflect.Field AL_SIZE;

    private static void ensureArrayListFields() {
        if (AL_ELEMENT_DATA != null) return;
        try {
            AL_ELEMENT_DATA = java.util.ArrayList.class.getDeclaredField("elementData");
            AL_SIZE         = java.util.ArrayList.class.getDeclaredField("size");
            AL_ELEMENT_DATA.setAccessible(true);
            AL_SIZE.setAccessible(true);
        } catch (NoSuchFieldException e) {
            throw new RuntimeException("ArrayList internals not found; add --add-opens=java.base/java.util=ALL-UNNAMED", e);
        }
    }

    /** Swap out the internal array and size so the list appears empty to JOL (O(1), no copying). */
    private static ArrayListSnapshot stealAndEmptyArrayList(java.util.ArrayList<?> list) {
        ensureArrayListFields();
        try {
            Object[] current = (Object[]) AL_ELEMENT_DATA.get(list);
            int size = (int) AL_SIZE.get(list);
            AL_ELEMENT_DATA.set(list, new Object[0]);
            AL_SIZE.set(list, 0);
            return new ArrayListSnapshot(current, size);
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Failed to swap ArrayList internals; add --add-opens=java.base/java.util=ALL-UNNAMED", e);
        }
    }

    /** Restore a list that was emptied by stealAndEmptyArrayList. */
    private static void restoreArrayList(java.util.ArrayList<?> list, ArrayListSnapshot snapshot) {
        if (snapshot == null) return;
        ensureArrayListFields();
        try {
            AL_ELEMENT_DATA.set(list, snapshot.elementData);
            AL_SIZE.set(list, snapshot.size);
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Failed to restore ArrayList internals", e);
        }
    }
    /**
     * EXACT with Java Object Layout (JOL): sum of all membership filters across all trees.
     * This is fast because the number of objects is O(number of levels), not O(window size).
     */
    public long jolExactMembershipBytes() {
        java.util.ArrayList<Object> roots = new java.util.ArrayList<>();
        for (tree.ImplicitTree<membership.Membership> t : this.trees) {
            int levels = t.totalFilters();               // pass-through to LevelDirectory.depth()
            for (int L = 0; L < levels; L++) {
                membership.Membership m = t.filterAtLevel(L);
                if (m != null) roots.add(m);
            }
        }
        if (roots.isEmpty()) return 0L;
        return GraphLayout.parseInstance(roots.toArray()).totalSize();
    }
    /**
     * EXACT with Java Object Layout (JOL): retained size of the trees list,
     * after (1) detaching estimators and (2) making every StreamBuffer's ArrayList appear empty.
     * This excludes the heavy per-symbol payload but counts all index structures exactly.
     */
    public long jolExactTreesNoPayload() {
        // Detach estimators to avoid counting them here
        java.util.ArrayList<estimators.Estimator> saved = new java.util.ArrayList<>(this.trees.size());
        for (tree.ImplicitTree<membership.Membership> t : this.trees) {
            saved.add(t.estimator);
            t.estimator = null;
        }

        // Empty each buffer's ArrayList in O(1)
        java.util.ArrayList<ArrayListSnapshot> snaps = new java.util.ArrayList<>(this.trees.size());
        for (tree.ImplicitTree<membership.Membership> t : this.trees) {
            snaps.add(stealAndEmptyArrayList(t.buffer.data));
        }

        try {
            return GraphLayout.parseInstance(this.trees).totalSize();
        } finally {
            // Restore buffers and estimators
            for (int i = 0; i < this.trees.size(); i++) {
                restoreArrayList(this.trees.get(i).buffer.data, snaps.get(i));
            }
            for (int i = 0; i < this.trees.size(); i++) {
                if (this.trees.get(i).estimator == null) this.trees.get(i).estimator = saved.get(i);
            }
        }
    }
    /**
     * HYBRID report (fast, robust):
     *  - Membership filters: EXACT with JOL
     *  - Tree skeleton (LevelDirectory shells, codecs/layouts): APPROX (tiny)
     *  - StreamBuffer payload: APPROX (ArrayList shell + Object[] refs + per-Integer payload)
     *  - AlphabetMapper: EXACT with JOL
     *  - Estimators: EXACT with JOL
     *
     * Prints five lines: Total / TreeStructures / StreamBuffer / AlphabetMapper / Estimators
     */
    public String jolHybridMemoryReport(int approxIntegerBytes) {
        // Exact membership bytes
        long membershipBytes = jolExactMembershipBytes();

        // Approximate: the small tree structures around the filters
        final int ALIGN = VM.current().objectAlignment();
        final int REF   = inferReferenceSizeFromVmDetails();
        final int AL_LIST_OBJ_OVERHEAD = 24;
        final int OBJ_ARRAY_HDR        = 16;
        final int SMALLS_PER_TREE      = 128;

        long treeStructuresMinusMembership = 0L;
        long streamBufferBytes = 0L;

        for (tree.ImplicitTree<membership.Membership> t : this.trees) {
            int levels =  t.totalFilters();

            // LevelDirectory shells (ArrayList + Object[] of references)
            long levelList = alignUpSafe(AL_LIST_OBJ_OVERHEAD, ALIGN);
            long levelArr  = alignUpSafe(OBJ_ARRAY_HDR + (long) levels * REF, ALIGN);
            treeStructuresMinusMembership = safeAdd(treeStructuresMinusMembership, safeAdd(levelList, levelArr));

            // Small per-tree structures (codec/layout/etc.)
            treeStructuresMinusMembership = safeAdd(treeStructuresMinusMembership, SMALLS_PER_TREE);

            // StreamBuffer: full (structure + payload) so you can see data cost
            int size = t.buffer.length();
            long bufShell = alignUpSafe(AL_LIST_OBJ_OVERHEAD, ALIGN);
            long bufArr   = alignUpSafe(OBJ_ARRAY_HDR + (long) size * REF, ALIGN);
            long bufInts  = (size > 0) ? (long) size * (long) approxIntegerBytes : 0L;
            streamBufferBytes = safeAdd(streamBufferBytes, safeAdd(safeAdd(bufShell, bufArr), bufInts));
        }

        // Exact with JOL: AlphabetMapper and Estimators (both small graphs)
        long alphabetBytes = (this.alphabetMap != null)
                ? GraphLayout.parseInstance(this.alphabetMap).totalSize()
                : 0L;

        long estimatorsBytes = 0L;
        if (!this.trees.isEmpty()) {
            java.util.ArrayList<Object> roots = new java.util.ArrayList<>(this.trees.size());
            for (tree.ImplicitTree<membership.Membership> t : this.trees) {
                if (t.estimator != null) roots.add(t.estimator);
            }
            if (!roots.isEmpty()) estimatorsBytes = GraphLayout.parseInstance(roots.toArray()).totalSize();
        }

        long treeStructuresBytes = safeAdd(membershipBytes, treeStructuresMinusMembership);
        long total = safeAdd(treeStructuresBytes, safeAdd(streamBufferBytes, safeAdd(alphabetBytes, estimatorsBytes)));

        java.util.Locale L = java.util.Locale.ROOT;
        return String.format(L, "Total (hybrid): %d B (%.3f MiB)", total, total/(1024.0*1024.0)) + "\n"
                + String.format(L, "TreeStructures (exact membership + approx shells): %d B (%.3f MiB)", treeStructuresBytes, treeStructuresBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "  ├─ Membership (JOL exact): %d B (%.3f MiB)", membershipBytes, membershipBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "  └─ Non-membership shells (approx): %d B (%.3f MiB)", treeStructuresMinusMembership, treeStructuresMinusMembership/(1024.0*1024.0)) + "\n"
                + String.format(L, "StreamBuffer (approx): %d B (%.3f MiB)", streamBufferBytes, streamBufferBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "AlphabetMapper (JOL): %d B (%.3f MiB)", alphabetBytes, alphabetBytes/(1024.0*1024.0)) + "\n"
                + String.format(L, "Estimators (JOL): %d B (%.3f MiB)", estimatorsBytes, estimatorsBytes/(1024.0*1024.0));
    }




}
