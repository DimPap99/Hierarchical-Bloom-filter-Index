import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.RegexIndex;
import estimators.CostFunctionMarkov;
import estimators.Estimator;
import estimators.HashMapEstimator;
import membership.BloomFilter;
import membership.Membership;
import search.*;
import utilities.*;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.function.Supplier;

public class ProcessStream {


    public static <T> ExperimentRunResult run(String inputFilePath,
                                              String queriesFilePath,
                                              IPMIndexing index,
                                              Reader<T> reader,
                                              int Ngram,
                                              int windowLen,
                                              boolean verbose,
                                              boolean queryResults, String readerMode) throws IOException {
        // 1) Load queries once
        List<String> queries = loadQueries(queriesFilePath);

        // Average query length (in n-grams)
        int sumQueryLen = 0;
        for (String q : queries) {
            Pattern p = new Pattern(q, Ngram);
            sumQueryLen += p.nGramToLong.length;
        }
        double avgQueryLen = queries.isEmpty() ? 0.0 : (sumQueryLen * 1.0 / queries.size());

        // 2) Stream the dataset and schedule full query workloads
        RingBuffer<T> window = new RingBuffer<>(Ngram);
        long totalInsertNanos = 0L;
        long insertionEvents = 0L;   // number of index.insert calls (ngrams inserted)
        long tokensInserted  = 0L;   // same as insertionEvents in this pipeline

        long nextQueryAt = Math.max(1, windowLen); // schedule queries every window length
        long sumQueryLoadMs = 0L;
        long queryLoads = 0L;

        while (reader.hasNext()) {
            T token = reader.next();
            window.append(token);
            if (!window.isFilled()) {
                continue;
            }
            long t0 = System.nanoTime();
            String s = window.snapshot().toString();
            if(s.equals("tcp:smtp:c->s:L0:R:T0")) {
                int b=2;
                System.out.println(insertionEvents);
            }

            index.insert(s);
            long t1 = System.nanoTime();
            totalInsertNanos += (t1 - t0);
            insertionEvents++;
            tokensInserted++;

            if (tokensInserted >= nextQueryAt) {
                long qMs = runQueryLoad(index, queries, Ngram, verbose, queryResults, readerMode);
                queryLoads++;
                sumQueryLoadMs += qMs;
                nextQueryAt += windowLen;
            }
        }

        double totalInsertMs = totalInsertNanos / 1_000_000.0;
        double avgInsertMsPerSymbol = insertionEvents == 0 ? 0.0 : (totalInsertMs / insertionEvents);
        double totalQueryMs = sumQueryLoadMs;
        double avgQueryLoadMs = queryLoads == 0 ? 0.0 : (totalQueryMs / queryLoads);

        return new ExperimentRunResult(totalQueryMs, totalInsertMs, new ArrayList<>(), avgQueryLen, avgInsertMsPerSymbol, avgQueryLoadMs, null);
    }



    public static void main(String[] args) throws Exception {
        // Paths and parameters moved here from ConfidenceExperiment-style configuration
        // Pull defaults from ConfidenceExperiment style constants
        String dataFile    = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data/packets_tokens_only.txt";
        String queriesFile = "/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries/mawi/mawi_1.txt";
        int textSizePow2   = 21;                   // 2^20 window and tree
        int subWindowLen   = 21;
        int windowLen      = 1 << textSizePow2;    // window length in symbols
        int treeLen        = 1 << subWindowLen;    // tree capacity in symbols
        double fpRate      = 0.005;                // Bloom FP rate
        int ngrams         = 1;                    // n-gram size
        int alphabet       = 502;                   // base alphabet (single symbol); expanded below
        double conf        = 0.99;                 // pruning confidence for CF plan
        String readerMode  = "segments2";          // "segments" or "chars"

        // Allow simple overrides via args: --data, --queries, --window, --tree, --ngrams, --alpha, --fpr, --conf
        for (int i = 0; i < args.length - 1; i += 2) {
            switch (args[i]) {
                case "--data" -> dataFile = args[i + 1];
                case "--queries" -> queriesFile = args[i + 1];
                case "--window" -> windowLen = Integer.parseInt(args[i + 1]);
                case "--tree" -> treeLen = Integer.parseInt(args[i + 1]);
                case "--ngrams" -> ngrams = Integer.parseInt(args[i + 1]);
                case "--alpha" -> alphabet = Integer.parseInt(args[i + 1]);
                case "--fpr" -> fpRate = Double.parseDouble(args[i + 1]);
                case "--conf" -> conf = Double.parseDouble(args[i + 1]);
                case "--reader" -> readerMode = args[i + 1];
                default -> { /* ignore unknown */ }
            }
        }

        // Expand alphabet to n-gram space
        int expandedAlphabet = (int) Math.pow(alphabet, ngrams);

        // Build index
        HBI hbi = newHbi(windowLen, treeLen, expandedAlphabet, fpRate, conf, ngrams);
        hbi.stats().setCollecting(false);
        hbi.stats().setExperimentMode(false);

        // Pick the reader type requested
        ExperimentRunResult res;
        if ("segments".equalsIgnoreCase(readerMode)) {
            SegmentReader seg = new SegmentReader(
                    dataFile,
                    queriesFile,
                    '\n',   // delimiter per line
                    false,  // do not include delimiter in the returned segment
                    true    // skip empty segments
            );
            SegmentReaderAdapter adapter = new SegmentReaderAdapter(seg);
            res = run(dataFile, queriesFile, hbi, adapter, ngrams, windowLen, false, false, readerMode);
        } else {
            IPMIndexing ipm = new RegexIndex();
            DatasetReader ds = new DatasetReader(dataFile, queriesFile, false);
            DatasetReaderAdapter adapter = new DatasetReaderAdapter(ds);
            res = run(dataFile, queriesFile, ipm, adapter, ngrams, windowLen, false, false, readerMode);
        }

        System.out.printf(Locale.ROOT,
                "Stream complete. insertMs=%.3f, avgInsertMs=%.6f, totalQueryMs=%.3f, avgQueryLoadMs=%.3f, avgQueryLen=%.2f\n",
                res.totalInsertTimeMs(), res.avgInsertMsPerSymbol(), res.totalRunTimeMs(), res.avgQueryLoadMs(), res.avgQuerySize());
    }

    private static long runQueryLoad(IPMIndexing index, List<String> queries, int ngram, boolean verbose, boolean collectResults, String readerMode) {
        long start = System.currentTimeMillis();
        ArrayList<Integer> report;
        for (String q : queries) {
            if (q == null || q.isEmpty()) continue;
            Pattern pat;
            if (readerMode == "segments") {
                String[] splitArray = q.split(" ");

                pat = new Pattern(splitArray, ngram);
            }else{
                pat = new  Pattern(q, ngram);
            }

            report = index.report(pat);

            if (verbose) {
                if (report.size() < 20 && q.length() < 80) System.out.println(q + ":" + report);
            }
            if (collectResults) {
                index.getLatestStats(); // touch stats if experiment mode is enabled
            }
        }
        return System.currentTimeMillis() - start;
    }

    private static List<String> loadQueries(String queriesFile) throws IOException {
        List<String> lines = Files.readAllLines(Path.of(queriesFile), StandardCharsets.UTF_8);
        ArrayList<String> out = new ArrayList<>(lines.size());
        for (String l : lines) {
            if (l == null) continue;
            String s = l.strip();
            if (!s.isEmpty()) out.add(s);
        }
        return out;
    }

    // Helper to construct HBI similar to ConfidenceExperiment/HBIDatasetBenchmark
    private static HBI newHbi(int windowLen, int treeLen, int alphabet, double fpRate, double conf, int ngrams) {
        Supplier<Estimator> estFactory = () -> new HashMapEstimator(treeLen);
        Supplier<Membership> memFactory = BloomFilter::new;
        Supplier<PruningPlan> prFactory = () -> new MostFreqPruning(conf);
        Verifier verifier = new VerifierLinearLeafProbe();
        return new HBI(
                new BlockSearch(),
                windowLen,
                fpRate,
                alphabet,
                treeLen,
                estFactory,
                memFactory,
                prFactory,
                verifier,
                new CostFunctionMarkov(),
                conf,
                ngrams
        );
    }
}
