package utilities;

import PMIndex.IPMIndexing;
import search.Pattern;
import tree.ssws.SuffixTree;

import java.io.IOException;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Helpers to run benchmarks in "segments" mode (one token per line).
 * Extracted from HBIDatasetBenchmarkMulti for reuse and modularity.
 */
public final class SegmentModeRunner {

    private SegmentModeRunner() {}

    /**
     * Build the index in segments mode by sliding an n-gram window over tokenized lines.
     */
    public static MultiQueryExperiment.InsertStats insertDatasetSegments(String datasetPath,
                                                                        IPMIndexing index,
                                                                        int nGram) throws Exception {
        RingBuffer<String> window = new RingBuffer<>(nGram);
        long startWallMs = System.currentTimeMillis();
        long insertionEvents = 0L;

        try (SegmentReader segReader = new SegmentReader(
                datasetPath,
                datasetPath,
                '\n',
                false,   // don't include delimiter
                true     // skip empty lines
        )) {
            for (String tok : segReader) {
                if (tok == null || tok.isEmpty()) continue;
                window.append(tok);
                if (!window.isFilled()) continue;
                String gram = window.snapshot().toString();
                index.insert(gram);
                insertionEvents++;
            }
        }

        long durationMs = System.currentTimeMillis() - startWallMs;
        double avgInsertMsPerSymbol = insertionEvents == 0 ? 0.0 : (durationMs / (double) insertionEvents);
        return new MultiQueryExperiment.InsertStats(durationMs, avgInsertMsPerSymbol);
    }

    /**
     * Execute queries from files in segments mode.
     */
    public static MultiQueryExperiment.MultiRunResult runQueriesSegments(String datasetPath,
                                                                         List<MultiQueryExperiment.QueryWorkload> workloads,
                                                                         IPMIndexing index,
                                                                         int nGram,
                                                                         MultiQueryExperiment.InsertStats insertStats,
                                                                         boolean verbose,
                                                                         boolean collectResults) throws IOException {

        Map<MultiQueryExperiment.QueryWorkload, ExperimentRunResult> results = new LinkedHashMap<>();

        for (MultiQueryExperiment.QueryWorkload workload : workloads) {
            long startTimeMs = System.currentTimeMillis();
            ArrayList<PatternResult> patternResults = new ArrayList<>();
            int totalQueryLen = 0;
            int queryCount = 0;
            int queryIndex = 0;

            try (SegmentReader segReader = new SegmentReader(
                    datasetPath,
                    workload.queryFile().toString(),
                    '\n',
                    false,
                    true
            )) {
                segReader.setQueryMode();
                for (String line : segReader) {
                    if (line == null || line.isEmpty()) continue;
                    String[] toks = line.split(" +");
                    Pattern qPat = new Pattern(toks, nGram);
                    ArrayList<Integer> report = index.report(qPat);
                    totalQueryLen += qPat.nGramToLong.length;
                    queryCount++;
                    if (collectResults) patternResults.add(index.getLatestStats());
                    if (verbose) {
                        if (report.size() < 20) {
                            System.out.println((line.length() < 80) ? (line + ":" + report) : ("Query number " + queryIndex + ": " + report));
                        } else {
                            System.out.println(line + ":" + report.size());
                        }
                    }
                    queryIndex++;
                }
            } catch (Exception e) {
                throw new RuntimeException(e);
            }

            long durationMs = System.currentTimeMillis() - startTimeMs;
            double avgQueryLen = (queryCount == 0) ? 0.0 : (totalQueryLen / (double) queryCount);

            ExperimentRunResult runResult = new ExperimentRunResult(
                    durationMs,
                    insertStats.insertDurationMs(),
                    patternResults,
                    avgQueryLen,
                    insertStats.avgInsertMsPerSymbol(),
                    durationMs,
                    null);

            results.put(workload, runResult);
        }

        return new MultiQueryExperiment.MultiRunResult(results);
    }
}

