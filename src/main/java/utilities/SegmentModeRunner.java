package utilities;

import PMIndex.IPMIndexing;
import search.Pattern;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

// Helpers to run benchmarks in segments mode: one token (record) per line in the dataset, and
// space-separated tokens per line in the query files.
//
// Files are read before the timers start, so insert and query times cover the index work only:
// building and hashing n-grams, inserting, and matching.
public final class SegmentModeRunner {

    private SegmentModeRunner() {}

    // Build the index in segments mode by sliding an n-gram window over the tokens.
    public static MultiQueryExperiment.InsertStats insertDatasetSegments(String datasetPath,
                                                                        IPMIndexing index,
                                                                        int nGram) throws Exception {
        List<String> tokens = readTokens(Path.of(datasetPath));
        RingBuffer<String> window = new RingBuffer<>(nGram);
        long insertionEvents = 0L;

        long startNanos = System.nanoTime();
        for (String tok : tokens) {
            window.append(tok);
            if (!window.isFilled()) continue;
            index.insert(window.join(Pattern.SEGMENT_JOINER));
            insertionEvents++;
        }
        double durationMs = (System.nanoTime() - startNanos) / 1e6;
        double avgInsertMsPerSymbol = insertionEvents == 0 ? 0.0 : (durationMs / insertionEvents);
        MultiQueryExperiment.InsertStats stats = new MultiQueryExperiment.InsertStats(durationMs, avgInsertMsPerSymbol);
        // Append a single terminator token for the suffix tree only, after measuring insert time.
        if (index instanceof PMIndex.SuffixTreeIndex) {
            try {
                index.insert("$");
            } catch (Exception ignored) {
                // best-effort
            }
        }
        return stats;
    }

    // Execute queries from files in segments mode. datasetPath is kept for the callers' signature.
    public static MultiQueryExperiment.MultiRunResult runQueriesSegments(String datasetPath,
                                                                         List<MultiQueryExperiment.QueryWorkload> workloads,
                                                                         IPMIndexing index,
                                                                         int nGram,
                                                                         MultiQueryExperiment.InsertStats insertStats,
                                                                         boolean verbose,
                                                                         boolean collectResults) throws IOException {

        Map<MultiQueryExperiment.QueryWorkload, ExperimentRunResult> results = new LinkedHashMap<>();

        for (MultiQueryExperiment.QueryWorkload workload : workloads) {
            // HBI statistics below are per workload, not running totals.
            if (index instanceof PMIndex.HBI hbi) {
                hbi.stats().resetCounters();
            }
            List<String[]> queries = new ArrayList<>();
            for (String line : readTokens(workload.queryFile())) {
                queries.add(line.split(" +"));
            }
            ArrayList<PatternResult> patternResults = new ArrayList<>();
            List<ArrayList<Integer>> answers = ResultVerifier.enabled() ? new ArrayList<>(queries.size()) : null;
            int totalQueryLen = 0;
            int queryIndex = 0;

            long startNanos = System.nanoTime();
            for (String[] toks : queries) {
                Pattern qPat = new Pattern(toks, nGram);
                ArrayList<Integer> report = index.report(qPat);
                totalQueryLen += qPat.nGramToLong.length;
                if (answers != null) answers.add(report);
                if (collectResults) patternResults.add(index.getLatestStats());
                if (verbose) {
                    String line = String.join(" ", toks);
                    if (report.size() < 20) {
                        System.out.println((line.length() < 80) ? (line + ":" + report) : ("Query number " + queryIndex + ": " + report));
                    } else {
                        System.out.println(line + ":" + report.size());
                    }
                }
                queryIndex++;
            }
            double durationMs = (System.nanoTime() - startNanos) / 1e6;
            if (answers != null) {
                ResultVerifier.check(index, nGram, workload.queryFile(), answers);
            }

            double avgQueryLen = queries.isEmpty() ? 0.0 : (totalQueryLen / (double) queries.size());
            double avgLpMs = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageLpTimeMillis() : 0.0;
            double avgCfLpMs = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageMinCostLpTimeMillis() : 0.0;
            double avgLpChosen = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageChosenLp() : 0.0;
            double avgCfLpChosen = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageCfChosenLp() : 0.0;
            ExperimentRunResult runResult = new ExperimentRunResult(
                    durationMs,
                    insertStats.insertDurationMs(),
                    patternResults,
                    avgQueryLen,
                    insertStats.avgInsertMsPerSymbol(),
                    durationMs,
                    avgLpMs,
                    avgCfLpMs,
                    avgLpChosen,
                    avgCfLpChosen,
                    null);

            results.put(workload, runResult);
        }

        return new MultiQueryExperiment.MultiRunResult(results);
    }

    // Non-empty lines of a file; a trailing '\r' (CRLF files) is not part of the token.
    private static List<String> readTokens(Path file) throws IOException {
        List<String> out = new ArrayList<>();
        try (BufferedReader in = Files.newBufferedReader(file, StandardCharsets.UTF_8)) {
            String line;
            while ((line = in.readLine()) != null) {
                if (!line.isEmpty()) {
                    out.add(line);
                }
            }
        }
        return out;
    }
}
