package utilities;

import PMIndex.IPMIndexing;
import search.Pattern;
import utilities.PatternResult;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public final class MultiQueryExperiment {

    private MultiQueryExperiment() {
    }

    public static MultiRunResult run(String datasetPath,
                                     List<QueryWorkload> workloads,
                                     IPMIndexing index,
                                     int nGram,
                                     boolean verbose,
                                     boolean queryResults) throws IOException {
        // Use populateIndex to ensure any index-specific finalization (e.g., suffix tree terminator) occurs
        InsertStats insertStats = populateIndex(datasetPath, index, nGram);
        return runQueries(workloads, index, nGram, insertStats, verbose, queryResults);
    }

    public static InsertStats populateIndex(String datasetPath,
                                            IPMIndexing index,
                                            int nGram) throws IOException {
        InsertStats stats = insertDataset(datasetPath, index, nGram);
        // Append a single terminator token for online suffix tree only, after measuring insert time.
        if (index instanceof PMIndex.SuffixTreeIndex) {
            try {
                index.insert("$");
            } catch (Exception ignored) {
                // Best-effort: do not fail the run if terminator insert throws
            }
        }
        return stats;
    }

    public static MultiRunResult runQueries(List<QueryWorkload> workloads,
                                            IPMIndexing index,
                                            int nGram,
                                            InsertStats insertStats,
                                            boolean verbose,
                                            boolean queryResults) throws IOException {
        Map<QueryWorkload, ExperimentRunResult> results = new LinkedHashMap<>();

        for (QueryWorkload workload : workloads) {
            QueryStats queryStats = executeQueries(workload.queryFile(), index, nGram, verbose, queryResults);

            double avgLpMs = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageLpTimeMillis() : 0.0;
            double avgCfLpMs = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageMinCostLpTimeMillis() : 0.0;
            double avgLpChosen = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageChosenLp() : 0.0;
            double avgCfLpChosen = (index instanceof PMIndex.HBI hbi)
                    ? hbi.stats().averageCfChosenLp() : 0.0;

            ExperimentRunResult runResult = new ExperimentRunResult(
                    queryStats.queryDurationMs(),
                    insertStats.insertDurationMs(),
                    queryStats.patternResults(),
                    queryStats.avgQueryLength(),
                    insertStats.avgInsertMsPerSymbol(),
                    queryStats.queryDurationMs(),
                    avgLpMs,
                    avgCfLpMs,
                    avgLpChosen,
                    avgCfLpChosen,
                    null);

            results.put(workload, runResult);
        }

        return new MultiRunResult(results);
    }

    private static InsertStats insertDataset(String datasetPath,
                                             IPMIndexing index,
                                             int nGram) throws IOException {
        CharRingBuffer window = new CharRingBuffer(nGram);
        long startTime = System.currentTimeMillis();
        long insertionCount = 0L;
        int totalSymbols = 0;

        try (BufferedReader reader = Files.newBufferedReader(Path.of(datasetPath), StandardCharsets.UTF_8)) {
            int ch;
            while ((ch = reader.read()) != -1) {
                char cChar = (char) ch;
                if (cChar == '\n' || cChar == '\r') {
                    continue;
                }
                window.append(cChar);
                totalSymbols++;
                if (window.isFilled()) {
                    index.insert(window.snapshot());
                    insertionCount++;
                }
            }
        }

        long insertDuration = System.currentTimeMillis() - startTime;
        double avgInsertMsPerSymbol = insertionCount == 0 ? 0.0 : insertDuration / (double) insertionCount;

        return new InsertStats(insertDuration, avgInsertMsPerSymbol);
    }

    private static QueryStats executeQueries(Path queryFile,
                                             IPMIndexing index,
                                             int nGram,
                                             boolean verbose,
                                             boolean collectResults) throws IOException {
        List<String> queries = readQueries(queryFile);
        long startTime = System.currentTimeMillis();
        ArrayList<PatternResult> patternResults = new ArrayList<>();
        int totalQueryLength = 0;
        int queryIndex = 0;

        for (String query : queries) {
            Pattern qPat = new Pattern(query, nGram);
            ArrayList<Integer> report = index.report(qPat);
            totalQueryLength += qPat.nGramToLong.length;
            if (collectResults) {
                patternResults.add(index.getLatestStats());
            }

            if (verbose) {
                if (report.size() < 20) {
                    if (query.length() < 80) {
                        System.out.println(query + ":" + report);
                    } else {
                        System.out.println("Query number " + queryIndex + ": " + report);
                    }
                } else {
                    System.out.println(query + ":" + report.size());
                }
            }
            queryIndex++;
        }

        double avgQueryLength = queries.isEmpty() ? 0.0 : totalQueryLength / (double) queries.size();

        long duration = System.currentTimeMillis() - startTime;
        return new QueryStats(duration, avgQueryLength, patternResults);
    }

    private static List<String> readQueries(Path queryFile) throws IOException {
        List<String> queries = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(queryFile, StandardCharsets.UTF_8)) {
            String line;
            while ((line = reader.readLine()) != null) {
                queries.add(line);
            }
        }
        return queries;
    }

    public record InsertStats(long insertDurationMs,
                              double avgInsertMsPerSymbol) {
    }

    public record QueryStats(long queryDurationMs,
                             double avgQueryLength,
                             ArrayList<PatternResult> patternResults) {
    }

    public record QueryWorkload(int patternLength, Path queryFile) {
    }

    public record MultiRunResult(Map<QueryWorkload, ExperimentRunResult> results) {
    }
}
