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
            // HBI statistics below are per workload, not running totals.
            if (index instanceof PMIndex.HBI hbi) {
                hbi.stats().resetCounters();
            }
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

    // The window is read into memory first, so the timer covers building n-grams and inserting them,
    // not reading the file.
    private static InsertStats insertDataset(String datasetPath,
                                             IPMIndexing index,
                                             int nGram) throws IOException {
        char[] text = readText(datasetPath);
        CharRingBuffer window = new CharRingBuffer(nGram);
        long insertionCount = 0L;

        long startNanos = System.nanoTime();
        for (char c : text) {
            window.append(c);
            if (window.isFilled()) {
                index.insert(window.snapshot());
                insertionCount++;
            }
        }
        double insertDurationMs = (System.nanoTime() - startNanos) / 1e6;
        double avgInsertMsPerSymbol = insertionCount == 0 ? 0.0 : insertDurationMs / insertionCount;

        return new InsertStats(insertDurationMs, avgInsertMsPerSymbol);
    }

    /** The dataset as the stream of characters the experiments insert: line breaks are dropped. */
    private static char[] readText(String datasetPath) throws IOException {
        String raw = Files.readString(Path.of(datasetPath), StandardCharsets.UTF_8);
        char[] text = new char[raw.length()];
        int n = 0;
        for (int i = 0; i < raw.length(); i++) {
            char c = raw.charAt(i);
            if (c != '\n' && c != '\r') {
                text[n++] = c;
            }
        }
        return (n == text.length) ? text : java.util.Arrays.copyOf(text, n);
    }

    private static QueryStats executeQueries(Path queryFile,
                                             IPMIndexing index,
                                             int nGram,
                                             boolean verbose,
                                             boolean collectResults) throws IOException {
        List<String> queries = readQueries(queryFile);
        ArrayList<PatternResult> patternResults = new ArrayList<>();
        List<ArrayList<Integer>> answers = ResultVerifier.enabled() ? new ArrayList<>(queries.size()) : null;
        int totalQueryLength = 0;
        int queryIndex = 0;

        long startNanos = System.nanoTime();
        for (String query : queries) {
            Pattern qPat = new Pattern(query, nGram);
            ArrayList<Integer> report = index.report(qPat);
            totalQueryLength += qPat.nGramToLong.length;
            if (answers != null) {
                answers.add(report);
            }
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

        double durationMs = (System.nanoTime() - startNanos) / 1e6;
        double avgQueryLength = queries.isEmpty() ? 0.0 : totalQueryLength / (double) queries.size();
        if (answers != null) {
            ResultVerifier.check(index, nGram, queryFile, answers);
        }
        return new QueryStats(durationMs, avgQueryLength, patternResults);
    }

    private static List<String> readQueries(Path queryFile) throws IOException {
        List<String> queries = new ArrayList<>();
        try (BufferedReader reader = Files.newBufferedReader(queryFile, StandardCharsets.UTF_8)) {
            String line;
            while ((line = reader.readLine()) != null) {
                if (!line.isEmpty()) {
                    queries.add(line);
                }
            }
        }
        return queries;
    }

    public record InsertStats(double insertDurationMs,
                              double avgInsertMsPerSymbol) {
    }

    public record QueryStats(double queryDurationMs,
                             double avgQueryLength,
                             ArrayList<PatternResult> patternResults) {
    }

    public record QueryWorkload(int patternLength, Path queryFile) {
    }

    public record MultiRunResult(Map<QueryWorkload, ExperimentRunResult> results) {
    }
}
