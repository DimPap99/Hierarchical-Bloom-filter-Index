import PMIndex.HBI;
import PMIndex.IPMIndexing;
import search.Pattern;
import utilities.CharRingBuffer;
import utilities.DatasetReader;
import utilities.ExperimentRunResult;
import utilities.PatternResult;
import utilities.RingBuffer;

import java.io.IOException;
import java.util.ArrayList;

public class Experiment {

    public static ExperimentRunResult run(String inputFilePath,
                                          String queriesFilePath,
                                          IPMIndexing index,
                                          int Ngram,
                                          boolean verbose,
                                          boolean queryResults) throws IOException {
        return run(inputFilePath, queriesFilePath, index, Ngram, verbose, queryResults, true);
    }

    public static ExperimentRunResult run(String inputFilePath,
                                          String queriesFilePath,
                                          IPMIndexing index,
                                          int Ngram,
                                          boolean verbose,
                                          boolean queryResults,
                                          boolean runQueries) throws IOException {
        ArrayList<PatternResult> queryResultsList = new ArrayList<>();
        long startTime = System.currentTimeMillis();
        RingBuffer<Character> window = new CharRingBuffer(Ngram);
        long insertDuration = 0L;
        long insertionCount = 0L;
        ArrayList<String> queries = new ArrayList<>();
        int totalSymbols = 0;
        int _i = 0;

        try (DatasetReader reader = new DatasetReader(inputFilePath, queriesFilePath, false)) {
            for (char cChar : reader) {
                window.append(cChar);
                totalSymbols++;

                if (window.isFilled()) {
                    String s = window.snapshot().toString();
                    if (_i == 50) {
                        int b = 2;
                    }
                    index.insert(s);
                    insertionCount++;
                }
                _i++;
            }

            insertDuration = System.currentTimeMillis() - startTime;
            double avgMsPerSymbol = totalSymbols == 0 ? 0.0 : insertDuration / (double) totalSymbols;
            System.out.println("Avg Insertion Time: " + avgMsPerSymbol);

            if (runQueries) {
                reader.setQueryMode();
                StringBuilder currentQuery = new StringBuilder();
                for (char ch : reader) {
                    if (ch == '\r') {
                        continue;
                    }
                    if (ch == '\n') {
                        queries.add(currentQuery.toString());
                        currentQuery.setLength(0);
                    } else {
                        currentQuery.append(ch);
                    }
                }
                if (currentQuery.length() > 0) {
                    queries.add(currentQuery.toString());
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        int avgQueryLength = 0;
        ArrayList<ArrayList<Integer>> matchRes = new ArrayList<>();
        long queryDuration = 0L;

        if (runQueries) {
            long queryStart = System.currentTimeMillis();
            int i = 0;
            for (String query : queries) {
                Pattern qPat = new Pattern(query, Ngram);
                ArrayList<Integer> report = index.report(qPat);
                avgQueryLength += qPat.nGramToLong.length;
                if (queryResults) {
                    PatternResult rr = index.getLatestStats();
                    queryResultsList.add(rr);
                }
                if (verbose) {
                    if (report.size() < 20) {
                        if (query.length() < 80) {
                            System.out.println(query + ":" + report);
                        } else {
                            System.out.println("Query number " + i + ": " + report);
                        }
                    } else {
                        System.out.println(query + ":" + report.size());
                    }
                }
                matchRes.add(report);
                i++;
            }
            if (!queries.isEmpty()) {
                avgQueryLength = avgQueryLength / queries.size();
            }
            queryDuration = System.currentTimeMillis() - queryStart;
            if (verbose) {
                System.out.println("Querying duration: " + queryDuration + " ms");
            }
        }

        double avgInsertMsPerSymbol = insertionCount == 0 ? 0.0 : (insertDuration / (double) insertionCount);
        double avgQueryLoadMs = runQueries ? (double) queryDuration : 0.0;
        double avgLpMs = (index instanceof HBI hbi) ? hbi.stats().averageLpTimeMillis() : 0.0;
        double avgCfLpMs = (index instanceof HBI hbi) ? hbi.stats().averageMinCostLpTimeMillis() : 0.0;
        return new ExperimentRunResult(queryDuration, insertDuration, queryResultsList, avgQueryLength, avgInsertMsPerSymbol, avgQueryLoadMs, avgLpMs, avgCfLpMs, matchRes);
    }
}
