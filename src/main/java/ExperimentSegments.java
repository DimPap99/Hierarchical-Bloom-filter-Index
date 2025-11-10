
import PMIndex.HBI;
import PMIndex.IPMIndexing;
import search.Pattern;
import utilities.*;

import java.io.IOException;
import java.util.ArrayList;

public class ExperimentSegments {

    /**
     * Run an experiment identical in outcome to Experiment.run, but using SegmentReader.
     *
     * @param inputFilePath     path to the dataset file
     * @param queriesFilePath   path to the queries file
     * @param index             the index implementing IPMIndexing
     * @param Ngram             the N-gram length (window size for emitted tokens)
     * @param verbose           whether to print per-query debug outputs
     * @param queryResults      whether to collect PatternResult statistics from the index
     * @return ExperimentRunResult containing query and insert durations, optional per-query results, and average query length
     * @throws IOException on I/O problems opening or reading the files
     */
    public static ExperimentRunResult run(String inputFilePath,
                                          String queriesFilePath,
                                          IPMIndexing index,
                                          int Ngram,
                                          boolean verbose,
                                          boolean queryResults) throws IOException {

        ArrayList<PatternResult> queryResultsList = new ArrayList<>();
        long startTime = System.currentTimeMillis();

        // Sliding window over characters to emit N-grams
        RingBuffer<String> window = new RingBuffer<String>(Ngram);

        long insertDuration;
        long insertionCount = 0L;
        ArrayList<String> queries = new ArrayList<>();
        int totalSymbols = 0;

        // -----------------------
        // Dataset ingestion phase
        // -----------------------
        try (SegmentReader reader = new SegmentReader(
                inputFilePath,
                queriesFilePath,
                '\n',   // delimiter
                false,  // do not include '\n' in returned segments
                true    // skip empty segments
        )) {
            // Stream the dataset in newline-delimited segments; push each character
            // into the ring buffer to produce N-grams exactly as in the char-by-char loop.
            for (String segment : reader) {

                window.append(segment);
                totalSymbols++;

                // Only once we have k chars -> emit the N-gram
                if (window.isFilled()) {
                    index.insert(window.snapshot().toString()); // the whole k-gram
                    insertionCount++;
                }


                // If you want to treat each newline as a real character in the stream,
                // append it here; by default we skip it to match the prior behavior
                // where line breaks were not included during dataset ingestion.
                // Example (disabled by default):
                // char newline = '\n';
                // window.append(newline);
                // totalSymbols++;
                // if (window.isFilled()) {
                //     index.insert(window.snapshot().toString());
                // }
            }

            insertDuration = System.currentTimeMillis() - startTime;
            double avgMsPerSymbol = (totalSymbols == 0) ? 0.0 : insertDuration / (double) totalSymbols;

            // ----------------
            // Querying phase
            // ----------------
            // Switch the same reader instance to the queries file; it will reopen internally.
            reader.setQueryMode();

            // Since SegmentReader normalizes CRLF when delimiter is '\n', there is no need
            // to manually strip '\r' or assemble lines; each segment is already a logical line.
            for (String q : reader) {
                // Skip empty lines already handled by skipEmptySegments=true.
                // If you need to keep empty queries, construct SegmentReader with skipEmptySegments=false.
                queries.add(q);
            }

            // avgMsPerSymbol is computed but not used in the return value here,
            // mirroring the original Experiment class behavior (you can log it if desired).
            // System.out.println("Avg ms per symbol: " + avgMsPerSymbol);

        } catch (Exception e) {
            throw new RuntimeException(e);
        }

        // ------------------------
        // Execute queries and time
        // ------------------------
        startTime = System.currentTimeMillis();
        int avgQueryLength = 0;
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
            i++;
        }

        if (!queries.isEmpty()) {
            avgQueryLength = avgQueryLength / queries.size();
        }

        long endTime = System.currentTimeMillis();
        long queryDuration = endTime - startTime;

        if (verbose) {
            System.out.println("Querying duration: " + queryDuration + " ms");
        }

        double avgInsertMsPerSymbol = insertionCount == 0 ? 0.0 : (insertDuration / (double) insertionCount);
        double avgQueryLoadMs = (double) queryDuration; // single load
        double avgLpMs = (index instanceof HBI hbi) ? hbi.stats().averageLpTimeMillis() : 0.0;
        double avgCfLpMs = (index instanceof HBI hbi) ? hbi.stats().averageMinCostLpTimeMillis() : 0.0;
        return new ExperimentRunResult(queryDuration, insertDuration, queryResultsList, avgQueryLength, avgInsertMsPerSymbol, avgQueryLoadMs, avgLpMs, avgCfLpMs, null);
    }
}
