import PMIndex.IPMIndexing;
import search.Pattern;
import utilities.*;

import java.io.IOException;
import java.util.ArrayList;

public class Experiment {

    public static ExperimentRunResult run(String inputFilePath, String queriesFilePath, IPMIndexing index, int Ngram, boolean verbose, boolean queryResults) throws IOException {
//        HBILogger.info("Running experiment for Index: " + index.getClass().getSimpleName());
        ArrayList<PatternResult> queryResultsList = new ArrayList<>();
        long startTime = System.currentTimeMillis();
        RingBuffer<Character> window = new CharRingBuffer(Ngram);
        long insertDuration = 0L;
        long insertionCount = 0L;
        ArrayList<String> queries = new ArrayList<>();
        int totalSymbols = 0;
        try (DatasetReader reader = new DatasetReader(inputFilePath, queriesFilePath, false)) {
            for (char cChar : reader) {
                window.append(cChar);
                totalSymbols++;

                /* Only once we have k chars -> emit the N-gram */
                if (window.isFilled()) {
                    index.insert(window.snapshot().toString());   // the whole k-gram
                    insertionCount++;
                }
            }

            insertDuration = System.currentTimeMillis() - startTime;
            double avgMsPerSymbol = totalSymbols == 0 ? 0.0 : insertDuration / (double) totalSymbols;
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
        } catch (Exception e) {
            throw new RuntimeException(e);
        }

//        HBILogger.info("Will start querying...");
        startTime = System.currentTimeMillis();
        int avgQueryLength = 0;
        int i = 0;
        for(String query : queries){
            Pattern qPat = new Pattern(query, Ngram);

            //            HBILogger.info("Query: " + query);
            ArrayList<Integer> report = index.report(qPat);
            avgQueryLength += qPat.nGramToInt.length;
            if(queryResults){
                PatternResult rr = index.getLatestStats();
                queryResultsList.add(rr);
            }
//            ArrayList<Long> s = index.getAvgTimes(new Pattern(query, Ngram));
            if(verbose){
                if(report.size() < 20){
                    if(query.length() < 80){System.out.println(query + ":" + report);}
                    else{ System.out.println("Query number " + i + ": " + report);
                    }
                }else{
                    System.out.println(query + ":" + report.size());
                }
            }
            i++;
        }
        if(!queries.isEmpty()){
            avgQueryLength = avgQueryLength / queries.size();
        }
        long endTime = System.currentTimeMillis();
        long queryDuration = endTime - startTime;
//        HBILogger.info("Report: " + report.toString());
        if(verbose){
            System.out.println("Querying duration: " +  queryDuration + " ms");

        }
        double avgInsertMsPerSymbol = insertionCount == 0 ? 0.0 : (insertDuration / (double) insertionCount);
        double avgQueryLoadMs = (double) queryDuration; // single load
        return new ExperimentRunResult(queryDuration, insertDuration, queryResultsList, avgQueryLength, avgInsertMsPerSymbol, avgQueryLoadMs);
    }
}
