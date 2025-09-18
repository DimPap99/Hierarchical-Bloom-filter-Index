import PMIndex.IPMIndexing;
import datagenerators.Generator;
import search.Pattern;
import utilities.*;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Deque;

public class Experiment {

    public static ExperimentRunResult run(String inputFilePath, String queriesFilePath, IPMIndexing index, int Ngram, boolean verbose, boolean queryResults) throws IOException {
//        HBILogger.info("Running experiment for Index: " + index.getClass().getSimpleName());
        int read;
        ArrayList<PatternResult> queryResultsList = new ArrayList<>();
        // Read characters

//        HBILogger.info("Reading input file...");
//        String data = Generator.generateUniform(1 << 16, 48, 122);
//                Files.write(Paths.get("uniform_text_big_16.txt"),
//                data.getBytes(StandardCharsets.UTF_8));   // no newline
        long startTime = System.currentTimeMillis();
        int c = 0;
        RingBuffer<Character> window = new CharRingBuffer(Ngram);
// ------------------------------------------------------------------

        try (FileReader fr = new FileReader(inputFilePath)) {
            int ch;
            while ((ch = fr.read()) != -1) {
                char cChar = (char) ch;
                window.append(cChar);

                /* Only once we have k chars -> emit the N-gram */
                if (window.isFilled()) {
                    index.insert(window.snapshot().toString());   // the whole k-gram
                }
            }
        }
        catch (IOException e){
            e.printStackTrace();
        }


//        for(char cc : data.toCharArray()){
//            index.insert(cc);
//        }
        long endTime = System.currentTimeMillis();
        long insertDuration = endTime - startTime;
        //timings.add(duration);
//        HBILogger.info("Report: " + report.toString());
//        HBILogger.info("Time taken: " + duration + " ms");
//
//        HBILogger.info("Reading queries file...");
        ArrayList<String> queries = new ArrayList<>();
        try{
            BufferedReader reader = new BufferedReader(new FileReader(queriesFilePath));
            String line;
            while((line = reader.readLine()) != null){
                queries.add(line);
            }
        }catch (IOException e){
            e.printStackTrace();
        }

//        HBILogger.info("Will start querying...");
        startTime = System.currentTimeMillis();
        int avgQueryLength = 0;
        int i =0;
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
        avgQueryLength = avgQueryLength / queries.size();
        endTime = System.currentTimeMillis();
        long queryDuration = endTime - startTime;
//        HBILogger.info("Report: " + report.toString());
        if(verbose){
            System.out.println("Querying duration: " +  queryDuration + " ms");

        }
        return new ExperimentRunResult(queryDuration,insertDuration, queryResultsList, avgQueryLength);
    }
}
