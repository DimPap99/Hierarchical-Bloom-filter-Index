import PMIndex.IPMIndexing;
import datagenerators.Generator;
import search.Pattern;
import utilities.CharRingBuffer;
import utilities.HBILogger;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Deque;

public class Experiment {

    public static double run(String inputFilePath, String queriesFilePath, IPMIndexing index, int Ngram, boolean verbose) throws IOException {
//        HBILogger.info("Running experiment for Index: " + index.getClass().getSimpleName());
        int read;
        // Read characters
//        HBILogger.info("Reading input file...");
//        String data = Generator.generateZipf(131072, 48, 122, 1.35);
//                Files.write(Paths.get("zipf_text2.txt"),
//                data.getBytes(StandardCharsets.UTF_8));   // no newline
        long startTime = System.currentTimeMillis();
        int c = 0;
        CharRingBuffer window = new CharRingBuffer(Ngram);
// ------------------------------------------------------------------

        try (FileReader fr = new FileReader(inputFilePath)) {
            int ch;
            while ((ch = fr.read()) != -1) {
                char cChar = (char) ch;
                window.append(cChar);

                /* Only once we have k chars -> emit the N-gram */
                if (window.isFilled()) {
                    index.insert(window.view());   // the whole k-gram
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
        long duration = endTime - startTime;
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
        for(String query : queries){
//            HBILogger.info("Query: " + query);
            ArrayList<Integer> report = index.report(new Pattern(query, Ngram));
            if(verbose){
                if(report.size() < 20){
                    System.out.println(query + ":" + report);
                }else{
                    System.out.println(query + ":" + report.size());
                }
            }


        }
        endTime = System.currentTimeMillis();
        duration = endTime - startTime;
//        HBILogger.info("Report: " + report.toString());
        if(verbose){
            System.out.println("Querying duration: " +  duration + " ms");

        }
        return duration;
    }
}
