import PMIndex.IPMIndexing;
import utilities.HBILogger;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Experiment {

    public static void run(String inputFilePath, String queriesFilePath, IPMIndexing index) {
        HBILogger.info("Running experiment for Index: " + index.getClass().getSimpleName());
        int read;
        // Read characters
        HBILogger.info("Reading input file...");
        long startTime = System.currentTimeMillis();
        int c = 0;
        try{
            FileReader fileReader = new FileReader(inputFilePath);
            int ch;
            while ((ch = fileReader.read()) != -1) {
                System.out.println(c);
                index.insert(Character.toString((char) ch));
                c++;
            }
            fileReader.close();
        HBILogger.info("Done reading input file.");
        }catch (IOException e){
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis();
        long duration = endTime - startTime;
//        HBILogger.info("Report: " + report.toString());
        HBILogger.info("Time taken: " + duration + " ms");

        HBILogger.info("Reading queries file...");
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

        HBILogger.info("Will start querying...");
        startTime = System.currentTimeMillis();
        for(String query : queries){
//            HBILogger.info("Query: " + query);
            ArrayList<Integer> report = index.report(query);
            if(report.size() < 20){
                System.out.println(report);
            }else{
                System.out.println(query + ":" + report.size());
            }

        }
        endTime = System.currentTimeMillis();
        duration = endTime - startTime;
//        HBILogger.info("Report: " + report.toString());
        HBILogger.info("Time taken: " + duration + " ms");
        HBILogger.info("---------------------------------------");
        HBILogger.info("---------------------------------------");
    }
}
