package utilities;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.util.ArrayList;

public class FileHelper {

    public void CreateCsvFromRunResult(String headers, ArrayList<RunResult> results, String filePath) throws IOException {
        StringBuilder builder = new StringBuilder();
        builder.append(headers);
        for(RunResult r : results) builder.append(r.toString());

        Files.write(Paths.get(filePath),
                builder.toString().getBytes(), StandardOpenOption.CREATE, StandardOpenOption.APPEND);
    }
}
