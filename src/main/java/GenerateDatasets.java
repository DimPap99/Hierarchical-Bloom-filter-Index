import datagenerators.Generator;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

// Generate Zipf-distributed synthetic datasets under data/<family>/<i>/<i>_<family>.txt.
public class GenerateDatasets {

    public static void main(String[] args) throws IOException {
        System.out.println("Generating datasets...");

        // Configuration
        int numDatasets = 1;
        int alphabetSize = 10;
        double exponent = 1;

        // Use a family name without a trailing slash. This becomes both a directory and part of the file name.
        String family = "wzipf_21_adv" + alphabetSize;

        // Root "data" directory of your repo. Keep it absolute if you prefer, or make it relative.
        Path dataRoot = Paths.get("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/data");
        Path queryRoot = Paths.get("/home/dimpap/Desktop/GraduationProject/Hierarchical-Bloom-filter-Index/Hierarchical-Bloom-filter-Index/queries");

        // Generate one Zipf string and reuse it for all N datasets. If you want different content per dataset,
        // move this call inside the loop and vary the seed or parameters.
        int length = 1 << 21;             // total number of characters
        int alphabetStart = 40;           // inclusive start character code
        int alphabetEnd = alphabetStart + alphabetSize; // generator-specific convention; keep as in your original

        for (int i = 1; i <= numDatasets; i++) {
            String data = Generator.generateZipf(length, alphabetStart, alphabetEnd, exponent, i);

            // Directory: <dataRoot>/<family>/<i>
            Path dir = dataRoot.resolve(family).resolve(Integer.toString(i));
            Path qdir = queryRoot.resolve(family).resolve(Integer.toString(i));
            // File name: <i>_<family>.txt
            String fileName = i + "_" + family + ".txt";

            // Full output path
            Path outputPath = dir.resolve(fileName);

            // Ensure the directory exists
            Files.createDirectories(dir);
            Files.createDirectories(qdir);

            // Write the data
            Files.writeString(outputPath, data, StandardCharsets.UTF_8);

            System.out.println("Wrote dataset " + i + " to: " + outputPath);
        }

        System.out.println("Done.");
    }
}
