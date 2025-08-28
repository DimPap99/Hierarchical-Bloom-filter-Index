import datagenerators.Generator;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;

public class GenerateDatasets {
    public static void main(String[] args) throws IOException {
        System.out.println("Generating datasets...");
        //uniform distributions
//        String data = Generator.generateUniform(1 << 16, 48, 122);
//        Files.write(Paths.get("data/uniform_text_16_experiment_small.txt"),
//        data.getBytes(StandardCharsets.UTF_8));
//
//        data = Generator.generateUniform(1 << 21, 48, 122);
//        Files.write(Paths.get("data/uniform_text_21_experiment.txt"),
//        data.getBytes(StandardCharsets.UTF_8));
//        //zipf distributions
//        data = Generator.generateZipf(1 << 16, 48, 122, 1);
//        Files.write(Paths.get("data/zipf_16_1.txt"),
//                data.getBytes(StandardCharsets.UTF_8));
//
//        data = Generator.generateZipf(1 << 21, 48, 122, 1);
//        Files.write(Paths.get("data/zipf_21_1.txt"),
//        data.getBytes(StandardCharsets.UTF_8));
//
//        data = Generator.generateZipf(1 << 21, 48, 122, 1.25);
//        Files.write(Paths.get("data/zipf_21_1_25.txt"),
//        data.getBytes(StandardCharsets.UTF_8));
//
//        data = Generator.generateZipf(1 << 21, 48, 122, 0.75);
//        Files.write(Paths.get("data/zipf_21_0_75.txt"),
//        data.getBytes(StandardCharsets.UTF_8));
    }
}
