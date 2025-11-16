import datagenerators.AdversarialGenerators;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Command-line helper that materialises the adversarial workloads discussed in the
 * evaluation doc. Running this class recreates both datasets under data/adversarial.
 */
public final class GenerateAdversarialDatasets {

    private static final int WINDOW_POWER = 21;
    private static final int WINDOW_LENGTH = 1 << WINDOW_POWER; // derived window size
    private static final String WINDOW_TOKEN = "w" + WINDOW_POWER;

    private static final int BLOCK_LEVEL = 12;
    private static final int BLOCK_LENGTH = 1 << BLOCK_LEVEL;
    private static final String BLOCK_ALPHABET = "ABCDEFGHIKLMNOP";
    private static final String BLOCK_CATEGORY = "blocks";
    private static final String BLOCK_BASE_NAME = "blocks_" + WINDOW_TOKEN;
    private static final String BLOCK_SINGLE_FILE = BLOCK_BASE_NAME + ".txt";
    private static final String BLOCK_FAMILY = BLOCK_BASE_NAME;

    private static final int DEBRUIJN_ORDER = 6;
    private static final int DEBRUIJN_BASE = 48; // treat alphabet base as raw code point (0-255-friendly)
    private static final int MAX_CHAR = 68;
    private static final int DEBRUIJN_SIGMA = 4;//MAX_CHAR - DEBRUIJN_BASE + 1;
    private static final char[] DEBRUIJN_ALPHABET = buildAlphabet();
    private static final String DEBRUIJN_CATEGORY = "debruijn";
    private static final String DEBRUIJN_BASE_NAME = "debruijn_k" + DEBRUIJN_ORDER + "_sigma" + DEBRUIJN_SIGMA + "_" + WINDOW_TOKEN;
    private static final String DEBRUIJN_SINGLE_FILE = DEBRUIJN_BASE_NAME + ".txt";
    private static final String DEBRUIJN_FAMILY = DEBRUIJN_BASE_NAME;
    private static final String DATASET_INDEX = "1";

    private GenerateAdversarialDatasets() {
    }

    public static void main(String[] args) throws IOException {
        Path adversarialRoot = Paths.get("data", "adversarial");
        Path adversarialQueryRoot = Paths.get("queries", "adversarial");
        Path blockDataRoot = adversarialRoot.resolve(BLOCK_CATEGORY);
        Path debruijnDataRoot = adversarialRoot.resolve(DEBRUIJN_CATEGORY);
        Path blockQueryRoot = adversarialQueryRoot.resolve(BLOCK_CATEGORY);
        Path debruijnQueryRoot = adversarialQueryRoot.resolve(DEBRUIJN_CATEGORY);

        Files.createDirectories(blockDataRoot);
        Files.createDirectories(debruijnDataRoot);
        Files.createDirectories(blockQueryRoot);
        Files.createDirectories(debruijnQueryRoot);

        writeBlockDataset(blockDataRoot, blockQueryRoot);
        writeDeBruijnDataset(debruijnDataRoot, debruijnQueryRoot);
    }

    private static void writeBlockDataset(Path dataRoot, Path queryRoot) throws IOException {
        char[] alphabet = BLOCK_ALPHABET.toCharArray();
        String data = AdversarialGenerators.generateAlternatingBlocks(WINDOW_LENGTH, BLOCK_LENGTH, alphabet);
        Path single = dataRoot.resolve(BLOCK_SINGLE_FILE);
        Files.writeString(single, data, StandardCharsets.UTF_8);
        System.out.printf("Wrote block-occupancy dataset (%d chars, %d-char blocks) to %s%n",
                WINDOW_LENGTH, BLOCK_LENGTH, single.toAbsolutePath());
        writeBenchmarkDataset(dataRoot, queryRoot, BLOCK_FAMILY, data,
                "Block-occupancy dataset (HBIDatasetBenchmarkMulti layout)");
    }

    private static void writeDeBruijnDataset(Path dataRoot, Path queryRoot) throws IOException {
        String data = AdversarialGenerators.generateDeBruijnSequence(DEBRUIJN_ALPHABET,
                DEBRUIJN_ORDER,
                WINDOW_LENGTH);
        Path single = dataRoot.resolve(DEBRUIJN_SINGLE_FILE);
        Files.writeString(single, data, StandardCharsets.UTF_8);
        System.out.printf("Wrote de Bruijn dataset (order=%d, sigma=%d) to %s%n",
                DEBRUIJN_ORDER, DEBRUIJN_ALPHABET.length, single.toAbsolutePath());
        writeBenchmarkDataset(dataRoot, queryRoot, DEBRUIJN_FAMILY, data,
                "De Bruijn dataset (HBIDatasetBenchmarkMulti layout)");
    }

    private static char[] buildAlphabet() {
        char[] chars = new char[DEBRUIJN_SIGMA];
        for (int i = 0; i < DEBRUIJN_SIGMA; i++) {
            chars[i] = (char) (DEBRUIJN_BASE + i);
        }
        return chars;
    }

    private static void writeBenchmarkDataset(Path dataRoot,
                                              Path queryRoot,
                                              String family,
                                              String data,
                                              String logLabel) throws IOException {
        Path datasetDir = dataRoot.resolve(family).resolve(DATASET_INDEX);
        Files.createDirectories(datasetDir);
        Path datasetFile = datasetDir.resolve(DATASET_INDEX + "_" + family + ".txt");
        Files.writeString(datasetFile, data, StandardCharsets.UTF_8);
        // Ensure the matching query directory exists so HBIDatasetBenchmarkMulti picks up workloads later.
        Path queryDir = queryRoot.resolve(family).resolve(DATASET_INDEX);
        Files.createDirectories(queryDir);
        System.out.printf("%s -> %s%n", logLabel, datasetFile.toAbsolutePath());
    }
}
