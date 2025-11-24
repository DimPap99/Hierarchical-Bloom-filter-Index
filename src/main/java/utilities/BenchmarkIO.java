package utilities;

import utilities.BenchmarkEnums.QueryType;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

// IO helpers for datasets, queries, and pattern lengths.
public final class BenchmarkIO {
    private BenchmarkIO() {}

    public static List<Path> listDatasetDirectories(Path root) throws IOException {
        try (Stream<Path> stream = Files.list(root)) {
            return stream
                    .filter(Files::isDirectory)
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), BenchmarkIO::compareNaturally))
                    .collect(Collectors.toList());
        }
    }

    public static Path findSingleDatasetFile(Path datasetDir) throws IOException {
        try (Stream<Path> stream = Files.list(datasetDir)) {
            return stream
                    .filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().endsWith(".txt"))
                    .findFirst()
                    .orElseThrow(() -> new IllegalStateException(
                            "No dataset .txt file found under " + datasetDir));
        }
    }

    public static List<Path> findQueryFiles(Path queryDir, QueryType type) throws IOException {
        String suffix = "." + type.fileToken() + ".txt";
        try (Stream<Path> stream = Files.list(queryDir)) {
            return stream
                    .filter(Files::isRegularFile)
                    .filter(path -> path.getFileName().toString().endsWith(suffix))
                    .sorted(Comparator.comparing(path -> path.getFileName().toString(), BenchmarkIO::compareNaturally))
                    .collect(Collectors.toList());
        }
    }

    public static int parsePatternLength(String fileName) {
        int dot = fileName.indexOf('.');
        if (dot < 0) throw new IllegalArgumentException("Cannot determine pattern length from " + fileName);
        String numeric = fileName.substring(0, dot);
        try {
            return Integer.parseInt(numeric);
        } catch (NumberFormatException ex) {
            throw new IllegalArgumentException("Invalid pattern length in " + fileName, ex);
        }
    }

    public static int compareNaturally(String left, String right) {
        try {
            int leftInt = Integer.parseInt(left);
            int rightInt = Integer.parseInt(right);
            return Integer.compare(leftInt, rightInt);
        } catch (NumberFormatException ignored) {
            return left.compareTo(right);
        }
    }
}
