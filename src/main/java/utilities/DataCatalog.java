package utilities;

import utilities.BenchmarkEnums.QueryType;
import utilities.Workload.ReaderMode;

import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

// Finds datasets and query files by short name, so a test never needs a full path.
//
// The data root is, in order: the path given to from(...), the HBI_DATA_ROOT environment variable, or
// the first "data" folder found in the working directory or one of its parents. The query root is the
// "queries" folder next to the data root unless given. Two layouts are understood:
//
//   flat     data/<name>.txt                + queries/<name>/<length>.txt          (pg2701, 1023)
//   sampled  data/<name>/<part>/<any>.txt   + queries/<name>/<part>/<length>.<type>.txt
//            (w19, w20, w21, wzipf_*, mawi*: <part> is one of the sampled windows, 1 to 5)
//
// Anything else (an adversarial file, a dataset whose query folder has another name) goes through
// files(...), which resolves relative paths under the two roots.
public final class DataCatalog {

    public static final String DATA_ROOT_ENV = "HBI_DATA_ROOT";

    private final Path dataRoot;
    private final Path queryRoot;

    private DataCatalog(Path dataRoot, Path queryRoot) {
        this.dataRoot = dataRoot.toAbsolutePath().normalize();
        this.queryRoot = queryRoot.toAbsolutePath().normalize();
    }

    // Default roots (HBI_DATA_ROOT, else the nearest data/ folder).
    public static DataCatalog fromDefaults() {
        return from(null, null);
    }

    // dataRoot and queryRoot may be null to use the defaults described above.
    public static DataCatalog from(Path dataRoot, Path queryRoot) {
        Path data = (dataRoot != null) ? dataRoot : defaultDataRoot();
        if (!Files.isDirectory(data)) {
            throw new IllegalArgumentException("Data root not found: " + data.toAbsolutePath()
                    + " (pass --data-root or set " + DATA_ROOT_ENV + ")");
        }
        Path queries = (queryRoot != null) ? queryRoot : data.toAbsolutePath().normalize().resolveSibling("queries");
        return new DataCatalog(data, queries);
    }

    public Path dataRoot() { return dataRoot; }
    public Path queryRoot() { return queryRoot; }

    // ---------------------------------------------------------------- layouts

    // data/<name>.txt with queries/<name>/<patternLength>.txt
    public Workload flat(String name, int patternLength, ReaderMode mode) {
        Path data = requireFile(dataRoot.resolve(name + ".txt"));
        Path queries = requireFile(queryRoot.resolve(name).resolve(patternLength + ".txt"));
        return new Workload(name + "/" + patternLength, data, queries, mode);
    }

    // data/<name>/<part>/*.txt with queries/<name>/<part>/<patternLength>.<type>.txt
    public Workload sampled(String name, int part, int patternLength, QueryType type, ReaderMode mode) {
        Path data = sampledDataFile(name, part);
        Path queries = requireFile(queryRoot.resolve(name).resolve(Integer.toString(part))
                .resolve(patternLength + "." + type.fileToken() + ".txt"));
        String label = name + "/" + part + "/" + patternLength + "." + type.fileToken();
        return new Workload(label, data, queries, mode);
    }

    // Explicit files. Relative paths resolve under the data root and the query root; absolute paths
    // are used as given. queriesPath may be null for a stream-only workload.
    public Workload files(String dataPath, String queriesPath, ReaderMode mode) {
        Path data = requireFile(dataRoot.resolve(dataPath));
        Path queries = (queriesPath == null) ? null : requireFile(queryRoot.resolve(queriesPath));
        return new Workload(null, data, queries, mode);
    }

    // Picks the layout from what is on disk: sampled when data/<name>/ is a folder, flat when
    // data/<name>.txt exists. part and type are only used by the sampled layout.
    public Workload byName(String name, int part, int patternLength, QueryType type, ReaderMode mode) {
        if (Files.isDirectory(dataRoot.resolve(name))) {
            return sampled(name, part, patternLength, type, mode);
        }
        if (Files.isRegularFile(dataRoot.resolve(name + ".txt"))) {
            return flat(name, patternLength, mode);
        }
        throw new IllegalArgumentException("No dataset '" + name + "' under " + dataRoot
                + " (expected " + name + ".txt or a folder " + name + "/)");
    }

    // ---------------------------------------------------------------- listings

    // The numbered parts of a sampled dataset, in numeric order.
    public List<Integer> parts(String name) {
        List<Integer> out = new ArrayList<>();
        try {
            for (Path dir : BenchmarkIO.listDatasetDirectories(requireDir(dataRoot.resolve(name)))) {
                try {
                    out.add(Integer.parseInt(dir.getFileName().toString()));
                } catch (NumberFormatException ignored) {
                    // not a part folder
                }
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        return out;
    }

    // Pattern lengths that have a query file of the given type for this part, ascending.
    public List<Integer> patternLengths(String name, int part, QueryType type) {
        Path dir = requireDir(queryRoot.resolve(name).resolve(Integer.toString(part)));
        List<Integer> out = new ArrayList<>();
        try {
            for (Path file : BenchmarkIO.findQueryFiles(dir, type)) {
                out.add(BenchmarkIO.parsePatternLength(file.getFileName().toString()));
            }
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
        return out;
    }

    // ---------------------------------------------------------------- command line

    // The data-selection flags shared by the experiment drivers. A driver passes every --key value pair
    // to accept(...) first and handles the keys it returns false for itself, then calls resolve().
    //
    //   --dataset <name>           w19, mawi21, pg2701, ... (default: the driver's)
    //   --part <n>                 sampled window 1..5 of a sampled dataset (default 1)
    //   --pattern-length <n>       which query file (default 10)
    //   --query-type <t>           uniform | rare | missing, sampled datasets only (default uniform)
    //   --reader chars|segments    how to read the files (default: segments for mawi*, chars otherwise);
    //                              --mode is accepted as the old name
    //   --data <path>, --queries <path>   explicit files instead of --dataset; relative to the roots
    //   --data-root <dir>, --query-root <dir>
    public static final class Selection {
        private String dataset;
        private int part = 1;
        private int patternLength = 10;
        private QueryType queryType = QueryType.UNIFORM;
        private ReaderMode mode;
        private String dataPath;
        private String queriesPath;
        private Path dataRoot;
        private Path queryRoot;

        public Selection(String defaultDataset) {
            this.dataset = defaultDataset;
        }

        // Driver defaults, overridden by the flags.
        public Selection part(int defaultPart) { this.part = defaultPart; return this; }
        public Selection patternLength(int defaultLength) { this.patternLength = defaultLength; return this; }
        public Selection queryType(QueryType defaultType) { this.queryType = defaultType; return this; }
        public Selection reader(ReaderMode defaultMode) { this.mode = defaultMode; return this; }

        public boolean accept(String key, String value) {
            switch (key) {
                case "dataset" -> dataset = value;
                case "part" -> part = Integer.parseInt(value);
                case "pattern-length" -> patternLength = Integer.parseInt(value);
                case "query-type" -> queryType = QueryType.fromString(value);
                case "reader", "mode" -> mode = ReaderMode.parse(value);
                case "data" -> dataPath = value;
                case "queries" -> queriesPath = value;
                case "data-root" -> dataRoot = Path.of(value);
                case "query-root" -> queryRoot = Path.of(value);
                default -> { return false; }
            }
            return true;
        }

        public String dataset() {
            return (dataPath != null) ? Path.of(dataPath).getFileName().toString() : dataset;
        }

        public ReaderMode mode() {
            if (mode != null) return mode;
            String name = dataset();
            return (name != null && name.toLowerCase(Locale.ROOT).startsWith("mawi"))
                    ? ReaderMode.SEGMENTS : ReaderMode.CHARS;
        }

        // Window length suggested by the dataset name ("w19" and "mawi19" give 2^19), else def.
        public int windowLength(int def) {
            return MultiBenchmarkOptions.deriveWindowLength(dataset(), def);
        }

        public Workload resolve() {
            DataCatalog catalog = DataCatalog.from(dataRoot, queryRoot);
            if (dataPath != null) {
                return catalog.files(dataPath, queriesPath, mode());
            }
            if (queriesPath != null) {
                throw new IllegalArgumentException("--queries needs --data");
            }
            return catalog.byName(dataset, part, patternLength, queryType, mode());
        }
    }

    // ---------------------------------------------------------------- helpers

    private Path sampledDataFile(String name, int part) {
        Path dir = requireDir(dataRoot.resolve(name).resolve(Integer.toString(part)));
        try {
            return BenchmarkIO.findSingleDatasetFile(dir);
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    private static Path defaultDataRoot() {
        String env = System.getenv(DATA_ROOT_ENV);
        if (env != null && !env.isBlank()) {
            return Path.of(env);
        }
        Path dir = Path.of("").toAbsolutePath();
        while (dir != null) {
            Path candidate = dir.resolve("data");
            if (Files.isDirectory(candidate)) {
                return candidate;
            }
            dir = dir.getParent();
        }
        return Path.of("data");
    }

    private static Path requireFile(Path path) {
        if (!Files.isRegularFile(path)) {
            throw new IllegalArgumentException("File not found: " + path);
        }
        return path;
    }

    private static Path requireDir(Path path) {
        if (!Files.isDirectory(path)) {
            throw new IllegalArgumentException("Folder not found: " + path);
        }
        return path;
    }

    @Override
    public String toString() {
        return String.format(Locale.ROOT, "data=%s queries=%s", dataRoot, queryRoot);
    }
}
