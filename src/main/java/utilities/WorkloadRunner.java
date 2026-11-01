package utilities;

import PMIndex.HBI;
import PMIndex.IPMIndexing;
import PMIndex.SuffixTreeIndex;
import search.Pattern;
import utilities.Workload.ReaderMode;

import java.nio.file.Path;
import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

// Runs a Workload against an index. This is the only place that inserts a dataset and runs queries,
// so the experiment classes only choose indexes and parameters.
//
//   insert(...)      fill an index with the whole dataset                     -> InsertStats
//   query(...)       run the workload's query file once                       -> QueryStats
//   run(...)         insert, then query once                                  -> ExperimentRunResult
//   runQueries(...)  query several files against one filled index            -> MultiRunResult
//   stream(...)      insert, and run every query after each window of symbols  -> ExperimentRunResult
//
// Files are read into memory before the timers start, so insert and query times cover the index
// work only: building n-gram keys, inserting, building patterns and matching.
public final class WorkloadRunner {

    private WorkloadRunner() {}

    // The queries of one pattern length (one query file).
    public record QuerySet(int patternLength, Path queryFile) {
    }

    public record InsertStats(double insertDurationMs, double avgInsertMsPerSymbol) {
    }

    public record QueryStats(double queryDurationMs,
                             double avgQueryLength,
                             ArrayList<PatternResult> patternResults,
                             ArrayList<ArrayList<Integer>> matches) {
    }

    public record MultiRunResult(Map<QuerySet, ExperimentRunResult> results) {
    }

    // What to print and keep while querying. Matches are needed to compare indexes; per-pattern stats
    // (index.getLatestStats()) for the estimator experiments.
    public record Options(boolean verbose, boolean collectStats, boolean collectMatches) {
        public static final Options QUIET = new Options(false, false, false);

        public Options verbose(boolean v) { return new Options(v, collectStats, collectMatches); }
        public Options collectStats(boolean v) { return new Options(verbose, v, collectMatches); }
        public Options collectMatches(boolean v) { return new Options(verbose, collectStats, v); }
    }

    // ---------------------------------------------------------------- insert

    // Inserts every n-gram of the dataset. The online suffix tree gets its "$" terminator afterwards,
    // outside the timer.
    public static InsertStats insert(Workload workload, IPMIndexing index, int ngram) {
        InsertStats stats = (workload.mode() == ReaderMode.CHARS)
                ? insertChars(workload.readChars(), index, ngram)
                : insertTokens(workload.readTokens(), index, ngram);
        if (index instanceof SuffixTreeIndex) {
            try {
                index.insert("$");
            } catch (Exception ignored) {
                // best effort: a missing terminator only affects the suffix tree baseline
            }
        }
        return stats;
    }

    private static InsertStats insertChars(char[] text, IPMIndexing index, int ngram) {
        CharRingBuffer window = new CharRingBuffer(ngram);
        long insertions = 0L;
        long startNanos = System.nanoTime();
        for (char c : text) {
            window.append(c);
            if (window.isFilled()) {
                index.insert(window.snapshot());
                insertions++;
            }
        }
        return insertStats(System.nanoTime() - startNanos, insertions);
    }

    private static InsertStats insertTokens(List<String> tokens, IPMIndexing index, int ngram) {
        RingBuffer<String> window = new RingBuffer<>(ngram);
        long insertions = 0L;
        long startNanos = System.nanoTime();
        for (String tok : tokens) {
            window.append(tok);
            if (window.isFilled()) {
                index.insert(window.join(Pattern.SEGMENT_JOINER));
                insertions++;
            }
        }
        return insertStats(System.nanoTime() - startNanos, insertions);
    }

    private static InsertStats insertStats(long nanos, long insertions) {
        double ms = nanos / 1e6;
        return new InsertStats(ms, insertions == 0 ? 0.0 : ms / insertions);
    }

    // ---------------------------------------------------------------- query

    // Runs every query of the workload once. With --verify on, the answers are checked against the
    // other indexes' answers for the same query file.
    public static QueryStats query(Workload workload, IPMIndexing index, int ngram, Options options) {
        List<Workload.Query> queries = workload.parsedQueries();
        ArrayList<PatternResult> patternResults = new ArrayList<>();
        boolean verify = ResultVerifier.enabled();
        ArrayList<ArrayList<Integer>> matches =
                (options.collectMatches() || verify) ? new ArrayList<>(queries.size()) : null;

        long[] loadNanos = {0L, 0L};
        runQueryLoad(queries, index, ngram, options, patternResults, matches, loadNanos);
        double durationMs = loadNanos[0] / 1e6;
        double avgQueryLength = averageQueryLength(queries, ngram);

        if (verify) {
            ResultVerifier.check(index, ngram, workload.queriesFile(), matches);
        }
        return new QueryStats(durationMs, avgQueryLength, patternResults,
                options.collectMatches() ? matches : null);
    }

    // Insert the whole dataset, then run the queries once.
    public static ExperimentRunResult run(Workload workload, IPMIndexing index, int ngram, Options options) {
        InsertStats insert = insert(workload, index, ngram);
        if (!workload.hasQueries()) {
            return result(index, insert, new QueryStats(0.0, 0.0, new ArrayList<>(), null), 0.0);
        }
        QueryStats q = query(workload, index, ngram, options);
        return result(index, insert, q, q.queryDurationMs());
    }

    // Runs each query set against an index already filled with dataset's data. HBI statistics are
    // reset before each set, so they describe that set only.
    public static MultiRunResult runQueries(Workload dataset,
                                            List<QuerySet> sets,
                                            IPMIndexing index,
                                            int ngram,
                                            InsertStats insert,
                                            Options options) {
        Map<QuerySet, ExperimentRunResult> results = new LinkedHashMap<>();
        for (QuerySet set : sets) {
            if (index instanceof HBI hbi) {
                hbi.stats().resetCounters();
            }
            QueryStats q = query(dataset.withQueries(set.queryFile()), index, ngram, options);
            results.put(set, result(index, insert, q, q.queryDurationMs()));
        }
        return new MultiRunResult(results);
    }

    // ---------------------------------------------------------------- stream

    // Streams the dataset into the index and runs the whole query file each time another
    // windowLength symbols have arrived, i.e. every time the window has been fully replaced.
    // totalRunTimeMs is the sum over all query loads and avgQueryLoadMs the mean per load.
    public static ExperimentRunResult stream(Workload workload,
                                             IPMIndexing index,
                                             int ngram,
                                             int windowLength,
                                             Options options) {
        List<Workload.Query> queries = workload.hasQueries() ? workload.parsedQueries() : List.of();
        ArrayList<PatternResult> patternResults = new ArrayList<>();
        ArrayList<ArrayList<Integer>> matches = options.collectMatches() ? new ArrayList<>() : null;
        long[] queryLoadNanos = {0L, 0L}; // total, count
        long period = Math.max(1, windowLength);

        RingBuffer<Object> window = new RingBuffer<>(ngram);
        long insertNanos = 0L;
        long insertions = 0L;
        long symbols = 0L;

        Iterable<?> stream = (workload.mode() == ReaderMode.CHARS)
                ? asCharacters(workload.readChars())
                : workload.readTokens();

        long chunkStart = System.nanoTime();
        for (Object symbol : stream) {
            window.append(symbol);
            symbols++;
            if (window.isFilled()) {
                index.insert(workload.key(window));
                insertions++;
            }
            if (!queries.isEmpty() && symbols % period == 0) {
                insertNanos += System.nanoTime() - chunkStart;
                runQueryLoad(queries, index, ngram, options, patternResults, matches, queryLoadNanos);
                chunkStart = System.nanoTime();
            }
        }
        insertNanos += System.nanoTime() - chunkStart;

        InsertStats insert = insertStats(insertNanos, insertions);
        double totalQueryMs = queryLoadNanos[0] / 1e6;
        double avgLoadMs = queryLoadNanos[1] == 0 ? 0.0 : totalQueryMs / queryLoadNanos[1];
        QueryStats q = new QueryStats(totalQueryMs, averageQueryLength(queries, ngram),
                patternResults, matches);
        return result(index, insert, q, avgLoadMs);
    }

    // ---------------------------------------------------------------- helpers

    // One pass over all queries. Adds the elapsed time to loadNanos[0] and counts the load in
    // loadNanos[1].
    private static void runQueryLoad(List<Workload.Query> queries,
                                     IPMIndexing index,
                                     int ngram,
                                     Options options,
                                     List<PatternResult> patternResults,
                                     List<ArrayList<Integer>> matches,
                                     long[] loadNanos) {
        long start = System.nanoTime();
        int i = 0;
        for (Workload.Query query : queries) {
            ArrayList<Integer> report = index.report(query.toPattern(ngram));
            if (matches != null) {
                matches.add(report);
            }
            if (options.collectStats()) {
                PatternResult stats = index.getLatestStats();
                if (stats != null) {
                    patternResults.add(stats);
                }
            }
            if (options.verbose()) {
                printReport(query.text(), i, report);
            }
            i++;
        }
        loadNanos[0] += System.nanoTime() - start;
        loadNanos[1]++;
    }

    // Mean number of n-grams per query.
    private static double averageQueryLength(List<Workload.Query> queries, int ngram) {
        if (queries.isEmpty()) return 0.0;
        long total = 0;
        for (Workload.Query q : queries) {
            total += q.toPattern(ngram).nGramToLong.length;
        }
        return total / (double) queries.size();
    }

    private static void printReport(String query, int i, List<Integer> report) {
        if (report.size() < 20) {
            System.out.println(query.length() < 80 ? (query + ":" + report) : ("Query number " + i + ": " + report));
        } else {
            System.out.println(query + ":" + report.size());
        }
    }

    private static ExperimentRunResult result(IPMIndexing index, InsertStats insert, QueryStats q, double avgQueryLoadMs) {
        HBI hbi = (index instanceof HBI h) ? h : null;
        return new ExperimentRunResult(
                q.queryDurationMs(),
                insert.insertDurationMs(),
                q.patternResults(),
                q.avgQueryLength(),
                insert.avgInsertMsPerSymbol(),
                avgQueryLoadMs,
                hbi == null ? 0.0 : hbi.stats().averageLpTimeMillis(),
                hbi == null ? 0.0 : hbi.stats().averageMinCostLpTimeMillis(),
                hbi == null ? 0.0 : hbi.stats().averageChosenLp(),
                hbi == null ? 0.0 : hbi.stats().averageCfChosenLp(),
                q.matches());
    }

    private static Iterable<Character> asCharacters(char[] text) {
        return () -> new java.util.Iterator<>() {
            private int i = 0;
            @Override public boolean hasNext() { return i < text.length; }
            @Override public Character next() { return text[i++]; }
        };
    }
}
