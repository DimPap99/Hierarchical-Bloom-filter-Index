public final class HBIDatasetBenchmarkMulti {
    public static void main(String[] args) throws Exception {
        // Parse CLI options and delegate to the orchestrator that drives
        // FPR x ngram x dataset loops and reporting.
        utilities.MultiBenchmarkOptions options = utilities.MultiBenchmarkOptions.parse(args);
        utilities.BenchmarkOrchestrator.run(options);
    }
}
