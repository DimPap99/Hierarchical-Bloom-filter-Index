package utilities;

import java.util.ArrayList;

public record ExperimentRunResult (
        double totalRunTimeMs,              // total query time (ms) across the run (or sum of query loads)
        double totalInsertTimeMs,           // total insertion time (ms)
        ArrayList<PatternResult> patternResults,
        double avgQuerySize,                // average query length (ngrams)
        double avgInsertMsPerSymbol,        // average time per insertion symbol (ms)
        double avgQueryLoadMs               // average time per query load (ms); for single-load experiments equals totalRunTimeMs
) {
}
