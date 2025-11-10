package utilities;

import java.util.ArrayList;

public record ExperimentRunResult (
        double totalRunTimeMs,              // total query time (ms) across the run (or sum of query loads)
        double totalInsertTimeMs,           // total insertion time (ms)
        ArrayList<PatternResult> patternResults,
        double avgQuerySize,                // average query length (ngrams)
        double avgInsertMsPerSymbol,        // average time per insertion symbol (ms)
        double avgQueryLoadMs ,             // average time per query load (ms); for single-load experiments equals totalRunTimeMs
        double avgLpMs,                     // average Lp estimation time per pattern (ms); 0 if unavailable
        double avgCfLpMs,                   // average CF-only Lp time per pattern (ms); 0 if unavailable
        double avgLpChosen,                 // average chosen Lp (current algo); 0 if unavailable
        double avgCfLpChosen,               // average chosen CF Lp (min-cost); 0 if unavailable
        ArrayList<ArrayList<Integer>> matchRes

) {
}
