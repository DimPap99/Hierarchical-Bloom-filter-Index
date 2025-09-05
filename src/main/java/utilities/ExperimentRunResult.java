package utilities;

import java.util.ArrayList;

public record ExperimentRunResult (double totalRunTimeMs, double totalInsertTimeMs, ArrayList<PatternResult> patternResults) {
}