package utilities;

import java.util.ArrayList;

public record ExperimentRunResultWithMatches(
        ExperimentRunResult result,
        ArrayList<String> queries,
        ArrayList<ArrayList<Integer>> matchPositions
) {}
