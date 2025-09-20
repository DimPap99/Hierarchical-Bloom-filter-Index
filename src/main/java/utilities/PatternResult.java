package utilities;

import search.Pattern;

public record PatternResult (double totalRunTime, int probes, int Lp, Pattern p, int cfLp, double predictedCost, int leafProbes, int arbitraryConfLp) {
}
