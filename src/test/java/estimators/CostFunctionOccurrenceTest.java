package estimators;

import PMIndex.HBI;
import PMIndex.HbiConfiguration;
import membership.MockMembership;
import org.junit.jupiter.api.Test;
import search.BlockSearch;
import search.MostFreqPruning;
import search.Pattern;
import search.VerifierLinearLeafProbe;
import utilities.MathUtils;
import utilities.PatternResult;

import java.util.Random;

import static org.junit.jupiter.api.Assertions.*;

// Checks CostFunctionOccurrence against the membership tests Block Search really performs from every
// start level, with exact membership so only the cost model can be wrong.
class CostFunctionOccurrenceTest {

    private static final int LOG2_W = 14;
    private static final int W = 1 << LOG2_W;

    // Random text over 26 letters with 6-grams: almost every n-gram occurs once, as on packet traces.
    // The product model predicts about one test per start node here and misses the occurrence path.
    @Test
    void sparseAlphabetPicksNearOptimalLevelAndPredictsTheOccurrencePath() {
        String text = randomText(W + 5, 26, 1L);
        Evaluation e = evaluate(text, 6, 24, 30, 2L, false);
        assertTrue(e.meanRegret() <= 1.15, "mean regret " + e.meanRegret());
        assertTrue(e.withinOne >= 0.9 * e.patterns, "within one level: " + e.withinOne + "/" + e.patterns);
        assertTrue(e.meanErrorAtPick() <= 0.35, "mean relative error at pick " + e.meanErrorAtPick());
    }

    // Dense iid text: nearly every symbol is in every large node, so the background term dominates
    // and the model behaves like the iid estimator of the thesis.
    @Test
    void denseAlphabetStillPicksNearOptimalLevel() {
        String text = randomText(W, 4, 3L);
        Evaluation e = evaluate(text, 1, 12, 30, 4L, false);
        assertTrue(e.meanRegret() <= 1.10, "mean regret " + e.meanRegret());
        assertTrue(e.withinOne >= 0.9 * e.patterns, "within one level: " + e.withinOne + "/" + e.patterns);
    }

    // Patterns that do not occur: a changed character turns the n-grams around it into absent keys,
    // so the search stops at the root and the model must not invent an occurrence path.
    @Test
    void missingPatternsStayCheap() {
        String text = randomText(W + 5, 26, 5L);
        Evaluation e = evaluate(text, 6, 24, 30, 6L, true);
        assertTrue(e.meanRegret() <= 1.15, "mean regret " + e.meanRegret());
    }

    private static Evaluation evaluate(String text, int n, int patternLength, int patterns, long seed,
                                       boolean mutate) {
        HBI hbi = newHbi(n);
        for (int i = n - 1; i < text.length(); i++) {
            hbi.insert(text.substring(i - n + 1, i + 1));
        }
        Random rnd = new Random(seed);
        int deepest = MathUtils.deepestVisitedLevel(W, patternLength);
        Evaluation e = new Evaluation();
        for (int k = 0; k < patterns; k++) {
            int start = rnd.nextInt(W - patternLength);
            char[] chars = text.substring(start, start + patternLength).toCharArray();
            if (mutate) {
                int at = 1 + rnd.nextInt(patternLength - 1);
                chars[at] = chars[at] == '#' ? '%' : '#';  // never in the text
            }
            String query = new String(chars);
            int[] actual = new int[deepest + 1];
            double[] predicted = new double[deepest + 1];
            int pick = -1;
            for (int lp = 0; lp <= deepest; lp++) {
                hbi.setLpOverride(lp);
                hbi.report(new Pattern(query, n));
                PatternResult pr = hbi.getLatestStats();
                actual[lp] = pr.probes();
                predicted[lp] = pr.predictedCost();
                pick = pr.cfLp();
            }
            int best = 0;
            for (int lp = 1; lp <= deepest; lp++) {
                if (actual[lp] < actual[best]) best = lp;
            }
            e.patterns++;
            e.regret += actual[pick] / (double) actual[best];
            e.errorAtPick += Math.abs(predicted[pick] - actual[pick]) / actual[pick];
            if (Math.abs(pick - best) <= 1) e.withinOne++;
        }
        return e;
    }

    private static HBI newHbi(int n) {
        HBI hbi = new HBI(HbiConfiguration.builder()
                .searchAlgorithm(new BlockSearch())
                .windowLength(W)
                .treeLength(W)
                .fpRate(0.0)
                .alphabetSize(1 << 20)
                .estimatorSupplier(() -> new HashMapEstimator(W))
                .membershipSupplier(MockMembership::new)
                .pruningPlanSupplier(() -> new MostFreqPruning(0.99, 0.0))
                .verifier(new VerifierLinearLeafProbe())
                .costFunction(new CostFunctionOccurrence())
                .confidence(0.99)
                .nGram(n)
                .experimentMode(true)
                .build());
        hbi.strides = false;
        hbi.stats().setCollecting(true);
        hbi.stats().setExperimentMode(true);
        return hbi;
    }

    private static String randomText(int length, int letters, long seed) {
        Random rnd = new Random(seed);
        StringBuilder sb = new StringBuilder(length);
        for (int i = 0; i < length; i++) {
            sb.append((char) ('a' + rnd.nextInt(letters)));
        }
        return sb.toString();
    }

    private static final class Evaluation {
        int patterns;
        int withinOne;
        double regret;
        double errorAtPick;

        double meanRegret() {
            return regret / patterns;
        }

        double meanErrorAtPick() {
            return errorAtPick / patterns;
        }
    }
}
