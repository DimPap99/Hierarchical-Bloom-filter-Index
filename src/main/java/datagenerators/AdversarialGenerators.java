package datagenerators;

import java.util.Arrays;
import java.util.Objects;

/**
 * Helpers that synthesise adversarial workloads meant to stress the HBI pruning logic.
 * Keeping the generators in Java allows us to reproduce datasets without depending on
 * external scripts and keeps the CLI integration simple.
 */
public final class AdversarialGenerators {

    private AdversarialGenerators() {
    }

    /**
     * Builds a stream made of alternating mono-character blocks. Block boundaries are aligned
     * so that each Bloom interval of size {@code blockLength} or smaller sees a single token.
     *
     * @param totalLength total number of characters to emit
     * @param blockLength length of each mono-character block (use powers of two to align with tree intervals)
     * @param alphabet    sequence of characters to cycle through for each block
     * @return dataset contents as a {@link String}
     */
    public static String generateAlternatingBlocks(int totalLength,
                                                   int blockLength,
                                                   char[] alphabet) {
        if (totalLength <= 0) {
            throw new IllegalArgumentException("totalLength must be positive");
        }
        if (blockLength <= 0) {
            throw new IllegalArgumentException("blockLength must be positive");
        }
        Objects.requireNonNull(alphabet, "alphabet");
        if (alphabet.length == 0) {
            throw new IllegalArgumentException("alphabet must contain at least one character");
        }

        StringBuilder sb = new StringBuilder(totalLength);
        int produced = 0;
        int symbolIndex = 0;
        while (produced < totalLength) {
            char symbol = alphabet[symbolIndex % alphabet.length];
            int runLength = Math.min(blockLength, totalLength - produced);
            sb.append(String.valueOf(symbol).repeat(runLength));
            produced += runLength;
            symbolIndex++;
        }
        return sb.toString();
    }

    /**
     * Emits a de Bruijn sequence for the specified alphabet and order. The full cycle is repeated
     * until {@code totalLength} characters are produced.
     *
     * @param alphabetSize number of symbols in the alphabet
     * @param order        de Bruijn order (equals the maximum n-gram in our workloads)
     * @param totalLength  desired output length
     * @param alphabetBase first code point in the alphabet; the alphabet is the contiguous range
     *                     {@code [alphabetBase, alphabetBase + alphabetSize)} (0-255 works just fine)
     * @return dataset contents as a {@link String}
     */
    public static String generateDeBruijnSequence(int alphabetSize,
                                                  int order,
                                                  int totalLength,
                                                  int alphabetBase) {
        if (alphabetBase < Character.MIN_VALUE) {
            throw new IllegalArgumentException("alphabetBase must be >= " + (int) Character.MIN_VALUE);
        }
        if (alphabetBase + alphabetSize > Character.MAX_VALUE + 1) {
            throw new IllegalArgumentException("alphabet exceeds character range");
        }
        char[] alphabet = new char[alphabetSize];
        for (int i = 0; i < alphabetSize; i++) {
            alphabet[i] = (char) (alphabetBase + i);
        }
        return generateDeBruijnSequence(alphabet, order, totalLength);
    }

    /**
     * Emits a de Bruijn sequence using an explicit alphabet.
     *
     * @param alphabet    array of symbols to use (copied directly, so any subset of 0-65535 is valid)
     * @param order       de Bruijn order / maximum n-gram length
     * @param totalLength desired output length
     * @return dataset contents as a {@link String}
     */
    public static String generateDeBruijnSequence(char[] alphabet,
                                                  int order,
                                                  int totalLength) {
        Objects.requireNonNull(alphabet, "alphabet");
        if (alphabet.length < 2) {
            throw new IllegalArgumentException("alphabet must contain at least two symbols");
        }
        if (order < 1) {
            throw new IllegalArgumentException("order must be positive");
        }
        if (totalLength <= 0) {
            throw new IllegalArgumentException("totalLength must be positive");
        }

        StringBuilder sb = new StringBuilder(totalLength);
        int k = alphabet.length;
        long workSize = (long) k * order;
        if (workSize > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("alphabetSize * order exceeds supported range");
        }
        int[] work = new int[(int) Math.max(1L, workSize)];

        while (sb.length() < totalLength) {
            Arrays.fill(work, 0);
            boolean completed = emitDeBruijnCycle(1, 1, k, order, work, alphabet, sb, totalLength);
            if (!completed) {
                break;
            }
        }
        return sb.toString();
    }

    // Emits a single de Bruijn cycle using Lyndon words; stops early when we reach the target length.
    private static boolean emitDeBruijnCycle(int t,
                                             int p,
                                             int k,
                                             int n,
                                             int[] a,
                                             char[] alphabet,
                                             StringBuilder out,
                                             int limit) {
        if (t > n) {
            if (n % p == 0) {
                for (int i = 1; i <= p; i++) {
                    out.append(alphabet[a[i]]);
                    if (out.length() >= limit) {
                        return false;
                    }
                }
            }
        } else {
            a[t] = a[t - p];
            if (!emitDeBruijnCycle(t + 1, p, k, n, a, alphabet, out, limit)) {
                return false;
            }
            for (int j = a[t - p] + 1; j < k; j++) {
                a[t] = j;
                if (!emitDeBruijnCycle(t + 1, t, k, n, a, alphabet, out, limit)) {
                    return false;
                }
            }
        }
        return true;
    }
}
