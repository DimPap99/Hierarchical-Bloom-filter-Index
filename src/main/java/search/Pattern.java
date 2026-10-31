package search;

import java.util.ArrayList;

public class Pattern{
    // Joins the tokens of one n-gram in segments mode. Segments are lines and query tokens are split on
    // spaces, so a line break can never occur inside a token and distinct token sequences stay distinct.
    public static final String SEGMENT_JOINER = "\n";

    public char[] text;
    public String patternTxt;
    public int originalSz;

    public String[] nGramArr;
    public long[] nGramToLong;      // filled by the index via StringKeyMapper, or by setTokens
    public long[] effectiveNgramArr; // every n-th token for stride mode, see buildStrideTokens
    public ArrayList<Integer> charStartLp;
    public int nGram = 1;
    public int size;
    // True once nGramToLong holds the final numeric tokens, so the index must not remap nGramArr.
    public boolean isConverted;

    public Pattern(String s, int nGram) {
        this.text = s.toCharArray();
        size = text.length;
        this.originalSz = this.text.length;
        this.patternTxt = s;
        this.nGram = nGram;
        int len = s.length() - nGram + 1;
        boolean shortPattern = len <= 0;
        if (shortPattern) {
            len = 1;
        }

        this.nGramArr   = new String[len];
        this.nGramToLong = new long[len];
        int addition =  originalSz%nGram == 0 ? 0 : 1;
        this.effectiveNgramArr = new long[originalSz/nGram + addition];
        if (shortPattern) {
            nGramArr[0] = s;
        } else {
            for (int i = 0; i < len; i++) {
                nGramArr[i] = s.substring(i, i + nGram);   // sliding window
            }
        }

    }

    // Segments mode: each element of s is one stream token.
    public Pattern(String[] s, int nGram) {
        this.size = s.length;
        this.originalSz = s.length;
        this.nGram = nGram;

        int len = size - nGram + 1;
        boolean shortPattern = len <= 0;
        if (shortPattern) {
            len = 1;
        }

        // sliding-window n-grams
        this.nGramArr = new String[len];
        this.nGramToLong = new long[len]; // stays allocated; encoding filled elsewhere
        int addition = (originalSz % nGram == 0) ? 0 : 1;
        this.effectiveNgramArr = new long[(originalSz / nGram) + addition];

        if (shortPattern) {
            nGramArr[0] = String.join(SEGMENT_JOINER, s);
        } else {
            for (int i = 0; i < len; i++) {
                nGramArr[i] = joinTokens(s, i, nGram);
            }
        }
        this.patternTxt = String.join(" ", s);
        this.text = this.patternTxt.toCharArray();
    }

    // Builds the n-gram key used for count consecutive tokens starting at from. Must match how the
    // stream side builds its n-grams in segments mode.
    public static String joinTokens(String[] tokens, int from, int count) {
        if (count == 1) {
            return tokens[from];
        }
        StringBuilder sb = new StringBuilder();
        for (int j = 0; j < count; j++) {
            if (j > 0) {
                sb.append(SEGMENT_JOINER);
            }
            sb.append(tokens[from + j]);
        }
        return sb.toString();
    }

    // Use already-mapped numeric tokens (one per overlapping n-gram) instead of mapping nGramArr.
    public void setTokens(long[] tokens) {
        if (tokens.length != nGramToLong.length) {
            throw new IllegalArgumentException(
                    "expected " + nGramToLong.length + " tokens but got " + tokens.length);
        }
        System.arraycopy(tokens, 0, nGramToLong, 0, tokens.length);
        this.isConverted = true;
    }

    // Fill effectiveNgramArr with the n-grams at offsets 0, n, 2n, ... (Sec. 4.3.1). When the
    // length is not a multiple of n, the last entry is the final overlapping n-gram so the tail
    // of the pattern is still covered.
    public void buildStrideTokens() {
        for (int nIdx = 0; nIdx < nGramToLong.length; nIdx += nGram) {
            effectiveNgramArr[nIdx / nGram] = nGramToLong[nIdx];
        }
        if (effectiveNgramArr.length > 0 && originalSz % nGram != 0) {
            effectiveNgramArr[effectiveNgramArr.length - 1] = nGramToLong[nGramToLong.length - 1];
        }
    }
}
