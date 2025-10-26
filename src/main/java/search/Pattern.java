package search;

import tree.ImplicitTree;

import java.util.ArrayList;
import java.util.Arrays;

public class Pattern{
    public char[] text;
    public String patternTxt;
    public int originalSz;
    public int[] pi;

    public String[] nGramArr;
    public long[] nGramToLong;      // filled by the index via StringKeyMapper
    public long[] effectiveNgramArr;
    public ArrayList<Integer> charStartLp;
    public int nGram = 1;
    public int size;
    public boolean isConverted;

    public Pattern(String s, int nGram) {
        this.text = s.toCharArray();
        size = text.length;
        this.originalSz = this.text.length;
        this.patternTxt = s;
        this.nGram = nGram;
        this.prefixFunction();
        int len = s.length() - nGram + 1;

        this.nGramArr   = new String[len];
        this.nGramToLong = new long[len];
        int addition =  originalSz%nGram == 0 ? 0 : 1;
        this.effectiveNgramArr = new long[originalSz/nGram + addition];
        for (int i = 0; i < len; i++) {
            nGramArr[i] = s.substring(i, i + nGram);   // sliding window
        }

    }

    public Pattern(String[] s, int nGram) {
        // this.text = s.toCharArray(); // not used in segments mode
        this.size = s.length;
        this.originalSz = s.length;
        // this.patternTxt = s;         // also unused/commented in your version
        this.nGram = nGram;

        int len = size - nGram + 1;

        // sliding-window n-grams
        this.nGramArr = new String[len];
        this.nGramToLong = new long[len]; // stays allocated; encoding filled elsewhere
        this.patternTxt = "";
        int addition = (originalSz % nGram == 0) ? 0 : 1;
        this.effectiveNgramArr = new long[(originalSz / nGram) + addition];

        for (int i = 0; i < len; i++) {
            StringBuilder sb = new StringBuilder();
            for (int j = 0; j < nGram; j++) {
                sb.append(s[i + j]);
            }
            nGramArr[i] = sb.toString();
            // nGramToLong[i] will get populated by the existing logic elsewhere
        }
        StringBuilder human = new StringBuilder();
        human.append(s[0]);
        for (int i = 1; i < s.length; i++) {
            human.append(' ');
            human.append(s[i]);
        }
        this.patternTxt = human.toString();
        // effectiveNgramArr is still allocated here and (per your note)
        // is populated later in the pipeline, so we leave it alone.
    }



//    public Pattern(String s, int nGram) {
//        this.text = s.toCharArray();
//        this.originalSz = this.text.length;
//        this.patternTxt = s;
//        this.nGram = nGram;
//        this.prefixFunction();
//        int len = s.length() - nGram + 1;
//
//        this.nGramArr   = new String[len];
//        this.nGramToInt = new int[len];        // or whatever you store there
//        int addition =  originalSz%nGram == 0 ? 0 : 1;
//        this.effectiveNgramArr = new int[originalSz/nGram + addition];
//        for (int i = 0; i < len; i++) {
//            nGramArr[i] = s.substring(i, i + nGram);   // sliding window
//            // nGramToInt[i] = … convert nGramArr[i] to int if needed
//        }
//
//    }




    /* ------- private copy of your helper so callers don't see it ------- */
    private int[] prefixFunction() {
        int[] pi = new int[this.text.length];
        int k = 0;
        for (int i = 1; i < this.text.length; ++i) {
            while (k > 0 && this.text[k] != this.text[i]) k = pi[k - 1];
            if (this.text[k] == this.text[i]) ++k;
            pi[i] = k;
        }
        return pi;
    }
}
