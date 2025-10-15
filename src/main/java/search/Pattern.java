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
//        this.text = s.toCharArray();
        size = s.length;
        this.originalSz = s.length;
//        this.patternTxt = s;
        this.nGram = nGram;

        int len = size - nGram + 1;

        this.nGramArr   = new String[len];
        Arrays.fill(this.nGramArr, "");
        this.nGramToLong = new long[len];
        int addition =  originalSz%nGram == 0 ? 0 : 1;
        this.effectiveNgramArr = new long[originalSz/nGram + addition];
        for (int i = 0; i < len; i++) {
            for (int j = i; j < i + nGram; j++ ) {
                nGramArr[i] += s[j];   // sliding window
            }
        }

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
