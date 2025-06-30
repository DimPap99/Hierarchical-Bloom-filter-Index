package search;

import java.util.ArrayList;

/** Immutable wrapper around the query plus its KMP π-table. */
public class Pattern{
    public char[] text;
    public String patternTxt;
    public int[] pi;

    public String[] nGramArr;
    public int[] nGramToInt;
    int nGram = 1;
    /** Factory that pre-computes π[] once. */
    public Pattern(String s, int nGram) {
        this.text = s.toCharArray();
        this.patternTxt = s;
        this.prefixFunction();
        int sz = (int)Math.ceil(s.length()/nGram);
        this.nGramArr = new String[sz];
        this.nGramToInt = new int[sz];
        for(int i = nGram; i <= s.length(); i++){
            this.nGramArr[i - 1] = s.substring(i-nGram, i);
        }
    }

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

