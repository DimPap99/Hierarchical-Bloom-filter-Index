package search;

import tree.ImplicitTree;

import java.util.ArrayList;

/** Immutable wrapper around the query plus its KMP π-table. */
public class Pattern{
    public char[] text;
    public String patternTxt;
    public int[] pi;

    public String[] nGramArr;
    public int[] nGramToInt;

    public ArrayList<Integer> charStartLp;
    int nGram = 1;

    public Pattern(String s, int nGram) {
        this.text = s.toCharArray();
        this.patternTxt = s;
        this.prefixFunction();
        int len = s.length() - nGram + 1;

        this.nGramArr   = new String[len];
        this.nGramToInt = new int[len];        // or whatever you store there
        for (int i = 0; i < len; i++) {
            nGramArr[i] = s.substring(i, i + nGram);   // sliding window
            // nGramToInt[i] = … convert nGramArr[i] to int if needed
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

