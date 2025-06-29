package search;

import java.util.ArrayList;

/** Immutable wrapper around the query plus its KMP π-table. */
public class Pattern{
    public char[] text;
    public int[] pi;

    boolean isMultiLevel = false;

    ArrayList<Integer> pruneLevels = new ArrayList<>();

    /** Factory that pre-computes π[] once. */
    public Pattern(String s, boolean isMultiLevel) {
        this.text = s.toCharArray();
        this.prefixFunction();
        this.isMultiLevel = isMultiLevel;

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

