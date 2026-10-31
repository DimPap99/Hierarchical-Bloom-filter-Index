package tree.ssws;

/**
 * The construction the paper uses: rank reduction ({@link RankReducer}), then the suffix tree of the ranked text
 * by {@link FarachColton}, written into a {@link StaticSuffixTree}. Gives the same trees as
 * {@link Dc3SuffixTreeBuilder}, but the suffix and LCP arrays take 2-3 times as long from 2^16 tokens up: every
 * level builds two tries, over-merges them and corrects the merge, where DC3 only radix-sorts.
 */
public final class FarachColtonSuffixTreeBuilder implements SuffixTreeBuilder {

    @Override
    public StaticSuffixTree build(long[] text, int offset, int length, ChildDictionary.Kind dictionary) {
        int[] s = RankReducer.rankWithTerminator(text, offset, length, 1);    // s[length] = s[length + 1] = 0
        int sigma = 0;
        for (int i = 0; i < length; i++) {
            sigma = Math.max(sigma, s[i]);
        }
        int[] sa = new int[length + 1];
        int[] lcp = new int[length];
        FarachColton.build(s, length, sigma, sa, lcp);
        return StaticSuffixTree.fromSuffixArray(text, offset, length, sa, lcp, dictionary);
    }
}
