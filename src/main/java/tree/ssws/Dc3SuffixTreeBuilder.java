package tree.ssws;

/**
 * The pipeline of the legacy {@link SuffixTree}, building a {@link StaticSuffixTree} instead: rank reduction
 * ({@link RankReducer}), the suffix array by DC3 ({@link SuffixArrayBuilder#buildSuffixArray}), the LCP array by
 * Kasai et al. ({@link SuffixArrayBuilder#buildLcpArray}), then the stack build from both arrays.
 */
public final class Dc3SuffixTreeBuilder implements SuffixTreeBuilder {

    @Override
    public StaticSuffixTree build(long[] text, int offset, int length, ChildDictionary.Kind dictionary) {
        int[] ranks = RankReducer.rankWithTerminator(text, offset, length);
        int[] sa = SuffixArrayBuilder.buildSuffixArray(ranks);
        int[] lcp = SuffixArrayBuilder.buildLcpArray(ranks, sa);
        return StaticSuffixTree.fromSuffixArray(text, offset, length, sa, lcp, dictionary);
    }
}
