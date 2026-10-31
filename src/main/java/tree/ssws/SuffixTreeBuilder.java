package tree.ssws;

/** Builds a {@link StaticSuffixTree}. Implementations must give identical trees for the same text. */
public interface SuffixTreeBuilder {

    /**
     * Builds the suffix tree of {@code text[offset, offset + length)} followed by a unique terminator. The tree keeps a
     * reference to {@code text}, which must not change afterwards.
     */
    StaticSuffixTree build(long[] text, int offset, int length, ChildDictionary.Kind dictionary);
}
