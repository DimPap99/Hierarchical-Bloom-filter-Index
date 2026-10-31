package tree.ssws;

import java.util.Locale;

/** Builds a {@link StaticSuffixTree}. Implementations must give identical trees for the same text. */
public interface SuffixTreeBuilder {

    /**
     * Builds the suffix tree of {@code text[offset, offset + length)} followed by a unique terminator. The tree keeps a
     * reference to {@code text} and reads those entries whenever it is searched, so they must not change while the
     * tree is in use; the rest of {@code text} may.
     */
    StaticSuffixTree build(long[] text, int offset, int length, ChildDictionary.Kind dictionary);

    enum Kind {
        /** Farach-Colton ({@link FarachColtonSuffixTreeBuilder}), the construction the paper uses. */
        FC,
        /** DC3 and Kasai ({@link Dc3SuffixTreeBuilder}), the pipeline of the legacy index; the default. */
        DC3;

        public SuffixTreeBuilder create() {
            return (this == FC) ? new FarachColtonSuffixTreeBuilder() : new Dc3SuffixTreeBuilder();
        }

        /** Parses {@code fc} or {@code dc3}, ignoring case. */
        public static Kind parse(String name) {
            return valueOf(name.trim().toUpperCase(Locale.ROOT));
        }
    }
}
