package tree.ssws;

/**
 * Static map from (internal node, first symbol of an outgoing edge) to the child, built once per tree.
 *
 * <p>Edges whose first symbol is the terminator are left out: no pattern symbol can match the terminator.
 * The trees never change after construction, so a static dictionary is enough (Bille et al., Sec. 2).
 */
public interface ChildDictionary {

    /** The child of {@code node} whose edge starts with {@code symbol}, or {@link StaticSuffixTree#NONE}. */
    int child(int node, long symbol);

    enum Kind {
        /**
         * Up to two children per node compared inline, more in a per-node fastutil {@code Long2IntOpenHashMap}.
         * The thesis's documented choice and the default.
         */
        FASTUTIL,
        /** One two-level perfect hash (Fredman, Komlós, Szemerédi) per tree over all (node, first symbol) pairs. */
        FKS;

        ChildDictionary build(StaticSuffixTree tree) {
            return switch (this) {
                case FASTUTIL -> new FastutilChildDictionary(tree);
                case FKS -> new FksChildDictionary(tree);
            };
        }

        /** Parses {@code fastutil} or {@code fks}, ignoring case. */
        public static Kind parse(String name) {
            return valueOf(name.trim().toUpperCase(java.util.Locale.ROOT));
        }
    }
}
