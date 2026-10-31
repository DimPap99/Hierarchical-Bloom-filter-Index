package utilities;

import tree.ssws.ChildDictionary;
import tree.ssws.SuffixTreeBuilder;

import java.util.Locale;

/**
 * Which SSWSI the benchmarks build, from {@code --sswsi-mode}, {@code --sswsi-builder} and {@code --sswsi-dict}.
 *
 * <p>{@code amortized} is the paper's amortized structure ({@link PMIndex.SlidingSuffixIndex}), with a suffix tree
 * builder and a child dictionary; {@code legacy} is the thesis-era structure with Phase 0's fixes
 * ({@link PMIndex.StreamingSlidingWindowIndex}), which has neither choice. The default is amortized with DC3 and
 * FKS: FKS is what the paper uses, and Farach-Colton, the paper's construction, builds the same trees 1.1-1.4 times
 * slower (docs/sswsi-results.md, Phase 2).
 */
public record SswsiSettings(Mode mode, SuffixTreeBuilder.Kind builder, ChildDictionary.Kind dictionary) {

    public enum Mode {
        LEGACY, AMORTIZED;

        public static Mode parse(String name) {
            return valueOf(name.trim().toUpperCase(Locale.ROOT));
        }
    }

    public static final SswsiSettings DEFAULT =
            new SswsiSettings(Mode.AMORTIZED, SuffixTreeBuilder.Kind.DC3, ChildDictionary.Kind.FKS);

    public static final SswsiSettings LEGACY = new SswsiSettings(Mode.LEGACY, null, null);

    public SswsiSettings {
        if (mode == null) {
            throw new IllegalArgumentException("SSWSI mode missing");
        }
        if ((mode == Mode.LEGACY) != (builder == null) || (mode == Mode.LEGACY) != (dictionary == null)) {
            throw new IllegalArgumentException("A builder and a dictionary are set exactly for amortized SSWSI");
        }
    }

    /**
     * Settings from the three flags; a null argument is a flag that was not given. A builder or dictionary with
     * {@code legacy} is an error, since the legacy index has neither choice.
     */
    public static SswsiSettings parse(String mode, String builder, String dictionary) {
        Mode m;
        try {
            m = (mode == null) ? DEFAULT.mode() : Mode.parse(mode);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("Unknown --sswsi-mode '" + mode + "' (expected legacy or amortized)");
        }
        if (m == Mode.LEGACY) {
            if (builder != null || dictionary != null) {
                throw new IllegalArgumentException(
                        "--sswsi-builder and --sswsi-dict apply only to --sswsi-mode amortized");
            }
            return LEGACY;
        }
        SuffixTreeBuilder.Kind b;
        try {
            b = (builder == null) ? DEFAULT.builder() : SuffixTreeBuilder.Kind.parse(builder);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("Unknown --sswsi-builder '" + builder + "' (expected fc or dc3)");
        }
        ChildDictionary.Kind d;
        try {
            d = (dictionary == null) ? DEFAULT.dictionary() : ChildDictionary.Kind.parse(dictionary);
        } catch (IllegalArgumentException e) {
            throw new IllegalArgumentException("Unknown --sswsi-dict '" + dictionary + "' (expected fastutil or fks)");
        }
        return new SswsiSettings(m, b, d);
    }

    /** {@code legacy}, or {@code amortized-<builder>-<dictionary>}, for example {@code amortized-dc3-fks}. */
    public String label() {
        if (mode == Mode.LEGACY) {
            return "legacy";
        }
        return "amortized-" + builder.name().toLowerCase(Locale.ROOT) + "-" + dictionary.name().toLowerCase(Locale.ROOT);
    }

    /** The amortized index needs a power-of-two window; fail before any data is read rather than on the first build. */
    public void requireValidWindow(int windowLength) {
        if (mode == Mode.AMORTIZED && (windowLength <= 0 || Integer.bitCount(windowLength) != 1)) {
            throw new IllegalArgumentException("--sswsi-mode amortized needs a power-of-two window length, got "
                    + windowLength + " (use --window-length, or --sswsi-mode legacy)");
        }
    }

    @Override
    public String toString() {
        return label();
    }
}
