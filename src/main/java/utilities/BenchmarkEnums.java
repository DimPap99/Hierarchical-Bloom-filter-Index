package utilities;

import java.util.EnumSet;

public final class BenchmarkEnums {
    private BenchmarkEnums() {}

    public enum IndexType {
        HBI("HBI", "hbi"),
        SUFFIX("SuffixTree", "suffix"),
        REGEX("Regex", "regex");

        private final String displayName;
        private final String csvLabel;

        IndexType(String displayName, String csvLabel) {
            this.displayName = displayName;
            this.csvLabel = csvLabel;
        }

        public String displayName() { return displayName; }
        public String csvLabel() { return csvLabel; }
    }

    public enum QueryType {
        MISSING("missing"),
        RARE("rare"),
        UNIFORM("uniform");

        private final String fileToken;
        QueryType(String token) { this.fileToken = token; }
        public String fileToken() { return fileToken; }

        public static QueryType fromString(String value) {
            return EnumSet.allOf(QueryType.class).stream()
                    .filter(type -> type.fileToken.equalsIgnoreCase(value))
                    .findFirst()
                    .orElseThrow(() -> new IllegalArgumentException("Unknown query type: " + value));
        }
    }
}

