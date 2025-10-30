package PMIndex;

/**
 * Marker interface for streaming sliding window string indices that follow the
 * SSWSI model from "Faster Sliding Window String Indexing in Streams" (CPM 2024).
 *
 * All concrete implementations also expose the operations mandated by
 * {@link IPMIndexing} so that they can participate in the existing benchmarking
 * harness that expects {@code IPMIndexing} instances.
 */
public interface SSWSI extends IPMIndexing {
}
