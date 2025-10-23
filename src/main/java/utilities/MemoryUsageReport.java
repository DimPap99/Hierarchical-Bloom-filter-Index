package utilities;

/** Simple container for a human-readable JOL report string and the total MiB value. */
public record MemoryUsageReport(String report, double totalMiB) {}

