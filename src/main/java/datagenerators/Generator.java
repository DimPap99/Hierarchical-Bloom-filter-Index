package datagenerators;

import java.util.Random;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.distribution.ZipfDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.Well19937c;

public class Generator {
    public static final char[] EXTRA_CHARS = (

            // Cyrillic lower-case (replace or extend as you wish)
                    // Greek lower-case
                    "αβγδεζηθικλμνξοπρστυφχψω" +
                    // Greek lower-case with tonos
                    "άέήίόύώ" +
                    // Greek upper-case
                    "ΑΒΓΔΕΖΗΘΙΚΛΜΝΞΟΠΡΣΤΥΦΧΨΩ" +
                    // Greek upper-case with tonos
                    "ΆΈΉΊΌΎΏ" +
                    // English upper-case
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZ" +
                    // English lower-case
                    "abcdefghijklmnopqrstuvwxyz"
    ).toCharArray();

    public static String generateUniform(int length, int min_domain, int max_domain) {
        if (min_domain < Character.MIN_VALUE) throw new IllegalArgumentException("min_domain < 0");
        if (max_domain > Character.MAX_VALUE) throw new IllegalArgumentException("max_domain > 0");

        Random rand = new Random();
        char[] chars = new char[length];
        for (int i = 0; i < length; i++) {
            // Character.MAX_VALUE is in the BMP (uses single char)
            chars[i] = (char) rand.nextInt(min_domain, max_domain);
        }
        return new String(chars);
    }

    public static String generateZipf(int length,
                                      int min_domain,
                                      int max_domain,
                                      double exponent,
                                      long seed) {
        if (min_domain < Character.MIN_VALUE) {
            throw new IllegalArgumentException("min_domain < Character.MIN_VALUE");
        }
        // We want codes in the half open interval [min_domain, max_domain)
        // so the number of distinct symbols is:
        int alphabetSize = max_domain - min_domain;
        if (alphabetSize <= 0) {
            throw new IllegalArgumentException("max_domain must be greater than min_domain");
        }
        if (max_domain - 1 > Character.MAX_VALUE) {
            throw new IllegalArgumentException("max_domain - 1 exceeds Character.MAX_VALUE");
        }

        // Seeded random number generator for reproducibility
        RandomGenerator rng = new Well19937c(seed);

        // ZipfDistribution samples integers in the closed interval [1, alphabetSize]
        ZipfDistribution dist = new ZipfDistribution(rng, alphabetSize, exponent);

        char[] chars = new char[length];
        for (int i = 0; i < length; i++) {
            int rank = dist.sample();                 // 1 .. alphabetSize
            int code = min_domain + (rank - 1);       // map to [min_domain, max_domain)
            chars[i] = (char) code;                   // Basic Multilingual Plane only
        }
        return new String(chars);
    }

}
