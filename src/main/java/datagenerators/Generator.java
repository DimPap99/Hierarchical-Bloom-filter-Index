package datagenerators;

import java.util.Random;
import org.apache.commons.math3.distribution.ZipfDistribution;

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

    public static String generateZipf(int length, int min_domain, int max_domain, double exponent) {
        ZipfDistribution dist = new ZipfDistribution(max_domain - min_domain, exponent);
        char[] chars = new char[length];
        for (int i = 0; i < length; i++) {
            chars[i] = (char) (dist.sample() - 1 + min_domain);
        }
        return new String(chars);
    }

}
