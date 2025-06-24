package membership;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Random;
import org.apache.commons.codec.digest.MurmurHash3;

public class BloomFilter implements Membership {
    BitSet filter;
    int[] seeds;
    private final byte[] leByteArr = new byte[8];
    private final int btLength = leByteArr.length;
    int n;
    double p;
    int m;
    int k;

    private static final long SALT1 = 0xC6BC279692B5CC83L;
    private static final long SALT2 = 0xD2B74407B1CE6E93L;

    @Override public String toString() {
        return String.format("n=%d, m=%d, k=%d, p(actual)=%.4g",
                n, m, k, p);
    }

    public BloomFilter(){};

    public void init(int n, double p){
        if (n <= 0)                  throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1)        throw new IllegalArgumentException("p must be in (0,1)");

        double ln2 = Math.log(2);
        //distinct elements
        this.n = n;
        //num of required bits
        this.m = (int) Math.ceil(-(n * Math.log(p)) / (ln2 * ln2));
        //num of haxsh funcs
        this.k = (int) Math.max(1, Math.round((m / (double) n) * ln2));
        //fp rate
        this.p    = Math.pow(1 - Math.exp(-k / (double) m * n), k);

        this.filter = new BitSet(m);
        initSeeds();
    };

    /**
     * https://en.wikipedia.org/wiki/Bloom_filter
     *
     * @param n
     * @param p
     */
    public BloomFilter(int n, double p) {
        if (n <= 0)                  throw new IllegalArgumentException("n must be > 0");
        if (p <= 0 || p >= 1)        throw new IllegalArgumentException("p must be in (0,1)");

        double ln2 = Math.log(2);
        this.n = n;
        this.m = (int) Math.ceil(-(n * Math.log(p)) / (ln2 * ln2));
        this.k = (int) Math.max(1, Math.round((m / (double) n) * ln2));
        this.p    = Math.pow(1 - Math.exp(-k / (double) m * n), k);

        this.filter = new BitSet(m);
        this.
        initSeeds();
    }

    private void initSeeds() {
        //Use only 2hf for Kirsch and Mitzenmacher Optimization
        this.seeds = new int[2];
        HashSet<Integer> seedHashset = new HashSet<Integer>();
        int seedsFound = 0;
        Random rand = new Random();
        while(seedsFound < seeds.length){
            int seed = rand.nextInt();
            if(!seedHashset.contains(seed)){
                seedHashset.add(seed);
                this.seeds[seedsFound] = seed;
                seedsFound++;
            }
        }
    }
    /* small helper: int → little-endian bytes */
    private void intToLE(int v, byte[] buf) {
        buf[0] = (byte)  v;
        buf[1] = (byte) (v >>> 8);
        buf[2] = (byte) (v >>> 16);
        buf[3] = (byte) (v >>> 24);
    }

    private void longToLE(long v, byte[] buf) {
        buf[0] = (byte)  v;
        buf[1] = (byte) (v >>> 8);
        buf[2] = (byte) (v >>> 16);
        buf[3] = (byte) (v >>> 24);
        buf[4] = (byte) (v >>> 32);
        buf[5] = (byte) (v >>> 40);
        buf[6] = (byte) (v >>> 48);
        buf[7] = (byte) (v >>> 56);
    }
    private static long mix64(long z) {
        z += 0x9E3779B97F4A7C15L;
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

//    public void insert(long key) {
//
//        longToLE(key, leByteArr);   // fill the 4-byte buffer
//
//        int h1 = Math.floorMod(
//                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[0]), m);
//        int h2 = Math.floorMod(
//                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[1]), m);
//
//        int idx = h1;
//        for (int i = 0; i < k; i++) {
//            filter.set(idx);
//            idx += h2;
//            if (idx >= m) idx -= m;
//        }
//    }
//
//    /* ================================================================== */
//    /*  CONTAINS (int key)                                                */
//    /* ================================================================== */
//    public boolean contains(long key) {
//
//        longToLE(key, leByteArr);   // reuse same buffer
//
//        int h1 = Math.floorMod(
//                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[0]), m);
//        int h2 = Math.floorMod(
//                MurmurHash3.hash32x86(leByteArr, 0, btLength, seeds[1]), m);
//
//        int idx = h1;
//        for (int i = 0; i < k; i++) {
//            if (!filter.get(idx))
//                return false;
//            idx += h2;
//            if (idx >= m) idx -= m;
//        }
//        return true;                // could still be FP
//    }
    public void insert(long key) {
        long h1 = mix64(key ^ SALT1);           // first hash
        long h2 = mix64(key ^ SALT2);           // second, independent hash

        int idx = (int) Math.floorMod(h1, m);   // h1 mod m  (m is int)
        int inc = (int) Math.floorMod(h2, m);   // h2 mod m  (step size)

        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += inc;
            if (idx >= m) idx -= m;             // 1× wrap
        }
    }

        public boolean contains(long key) {
            long h1 = mix64(key ^ SALT1);
            long h2 = mix64(key ^ SALT2);

            int idx = (int) Math.floorMod(h1, m);
            int inc = (int) Math.floorMod(h2, m);

            for (int i = 0; i < k; i++) {
                if (!filter.get(idx)) return false;
                idx += inc;
                if (idx >= m) idx -= m;
            }
            return true;                            // maybe FP
        }


}