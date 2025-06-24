package membership;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Random;
import org.apache.commons.codec.digest.MurmurHash3;

public class BloomFilter implements Membership {
    BitSet filter;
    int[] seeds;
    private final byte[] le4 = new byte[4];

    int n;
    double p;
    int m;
    int k;


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

    public void insert(int key) {

        intToLE(key, le4);   // fill the 4-byte buffer

        int h1 = Math.floorMod(
                MurmurHash3.hash32x86(le4, 0, 4, seeds[0]), m);
        int h2 = Math.floorMod(
                MurmurHash3.hash32x86(le4, 0, 4, seeds[1]), m);

        int idx = h1;
        for (int i = 0; i < k; i++) {
            filter.set(idx);
            idx += h2;
            if (idx >= m) idx -= m;
        }
    }

    /* ================================================================== */
    /*  CONTAINS (int key)                                                */
    /* ================================================================== */
    public boolean contains(int key) {

        intToLE(key, le4);   // reuse same buffer

        int h1 = Math.floorMod(
                MurmurHash3.hash32x86(le4, 0, 4, seeds[0]), m);
        int h2 = Math.floorMod(
                MurmurHash3.hash32x86(le4, 0, 4, seeds[1]), m);

        int idx = h1;
        for (int i = 0; i < k; i++) {
            if (!filter.get(idx))
                return false;
            idx += h2;
            if (idx >= m) idx -= m;
        }
        return true;                // could still be FP
    }


}