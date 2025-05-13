package membership;

import java.util.BitSet;
import java.util.Random;
import org.apache.commons.codec.digest.MurmurHash3;

public class BloomFilter implements Membership {
    BitSet filter;
    int[] seeds;

    int n;
    double p;
    int m;
    int k;

    @Override
    public String toString() {
        return "n: " + n + ", p: " + p + ", m: " + m + ", k: " + k;
    }

    /**
     * https://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=8586915
     *
     * @param n Expected inserted elements
     * @param m Length of bit array
     */
    public BloomFilter(int n, int m) {
        this.n = n;
        this.m = m;
        this.k = (int) Math.max(1, Math.round(((double) m / n) * Math.log(2)));
        this.p = Math.pow(0.5, this.k);
        initBloomFilter();
    }

    /**
     * https://hur.st/bloomfilter
     * NOT IMPLEMENTED
     *
     * @param n Expected inserted elements
     * @param p False positive probability
     * @param m Length of bit array
     * @param k Number of hash functions
     */
    public BloomFilter(int n, double p, int m, int k) {

    }

    /**
     * https://en.wikipedia.org/wiki/Bloom_filter
     * BROKEN FORMULAS
     *
     * @param n
     * @param p
     */
    public BloomFilter(int n, double p) {
        int m = (int) Math.round(-((n * Math.log(Math.exp(1)) / Math.pow(Math.log(2), 2))));
        int k = (int) Math.round(-(Math.log(p) / Math.log(2)));
        this.n = n;
        this.p = p;
        this.m = m;
        this.k = k;
        initBloomFilter();
    }

    private void initBloomFilter() {
        this.filter = new BitSet(m);
        this.seeds = new int[k];
        Random rand = new Random();
        for (int i = 0; i < seeds.length; i++) {
            this.seeds[i] = rand.nextInt();
        }
    }

    public void insert(String key) {
        for (int i = 0; i < seeds.length; i++) {
            this.filter.set(Math.floorMod(hash(key, this.seeds[i]), this.m));
        }
    }

    public boolean contains(String key) {
        for (int i = 0; i < seeds.length; i++) {
            if (this.filter.get(Math.floorMod(hash(key, this.seeds[i]), this.m)) == false) return false;
        }
        return true;
    }

    private int hash(String word, int seed) {
        byte[] bytes = word.getBytes();
        return MurmurHash3.hash32x86(bytes, 0, bytes.length, seed);
    }
}