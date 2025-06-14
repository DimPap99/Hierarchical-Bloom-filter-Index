package membership;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Random;
import org.apache.commons.codec.digest.MurmurHash3;

public class BloomFilter implements Membership {
    BitSet filter;
    int[] seeds;

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
        this.seeds = new int[k];
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

    public void insert(String key) {
        byte[] bytes = key.getBytes();
        for (int seed : seeds)
            filter.set(Math.floorMod(MurmurHash3.hash32x86(bytes, 0, bytes.length, seed), this.m));
    }

    public boolean contains(String key) {
        byte[] bytes = key.getBytes();
        for (int seed : seeds)
            if (!filter.get(Math.floorMod(MurmurHash3.hash32x86(bytes, 0, bytes.length, seed), this.m)))
                return false;
        return true;
    }
}