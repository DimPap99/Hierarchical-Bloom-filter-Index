package utilities;

import it.unimi.dsi.fastutil.longs.Long2IntOpenHashMap;
import org.apache.commons.codec.digest.MurmurHash3;

import java.nio.charset.StandardCharsets;

/**
 * Maps string tokens to dense integer ids via 64-bit Murmur3 hashing. We keep only
 * the hash-to-id table (no original strings), so collisions are vanishingly unlikely
 * while the downstream code can continue to operate on compact ints.
 */
public final class HashedStringMapper {

    private final Long2IntOpenHashMap hashToId;
    private int nextId = 1; // 0 reserved for sentinel

    public HashedStringMapper(int expectedDistinct) {
        int capacity = Math.max(16, expectedDistinct);
        this.hashToId = new Long2IntOpenHashMap(capacity);
        this.hashToId.defaultReturnValue(-1);
    }

    public int getId(String token) {
        if (token == null) {
            return 0;
        }
        byte[] bytes = token.getBytes(StandardCharsets.UTF_8);
        long[] h = MurmurHash3.hash128x64(bytes);
        long hash = (h[0] ^ h[1]) & 0x7fffffffffffffffL;
        int id = hashToId.get(hash);
        if (id == -1) {
            id = nextId++;
            hashToId.put(hash, id);
        }
        return id;
    }

    public int size() {
        return hashToId.size();
    }

    public void clear() {
        hashToId.clear();
        nextId = 1;
    }
}

