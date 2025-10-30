package utilities;

import it.unimi.dsi.fastutil.objects.Object2IntOpenHashMap;

public class AlphabetMapper<T> {
    // 0 is reserved as a global sentinel for suffix trees
    int nextId = 1;
    float loadFactor = 0.75f;

    // Primitive map to avoid boxing and reduce overhead; public for incidental tooling
    public Object2IntOpenHashMap<T> wordToId;

    int capacity;

    public AlphabetMapper(int capacity) {
        this.nextId = 1; // id 0 is reserved as a sentinel for the suffix tree DC3
        this.capacity = Math.max(1, capacity);

        // Pre-size to the expected alphabet size to avoid rehashing.
        this.wordToId = new Object2IntOpenHashMap<>(this.capacity, loadFactor);
        // Use -1 as the default return to distinguish from valid ids (>=1) and sentinel 0.
        this.wordToId.defaultReturnValue(-1);
    }

    public int getSize() {
        return wordToId.size();
    }

    public int getCapacity() {
        return capacity;
    }

    public int insert(T item) {
        int id = wordToId.getInt(item);
        if (id == -1) {
            id = nextId;
            nextId++;
            wordToId.put(item, id);
        }
        return id;
    }

    // Insert-on-miss mapping
    public int getId(T item) {
        int id = wordToId.getInt(item);
        if (id == -1) {
            id = insert(item);
        }
        return id;
    }

    public void clear() {
        wordToId.clear();
        nextId = 1; // keep 0 reserved for sentinel
    }
}
