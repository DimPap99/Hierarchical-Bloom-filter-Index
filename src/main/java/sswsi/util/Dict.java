package sswsi.util;

import java.util.Map;

// Abstraction for dynamic dictionaries D(v) mapping first symbol to child edge.
public interface Dict<K, V> {

    V get(K key);

    void put(K key, V value);

    boolean containsKey(K key);

    Iterable<Map.Entry<K, V>> entries();
}
