package sswsi.util;

import java.util.Map;

/**
 * Abstraction for the dynamic dictionaries D(v) that map the first symbol on an
 * outgoing edge to the child edge within the online suffix tree. The paper's
 * deterministic and randomized implementations differ only by the underlying
 * dictionary, so we keep the interface small.
 */
public interface Dict<K, V> {

    V get(K key);

    void put(K key, V value);

    boolean containsKey(K key);

    Iterable<Map.Entry<K, V>> entries();
}
