package sswsi.util;

import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.objects.ObjectIterator;

import java.util.AbstractMap;
import java.util.Iterator;
import java.util.Map;

// Randomized dictionary backed by FastUtil's Int2ObjectOpenHashMap.
public final class DictHash<V> implements Dict<Integer, V> {

    private final Int2ObjectOpenHashMap<V> map;

    public DictHash() {
        this.map = new Int2ObjectOpenHashMap<>();
        this.map.defaultReturnValue(null);
    }

    @Override
    public V get(Integer key) {
        return map.get(key);
    }

    @Override
    public void put(Integer key, V value) {
        map.put(key, value);
    }

    @Override
    public boolean containsKey(Integer key) {
        return map.containsKey(key);
    }

    @Override
    public Iterable<Map.Entry<Integer, V>> entries() {
        return () -> new Iterator<Map.Entry<Integer, V>>() {
            final ObjectIterator<Int2ObjectMap.Entry<V>> it = map.int2ObjectEntrySet().fastIterator();

            @Override
            public boolean hasNext() {
                return it.hasNext();
            }

            @Override
            public Map.Entry<Integer, V> next() {
                Int2ObjectMap.Entry<V> e = it.next();
                return new AbstractMap.SimpleImmutableEntry<>(e.getIntKey(), e.getValue());
            }
        };
    }
}
