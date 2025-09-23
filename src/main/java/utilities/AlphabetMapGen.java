package utilities;

import java.util.HashMap;
import java.util.List;

public class AlphabetMapGen<K> {

    // Maps each n-gram (as a String) to a unique, dense id.
    public final HashMap<String,Integer> alphabetMap = new HashMap<>();

    private final int nGram;          // length of the n-gram
    private final List<K> symbols;    // the “alphabet” to combine
    private int nextId = 0;           // running id assigner

    public AlphabetMapGen(int nGram, List<K> words) {
        if (nGram <= 0) throw new IllegalArgumentException("nGram must be ≥1");
        if (words.isEmpty()) throw new IllegalArgumentException("alphabet empty");

        this.nGram   = nGram;
        this.symbols = words;
        buildMap();                   // populate alphabetMap eagerly
    }


    // Recursively builds all n-grams and puts them into  alphabetMap}.
    private void buildMap() {
        StringBuilder sb = new StringBuilder(nGram);
        dfs(0, sb);
    }

    // depth-first Cartesian product
    private void dfs(int depth, StringBuilder sb) {
        if (depth == nGram) {                     // full length reached
            alphabetMap.put(sb.toString(), nextId++);
            return;
        }
        for (K sym : symbols) {
            sb.append(sym.toString());            // add one symbol
            dfs(depth + 1, sb);                   // recurse deeper
            sb.setLength(sb.length() - 1);        // pop – reuse builder
        }
    }
}

