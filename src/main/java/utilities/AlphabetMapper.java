package utilities;

import java.util.HashMap;

public class AlphabetMapper <T>{
    int nextId;//init id is -1

    float loadFactor = 0.75f;
    public HashMap<T, Integer> wordToId;
    public HashMap<Integer, T> idToWord;

    //TODO: for hashmap space allocation maybe follow a strategy as cpp vector. Start with init capacity and double once we reach it
    // ---> with a lot of ngrams still explodes memory
    int capacity;
    public AlphabetMapper(int capacity) {
        this.nextId = 0;
//        if(capacity > 0){
            this.capacity = capacity;
            //pre assign enough space to hashmap so that we avoid resizes
//            int internalCapacity = Math.max(16, (int) Math.ceil(capacity / loadFactor));
            this.wordToId = new HashMap<>();
            this.idToWord = new HashMap<>();
//        }else{
//            throw new IllegalArgumentException("Negative or zero capacity");
//        }
    }

    public int getSize() {
        return wordToId.size();
    }
    public int getCapacity() {
        return capacity;
    }
    public int insert(T item){
        int id = wordToId.getOrDefault(item, -1);

        if(id == -1){
            id = nextId;
            this.nextId++;
//            if(id >= this.capacity){ throw new IllegalStateException("Exceeded capacity.");}
            wordToId.put(item, id);
            idToWord.put(id, item);
        }
        return id;
    }
    
    //We ensure programmatically wherever this is used that the key is in fact inside the map
    public int getId(T item){
        int retVal =  wordToId.getOrDefault(item, -999);
        if(retVal == -999) retVal =  insert(item);
        return retVal;
    }

    public void clear() {
        wordToId.clear();
        idToWord.clear();
        nextId = 0;
    }


}
