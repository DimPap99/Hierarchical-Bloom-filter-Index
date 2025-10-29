package utilities;

import java.util.HashMap;

public class AlphabetMapper <T>{
    int nextId;//id 0 is reserved as a sentinel for the suffix tree DC3
    float loadFactor = 0.75f;
    public HashMap<T, Integer> wordToId;
//    public HashMap<Integer, T> idToWord;


    int capacity;
    public AlphabetMapper(int capacity) {
        this.nextId = 1;//id 0 is reserved as a sentinel for the suffix tree DC3
            this.capacity = capacity;

            this.wordToId = new HashMap<>();
//            this.idToWord = new HashMap<>();

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
//            idToWord.put(id, item);
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
//        idToWord.clear();
        nextId = 0;
    }


}
