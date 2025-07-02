package utilities;

public class SlidingStr {

    int k;
    StringBuilder str = new StringBuilder();
    public SlidingStr(int k) {
        this.k = k;
    }



    public boolean isFull(){
        return str.length() >= this.k;
    }
    public void append(char c){
        str.append(c);
        if(str.length()>this.k) str.deleteCharAt(0);
    }


    public String getStr(){
        return this.str.toString();
    }
}
