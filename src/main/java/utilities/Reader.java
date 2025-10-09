package utilities;

public interface Reader<T> {
    boolean hasNext();
    T next();
    void setQueryMode();

}
