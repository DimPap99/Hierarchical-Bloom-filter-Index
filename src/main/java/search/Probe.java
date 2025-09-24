// search/Probe.java
package search;

public final class Probe {
    private final int  consumed;
    private final boolean complete;
    private final boolean testedFirst;

    public Probe(int consumed, boolean complete) {
        this(consumed, complete, true); // preserve old ctor semantics
    }
    public Probe(int consumed, boolean complete, boolean testedFirst) {
        this.consumed = consumed;
        this.complete = complete;
        this.testedFirst = testedFirst;
    }
    public int consumed()      { return consumed; }
    public boolean complete()  { return complete; }
    public boolean testedFirst(){ return testedFirst; }
}
