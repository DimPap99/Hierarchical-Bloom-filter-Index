package search;

public final class Probe {
    private int consumed;
    private boolean complete;
    private boolean testedFirst;

    public Probe() {}

    public Probe(int consumed, boolean complete) {
        this(consumed, complete, true);
    }

    public Probe(int consumed, boolean complete, boolean testedFirst) {
        this.consumed = consumed;
        this.complete = complete;
        this.testedFirst = testedFirst;
    }

    public void set(int consumed, boolean complete, boolean testedFirst) {
        this.consumed = consumed;
        this.complete = complete;
        this.testedFirst = testedFirst;
    }

    public int consumed()       { return consumed; }
    public boolean complete()   { return complete; }
    public boolean testedFirst(){ return testedFirst; }
}
