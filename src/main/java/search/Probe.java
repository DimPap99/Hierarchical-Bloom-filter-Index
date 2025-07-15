package search;

record Probe(int consumed, boolean complete) {
    public int consumed() { return consumed; }
    public boolean complete() { return complete; }
} // complete==true if no miss in this node

