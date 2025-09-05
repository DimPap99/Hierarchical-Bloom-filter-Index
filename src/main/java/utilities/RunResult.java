package utilities;

import java.util.ArrayList;

public record RunResult(double alpha, double runAvgMs, double avpLp, int probes){


    public void print(){
        System.out.println(alpha +","+runAvgMs+","+avpLp+ "," +probes);
    }
}
