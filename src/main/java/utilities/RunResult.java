package utilities;

import java.util.ArrayList;

public record RunResult(double alpha, double runAvgMs, double avpLp){


    public void print(){
        System.out.println(alpha +","+runAvgMs+","+avpLp);
    }
}
