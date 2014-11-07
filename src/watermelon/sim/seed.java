package watermelon.sim;

import java.util.ArrayList;


public class seed {
    public double x;
    public double y;
    public boolean tetraploid;
    public double score;
    public boolean marked;
    
    public seed() { x = 0.0; y = 0.0; tetraploid = false; marked = false;}

    public seed(double xx, double yy, boolean tetra) {
        x = xx;
        y = yy;
        tetraploid = tetra;
        marked = false;
    }
    public String toString() {
      return "(" + x + ", " + y + ")";
    }
}
