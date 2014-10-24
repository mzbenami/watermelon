package watermelon.sim;

import watermelon.sim.Pair;

public class Pair {
    public double x;
    public double y;

    public Pair() { x = 0.0; y = 0.0; }

    public Pair(double xx, double yy) {
        x = xx;
        y = yy;
    }

    public Pair(Pair o) {
        this.x = o.x;
        this.y = o.y;
    }

    public boolean equals(Pair o) {
        return o.x == x && o.y == y;
    }
}