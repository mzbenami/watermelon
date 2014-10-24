package watermelon.sim;

import java.util.ArrayList;

import watermelon.sim.seed;

public abstract class Player {

	public Player() {
	}

	public abstract void init();

	public abstract ArrayList<seed> move(ArrayList<Pair> treelist, double w,
			double l, double s); // positions of all the outpost, playerid

}
