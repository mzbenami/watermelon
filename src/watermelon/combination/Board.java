package watermelon.combination;

import java.util.*;
import watermelon.sim.Pair;
import watermelon.sim.seed;

import java.awt.Point;

public class Board {

	public ArrayList<seed> seedlist;
	public ArrayList<seed> ghostlist;
	public Hashtable<Point,ArrayList<seed>> grid;

	public double row_space;
	public double column_space;

	public Board(ArrayList<seed> seedlist, ArrayList<seed> ghostlist, Hashtable<Point, ArrayList<seed>> grid) {
		this.seedlist = seedlist;
		this.ghostlist = ghostlist;
		this.grid = grid;
	}

	public Board(ArrayList<seed> seedlist) {
		this.seedlist = seedlist;
	}

	public Board(ArrayList<seed> seedlist, ArrayList<seed> ghostlist) {
		this.seedlist = seedlist;
		this.ghostlist = ghostlist;
	}
}