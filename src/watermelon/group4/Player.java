package watermelon.group4;

import java.util.*;

import watermelon.sim.Pair;
import watermelon.sim.Point;
import watermelon.sim.seed;

public class Player extends watermelon.sim.Player {
	static double distowall = 1.00;
	static double distotree = 1.00;
	static double distoseed = 2.00;

	public void init() {

	}

	static double distance(seed tmp, Pair pair) {
		return Math.sqrt((tmp.x - pair.x) * (tmp.x - pair.x) + (tmp.y - pair.y) * (tmp.y - pair.y));
	}

	// Return: the next position
	// my position: dogs[id-1]
	static double distance(seed a, Point b) {
		return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
	}

	@Override
	public ArrayList<seed> move(ArrayList<Pair> treelist, double width, double length, double s) {
		// TODO Auto-generated method stub

		boolean lastime = false;

		ArrayList<seed> seedlist = new ArrayList<seed>();
		for (double i = distowall; i < width - distowall; i = i + distoseed) {
			int rowCounter = 0;
			for (double j = distowall; j < length - distowall; j = j + distoseed) {
				seed tmp;
				tmp = new seed(i, j, lastime);
				boolean add = true;
				for (int f = 0; f < treelist.size(); f++) {
					if (distance(tmp, treelist.get(f)) < distotree) {
						add = false;
						break;
					}
				}
				if (add) {
					seedlist.add(tmp);
				}
				
				lastime = !lastime;
				rowCounter++;
			}
			if (rowCounter % 2 == 0) {
				lastime = !lastime;
			}

		}
		System.out.printf("seedlist size is %d", seedlist.size());
		return seedlist;
	}
}