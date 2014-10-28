package watermelon.group4;

import java.util.*;
import java.lang.Math;

import watermelon.sim.Pair;
import watermelon.sim.Point;
import watermelon.sim.seed;

public class Player extends watermelon.sim.Player {
	static double distowall = 1.0000000001;
	static double distotree = 2.0000000001;
	static double distoseed = 2.000000001;

	static final int ALT_GRID_MOVE = 0;
	static final int ALT_GRID_STAG_MOVE = 1;
	static final int CHOOSE_ALT_GRID = 2;

	private boolean hasInitialized = false;
	private int ourMove;
	private double s;

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


	static double distanceseed(seed a, seed b) {
		return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
	}

	@Override
	public ArrayList<seed> move(ArrayList<Pair> treelist, double width, double length, double s) {
		
		if (!hasInitialized) {
			this.s = s;
			hasInitialized = true;
		}

		ArrayList<ArrayList<seed>> solutionList = new ArrayList<ArrayList<seed>>();

		solutionList.add(altGridMove(treelist, width, length, s));
		solutionList.add(staggeredMove(treelist, width, length, s));
		return chooseAltGrid(solutionList);

	}


	public ArrayList<seed> altGridMove(ArrayList<Pair> treelist, double width, double length, double s) {
		// TODO Auto-generated method stub

		boolean lastime = false;

		ArrayList<seed> seedlist = new ArrayList<seed>();
		for (double i = distowall; i < width - distowall; i = i + distoseed) {
			int rowCounter = 0;
			for (double j = distowall; j < length - distowall; j = j + distoseed + .001) {
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
		return seedlist;
	}

	public ArrayList<seed> staggeredMove(ArrayList<Pair> treelist, double width, double length, double s) {
		// TODO Auto-generated method stub

		boolean lastime = false;

		ArrayList<seed> seedlist = new ArrayList<seed>();
		int rowCounter = 0;
		for (double j = distowall; j < length - distowall; j = j + Math.tan(Math.toRadians(60.00001))) {
			for (double i = distowall; i < width - distowall; i = i + distoseed) {
		
				seed tmp;
				
				double stag_i = i;
				if 	(rowCounter % 2 == 1) {
					if (i + 1 < width - distowall) {
						stag_i = i + 1;
					} else {
						continue;
					}
				}
				
				tmp = new seed(stag_i, j, lastime);
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
			}

			rowCounter ++;
			if (rowCounter % 2 == 0) {
				lastime = !lastime;
			}

		}

		return seedlist;
	}

	public ArrayList<seed> chooseAltGrid(ArrayList<ArrayList<seed>> solutionList) {
		double highScore = 0.0;
		double temp = 0.0;
		ArrayList<seed> finalList = null;
		for (ArrayList<seed> solution : solutionList){
			temp = calculatescore(solution);
			if (temp > highScore){
				highScore = temp;
				finalList = solution;
			}
		}
		
		return finalList;
	}

	private double calculatescore(ArrayList<seed> seedlist) {
		double total = 0.0;

		for (int i = 0; i < seedlist.size(); i++) {
			double score;
			double chance = 0.0;
			double totaldis = 0.0;
			double difdis = 0.0;
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i) {
					totaldis = totaldis
							+ Math.pow(
									distanceseed(seedlist.get(i),
											seedlist.get(j)), -2);
				}
			}
			for (int j = 0; j < seedlist.size(); j++) {
				if (j != i
						&& ((seedlist.get(i).tetraploid && !seedlist.get(j).tetraploid) || (!seedlist
								.get(i).tetraploid && seedlist.get(j).tetraploid))) {
					difdis = difdis
							+ Math.pow(
									distanceseed(seedlist.get(i),
											seedlist.get(j)), -2);
				}
			}
			//System.out.println(totaldis);
			//System.out.println(difdis);
			chance = difdis / totaldis;
			score = chance + (1 - chance) * s;
			total = total + score;
		}
		return total;
	}
}