package watermelon.michael;

import java.util.*;
import java.lang.Math;

import watermelon.sim.Pair;
import watermelon.sim.Point;
import watermelon.sim.seed;

public class Player extends watermelon.sim.Player {
	static double distowall = 1.00;
	static double distotree = 2.00;
	static double distoseed = 2.00;

	static final int ALT_GRID_MOVE = 0;
	static final int ALT_GRID_STAG_MOVE = 1;
	static final int CHOOSE_ALT_GRID = 2;

	private boolean hasInitialized = false;
	private int ourMove;
	private double s;
	private double width;
	private double length;
	private ArrayList<Pair> treelist;

	public void init() {
	}

	static double distance(seed tmp, Pair pair) {
		return Math.sqrt((tmp.x - pair.x) * (tmp.x - pair.x) + (tmp.y - pair.y) * (tmp.y - pair.y));
	}

	static double distance(Pair tmp, Pair pair) {
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
			this.width = width;
			this.length = length;
			this.treelist = treelist;
			hasInitialized = true;
		}

		ArrayList<ArrayList<seed>> solutionList = new ArrayList<ArrayList<seed>>();

		solutionList.add(altGridMove(treelist, width, length, s));
		solutionList.add(staggeredMove(treelist, width, length, s));
		solutionList.addAll(treeLayouts());
		return chooseAltGrid(solutionList);

	}


	public ArrayList<seed> altGridMove(ArrayList<Pair> treelist, double width, double length, double s) {
		// TODO Auto-generated method stub

		boolean lastime = false;

		ArrayList<seed> seedlist = new ArrayList<seed>();
		int rowCounter = 0;
		for (double j = distowall; j <= length - distowall; j = j + distoseed) {
			for (double i = distowall; i <= width - distowall; i = i + distoseed) {

				seed tmp;
				tmp = new seed(i, j, lastime);
				boolean add = true;
				for (int f = 0; f < treelist.size(); f++) {
					if (distance(tmp, treelist.get(f)) <= distotree) {
						add = false;
						break;
					}
				}
				if (add) {
					seedlist.add(tmp);
				}
				
				lastime = !lastime;
			}
			lastime = !lastime;
		}
		return seedlist;
	}

	public ArrayList<seed> staggeredMove(ArrayList<Pair> treelist, double width, double length, double s) {
		// TODO Auto-generated method stub

		boolean lastime = false;

		ArrayList<seed> seedlist = new ArrayList<seed>();
		int rowCounter = 0;
		for (double j = distowall; j <= length - distowall; j = j + Math.tan(Math.toRadians(60.00))) {
			for (double i = distowall; i <= width - distowall; i = i + distoseed) {
		
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
					if (distance(tmp, treelist.get(f)) <= distotree) {
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

	private ArrayList<ArrayList<seed>> treeLayouts() {


		ArrayList<ArrayList<seed>> treeLayouts = new ArrayList<ArrayList<seed>>();

		ArrayList<Pair> treelist = filterTreeList(this.treelist);

		if (treelist.size() < 2) {
			if (treelist.size() == 1) {
				treeLayouts.add(boardFromPointAndAngle(treelist.get(0), 0.0));
			} else {
				treeLayouts.add(boardFromPointAndAngle(new Pair(1.0,1.0), 0.0));
			}
			
			return treeLayouts;
		}

		for (int i = 0; i < treelist.size() - 1; i++) {
			for (int j = i + 1; j < treelist.size(); j++) {
				double angleClass = angleClassBetweenPoints(treelist.get(i), treelist.get(j));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(i), angleClass));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(i), angleClass + 15));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(i), angleClass + 30));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(i), angleClass + 45));

				treeLayouts.add(boardFromPointAndAngle(treelist.get(j), angleClass));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(j), angleClass + 15));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(j), angleClass + 30));
				treeLayouts.add(boardFromPointAndAngle(treelist.get(j), angleClass + 45));
			}
		}

		return treeLayouts;

	}

	private ArrayList<seed> boardFromPointAndAngle(Pair point, double angle) {
		
		ArrayList<seed> boardSeedList = new ArrayList<seed>();

		ArrayList<seed> lineSeedList = lineFromPointAndAngle(point, angle);
		seed refSeed = null;

		if (lineSeedList.size() >= 1) {
			refSeed = lineSeedList.get(0);
			boardSeedList.addAll(lineSeedList);
		} else {
			return boardSeedList;
		}

		double x_unit = distoseed * Math.cos(Math.toRadians(angle + 60.0));
		double y_unit = distoseed * Math.sin(Math.toRadians(angle + 60.0));

		while (lineSeedList.size() >= 1) {
			seed startPoint = lineSeedList.get(0);
			lineSeedList = lineFromPointAndAngle(new Pair(startPoint.x + x_unit, startPoint.y + y_unit), angle);
			boardSeedList.addAll(lineSeedList);
		}

		lineSeedList = lineFromPointAndAngle(new Pair(refSeed.x - x_unit, refSeed.y - y_unit), angle);

		while (lineSeedList.size() >= 1) {
			boardSeedList.addAll(lineSeedList);
			seed startPoint = lineSeedList.get(0);
			lineSeedList = lineFromPointAndAngle(new Pair(startPoint.x - x_unit, startPoint.y - y_unit), angle);
		}

		return boardSeedList;
	}

	private ArrayList<seed> lineFromPointAndAngle(Pair point, double angle) {
		// TODO Auto-generated method stub

		ArrayList<seed> seedlist = new ArrayList<seed>();

		Pair startPoint = cornersCloserTo(point, angle, true);
		
		if (startPoint == null) {
			startPoint = cornersCloserTo(point, angle, false);
		}

		if (startPoint == null) {
			return seedlist;
		}

		double i = startPoint.x;
		double j = startPoint.y;

		double x_unit = distoseed * Math.cos(Math.toRadians(angle));
		double y_unit = distoseed * Math.sin(Math.toRadians(angle));

		boolean lastime = true;
		while (i + distowall <= width && j + distowall <= length && i - distowall >= 0.0 && j - distowall >= 0.0) {

				seed tmp = new seed(i, j, lastime);
				boolean add = true;
				for (int f = 0; f < treelist.size(); f++) {
					if (distance(tmp, treelist.get(f)) <= distotree) {
						add = false;
						break;
					}
				}
				
				if (add) {
					seedlist.add(tmp);
				}

				lastime = !lastime;
				i += x_unit;
				j += y_unit;

		}

		i = startPoint.x - x_unit;
		j = startPoint.y - y_unit;

		lastime = true;
		while (i + distowall <= width && j + distowall <= length && i - distowall >= 0.0 && j - distowall >= 0.0) {

				seed tmp = new seed(i, j, lastime);
				boolean add = true;
				for (int f = 0; f < treelist.size(); f++) {
					if (distance(tmp, treelist.get(f)) <= distotree) {
						add = false;
						break;
					}
				}
				
				if (add) {
					seedlist.add(tmp);
				}

				lastime = !lastime;
				i -= x_unit;
				j -= y_unit;

		}		
		return seedlist;
	}

	private Pair cornersCloserTo(Pair a, double angle, boolean positive) {
		
		if (isInBounds(a)) {
			return a;
		}

		double x_unit = distoseed * Math.cos(Math.toRadians(angle));
		double y_unit = distoseed * Math.sin(Math.toRadians(angle));

		if (!positive) {
			x_unit = (-1) * x_unit;
			y_unit = (-1) * y_unit;
		}
		
		Pair b = new Pair(a.x + x_unit, a.y + y_unit);

		double aToNw = distance(a, new Pair(0,0));
		double aToNe = distance(a, new Pair(width, 0));
		double aToSw = distance(a, new Pair(0, length));
		double aToSe = distance(a, new Pair(width, length));

		double bToNw = distance(b, new Pair(0,0));
		double bToNe = distance(b, new Pair(width, 0));
		double bToSw = distance(b, new Pair(0, length));
		double bToSe = distance(b, new Pair(width, length));

		if (isInBounds(b)) {
			return b;
		}

		if (bToNw < aToNw || bToNe < aToNe || bToSw < aToSw || bToSe < aToSe) {
			return cornersCloserTo(b, angle, positive);		
		}

		return null;
	}

	private ArrayList<Pair> filterTreeList(ArrayList<Pair> treelist) {

		ArrayList<Pair> filterList = new ArrayList<Pair>();

		for (Pair tree : treelist) {
			if (isInBounds(tree)) {
				filterList.add(tree);
			}
		}

		return filterList;
	}

	private boolean isInBounds(Pair point) {
		if (point.x + distowall <= width && point.y + distowall <= length && point.x - distowall >= 0.0 && point.y - distowall >= 0.0) {
			return true;
		}

		return false;
	}

	private double angleClassBetweenPoints(Pair a, Pair b) {
		double angle = Math.toDegrees(Math.atan((b.y - a.y)  / (b.x - a.x)));

		// while (angle < 0) {
		// 	angle += 360.0;
		// }

		return angle;
	}

	private ArrayList<seed> iterativeColoring(ArrayList<seed> originalBoard) {

		System.out.println("Score before coloring: " + calculatescore(originalBoard));

		int iterations = 20;
		for (int i = 0; i < iterations; i++) {
			Collections.shuffle(originalBoard);
			for (int j = 0; j < originalBoard.size(); j++) {
				double scoreBeforeFlip = calculatescoreForSeed(j, originalBoard);
				originalBoard.get(j).tetraploid = !originalBoard.get(j).tetraploid;
				double scoreAfterFlip = calculatescoreForSeed(j, originalBoard);

				if (scoreAfterFlip < scoreBeforeFlip) {
					originalBoard.get(j).tetraploid = !originalBoard.get(j).tetraploid;
				}
			}
		}

		return originalBoard;
	}


	public ArrayList<seed> chooseAltGrid(ArrayList<ArrayList<seed>> solutionList) {
		double highScore = 0.0;
		double temp = 0.0;
		ArrayList<seed> finalList = null;
		for (ArrayList<seed> solution : solutionList){
			solution = iterativeColoring(solution);
			temp = calculatescore(solution);
			System.out.println("Score for board: " + temp);
			if (temp >= highScore){
				highScore = temp;
				finalList = solution;
			}
			System.out.println("High score: " + highScore);
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

	private double calculatescoreForSeed(int i, ArrayList<seed> seedlist) {
		double total = 0.0;
	
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

		return total;
	}
}