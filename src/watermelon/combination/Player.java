package watermelon.combination;

import java.util.*;
import java.awt.Point;
import java.lang.Math;

import watermelon.sim.Pair;
import watermelon.sim.seed;

public class Player extends watermelon.sim.Player {
  static double distowall = 1.0;
  static double distotree = 2.0;
  static double distoseed = 2.0;
  static ArrayList<Pair> trees = new ArrayList<Pair>();
  
  private Random rand = new Random();
  private double stag_distoseed = Math.tan(Math.toRadians(60.00));
  private boolean hasInitialized = false;
  private double s;
  private double w,l;

  final int MAX_ANGLE = 45; //Maximum rotation for tree packing
  final int ANGLE_INCREMENT = 15; //Iterator for rotation in tree packing
  final int SCORE_THRESHOLD = 20;
  private double angle;

  final int GRID_WIDTH = 2;

  @Override
  //Top-level move method
  public ArrayList<seed> move(ArrayList<Pair> treelist, double width, double length, double s) {

    if (!hasInitialized) {
      this.s = s;
      this.trees = treelist;
      this.w = width;
      this.l = length;
      hasInitialized = true;
    }

    ArrayList<ArrayList<seed>> solutionList = new ArrayList<ArrayList<seed>>();
    solutionList.add(altGridMove());
    solutionList.add(staggeredMove());
    solutionList.addAll(treeLayouts());
    System.out.println("Done packing");
    return chooseAltGrid(solutionList);

  }

  //Rectilinear packing
  public ArrayList<seed> altGridMove() {
    return tightPacking(false,distoseed);
  }

  //Hexagonal Packing
  public ArrayList<seed> staggeredMove() {

    return tightPacking(true,stag_distoseed);
  }

  //Generic Packing
  public ArrayList<seed> tightPacking(boolean isStaggered, double lengthIterator) {

    boolean lastime = false;
    double x=0.0;
    int rowCounter = 0;
    ArrayList<seed> seedlist = new ArrayList<seed>();
    ArrayList<seed> ghostseeds = new ArrayList<seed>();//Keep track of seeds we can't plant because of trees

    double j;
    for (j = distowall; j <= l - distowall; j = j + lengthIterator) {
      for (double i = distowall; i <= w - distowall; i = i + distoseed) {

        x = i;
        if (isStaggered && rowCounter % 2 == 1){
          if (i + 1 < w - distowall) {
            x = i + 1;
          } else {
            continue;
          }
        }

        lastime = addSeed(seedlist,ghostseeds,x,j,lastime);
      }
      rowCounter ++;
      if (rowCounter % 2 == 0) {
        lastime = !lastime;
      }
    }

    recolorTreeNeighbors(seedlist,ghostseeds);

    return seedlist;
  }

  //Tree packing
  private ArrayList<ArrayList<seed>> treeLayouts() {

    ArrayList<ArrayList<seed>> treeLayouts = new ArrayList<ArrayList<seed>>();
    Pair tree1,tree2;
    double search_angle;

    //trees that are fully inbounds
    ArrayList<Pair> treelist = new ArrayList<Pair>();
    for (Pair tree : trees) {
      if (isInBounds(tree)) {
        treelist.add(new Pair(tree.x,tree.y));
      }
    }


    if (treelist.size() < 2) {
      if (treelist.size() == 1) {
        //Just pack around one tree
        angle = 0.0;
        treeLayouts.add(boardFromPointAndAngle(treelist.get(0)));
      } 
      return treeLayouts;
    }

    int outer = 0;
    for (int i = 0; i < treelist.size() - 1; i++) {
      for (int j = i + 1; j < treelist.size(); j++) {
        tree1 = treelist.get(i);
        tree2 = treelist.get(j);
        search_angle = Math.toDegrees(Math.atan((tree2.y - tree1.y) / (tree2.x - tree1.x))) ;
        for (double k = search_angle ; k <= search_angle + MAX_ANGLE; k += ANGLE_INCREMENT){
          angle = k;
          treeLayouts.add(boardFromPointAndAngle(treelist.get(i)));
          treeLayouts.add(boardFromPointAndAngle(treelist.get(j)));
        }
        //System.out.println("Inner");
      }
      //System.out.println("Outer loop " + (outer++));
    }

    return treeLayouts;

  }

  private ArrayList<seed> boardFromPointAndAngle(Pair point) {

    //System.out.println("Starting board");
    ArrayList<seed> boardSeedList = new ArrayList<seed>();
    ArrayList<seed> boardGhostList = new ArrayList<seed>();

    seed refSeed = null;
    seed startPoint = null;

    double x_unit = distoseed * Math.cos(Math.toRadians(angle + 60.0));
    double y_unit = distoseed * Math.sin(Math.toRadians(angle + 60.0));

    ArrayList<seed> lineSeedList = lineFromPointAndAngle(point,boardGhostList);

    if (lineSeedList.size() >= 1) {
      refSeed = lineSeedList.get(0);
      boardSeedList.addAll(lineSeedList);
    } else {
      return boardSeedList;
    }

    while (lineSeedList.size() >= 1) {
      lineSeedList = addLineToBoard(lineSeedList,boardSeedList,boardGhostList,x_unit,y_unit);
    }

    lineSeedList = lineFromPointAndAngle(new Pair(refSeed.x - x_unit, refSeed.y - y_unit),boardGhostList);
    boardSeedList.addAll(lineSeedList);

    while (lineSeedList.size() >= 1) {
      lineSeedList = addLineToBoard(lineSeedList,boardSeedList,boardGhostList,(-1) * x_unit, (-1) * y_unit);
    }

    recolorTreeNeighbors(boardSeedList,boardGhostList);

    //System.out.println("Finished with board");
    return boardSeedList;
  }

  private ArrayList<seed> addLineToBoard(ArrayList<seed> line,
                                         ArrayList<seed> board, 
                                         ArrayList<seed> ghosts,
                                         double x_inc, 
                                         double y_inc){
    seed startPoint = line.get(0);
    line = lineFromPointAndAngle(new Pair(startPoint.x + x_inc, startPoint.y + y_inc),ghosts);
    board.addAll(line);
    return line;
  }

  private ArrayList<seed> lineFromPointAndAngle(Pair point, ArrayList<seed> ghosts) {
    // TODO Auto-generated method stub

    //System.out.println("Starting Line");
    ArrayList<seed> seedlist = new ArrayList<seed>();

    Pair startPoint = cornersCloserTo(point, true);

    if (startPoint == null) {
      startPoint = cornersCloserTo(point, false);
    }

    if (startPoint == null) {
      return seedlist;
    }

    double x_unit = distoseed * Math.cos(Math.toRadians(angle));
    double y_unit = distoseed * Math.sin(Math.toRadians(angle));

    treePack(seedlist,ghosts,startPoint,x_unit,y_unit);
    startPoint.x -= x_unit;
    startPoint.y -= y_unit;
    treePack(seedlist,ghosts,startPoint,(-1) * x_unit, (-1) * y_unit);

    //System.out.println("Finished with Line");
    return seedlist;
  }

  private void  treePack(ArrayList<seed> seedlist, ArrayList<seed> ghostseeds, Pair startPoint, double x_inc, double y_inc){
    boolean lastime = true;
    double i = startPoint.x;
    double j = startPoint.y;

    while (i + distowall <= w && j + distowall <= l && i - distowall >= 0.0 && j - distowall >= 0.0) {
      lastime = addSeed(seedlist, ghostseeds,i,j,lastime);
      i += x_inc;
      j += y_inc;

    }
  }

  private Pair cornersCloserTo(Pair a, boolean positive) {

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
    double aToNe = distance(a, new Pair(w, 0));
    double aToSw = distance(a, new Pair(0, l));
    double aToSe = distance(a, new Pair(w, l));

    double bToNw = distance(b, new Pair(0,0));
    double bToNe = distance(b, new Pair(w, 0));
    double bToSw = distance(b, new Pair(0, l));
    double bToSe = distance(b, new Pair(w, l));

    if (isInBounds(b)) {
      return b;
    }

    if (bToNw < aToNw || bToNe < aToNe || bToSw < aToSw || bToSe < aToSe) {
      return cornersCloserTo(b, positive);		
    }

    return null;
  }


  //Recolor seeds neighboring trees if it will increase our score
  public void recolorTreeNeighbors(ArrayList<seed> seeds, ArrayList<seed> ghosts){
    //System.out.println("About to generate Graph");
    Hashtable<seed,ArrayList<seed>> graph = generateGraph(seeds);
    LinkedList<seed> flippable = new LinkedList<seed>(); //seeds that we will consider flipping
    seed tmp;
    double highScore = 0.0;
    double temp = 0.0;


    //populate flippable seed list
    for (seed s : seeds){
      for (seed g : ghosts) {
          //northern and southern petals
          if (s.y - stag_distoseed == g.y || s.y + stag_distoseed == g.y){
            if (s.x - 1 == g.x || s.x + 1 == g.x){
              flippable.add(s);
            }
          }
          else if (s.y == g.y){
            if (s.x - distoseed == g.x || s.x + distoseed == g.x){ //W
              flippable.add(s);
            }
          }
        }
    }

    //Shuffle seeds so we don't just do top-down left-right every time
    Collections.shuffle(flippable,rand);

    highScore = calculatescore(seeds);
    while (flippable.size() > 0){
      
      //Try flipping first seed in list, adding its neighbors (and itself) back into
      //the list if our score was increased.

      tmp = flippable.poll();
      tmp.tetraploid = !tmp.tetraploid;
      temp = calculatescore(seeds);
      if (temp > highScore){
        highScore = temp;
        flippable.add(tmp);
        for (seed q : graph.get(tmp)){
          if (!flippable.contains(q))
            flippable.add(q);
        }
      }
      else{
        tmp.tetraploid = !tmp.tetraploid;
      }
    }
  }


  //Generate a seed graph using the radar technique
  public Hashtable<seed,ArrayList<seed>> generateGraph(ArrayList<seed> seedlist){
    Hashtable<Point,ArrayList<seed>> grid = new Hashtable<Point,ArrayList<seed>>();
    Hashtable<seed,ArrayList<seed>> neighbors = new Hashtable<seed,ArrayList<seed>>();
    
    int rotate;
    int box_x,box_y;
    double length;
    double x,y;
    double connect,withinReach;//distaces between seeds and end of "radar"
    boolean foundNeighbor;

    Pair endRadar = new Pair();
    Pair origRadar = new Pair();

    for (seed s : seedlist){
      neighbors.put(s, new ArrayList<seed>(seedlist.size()));
    }
    double biggestL = distoseed;

    //Create my grid
    //System.out.println("Looping");
    for (int i = 0; i < w; i+=GRID_WIDTH){
      for (int j = 0; j < l; j+=GRID_WIDTH){
        grid.put(new Point(i,j),new ArrayList<seed>());
      }
    }

    //System.out.println("Making grid");
    for (seed s : seedlist) {
      ArrayList<seed> cellList;
      Point containingCell = new Point((int)(s.x - (s.x % GRID_WIDTH)), (int) (s.y  - (s.y % GRID_WIDTH)));

      cellList = grid.get(containingCell);
      cellList.add(s);
      grid.put(containingCell,cellList);
      //Check the reach of the seed
      for (int d = 0; d < 360; d++){
        x = s.x + distowall * Math.cos(Math.toRadians(d));
        y = s.y + distowall * Math.sin(Math.toRadians(d));
        containingCell.x = (int) (x - (x % GRID_WIDTH));
        containingCell.y = (int) (y - (y % GRID_WIDTH));
        //Add new seed to grid cell
        if(isInBounds(containingCell) && !grid.get(containingCell).contains(s)){
          cellList = grid.get(containingCell);
          cellList.add(s);
          grid.put(containingCell,cellList);
        }
      }
    }

    //System.out.println("Done with grid");
      long lstart = System.currentTimeMillis();
    for (seed s : seedlist){
      //Adjacency list for this seed
      ArrayList<seed> adjacencies = neighbors.get(s);
      //Adjacency list for seeds that we find in our search
      ArrayList<seed> nAdjacencies = new ArrayList<seed>();;

      for (int degrees = 0; degrees < 360; degrees ++){
        length = distoseed;//Reset length of radar
        foundNeighbor = false;

        while (length < 45){

          //Get radar position
          x = length * Math.cos(Math.toRadians(degrees));
          y = length * Math.sin(Math.toRadians(degrees));
          origRadar = new Pair(s.x + x, s.y + y);
          box_x = (int)(s.x + x);
          box_y = (int)( s.y + y);
          Point bounding = new Point((int)(box_x - (box_x % GRID_WIDTH)), (int)(box_y - (box_y % GRID_WIDTH)));


          //If we're out of bounds, rotate radar
          if (!isInBounds(origRadar) || ! isInBounds(bounding)){
            break;
          }

          
          //System.out.println("Bounding " + bounding + " has a total number of seeds: " + grid.get(bounding).size());
          for (seed r : grid.get(bounding)){
            //System.out.println("CHecking seed " + r + " with respect to seed " + s + " at degrees " + degrees);
            if (r.equals(s)){
              if (grid.get(bounding).size() == 1)
                foundNeighbor = true;//Make sure we still break out of the while loop
              continue;
            }
            endRadar.x = origRadar.x;
            endRadar.y = origRadar.y;
            withinReach = distance(r,endRadar);
            connect = distanceseed(r,s);
            if (withinReach <= distowall){
              if (!adjacencies.contains(r)) {
                //Our radar's endpoint is within r's 1 meter radius
                adjacencies.add(r);
                foundNeighbor = true;
                if(!neighbors.get(r).contains(s))
                  nAdjacencies = neighbors.get(r);
                  nAdjacencies.add(s);
                  neighbors.put(r,nAdjacencies);
              }
              //Find the boundary of the seed that we found with the radar
              //System.out.println("Moving to boundary of seed");
              while (degrees < 360 && (radarToSeed(s,r,degrees,connect) <= distowall || !isInBounds(origRadar) || radarInBounding(origRadar,bounding))){
                degrees++;
                origRadar.x = length * Math.cos(Math.toRadians(degrees));
                origRadar.y = length * Math.sin(Math.toRadians(degrees));
              }
              //System.out.println("Moved forward " + degrees + " degrees for seed " + s);
              break;
            }
          }

          //Didn't find any seeds, so make radar longer
          if (!foundNeighbor){
            length += 0.05;
            if (length > biggestL)
              biggestL = length;
          }
          else{
            break;
          }
        }
      }
      //System.out.println("Biggest Length: " + biggestL);

      //Add our adjacency list to our graph
      neighbors.put(s,adjacencies);

    }
        long lend = System.currentTimeMillis();
    //System.out.println("Degree loop " + (lend - lstart));
    return neighbors;
  }

  public boolean radarInBounding(Pair radar, Point bounding){
    return (radar.x >= bounding.x && radar.y >= bounding.y && radar.x <= bounding.x + GRID_WIDTH && radar.y <= bounding.y + GRID_WIDTH);
  }

  //Return the distance between a seed's center and
  //a given radar
  public double radarToSeed (seed source, seed r, int angle, double connect){
    double x,y;
    Pair endRadar;
    x = connect * Math.cos(Math.toRadians(angle));
    y = connect * Math.sin(Math.toRadians(angle));
    endRadar = new Pair(source.x + x, source.y + y);
    return distance(r,endRadar);
  }

  //Add a seed unless a tree is blocking it
  public boolean addSeed(ArrayList<seed> slist, ArrayList<seed> ghosts, double x, double y, boolean ploidy){
    seed tmp;
    boolean add = true;
    tmp = new seed(x,y,ploidy);

    for (int f = 0; f < trees.size(); f++){
      if (distance(tmp, trees.get(f)) < distotree) {
        add = false;
        break;
      }
    }

    if (add) 
      slist.add(tmp);
    else
      ghosts.add(tmp);

    return !ploidy;//return ploidy of next seed to plant
  }





  //Given a rectilinear packing and a hexagonal packing, pick the one yielding
  //a higher score.
  public ArrayList<seed> chooseAltGrid(ArrayList<ArrayList<seed>> solutionList) {
    double highScore = 0.0;
    double origScore = highScore;
    double temp = 0.0;
    seed tmp;
    ArrayList<seed> finalList = null;
    ArrayList<seed> highestBoard = null;
    System.out.println(solutionList.size());
    int count = 0;
    double highestBoardScore = 0.0;
    //TODO: Redo this to the original second phase method (find highest board then recolor it)
    
    for (ArrayList<seed> solution : solutionList){
      temp = calculatescore(solution);
      if (highScore < temp){
        highScore = temp;
        finalList = solution;
      }
    }
    System.out.println("Score for board highest board is " + highScore);
    for (int i = 0; i < 200; i++){
      //Try flipping every seed
      for (seed s : finalList){
        tmp = s;
        tmp.tetraploid = !tmp.tetraploid;
        //If our change doesn't help our score, undo it.
        //If it does, randomly decide if we should undo it or keep it,
        //since the local choice might not be the best global choice.
        if (calculatescore(finalList) < highScore){
          tmp.tetraploid = !tmp.tetraploid;
        }
      }
      temp = calculatescore(finalList);
      //System.out.println("Score is " + temp);
      if (temp > highScore){
        highScore = temp;
      }
      else{
        break;
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
  
  public void init() {
  }

  static double distance(Pair p1, Pair p2){
    return distance(new seed(p1.x,p1.y,false),p2);
  }

  static double distance(seed tmp, Pair pair) {
    return Math.sqrt((tmp.x - pair.x) * (tmp.x - pair.x) + (tmp.y - pair.y) * (tmp.y - pair.y));
  }

  static double distanceseed(seed a, seed b) {
    return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
  }

  public boolean isInBounds (Pair p){
    return p.x >= 0.0 && p.y >= 0.0 && p.x <= w && p.y <= l;
  }
  public boolean isInBounds (Point p){
    return p.x >= 0.0 && p.y >= 0.0 && p.x < w && p.y < l;
  }

}