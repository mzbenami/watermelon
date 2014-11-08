package watermelon.ghostseed;

import java.util.*;
import java.lang.Math;

import watermelon.sim.Pair;
import watermelon.sim.Point;
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


  //Recolor seeds neighboring trees if it will increase our score
  public void recolorTreeNeighbors(ArrayList<seed> seeds, ArrayList<seed> ghosts){
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
    Hashtable<seed,ArrayList<seed>> neighbors = new Hashtable<seed,ArrayList<seed>>();
    
    int deg,rotate;
    double length;
    double x,y;
    double connect,withinReach;//distaces between seeds and end of "radar"
    boolean foundNeighbor;

    Pair endRadar = new Pair();
    Pair origRadar = new Pair();

    //Represents every degree we have to search per seed
    LinkedList<Integer> degrees = new LinkedList<Integer>();
    LinkedList<Integer> fullDegrees = new LinkedList<Integer>();

    for (int i = 0; i < 360; i++)
      fullDegrees.add(new Integer(i));

    for (seed s : seedlist){
      //Reset degrees list and length of radar
      degrees = (LinkedList<Integer>) fullDegrees.clone();
      length = distoseed;

      //Adjacency list for this seed
      ArrayList<seed> adjacencies = new ArrayList<seed>();

      //This determines when we increment our radar
      int counter = 0;

      while (!degrees.isEmpty() && length < 45){
        //Set new iteration values
        deg = degrees.getFirst();
        foundNeighbor = false;
        counter++;

        //Increase our radar length if we've made it through all our degrees once
        if (counter > degrees.size()){
          length += 0.05;
          counter = 0;
        }

        //Get radar position
        x = length * Math.cos(Math.toRadians(deg));
        y = length * Math.sin(Math.toRadians(deg));
        origRadar = new Pair(s.x + x, s.y + y);

        //If we're out of bounds, remove this degree from our search
        if (!isInBounds(origRadar)){
          degrees.remove(new Integer(deg));
          continue;
        }

        for (seed r : seedlist){
          if (r.equals(s))
            continue;
          endRadar.x = origRadar.x;
          endRadar.y = origRadar.y;
          withinReach = distance(r,endRadar);
          connect = distanceseed(r,s);
          if (withinReach <= distowall && !adjacencies.contains(r)) {
            //Our radar's endpoint is within r's 1 meter radius
            adjacencies.add(r);
            rotate = deg;

            //Find the boundary of the seed that we found with the radar
            while (radarToSeed(s,r,rotate,connect) <= distowall){
              rotate = (rotate + 1) % 360;
            }

            //Remove all degrees that are covered by this seed
            //by rotating the radar in the other direction
            degrees.remove(new Integer(rotate));
            deg = rotate - 1;
            while (radarToSeed(s,r,deg,connect) <= distowall){
              degrees.remove(new Integer(deg));
              deg = (deg - 1) % 360;
              if (deg < 0)
                deg += 360;
            }
            degrees.remove(new Integer(deg));
            counter = 0;//We removed degrees from our search, so reset to zero
            foundNeighbor = true;
            break;
          }
        }

        //Didn't find any seeds, so move this degree to the end of our list
        //so we'll search it again with a longer radar.
        if (!foundNeighbor){
          degrees.addLast(degrees.poll());
        }
      }

      //Add our adjacency list to our graph
      neighbors.put(s,adjacencies);

    }
    return neighbors;
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
    double temp = 0.0;
    ArrayList<seed> finalList = null;
    for (ArrayList<seed> solution : solutionList){
      temp = calculatescore(solution);
      if (temp > highScore){
        highScore = temp;
        finalList = solution;
      }
    }
    seed tmp;
    ArrayList<seed> actual = (ArrayList<seed>)finalList.clone();
    //Do this process 20 ttimes, picking the best result
    for (int i = 0; i < 200; i++){
      //Try flipping every seed
      for (seed s : actual){
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
      System.out.println("Score is " + temp);
      if (temp > highScore){
        highScore = temp;
      }
      else
        break;
    }

    return actual;
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

  static double distance(seed tmp, Pair pair) {
    return Math.sqrt((tmp.x - pair.x) * (tmp.x - pair.x) + (tmp.y - pair.y) * (tmp.y - pair.y));
  }

  static double distanceseed(seed a, seed b) {
    return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
  }

  public boolean isInBounds (Pair p){
    return p.x >= 0.0 && p.y >= 0.0 && p.x <= w && p.y <= l;
  }

}
