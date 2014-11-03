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

  private double stag_distoseed = Math.tan(Math.toRadians(60.00001));
  private boolean hasInitialized = false;
  private double s;

  public void init() {
  }

  static double distance(seed tmp, Pair pair) {
    return Math.sqrt((tmp.x - pair.x) * (tmp.x - pair.x) + (tmp.y - pair.y) * (tmp.y - pair.y));
  }

  static double distance(seed a, Point b) {
    return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
  }

  static double distanceseed(seed a, seed b) {
    return Math.sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y));
  }

  @Override
  //Top-level move method
  public ArrayList<seed> move(ArrayList<Pair> treelist, double width, double length, double s) {

    if (!hasInitialized) {
      this.s = s;
      this.trees = treelist;
      hasInitialized = true;
    }

    ArrayList<ArrayList<seed>> solutionList = new ArrayList<ArrayList<seed>>();

    solutionList.add(altGridMove(width, length, s));
    solutionList.add(staggeredMove(width, length, s));
    return chooseAltGrid(solutionList);

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


  //Rectilinear packing
  public ArrayList<seed> altGridMove(double width, double length, double s) {

    boolean lastime = false;
    int rowCounter;
    ArrayList<seed> seedlist = new ArrayList<seed>();
    ArrayList<seed> ghostseeds = new ArrayList<seed>();//Keep track of seeds we can't plant because of trees

    for (double i = distowall; i <= width - distowall; i = i + distoseed) {
      rowCounter = 0;
      for (double j = distowall; j <= length - distowall; j = j + distoseed) {
        lastime = addSeed(seedlist,ghostseeds,i,j,lastime);
        rowCounter++;
      }

      if (rowCounter % 2 == 0) {
        lastime = !lastime;
      }

    }
    recolorTreeNeighbors(seedlist,ghostseeds);
    return seedlist;
  }


  //Hexagonal Packing
  public ArrayList<seed> staggeredMove(double width, double length, double s) {

    boolean lastime = false;
    double stag_i;
    int rowCounter = 0;
    ArrayList<seed> seedlist = new ArrayList<seed>();
    ArrayList<seed> ghostseeds = new ArrayList<seed>();//Keep track of seeds we can't plant because of trees

    for (double j = distowall; j <= length - distowall; j = j + stag_distoseed) {
      for (double i = distowall; i <= width - distowall; i = i + distoseed) {
        stag_i = i;
        if (rowCounter % 2 == 1) {
          if (i + 1 < width - distowall) {
            stag_i = i + 1;
          } else {
            continue;
          }
        }

        lastime = addSeed(seedlist,ghostseeds,stag_i,j,lastime);
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
    Collections.shuffle(flippable);

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
        addNeighbors(tmp,seeds,flippable);
      }
      else{
        tmp.tetraploid = !tmp.tetraploid;
      }
    }
  }

  //TODO: Do we want to add "neighbors" across trees as well?
  public void addNeighbors(seed s, ArrayList<seed> seeds, LinkedList<seed> flippable){
    for (seed sd : seeds){
        if (s.y - stag_distoseed == sd.y || s.y + stag_distoseed == sd.y){
          if (s.x - 1 == sd.x || s.x + 1 == sd.x){
            flippable.add(sd);
          }
        }
        else if (s.y == sd.y){
          if (s.x - distoseed == sd.x || s.x + distoseed == sd.x){ 
            flippable.add(sd);
          }
        }
    }
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

    //TODO: Try some random flipping as that seems to increase our score a bit
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
