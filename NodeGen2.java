import java.util.*;
import java.util.function.Function;
import java.awt.Color;

public class NodeGen {

  // stores the location of each player and each pf their colors
  public static HashMap<Point, Color> playerColors = null;

  public static void main(String[] args) {
    int size = 10;
    if (args.length != 0) {
      try {
        size = Integer.parseInt(args[0]);
      } catch (Exception e) {
        Core.exit("Please enter a valid size as number");
      }
    }

    // declare space
    int w = 10;
    int h = 10;

    // phase 1.1
    ArrayList<Point> ps = generatePoints(w, h, size);

    // draw the points
    int pixls = 80;
    StdDraw.setCanvasSize(w * pixls, h * pixls);
    StdDraw.setXscale(0, w);
    StdDraw.setYscale(0, h);

    // phase 1.2
    HashMap<Point, ArrayList<Point>> G = createGraph(ps);

    // phase 2.1
    int players = 4;
    HashMap<Point, ArrayList<Point>> graph = generatePlayerSpawns(G, players);

    // create array list of all player spawns
    ArrayList<Point> vlis = new ArrayList<>();
    for (Point p : playerColors.keySet()) {
      vlis.add(p);
    }
    // phase 2.2
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs = CoreGeom.voronoiGraph(vlis, graph);
    for (Point player : subgraphs.keySet()) {
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(player);
      drawEdges(subgraph, playerColors.get(player));
    }


    while (!acceptable(vlis, subgraphs)) {
      // 1) first set spawn points at centroids

      Core.freeze();
      // remove spawn points
      for (Point spawn : vlis) {
        ArrayList<Point> e = graph.get(spawn);
        if (e.size() != 2) {
          Core.exit("e.size() != 2");
        }
        Point left = e.get(0);
        Point right = e.get(1);
        // severing its connections between its edge's endpoints
        graph.get(spawn).remove(left);
        graph.get(left).remove(spawn);
        graph.get(spawn).remove(right);
        graph.get(right).remove(spawn);
        // removing spawn
        graph.remove(spawn);
        // connecting the endpoints to each other (as they originally were)
        graph.get(left).add(right);
        graph.get(right).add(left);
      }

      StdDraw.clear(StdDraw.WHITE);
      drawGraph(graph, StdDraw.BLACK);
      Core.freeze("removed spawn points");


      // Get the centroid of each graphs
      HashMap<Point, Point[]> centroids = new HashMap<>();
      for (Point player : subgraphs.keySet()) {
        HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(player);


        ArrayList<Point> e = subgraph.get(player);
        if (e.size() != 2) {
          Core.exit("e.size() != 2");
        }
        Point left = e.get(0);
        Point right = e.get(1);
        // severing its connections between its edge's endpoints
        subgraph.get(player).remove(left);
        subgraph.get(left).remove(player);
        subgraph.get(player).remove(right);
        subgraph.get(right).remove(player);
        // removing player
        subgraph.remove(player);
        // connecting the endpoints to each other (as they originally were)
        subgraph.get(left).add(right);
        subgraph.get(right).add(left);
        // now find the centroid

        Point[] centroi = CoreGeom.graphCentroid(subgraph);
        Point centroid = centroi[1];
        // draw centroid
        Color color = playerColors.get(player);
        StdDraw.filledSquare(centroid, 0.05, color);
        // look for which edge it belongs to
        // calc which edges of original graph have the centroids
        int i = 0;
        for (Point a : graph.keySet()) {
          for (Point b : graph.get(a)) {
            if (Core.epsilon(a.dist(centroid) + centroid.dist(b), a.dist(b))) {
              // then we have found the edge
              // save centroid
              centroids.put(player, new Point[]{a, centroid, b});
              i++;
              if (i > 2) {
                Core.freeze("I've got a bad feeling about this...");
              }
            }
          }
        }
      }

      // now we have the centroids and where they are
      Core.freeze("these are the centroids");

      // change spawn points to coordinates of centroid
      for (Point spawn : vlis) {
        // get player centroid
        Point[] centroi = centroids.get(spawn);
        Point left = centroi[0], centroid = centroi[1], right = centroi[2];
        // set spawn's coordinates
        spawn.setTo(centroid);
        // sever left from right
        graph.get(left).remove(right);
        graph.get(right).remove(left);
        // add spawn
        graph.put(spawn, new ArrayList<Point>());
        // connect left and right to centroid
        graph.get(spawn).add(left);
        graph.get(left).add(spawn);
        graph.get(spawn).add(right);
        graph.get(right).add(spawn);
      }

      StdDraw.clear(StdDraw.WHITE);
      drawGraph(graph, StdDraw.BLACK);
      for (Point player : vlis) {
        StdDraw.filledCircle(player, 0.05, playerColors.get(player));
      }
      Core.freeze("added new spawn points as centroids");

      // 2) now get the new subgraphs
      subgraphs = CoreGeom.voronoiGraph(vlis, graph);
      Core.log("-");
    }
  }

  public static boolean acceptable(
    ArrayList<Point> vlis,
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs) {

    // get the lengths of each graph.
    int players = vlis.size();
    double[] lengths = new double[players];
    for (int i = 0; i < players; i++) {
      Point player = vlis.get(i);
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(player);
      lengths[i] = CoreGeom.graphLength(subgraph);
    }

    // get stdev
    double stddev = StdStats.stddev(lengths);
    double mean = StdStats.mean(lengths);
    // check whether any are out of range
    for (double length : lengths) {
      if (Math.abs(length - mean) > stddev) {
        // report unnacceptable:
        Core.log("unnacceptable");
        System.out.println("lengths: ");
        for (int i = 0; i < players; i++) {
          System.out.println("\t"+StdDraw.colorName(playerColors.get(vlis.get(i)))+" = "+lengths[i]);
        }
        System.out.println();
        System.out.println("stddev = " + stddev);
        System.out.println("mean = " + mean);

        return false;
      }
    }

    Core.log("ACCEPTED Graph!!");
    System.out.println("lengths: ");
    for (int i = 0; i < players; i++) {
      System.out.println("\t"+StdDraw.colorName(playerColors.get(vlis.get(i)))+" = "+lengths[i]);
    }
    System.out.println();
    System.out.println("stddev = " + stddev);
    System.out.println("mean = " + mean);
    return true;
  }

  // phase 2.1
  public static HashMap<Point, ArrayList<Point>> generatePlayerSpawns(
    HashMap<Point, ArrayList<Point>> G, int players) {
    // now onto phase 2:
    // we must perform voronoi on the map
    // first choose random places to place players
    int amnt = players;
    // copy graph of G
    HashMap<Point, ArrayList<Point>> H = new HashMap<>();
    for (Point p : G.keySet()) {
      ArrayList<Point> adj = new ArrayList<>();
      for (Point a : G.get(p)) {
        adj.add(a);
      }
      H.put(p, adj);
    }
    Point[] playerSpawns = new Point[4];
    // make copy of graph
    HashMap<Point, ArrayList<Point>> graph = new HashMap<>();
    for (Point p : G.keySet()) {
      ArrayList<Point> adj = new ArrayList<>();
      for (Point a : G.get(p)) {
        adj.add(a);
      }
      graph.put(p, adj);
    }
    while (amnt > 0) {
      // choose random point in H
      Point rand = Core.randomKey(H);
      // choose a random edge
      Point adj = H.get(rand).get((int)(Math.random()*H.get(rand).size()));
      // add new spread point there
      double randDouble = Math.random();
      Core.log("randDouble = " + randDouble);
      // Point playerSpawn = new Point((rand.x + adj.x)*randDouble, (rand.y + adj.y)*randDouble);
      Point playerSpawn = new Point(rand.x*randDouble + adj.x*(1-randDouble), rand.y*randDouble + adj.y*(1-randDouble));
      // add to player spawns
      playerSpawns[--amnt] = playerSpawn;
      // remove them from H
      ArrayList<Point> adjToRand = new ArrayList<>();
      for (Point p : H.get(rand)) {
        adjToRand.add(p);
      }
      ArrayList<Point> adjToAdj = new ArrayList<>();
      for (Point p : H.get(adj)) {
        adjToAdj.add(p);
      }
      for (Point p : adjToRand) {
        H.get(rand).remove(p);
        H.get(p).remove(rand);
      }
      for (Point p : adjToAdj) {
        H.get(adj).remove(p);
        H.get(p).remove(adj);
      }
      H.remove(rand);
      H.remove(adj);

      // now make the edit to graph
      // add spawn to graph
      graph.put(playerSpawn, new ArrayList<Point>());
      // add connection to adj and rand
      graph.get(playerSpawn).add(rand);
      graph.get(playerSpawn).add(adj);
      graph.get(rand).add(playerSpawn);
      graph.get(adj).add(playerSpawn);
      // sever connection between rand and adj
      graph.get(rand).remove(adj);
      graph.get(adj).remove(rand);
    }

    StdDraw.clear(StdDraw.WHITE);
    drawGraph(graph);
    if (isBroken(graph) != null) {
      Core.exit("broke graph pre post_fix.png");
    }

    playerColors = new HashMap<>();
    playerColors.put(playerSpawns[0], StdDraw.RED);
    playerColors.put(playerSpawns[1], StdDraw.BLUE);
    playerColors.put(playerSpawns[2], StdDraw.GREEN);
    playerColors.put(playerSpawns[3], StdDraw.PURPLE);
    // draw all player spawns on G
    for (Point p : playerSpawns) {
      StdDraw.setPenColor(playerColors.get(p));
      StdDraw.filledCircle(p.x, p.y, 0.075);
    }
    StdDraw.save("players.png");
    // int iia = (new int[1])[1];

    // now we have player spawns, perform voronoi
    // now add the playerSpawns
    ArrayList<Point> vlis = new ArrayList<>();
    for (Point p : playerSpawns) {
      vlis.add(p);
    }

    return graph;
  }

  // phase 1.2
  public static HashMap<Point, ArrayList<Point>> createGraph(ArrayList<Point> ps) {

    // delaunize points
    HashMap<Point, ArrayList<Point>> del = Delaunay.delaunize(ps);
    StdDraw.clear(StdDraw.WHITE);
    drawGraph(del);
    StdDraw.save("del.png");


    // create graph G
    HashMap<Point, ArrayList<Point>> G = new HashMap<>();
    for (Point p : ps) {
      G.put(p, new ArrayList<Point>(0));
    }
    // setup lis and connections in G
    ArrayList<Point> lis = new ArrayList<>(0);
    for (Point p : ps) {
      // if degree of p in del is <= 3
      if (del.get(p).size() <= 3) {
        // create same connections into G
        for (Point adj : del.get(p)) {
          G.get(p).add(adj);
          G.get(adj).add(p);
        }
      } else {
        lis.add(p);
      }
    }

    loop:
    while (!lis.isEmpty()) {
      Point n = lis.get((int)(Math.random() * lis.size()));
      if (G.get(n).size() >= 3) {
        lis.remove(n);
        continue loop;
      }
      // options = list of all nodes adj to n in del
      ArrayList<Point> options = new ArrayList<Point>();
      for (Point p : del.get(n)) {
        if (G.get(p).size() < 3) {
          options.add(p);
        }
      }

      // remove all points adj to n in G
      for (Point p : G.get(n)) {
        if (options.contains(p)) {
          options.remove(p);
        }
      }

      // find the maximum degree (based on G)
      int maxDeg = -1;
      for (Point p : options) {
        maxDeg = Math.max(G.get(p).size(), maxDeg);
      }
      // remove points from options that don't have same degree (based on G)
      for (int i = 0; i < options.size(); i++) {
        Point p = options.get(i);
        if (G.get(p).size() < maxDeg) {
          options.remove(p);
          i--;
        }
      }

      // remove n if options are empty
      if (options.isEmpty()) {
        lis.remove(n);
        continue loop;
      }

      // select random node among options to connect to
      Point a = options.get((int)(Math.random() * options.size()));
      // setup connections in G
      G.get(n).add(a);
      G.get(a).add(n);
    }

    StdDraw.clear();
    drawGraph(G);
    if (isBroken(G) != null) {
      Core.exit("broke graph pre first.png");
    }

    StdDraw.save("first.png");

    // post fix
    post_fix:
    for (Point p : G.keySet()) {
      // if there is a node with no connection
      if (G.get(p).size() <= 1) {
        // and it is adj in del to 2 other vertices that are connected to each other
        // sever their connection to each other and connect to this vertex

        // get list of points adj to p in del
        ArrayList<Point> adj = new ArrayList<>();
        for (Point a : del.get(p)) {
          adj.add(a);
        }
        // check if any of them are connected to each other in G
        for (int i = 0; i < adj.size(); i++) {
          Point a = adj.get(i);
          for (int j = i+1; j < adj.size(); j++) {
            Point b = adj.get(j);

            if (G.get(a).contains(b) && G.get(b).contains(a)) {
              // sever their connection and connect to this
              G.get(a).remove(b);
              G.get(b).remove(a);
              // connect them to p
              G.get(p).add(a);
              G.get(p).add(b);
              G.get(a).add(p);
              G.get(b).add(p);

              continue post_fix;
            }
          }
        }
      }
    }

    StdDraw.clear(StdDraw.WHITE);
    drawGraph(G);
    if (isBroken(G) != null) {
      Core.exit("broke graph pre post_fix.png");
    }
    StdDraw.save("post_fix.png");

    return G;
  }

  // phase 1.1
  private static ArrayList<Point> generatePoints(int w, int h, int size) {
    double x0 = 0, x1 = w, y0 = 0, y1 = h;

    // create a list of random nodes in the space
    ArrayList<Point> ps = CoreGeom.randomPoints(w, h, size, 1);

    // find current box
    double left = ps.get(0).x, top = ps.get(0).y;
    double right = left, bot = top;
    for (Point p : ps) {
      left = Math.min(p.x, left);
      right = Math.max(p.x, right);
      bot = Math.min(p.y, bot);
      top = Math.max(p.y, top);
    }

    // find ideal box
    double left_d = x0 + (x1 - x0) * 0.05;
    double right_d = x1 - (x1 - x0) * 0.05;
    double bot_d = y0 + (y1 - y0) * 0.05;
    double top_d = y1 - (y1 - y0) * 0.05;

    // find multipliers
    double vertical_multiplier = (top_d - bot_d) / (top - bot);
    double horizontal_multiplier = (right_d - left_d) / (right - left);

    // find the centre of the current box
    double centre_x = left + (right - left) * 0.5;
    double centre_y = bot + (top - bot) * 0.5;

    // subtract centre from a;; points
    for (Point p : ps) {
      p.x -= centre_x;
      p.y -= centre_y;
    }

    // multiply all points by multiplier
    for (Point p : ps) {
      p.x *= horizontal_multiplier;
      p.y *= vertical_multiplier;
    }

    // find the centre of the true box
    double centre_x_d = left_d + (right_d - left_d) * 0.5;
    double centre_y_d = bot_d + (top_d - bot_d) * 0.5;

    // add true centre to all points
    for (Point p : ps) {
      p.x += centre_x_d;
      p.y += centre_y_d;
    }

    return ps;
  }

  private static void drawGraph(HashMap<Point, ArrayList<Point>> graph) {
    double radius = 0.1;
    StdDraw.setPenColor(StdDraw.BLACK);
    for (Point p : graph.keySet()) {
      StdDraw.filledCircle(p, radius);
      for (Point adj : graph.get(p)) {
        StdDraw.line(p, adj);
      }
    }
  }

  private static void drawGraph(HashMap<Point, ArrayList<Point>> graph, Color color) {
    double radius = 0.1;
    StdDraw.setPenColor(color);
    for (Point p : graph.keySet()) {
      StdDraw.filledCircle(p, radius);
      for (Point adj : graph.get(p)) {
        StdDraw.line(p, adj);
      }
    }
  }

  private static void drawEdges(HashMap<Point, ArrayList<Point>> graph, Color color) {
    double radius = 0.1;
    StdDraw.setPenColor(color);
    for (Point p : graph.keySet()) {
      for (Point adj : graph.get(p)) {
        StdDraw.line(p, adj);
      }
    }
  }


  private static Point[] isBroken(HashMap<Point, ArrayList<Point>> graph) {
    for (Point p : graph.keySet()) {
      for (Point adj : graph.get(p)) {
        if (!graph.get(adj).contains(p)) {
          return new Point[]{p, adj};
        }
      }
    }
    return null;
  }


}
