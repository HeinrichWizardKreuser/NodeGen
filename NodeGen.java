import java.util.*;
import java.util.function.Function;
import java.awt.Color;

public class NodeGen {

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

    // System.out.println("left = " + left);
    // System.out.println("right = " + right);
    // System.out.println("top = " + top);
    // System.out.println("bot = " + bot);

    // find ideal box
    double left_d = x0 + (x1 - x0) * 0.05;
    double right_d = x1 - (x1 - x0) * 0.05;
    double bot_d = y0 + (y1 - y0) * 0.05;
    double top_d = y1 - (y1 - y0) * 0.05;

    // System.out.println("left_d = " + left_d);
    // System.out.println("right_d = " + right_d);
    // System.out.println("top_d = " + top_d);
    // System.out.println("bot_d = " + bot_d);


    // find multipliers
    double vertical_multiplier = (top_d - bot_d) / (top - bot);
    double horizontal_multiplier = (right_d - left_d) / (right - left);
    // System.out.println("vertical_multiplier = " + vertical_multiplier);
    // System.out.println("horizontal_multiplier = " + horizontal_multiplier);


    // find the centre of the current box
    double centre_x = left + (right - left) * 0.5;
    double centre_y = bot + (top - bot) * 0.5;
    // System.out.println("centre_x = " + centre_x);
    // System.out.println("centre_y = " + centre_y);

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
    // System.out.println("centre_x_d = " + centre_x_d);
    // System.out.println("centre_y_d = " + centre_y_d);

    // add true centre to all points
    for (Point p : ps) {
      p.x += centre_x_d;
      p.y += centre_y_d;
    }




    // draw the points
    int pixls = 80;
    StdDraw.setCanvasSize(w * pixls, h * pixls);
    StdDraw.setXscale(0, w);
    StdDraw.setYscale(0, h);

    StdDraw.clear(StdDraw.WHITE);
    drawPoints(ps);

    // delaunize points
    HashMap<Point, ArrayList<Point>> del = Delaunay.delaunize(ps);
    StdDraw.clear(StdDraw.WHITE);
    drawPoints(ps);
    drawEdges(del);
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

    // Core.freeze("updating:");
    StdDraw.clear(StdDraw.WHITE);
    drawPoints(ps);
    drawEdges(G);

    // Core.freeze("about to begin:");

    loop:
    while (!lis.isEmpty()) {
      Point n = lis.get((int)(Math.random() * lis.size()));
      if (G.get(n).size() >= 3) {
        lis.remove(n);
        continue loop;
      }
      StdDraw.setPenColor(StdDraw.GREEN);
      StdDraw.filledCircle(n.x, n.y, 0.1);
      // Core.freeze();
      // options = list of all nodes adj to n in del
      ArrayList<Point> options = new ArrayList<Point>();
      for (Point p : del.get(n)) {
        if (G.get(p).size() < 3) {
          options.add(p);
        }
      }

      drawPoints(ps);
      drawPoints(options, StdDraw.RED);
      // Core.freeze("all adj in del with deg < 3 in G");

      // remove all points adj to n in G
      for (Point p : G.get(n)) {
        if (options.contains(p)) {
          options.remove(p);
        }
      }
      drawPoints(ps);
      drawPoints(options, StdDraw.RED);
      // Core.freeze("all not adj in G");

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

      for (Point p : options) {
        StdDraw.setPenColor(StdDraw.RED);
        StdDraw.filledCircle(p.x, p.y, 0.1);
      }
      // Core.freeze();

      // select random node among options to connect to
      Point a = options.get((int)(Math.random() * options.size()));
      // setup connections in G
      G.get(n).add(a);
      G.get(a).add(n);

      // draw updates
      StdDraw.clear(StdDraw.WHITE);
      drawPoints(ps);
      drawEdges(G);
      // Core.freeze();
    }

    StdDraw.save("first.png");

    //Core.freeze("applying post fix");

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

//__________________________________________________________________________________________________
// =================================================================================================
//__________________________________________________________________________________________________




    StdDraw.clear(StdDraw.WHITE);
    drawPoints(ps);
    drawEdges(G);

    StdDraw.save("post_fix.png");


    // now onto phase 2:
    // we must perform voronoi on the map
    // first choose random places to place players
    int players = 4;
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

    // // fix anything wrong in G
    // // just draw the graph for me real quick
    // for (Point p : G.keySet()) {
    //   for (Point adj : G.get(p)) {
    //     if (!G.get(adj).contains(p)) {
    //       StdDraw.filledCircle(p, 0.05, StdDraw.CYAN);
    //       StdDraw.filledCircle(adj, 0.05, StdDraw.YELLOW);
    //       Core.freeze();
    //     }
    //   }
    // }
    // Core.freeze();





    playerColors = new HashMap<>();
    playerColors.put(playerSpawns[0], StdDraw.RED);
    playerColors.put(playerSpawns[1], StdDraw.BLUE);
    playerColors.put(playerSpawns[2], StdDraw.GREEN);
    playerColors.put(playerSpawns[3], StdDraw.YELLOW);
    // draw all player spawns on G
    for (Point p : playerSpawns) {
      StdDraw.setPenColor(playerColors.get(p));
      StdDraw.filledCircle(p.x, p.y, 0.075);
    }


    // now we have player spawns, perform voronoi
    // now add the playerSpawns
    ArrayList<Point> vlis = new ArrayList<>();
    for (Point p : playerSpawns) {
      vlis.add(p);
    }
    voronoiGraph(vlis, graph);
  }

  public static HashMap<Point, Color> playerColors = null;

  public static void voronoiGraph(ArrayList<Point> vlis, HashMap<Point, ArrayList<Point>> graph) {
    // fix anything wrong in G
    // just draw the graph for me real quick
    Point[] e = isBroken(graph);
    if (e != null) {
      StdDraw.filledCircle(e[0], 0.05, StdDraw.CYAN);
      StdDraw.filledCircle(e[1], 0.05, StdDraw.YELLOW);
      int iia = (new int[1])[1];
    }

    // voro = map of distance between playerspawns and their distnaces to other points in the graph
    HashMap<Point, HashMap<Point, Double>> voro = new HashMap<Point, HashMap<Point, Double>>();
    double d = 0;
    for (Point a : vlis) {
      Core.log("____________");
      // a is the starting node
      // lis is the current nodes to look at
      ArrayList<Point> lis = new ArrayList<>();
      lis.add(a);
      // explored is all of the nodes that we have looked at
      HashMap<Point, Double> explored = new HashMap<>();
      explored.put(a, 0d);
      whileLisNotEmpty:
      while (!lis.isEmpty()) {

        //Point p = lis.remove(0);
        // remove smallest from list
        Point p = lis.get(0);
        for (int i = 1; i < lis.size(); i++) {
          if (explored.get(lis.get(i)) < explored.get(p)) {
            p = lis.get(i);
          }
        }
        lis.remove(p);


        StdDraw.clear();
        drawGraph(graph);
        StdDraw.filledCircle(p, 0.2, StdDraw.GREEN);
        StdDraw.text(p, "p", StdDraw.BLACK);
        for (Point adj : graph.get(p)) {
          StdDraw.filledCircle(adj, 0.2, StdDraw.BLUE);
        }
        for (Point exp : explored.keySet()) {
          StdDraw.filledCircle(exp, 0.15, StdDraw.YELLOW);
          StdDraw.text(exp, Math.round(explored.get(exp))+"", StdDraw.BLACK);
        }
        // draw all next on list
        for (Point next : lis) {
          StdDraw.filledCircle(next, 0.05, StdDraw.PURPLE);
        }
        // draw on adjacent
        Core.freeze("looking at p");
        // show p
        // show all that has been looked at and their values

        //StdDraw.filledCircle(p, 0.05, StdDraw.CYAN);
        for (Point adj : graph.get(p)) {
          //StdDraw.filledCircle(adj, 0.05, StdDraw.YELLOW);
          d = explored.get(p) + p.dist(adj);
          if (explored.containsKey(adj) && d >= explored.get(adj)) {
            continue whileLisNotEmpty;
          }
          //Core.freeze();
          explored.put(adj, d);
          // now we need to add adj to lis, but via binary insertion
          //binaryInsert(adj, d, lis, explored);
          lis.add(p);
          //StdDraw.text(adj, d+"", StdDraw.GREEN);

          /*
          // or alternatively...:
          int n = lis.size();
          int i = 0;
          while (i < n) {
            if (d < explored.get(lis.get(i))) {
              i++;
            } else {
              break;
            }
          }
          lis.add(i, adj);
          */

        }
      }
      voro.put(a, explored);
      //int iia = (new int[1])[1];

    }


    for (Point a : voro.keySet()) {
      StdDraw.setPenColor(playerColors.get(a));
      for (Point p : voro.get(a).keySet()) {
        StdDraw.filledCircle(p, 0.05);
      }
      StdDraw.text(a, "a", StdDraw.PURPLE);

      Core.freeze("a");
    }


    // closestTo = all points in G mapped to the point in voro they are closest to
    HashMap<Point, Point> closestTo = new HashMap<Point, Point>();
    ArrayList<Point[]> edges = new ArrayList<>();
    Point minDistP = null;
    for (Point p : graph.keySet()) {
      d = Double.POSITIVE_INFINITY;
      for (Point a : voro.keySet()) {

        if (!voro.get(a).containsKey(p)) {
          StdDraw.filledCircle(a, 0.1, StdDraw.CYAN);
          StdDraw.text(a, "a", StdDraw.BLACK);
          StdDraw.filledCircle(p, 0.1, StdDraw.RED);
          StdDraw.text(p, "p", StdDraw.WHITE);
          Core.freeze();
        } else {

        }

        if (voro.get(a).get(p) < d) {
          d = voro.get(a).get(p);
          minDistP = a;
        }
      }
      closestTo.put(p, minDistP);
      // add edge p - adj
      for (Point adj : graph.get(p)) {
        if (closestTo.containsKey(adj)) {
          if (closestTo.get(adj) != p) {
            edges.add(new Point[]{p, adj});
          }
        }
      }
    }


    // now draw yo
    // draw green for all that are closest


  }

  private static void binaryInsert(Point insertP, double insertCost,
    ArrayList<Point> lis, HashMap<Point, Double> explored) {

    int n = lis.size();
    if (n == 0) {
      lis.add(insertP);
      return;
    }

    Function<Integer, Double> cost = (index) -> {
      return explored.get(lis.get(index));
    };

    int lo = 0;
    int hi = n - 1;
    double oo = cost.apply(hi) - cost.apply(lo);
    if (oo == 0) {
      // assume min
      oo = 1;
    }

    // new element is before first
    if (oo * (insertCost - cost.apply(lo)) <= 0) {
      lis.add(0, insertP);
    }
    // new element is after last
    else if (oo * (insertCost - cost.apply(hi)) >= 0) {
      lis.add(insertP);
    }
    // is in middle
    else {
      while (lo + 1 < hi) {
        int mid = (lo + hi) / 2;
        double ox = (insertCost - cost.apply(mid)) * oo;
        // new element is after mid
        if (ox >= 0) {
          lo = mid;
        }
        // new element is before mid
        if (ox <= 0) {
          hi = mid;
        }
      }
      lis.add(hi, insertP);
    }
  }


  private static void drawPoints(ArrayList<Point> ps) {
    double radius = 0.1;
    StdDraw.setPenColor(StdDraw.BLACK);
    for (Point p : ps) {
      StdDraw.filledCircle(p.x, p.y, radius);
    }
  }

  private static void drawPoints(ArrayList<Point> ps, Color color) {
    double radius = 0.1;
    StdDraw.setPenColor(color);
    for (Point p : ps) {
      StdDraw.filledCircle(p.x, p.y, radius);
    }
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

  private static void drawEdges(HashMap<Point, ArrayList<Point>> del) {
    StdDraw.setPenColor(StdDraw.BLACK);
    for (Point p : del.keySet()) {
      for (Point adj : del.get(p)) {
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

/*
  // create a shortest spanning tree from the given points
  // return the graph of this tree
  private static HashMap<Point, ArrayList<Point>> nodeGen(ArrayList<Point> ps) {
    int n = ps.size();
    // get distances, also get max
    int max = 0;
    double[][] dists = new double[n][n];
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        dists[i][j] = ps.get(i).dist(ps.get(j));
        max = Math.max(max, dists[i][j]);
      }
    }
    max++;

    // get the minimum


    double radius = 0.1;
    StdDraw.setPenColor(StdDraw.BLACK);
    for (Point p : ps) {
      StdDraw.filledCircle(p.x, p.y, radius);
    }

  }

  private static Point[] mini(ArrayList<Point> ps, double[][] dists, int max) {
    Point[] mini = new Point(-1, -1); // the edge with minimum length
    int n = ps.size();
    int min = max;
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        if (dists[i][j] < min) {
          min = dists[i][j];
          mini.x = ps.x;
          mini.y = ps.y;
        }
      }
    }
    return mini;
  }
*/
}
