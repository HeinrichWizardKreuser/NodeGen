import java.util.*;
import java.util.function.Function;
import java.awt.Color;

public class NodeGen2 {

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

    ArrayList<Point> ps = generatePoints(w, h, size);


    // draw the points
    int pixls = 80;
    StdDraw.setCanvasSize(w * pixls, h * pixls);
    StdDraw.setXscale(0, w);
    StdDraw.setYscale(0, h);

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


//__________________________________________________________________________________________________
// =================================================================================================
//__________________________________________________________________________________________________

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



    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs = voronoiGraph(vlis, graph);

    // now draw the subgraphs
    StdDraw.clear(StdDraw.WHITE);
    for (Point player : subgraphs.keySet()) {
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(player);
      Color color = playerColors.get(player);
      drawGraph(subgraph, color);
    }

    Core.freeze();

    // get the lengths of each graph.
    // if the lengths are not within acceptable range, move the spawn points to
      // the centroids of these subgraphs

    // with subgraphs succesfully gotten, we must get the centroid of these graphs
    for (Point player : subgraphs.keySet()) {
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(player);
      Color color = playerColors.get(player);
      Point centroid = graphCentroid(subgraph)[1];
      StdDraw.filledSquare(centroid, 0.05, color);
    }

    // now we have the centroid...what now?
    // now we need to move the spawn zone to this point.
    // to do this, we must (for each spawn point):
      // remove the spawn point by
        // severing its connections between its edge's endpoints
        // connecting the endpoints to each other (as they originally were)
      // place spawn points at the coordinates of the centroids

    // re do voronoization and get new pathlengths
    // repeat until lengths are optimal
  }


  // get the centroid from the given graph
  public static Point[] graphCentroid(HashMap<Point, ArrayList<Point>> graph) {
    // map each node to the value of the longest distance to another node in the graph
    HashMap<Point, Double> furthestCosts = new HashMap<>();
    for (Point p : graph.keySet()) {
      // get the cost matrix
      HashMap<Point, Double> cost = CoreGeom.dijkstraGraph(p, graph);
      // get the point with the highest value
      Point furthest = Core.getMaxKey(cost);
      double furthestCost = cost.get(furthest);
      // add to mapping
      furthestCosts.put(p, furthestCost);
    }
    // the node with the minimum max cost will be part of the edge containing the centroid
    Point e0 = Core.getMinKey(furthestCosts);
    // now search among adjacent nodes for the one with the smallest furthest cost
    Point e1 = Core.randomKey(graph.get(e0));
    for (Point adj : graph.get(e0)) {
      if (furthestCosts.get(adj) < furthestCosts.get(e1)) {
        e1 = adj;
      }
    }
    // now we have the edge. Get weighted ave on that line based on the costs
    double e0Cost = furthestCosts.get(e0), e1Cost = furthestCosts.get(e1);

    double d = e0Cost/(e0Cost + e1Cost);

    Point ave = new Point(e0, e1, d);

    return new Point[]{e0, ave, e1};
  }

  public static HashMap<Point, Color> playerColors = null;

  public static HashMap<Point, HashMap<Point, ArrayList<Point>>> voronoiGraph(
    ArrayList<Point> vlis,
    HashMap<Point, ArrayList<Point>> graph) {

    // voro = map of distance between playerspawns and their distnaces to other points in the graph
    HashMap<Point, HashMap<Point, Double>> voro = new HashMap<Point, HashMap<Point, Double>>();

    // let's just do a breadth first search
    for (Point a : vlis) {

      ArrayList<Point> curr = new ArrayList<>();
      curr.add(a);
      HashMap<Point, Double> cost = new HashMap<>();
      cost.put(a, 0d);
      ArrayList<Point> explored = new ArrayList<>();
      while (!curr.isEmpty()) {
        // find point with lowest score in cost
        Point p = curr.get(0);
        for (Point p1 : curr) {
          if (cost.get(p1) < cost.get(p)) {
            p = p1;
          }
        }
        curr.remove(p);
        // explore point
        explored.add(p);
        // now go to all adjacent points
        for (Point adj : graph.get(p)) {
          double newDist = cost.get(p) + p.dist(adj);
          if (explored.contains(adj) && newDist > cost.get(adj)) {
            continue;
          }
          curr.add(adj);
          explored.add(adj);
          cost.put(adj, newDist);

          // StdDraw.filledCircle(adj, 0.05, StdDraw.CYAN);
          // StdDraw.text(adj, Math.round(cost.get(adj))+"", StdDraw.RED);
          // Core.freeze();
        }
      }
      // print out cost
      // System.out.println("____________________");
      // for (Point p : cost.keySet()) {
      //   System.out.println(p.toString() + " = " + cost.get(p));
      // }

      voro.put(a, cost);
    }

    // closestTo = all points in G mapped to the point in voro they are closest to
    HashMap<Point, Point> closestTo = new HashMap<Point, Point>();
    // edgeConflicts = list of all edges that have endpoints closer to opposing players
    ArrayList<Point[]> edgeConflicts = new ArrayList<>();
    // edgeFriendly = list of all edges of which both ends belong to the same player
    ArrayList<Point[]> edgeFriendly = new ArrayList<>();
    // find the player to which each node is closest to
    for (Point p : graph.keySet()) {
      // Core.println("____");
      // find which player p is closest to
      Point minDistP = Core.randomKey(voro);
      double minDist = voro.get(minDistP).get(p);
      for (Point a : voro.keySet()) {
        HashMap<Point, Double> cost = voro.get(a);
        // Core.log("comparing: minDist " + minDist + " < " + cost.get(p));
        if (cost.get(p) < minDist) {
          minDist = cost.get(p);
          minDistP = a;
        }
      }
      // if (minDistP == null) {
      //   Core.exit("minDistP == null");
      // }
      // Core.log("winner = " + minDist);
      closestTo.put(p, minDistP);

      // check if p is adj to an opposing player
      for (Point adj : graph.get(p)) {
        // if adj has also been treated
        if (closestTo.containsKey(adj)) {
          // are they closer to different players?
          if (closestTo.get(adj) != p) {
            edgeConflicts.add(new Point[]{p, adj});
          } else {
            edgeFriendly.add(new Point[]{p, adj});
          }
        }
      }
    }

    // define the edges on which there are conflicts
    HashMap<Point[], Point> borders = new HashMap<>();
    // HashMap<Point, Point> borderAccess = new HashMap<>();
    for (Point[] e : edgeConflicts) {
      double d = e[0].dist(e[1]);
      double len = (voro.get(closestTo.get(e[0])).get(e[0]) +
                    voro.get(closestTo.get(e[1])).get(e[1]) + d)/2 -
                    voro.get(closestTo.get(e[0])).get(e[0]);

      Point border = new Point(e[0], e[1], 1 - len/d);
      borders.put(e, border);

      // borderAccess.put(e[0], border);
      // borderAccess.put(e[1], border);
    }

    StdDraw.clear(StdDraw.WHITE);
    //drawGraph(graph);
    // draw the friendly edges
    for (Point[] e : edgeFriendly) {
      StdDraw.line(e[0], e[1], playerColors.get(e[0]));
    }
    // draw the border lines
    for (Point[] e : borders.keySet()) {
      StdDraw.line(e[0], borders.get(e), playerColors.get(closestTo.get(e[0])));
      StdDraw.line(e[1], borders.get(e), playerColors.get(closestTo.get(e[1])));
    }

    // color all nodes in the color of the player they are closest to
    for (Point p : graph.keySet()) {
      Point belongsTo = closestTo.get(p);
      StdDraw.filledCircle(p, 0.05, playerColors.get(belongsTo));
    }
    // make the specific spawns a bit bigger
    for (Point p : voro.keySet()) {
      StdDraw.filledCircle(p, 0.075, playerColors.get(p));
    }

    StdDraw.save("voronoi.png");
    Core.freeze();

    // __________________BALANCING______________________
    // create subgraph of each player
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs = new HashMap<>();
    for (Point p : voro.keySet()) {
      subgraphs.put(p, new HashMap<>());
    }
    // add each point
    for (Point p : graph.keySet()) {
      // check which player this belongs to
      Point belongsTo = closestTo.get(p);
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(belongsTo);
      subgraph.put(p, new ArrayList<Point>());
      for (Point adj : graph.get(p)) {
        Point adjBelongsTo = closestTo.get(adj);
        if (belongsTo == adjBelongsTo) {
          subgraph.get(p).add(adj); // we do not do this for adj to p since that will be done on its own
        } else {
          // find border that has this node and adj in edge
          for (Point[] e : borders.keySet()) {
            if ((e[0] == p && e[1] == adj) || (e[0] == adj && e[1] == p)) {
              Point border = borders.get(e);
              // now we must add a new node - from borderAccess
              subgraph.get(p).add(border);
              // however we will do the border to p for this case
              subgraph.put(border, new ArrayList<Point>());
              subgraph.get(border).add(p);
              break;
            }
          }
        }
      }
    }


    return subgraphs;
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
}
