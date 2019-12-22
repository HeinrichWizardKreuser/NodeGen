/*

*/


import java.util.*;
import java.util.function.Function;
import java.awt.Color;

public class NodeGen2 {

  // stores the location of each player and each pf their colors
  public static HashMap<Point, Color> playerColors = null;
  public static ArrayList<Color> colors = null;
  public static HashMap<Color, ArrayList<Double>> progression;

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
    // draw the points
    int pixls = 80;
    StdDraw.setCanvasSize(w * pixls, h * pixls);
    StdDraw.setXscale(0, w);
    StdDraw.setYscale(0, h);
    CoreDraw.spawnX += w * pixls;


    // phase 1.1
    ArrayList<Point> ps = generatePoints(w, h, size);
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

    CoreDraw stats = new CoreDraw(56, 20, 20, "Progression");
    // keep track of progression
    // mak colors to arraylist of doubles
    progression = new HashMap<>();
    for (Color c : colors) {
      progression.put(c, new ArrayList<>());
    }
    StdDraw.enableDoubleBuffering();
    // while lengths are unnacceptable, perfrom phase 2.2 again
    while (!acceptable(subgraphs)) {

      drawProgression(stats);
      StdDraw.show();


      // draw graph in background
      StdDraw.clear(StdDraw.WHITE);
      StdDraw.drawGraph(graph, 0.2, StdDraw.LIGHT_GRAY);
      // Draw the edges of each subgraph
      for (Point spawn : subgraphs.keySet()) {
        HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
        Color color = playerColors.get(spawn);
        StdDraw.drawEdges(subgraph, color);
        // draw the nodes of the current spawns
        StdDraw.filledCircle(spawn, 0.1, color);
      }
      // Core.freeze("subgraphs");



      HashMap<Point, Point[]> absoluteCentres = new HashMap<>();
      // check if stats are too similar so far
      if (slowProgression()) {
        Core.freeze("progression deemed slow, commencing extreme balance");
        // then apply genetic algorithm
        // make the absolute centres somethings else
        // find the losing player
        int N = progression.get(StdDraw.RED).size();
        double min = progression.get(StdDraw.RED).get(N - 1);
        Color minColor = StdDraw.RED;
        for (Color c : progression.keySet()) {
          double val = progression.get(c).get(N - 1);
          if (val < min) {
            min = val;
            minColor = c;
          }
        }

        // now we have determined the color of the player who is losing
        // get the spawn point of that player
        Point losingPlayer = null;
        for (Point spawn : playerColors.keySet()) {
          if (playerColors.get(spawn) == minColor) {
            losingPlayer = spawn;
            break;
          }
        }
        for (Point spawn : subgraphs.keySet()) {
          if (spawn == losingPlayer) {
            continue;
          }
          HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
          // get distance to losingplayer
          Point[] furthestNewSpawn = new Point[]{subgraph.get(spawn).get(0), spawn, subgraph.get(spawn).get(1)};
          System.out.println("1");
          double minLen = CoreGeom.lineLength(CoreGeom.aStar(graph, spawn, losingPlayer));
          // get each edge in the graph
          for (Point[] edge : CoreGeom.edgeList(subgraph)) {
            System.out.println("2");
            // check what the pathlength would be if it the spawn was on this edge instead
            // get edge points
            Point e0 = edge[0], e1 = edge[1];
            // create point in centre
            Point centre = new Point(e0, e1);
            // add centre into graph
            CoreGeom.addBetween(e0, e1, graph, centre);
            // check the shortest distance to the losing player if this was in the real graph
            ArrayList<Point> path = CoreGeom.aStar(graph, centre, losingPlayer);
            // get length
            double len = CoreGeom.lineLength(path);
            // save if the length is the lowest
            if (len < minLen) {
              minLen = len;
              furthestNewSpawn = new Point[]{e0, centre, e1};
            }
            // remove from graph
            CoreGeom.remove(centre, graph);
            System.out.println("3");
          }
          System.out.println("4");
          // so now we have the point where the spawn ought to be.
          // so we'll move it there later, save it for now
          absoluteCentres.put(spawn, furthestNewSpawn);
        }




      } else {
        // get absolute centres
        for (Point spawn : subgraphs.keySet()) {
          HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
          Point[] absoluteCentre = CoreGeom.absoluteCentre(subgraph);
          absoluteCentres.put(spawn, absoluteCentre);
        }
      }





      // Draw the absolute centres for each color
      for (Point spawn : absoluteCentres.keySet()) {
        Point[] e = absoluteCentres.get(spawn);
        Point e0 = e[0], e1 = e[2];
        Point absoluteCentre = e[1];
        Color color = playerColors.get(spawn);
        StdDraw.line(e0, e1, StdDraw.YELLOW);
        StdDraw.filledSquare(absoluteCentre, 0.1, color);
      }
      // Core.freeze("absolute centres");




      // add them to the graph
      for (Point spawn : absoluteCentres.keySet()) {
        Point[] e = absoluteCentres.get(spawn);
        Point e0 = e[0], e1 = e[2];
        Point absoluteCentre = e[1];
        // add it to the graph
        boolean containsE0 = graph.containsKey(e0), containsE1 = graph.containsKey(e1);
        if (containsE0 && containsE1) { // best case scenario, then both e0 and e1 are in the graph
          // add absoluteCentre between them e0 and e1
          CoreGeom.sever(e0, e1, graph);
          graph.put(absoluteCentre, new ArrayList<Point>(0));
          CoreGeom.connect(absoluteCentre, e0, graph);
          CoreGeom.connect(absoluteCentre, e1, graph);
        } else if (!containsE0 && !containsE1) {
          // worst case scenario. this is impossible. This shows that e0 and e1 are not real. Odd
          Core.exit("e0 and e1 does not exist within the graph!");
        } else {
          Point eLeaf = null, eReal = null;
          if (containsE0) {
            Core.log("e1 is a leaf!");
            eLeaf = e1;
            eReal = e0;
          } else if (containsE1) {
            Core.log("e0 is a leaf!");
            eLeaf = e0;
            eReal = e1;
          }

          // finds the leaf with identical coordinates but not same and then
          findAdjToOtherLeaf:
          if (true) {
            for (Point spawn2 : subgraphs.keySet()) {
              if (spawn2 == spawn) {
                continue;
              }
              HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn2);
              for (Point p : subgraph.keySet()) {
                // check if leaf node
                if (subgraph.get(p).size() == 1) {
                  // check if .equals(e0) but != e0
                  if (p.equals(eLeaf)) {
                    Point adjToLeaf = subgraph.get(p).get(0);

                    CoreGeom.sever(eReal, adjToLeaf, graph);
                    graph.put(absoluteCentre, new ArrayList<Point>(0));
                    CoreGeom.connect(absoluteCentre, adjToLeaf, graph);
                    CoreGeom.connect(absoluteCentre, eReal, graph);

                    break findAdjToOtherLeaf;
                  }
                }
              }
            }
            Core.log("errrrrrr");
            while (true) {}
          }

        }
      }




      // remove the spawns from graph
      for (Point spawn : subgraphs.keySet()) {
        // get points adjacent to spawn
        ArrayList<Point> adjToSpawn = graph.get(spawn);
        if (adjToSpawn.size() != 2) {
          Core.exit("adjToSpawn.size() = " + adjToSpawn.size());
        }
        // remove spawn
        graph.remove(spawn);
        // make those previously adjacent to spawn no longer adjacent
        Point e0 = adjToSpawn.get(0), e1 = adjToSpawn.get(1);
        graph.get(e0).remove(spawn);
        graph.get(e1).remove(spawn);
        // connect them to each other instead (as they originally were)
        CoreGeom.connect(e0, e1, graph);
      }



      // confirm that the graph no longer contains any trace of the spawn
      for (Point p : graph.keySet()) {
        if (vlis.contains(p)) {
          Core.exit("graph still contains a spawn as key");
        }
        for (Point adj : graph.get(p)) {
          if (vlis.contains(adj)) {
            Core.exit("there still exists a node adjacent to spawn");
          }
        }
      }



      // spawn is removed, we must set the spawns to be the absolute centres.
      for (Point spawn : absoluteCentres.keySet()) {
        Point absoluteCentre = absoluteCentres.get(spawn)[1];
        // make the absoluteCentres the new spawns
        // update the colors to reflect that the spawns no longer exist
        Color color = playerColors.remove(spawn);
        playerColors.put(absoluteCentre, color);
        // update vlis to show the new players
        vlis.remove(spawn);
        vlis.add(absoluteCentre);
      }



      // reset subgraphs
      subgraphs = CoreGeom.voronoiGraph(vlis, graph);

      StdDraw.show();

    }

    // now draw the final graph
    StdDraw.clear(StdDraw.WHITE);
    StdDraw.drawGraph(graph, 0.2, StdDraw.LIGHT_GRAY);
    // Draw the edges of each subgraph
    for (Point spawn : subgraphs.keySet()) {
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
      Color color = playerColors.get(spawn);
      StdDraw.drawEdges(subgraph, color);
      // draw the nodes of the current spawns
      StdDraw.filledCircle(spawn, 0.1, color);
    }

    drawProgression(stats);

  }

  // check if the latest stats show that there hasn't been much change
  public static boolean slowProgression() {
    if (progression.get(StdDraw.RED).size() < 2) {
      return false;
    }
    int lastI = progression.get(StdDraw.RED).size() - 1;
    int secondLastI = lastI - 1;
    for (Color c : progression.keySet()) {
      ArrayList<Double> stats = progression.get(c);
      // check if last few values are the same
      double lastVal = stats.get(lastI), secondLastVal = stats.get(stats.size() - 2);
      if (!Core.epsilon(stats.get(lastI), stats.get(secondLastI))) {
        return false;
      }
    }
    return true;
  }

  public static boolean acceptable(HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs) {
    // check length of each subgraph
    HashMap<Point, Double> graphLengths = new HashMap<>();
    for (Point spawn : subgraphs.keySet()) {
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
      // get path length
      double length = CoreGeom.graphLength(subgraph);
      graphLengths.put(spawn, length);
    }
    // get max and min
    Point minSpawn = Core.getMinKey(graphLengths), maxSpawn = Core.getMaxKey(graphLengths);
    double min = graphLengths.get(minSpawn), max = graphLengths.get(maxSpawn);
    // check if ratio is okay
    double ratio = (max - min) / max;
    if (ratio > 0.20) {
      Core.log("unnacceptable");
      System.out.println("lengths: ");

      inOrderOfColorLoop:
      for (Color color : colors) {
        // find player that belongs to color
        for (Point spawn : playerColors.keySet()) {
          if (color == playerColors.get(spawn)) {
            double length = graphLengths.get(spawn);
            System.out.println("\t"+StdDraw.colorName(playerColors.get(spawn))+" = "+length);
            // add to progression
            progression.get(color).add(length);
            continue inOrderOfColorLoop;
          }
        }
      }

      System.out.println();
      return false;
    }

    Core.log("ACCEPTED Graph!!");
    System.out.println("lengths: ");
    for (Point spawn : subgraphs.keySet()) {
      double length = graphLengths.get(spawn);
      System.out.println("\t"+StdDraw.colorName(playerColors.get(spawn))+" = "+length);
      // add to progression
      progression.get(playerColors.get(spawn)).add(length);
    }
    System.out.println();
    return true;
  }

  public static void drawProgression(CoreDraw stats) {
    // filter down to what can fit on the screen

    // draw up the progression
    int x = progression.get(StdDraw.RED).size();
    double min = progression.get(StdDraw.RED).get(0);
    double max = min;
    if (x > 1) {
      // first get the highest one and the lowest one
      for (Color c : colors) {
        for (Double d : progression.get(c)) {
          min = Math.min(min, d);
          max = Math.max(max, d);
        }
      }
      // now shape the stats graph so that the height is there
      double scalar = stats.frameH / (max - min);
      final double mini = min;
      Function<Double, Double> scale = (d) -> {
        return (d - mini) * scalar;
      };
      // now we have a way to scale all of the amounts
      // we must now draw the entire graph
      stats.clear(StdDraw.WHITE);
      for (Color c : colors) {
        stats.setPenColor(c);
        ArrayList<Double> lis = progression.get(c);
        // draw a line from the prev height to the current
        // if x is greater than frameW, start i and a different place
        if (x > stats.frameW) {
          for (int i = x-stats.frameW; i < x; i++) {
            Point prev = new Point(i-1 - (x-stats.frameW), scale.apply(lis.get(i-1)));
            Point next = new Point(i - (x-stats.frameW), scale.apply(lis.get(i)));
            stats.line(prev, next);
          }
        } else {
          for (int i = 1; i < x; i++) {
            Point prev = new Point(i-1, scale.apply(lis.get(i-1)));
            Point next = new Point(i, scale.apply(lis.get(i)));
            stats.line(prev, next);
          }
        }
      }
    }
  }
  // phase 2.1
  public static HashMap<Point, ArrayList<Point>> generatePlayerSpawns(HashMap<Point, ArrayList<Point>> G, int players) {
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
      if (H.isEmpty()) {
        Core.log("PHASE 2.1: H is Empty...trying again...");
        return generatePlayerSpawns(G, players);
      }
      // choose random point in H
      Point rand = Core.randomKey(H);
      if (H.get(rand).isEmpty()) {
        H.remove(rand);
        continue;
      }
      // choose a random edge
      Point adj = Core.randomKey(H.get(rand));
      // add playerspawn here
      Point playerSpawn = new Point(rand, adj, 0.2+(Math.random()*0.6));
      // add to player spawns
      playerSpawns[--amnt] = playerSpawn;

      //______________this is for H, for generating points______________
      // remove them from H
      ArrayList<Point> adjToRand = new ArrayList<>();
      for (Point p : H.get(rand)) {
        adjToRand.add(p);
      }
      for (Point p : adjToRand) {
        H.get(rand).remove(p);
        H.get(p).remove(rand);
      }
      ArrayList<Point> adjToAdj = new ArrayList<>();
      for (Point p : H.get(adj)) {
        adjToAdj.add(p);
      }
      for (Point p : adjToAdj) {
        H.get(adj).remove(p);
        H.get(p).remove(adj);
      }
      H.remove(rand);
      H.remove(adj);

      //___________this is for the actual graph__________________
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

    if (!CoreGeom.isUndirected(G) || !CoreGeom.isConnected(G)) {
      Core.log("PHASE 2.1: graph creation failed pre post_fix...trying again...");
      return generatePlayerSpawns(G, players);
    }

    playerColors = new HashMap<>();
    playerColors.put(playerSpawns[0], StdDraw.RED);
    playerColors.put(playerSpawns[1], StdDraw.BLUE);
    playerColors.put(playerSpawns[2], StdDraw.GREEN);
    playerColors.put(playerSpawns[3], StdDraw.PURPLE);

    colors = new ArrayList<Color>(Arrays.asList(new Color[]{
      StdDraw.RED,
      StdDraw.BLUE,
      StdDraw.GREEN,
      StdDraw.PURPLE
    }));

    return graph;
  }

  // phase 1.2
  public static HashMap<Point, ArrayList<Point>> createGraph(ArrayList<Point> ps) {

    // delaunize points
    HashMap<Point, ArrayList<Point>> del = Delaunay.delaunize(ps);

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
      Point n = Core.randomKey(lis);
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

    if (!CoreGeom.isUndirected(G)) {
      Core.log("undirected -_- phase 1.2 pre post_fix");
      Core.exit("");
    } else if (!CoreGeom.isConnected(G)) {
      Core.log("not connected -_- phase 1.2 pre post_fix");
      Core.exit("");
    }

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

    if (!CoreGeom.isUndirected(G) || !CoreGeom.isConnected(G)) {
      Core.log("PHASE 1.2: graph creation failed after post_fix...trying again...");
      return createGraph(ps);
    }

    return G;
  }

  // phase 1.1
  public static ArrayList<Point> generatePoints(int w, int h, int size) {
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
