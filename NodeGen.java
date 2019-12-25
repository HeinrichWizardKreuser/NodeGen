/*
  TODO:
    uograde movement of spawns to Genetic Algorithm approach
*/


import java.util.*;
import java.util.function.Function;
import java.awt.Color;

public class NodeGen {

  // stores the location of each player and each pf their colors
  public static ArrayList<Color> colors = null;
  public static HashMap<Color, ArrayList<Double>> progression;
  public static double mutate = 0.8;

  public static ArrayList<Point> _vlis;
  public static HashMap<Point, ArrayList<Point>> _graph;
  public static HashMap<Point, HashMap<Point, ArrayList<Point>>> _subgraphs;
  public static HashMap<Point, Color> _playerColors;
  public static HashMap<Point, Point[]> _spawnEdges;
  public static CoreDraw _gui;

  private static class Graph {

    public ArrayList<Point> vlis;
    public HashMap<Point, ArrayList<Point>> graph;
    public HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs;
    public HashMap<Point, Color> playerColors;
    // spawn -> e0-spawn-e1 (e0 & e1 definitely in graph, not leaves of subgraph)
    public HashMap<Point, Point[]> spawnEdges;

    public Graph(HashMap<Point, ArrayList<Point>> G, int players) {
      // make copy of graph
      this.graph = CoreGeom.copy(G);
      HashMap<Point, ArrayList<Point>> H = CoreGeom.copy(G);
      // now create random playerSpawns
      int amnt = players;
      Point[] playerSpawns = new Point[amnt];
      this.spawnEdges = new HashMap<Point, Point[]>();
      while (amnt > 0) {
        if (H.isEmpty()) {
          Core.log("PHASE 2.1: H is Empty...trying again...");
          Graph g = new Graph(G, players);
          this.vlis = g.vlis;
          this.graph = g.graph;
          this.subgraphs = g.subgraphs;
          this.playerColors = g.playerColors;
          this.spawnEdges = g.spawnEdges;
          return;
        }
        // choose random point in H
        Point rand = Core.randomKey(H);
        // if adjacent to nothing, remove from options and try again
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

        // add playerSpawn tp the graph
        CoreGeom.addBetween(rand, adj, graph, playerSpawn);

        this.spawnEdges.put(playerSpawn, new Point[]{rand, playerSpawn, adj});
      }

      if (!CoreGeom.isUndirected(G) || !CoreGeom.isConnected(G)) {
        Core.log("PHASE 2.1: graph creation failed pre post_fix...trying again...");
        Graph g = new Graph(G, players);
        this.vlis = g.vlis;
        this.graph = g.graph;
        this.subgraphs = g.subgraphs;
        this.playerColors = g.playerColors;
        this.spawnEdges = g.spawnEdges;
        return;
      }
      // create vlis
      this.vlis = new ArrayList<Point>(Arrays.asList(playerSpawns));
      // map to each color
      this.playerColors = new HashMap<Point, Color>();
      this.playerColors.put(playerSpawns[0], CoreDraw.RED);
      this.playerColors.put(playerSpawns[1], CoreDraw.BLUE);
      this.playerColors.put(playerSpawns[2], CoreDraw.GREEN);
      this.playerColors.put(playerSpawns[3], CoreDraw.PURPLE);
      // create the subgraph
      this.subgraphs = CoreGeom.voronoiGraph(this.vlis, this.graph);
    }

    // create a new graph that is a mutation of the given
    public Graph(Graph g) {
      this.vlis = Core.copy(g.vlis);
      this.graph = CoreGeom.copy(g.graph);
      this.subgraphs = new HashMap<Point, HashMap<Point, ArrayList<Point>>>();
      for (Point spawn : g.subgraphs.keySet()) {
        this.subgraphs.put(spawn, CoreGeom.copy(g.subgraphs.get(spawn)));
      }
      this.playerColors = Core.copy(g.playerColors);
      this.spawnEdges = Core.copy(g.spawnEdges);
      // choose one of the edges and then perform one of the stuff
      Point spawn = Core.randomKey(g.vlis);
      double rand = Math.random();
      // check what to do
      if (rand < 0.33) {
        shuffle(spawn);
      } else if (rand < 0.66) {
        move(spawn, -1);
      } else if (rand < 1) {
        move(spawn, 1);
      }
    }

    // moves the spawn point somewhere else
    public void shuffle(Point spawn) {
      // get the subgraph
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
      // search for new position in subgraph tree
      ArrayList<Point[]> edgeList = CoreGeom.edgeList(subgraph);
      Point[] newEdge = null;
      loop:
      while (true) {
        if (edgeList.isEmpty()) {
          Core.exit("edgeList.isEmpty()");
        }
        newEdge = Core.randomKey(edgeList);
        // if newEdge contains a leaf node, remove and try again
        if (!graph.containsKey(newEdge[0]) || !graph.containsKey(newEdge[1])) {
          // then one of them is a leaf
          edgeList.remove(newEdge);
          continue loop;
        }
        break loop;
      }
      // create new spawn between nodes
      Point newSpawn = new Point(newEdge[0], newEdge[1], 0.2 + Math.random()*0.6);
      // insert new
      CoreGeom.addBetween(newEdge[0], newEdge[1], graph, newSpawn);
      // remove original
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
      // reconnect adjToSpawn
      CoreGeom.connect(e0, e1, graph);

      // set spawn to newSpawn
      // update the colors to reflect that the spawns no longer exist
      Color color = playerColors.remove(spawn);
      playerColors.put(newSpawn, color);
      // update vlis to show the new players
      vlis.remove(spawn);
      vlis.add(newSpawn);
      // reset subgraphs
      subgraphs = CoreGeom.voronoiGraph(vlis, graph);

      // replace spawnEdge
      spawnEdges.remove(spawn);
      // add new spawn between whichever nodes it is in the graph
      adjToSpawn = graph.get(newSpawn);
      if (adjToSpawn.size() != 2) {
        Core.exit("_adjToSpawn.size() = " + adjToSpawn.size());
      }
      spawnEdges.put(newSpawn, new Point[]{adjToSpawn.get(0), newSpawn, adjToSpawn.get(1)});
    }

    // moves the spawn point somewhere else
    public void lloyd(Point spawn) {
      // get the subgraph
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
      // search for new position in subgraph tree
      Point[] e = CoreGeom.absoluteCentre(subgraph);
      // add absolute centre to the graph
      Point e0 = e[0], e1 = e[2];
      Point absoluteCentre = e[1];
      // add it to the graph
      boolean containsE0 = graph.containsKey(e0), containsE1 = graph.containsKey(e1);
      if (containsE0 && containsE1) { // best case scenario, then both e0 and e1 are in the graph
        // add absoluteCentre between them e0 and e1
        CoreGeom.addBetween(e0, e1, graph, absoluteCentre);
        // add to spawn edges
        spawnEdges.put(absoluteCentre, e);
      } else if (!containsE0 && !containsE1) {
        // worst case scenario. this is impossible. This shows that e0 and e1 are not real. Odd
        Core.exit("e0 and e1 does not exist within the graph!");
      } else {
        Point eLeaf = null, eReal = null;
        if (containsE0) {
          eLeaf = e1;
          eReal = e0;
        } else if (containsE1) {
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
            subgraph = subgraphs.get(spawn2);
            for (Point p : subgraph.keySet()) {
              // check if leaf node
              if (subgraph.get(p).size() == 1) {
                // check if .equals(e0) but != e0
                if (p.equals(eLeaf)) {
                  Point adjToLeaf = subgraph.get(p).get(0);
                  CoreGeom.addBetween(eReal, adjToLeaf, graph, absoluteCentre);
                  // add to spawn edges
                  spawnEdges.put(absoluteCentre, new Point[]{eReal, absoluteCentre, adjToLeaf});
                  break findAdjToOtherLeaf;
                }
              }
            }
          }
          Core.log("errrrrrr");
          while (true) {}
        }
      }

      // remove original from the graph
      // get points adjacent to spawn
      ArrayList<Point> adjToSpawn = graph.get(spawn);
      if (adjToSpawn.size() != 2) {
        Core.exit("adjToSpawn.size() = " + adjToSpawn.size());
      }
      // remove spawn
      graph.remove(spawn);
      // make those previously adjacent to spawn no longer adjacent
      e0 = adjToSpawn.get(0);
      e1 = adjToSpawn.get(1);
      graph.get(e0).remove(spawn);
      graph.get(e1).remove(spawn);
      // reconnect adjToSpawn
      CoreGeom.connect(e0, e1, graph);

      // set spawn to absolute centre
      // make the absoluteCentres the new spawns
      // update the colors to reflect that the spawns no longer exist
      Color color = playerColors.remove(spawn);
      playerColors.put(absoluteCentre, color);
      // update vlis to show the new players
      vlis.remove(spawn);
      vlis.add(absoluteCentre);
      // reset subgraphs
      subgraphs = CoreGeom.voronoiGraph(vlis, graph);
    }

    // moves the spawn point sligthly left or slightly right
    public void move(Point spawn, double direction) {
      // // get the subgraph of the spawn
      // HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
      // get the edge that it is on
      if (!spawnEdges.containsKey(spawn)) {
        Core.log("spawnEdges.size() = " + spawnEdges.size());
      }
      Point[] edge = spawnEdges.get(spawn);
      // move the spawn in the direction indicated by direction
      Point eO = edge[direction < 0 ? 0 : 2];
      // create new spawn
      Point newSpawn = new Point(eO, edge[1], Math.random());
      // confirm that eO does exist
      // check eO exists and if it does not
      if (!graph.containsKey(eO)) {
        Core.exit("!graph.contains(eO)");
      }
      // insert new spawn
      CoreGeom.addBetween(eO, edge[1], graph, newSpawn);
      // replace spawnEdge
      spawnEdges.remove(spawn);
      spawnEdges.put(newSpawn, new Point[]{eO, newSpawn, edge[1]});


      // remove original from the graph
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
      // reconnect adjToSpawn
      CoreGeom.connect(e0, e1, graph);

      // set spawn to newSpawn
      // update the colors to reflect that the spawns no longer exist
      Color color = playerColors.remove(spawn);
      playerColors.put(newSpawn, color);
      // update vlis to show the new players
      vlis.remove(spawn);
      vlis.add(newSpawn);
      // reset subgraphs
      subgraphs = CoreGeom.voronoiGraph(vlis, graph);

    }

    // calculates the cost of the graph (the ratio)
    public double cost() {
      // check length of each subgraph
      HashMap<Point, Double> graphLengths = new HashMap<>();
      for (Point spawn : this.subgraphs.keySet()) {
        // get path length
        double length = CoreGeom.graphLength(this.subgraphs.get(spawn));
        graphLengths.put(spawn, length);
      }
      // get max and min
      Point minSpawn = Core.getMinKey(graphLengths), maxSpawn = Core.getMaxKey(graphLengths);
      double min = graphLengths.get(minSpawn), max = graphLengths.get(maxSpawn);
      // check if ratio is okay
      return (max - min) / max;
    }

    // draws the subgraphs on the given gui
    public void draw(CoreDraw gui) {
      // clear background
      gui.clear(CoreDraw.WHITE);
      // draw the graph
      gui.setPenColor(CoreDraw.LIGHT_GRAY);
      for (Point p : this.graph.keySet()) {
        gui.filledCircle(p, 0.2);
        // for (Point adj : this.graph.get(p)) {
        //   gui.line(p, adj);
        // }
      }
      // Draw the edges of each subgraph
      for (Point spawn : this.subgraphs.keySet()) {
        HashMap<Point, ArrayList<Point>> subgraph = this.subgraphs.get(spawn);
        Color color = this.playerColors.get(spawn);
        gui.drawEdges(subgraph, color);
        // draw the nodes of the current spawns
        gui.filledCircle(spawn, 0.1, color);
      }
    }

    // performs loydd for a given number of iterations
    public void lloyd(int iterations) {

      HashMap<Color, ArrayList<Double>> progression = new HashMap<>();
      for (Color c : colors) {
        progression.put(c, new ArrayList<Double>());
      }

      loop:
      for (int i = 0; i < iterations; i++) {
        // check if stats are too similar so far
        HashMap<Point, Point[]> absoluteCentres = new HashMap<>();

        // get absolute centres
        for (Point spawn : subgraphs.keySet()) {
          HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);
          Point[] absoluteCentre = CoreGeom.absoluteCentre(subgraph);
          absoluteCentres.put(spawn, absoluteCentre);
        }

        // add them to the graph
        for (Point spawn : absoluteCentres.keySet()) {
          Point[] e = absoluteCentres.get(spawn);
          Point e0 = e[0], e1 = e[2];
          Point absoluteCentre = e[1];
          // add it to the graph
          boolean containsE0 = graph.containsKey(e0), containsE1 = graph.containsKey(e1);
          if (containsE0 && containsE1) { // best case scenario, then both e0 and e1 are in the graph
            // add absoluteCentre between them e0 and e1
            CoreGeom.addBetween(e0, e1, graph, absoluteCentre);
          } else if (!containsE0 && !containsE1) {
            // worst case scenario. this is impossible. This shows that e0 and e1 are not real. Odd
            Core.exit("e0 and e1 does not exist within the graph!");
          } else {
            Point eLeaf = null, eReal = null;
            if (containsE0) {
              // Core.log("e1 is a leaf!");
              eLeaf = e1;
              eReal = e0;
            } else if (containsE1) {
              // Core.log("e0 is a leaf!");
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
                      CoreGeom.addBetween(eReal, adjToLeaf, graph, absoluteCentre);
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

        if (!CoreGeom.isUndirected(graph)) {
          Core.exit("graph is directed");
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
          // recconnect adjToSpawn
          CoreGeom.connect(e0, e1, graph);
        }

        // confirm that the graph no longer contains any trace of the spawn
        checkIfBroken:
        if (true) {
          fixIt:
          if (true) {
            for (Point p : graph.keySet()) {
              if (vlis.contains(p)) {
                Core.exit("graph still contains a spawn as key");
              }
              for (Point adj : graph.get(p)) {
                if (vlis.contains(adj)) {
                  Core.log("there still exists a node adjacent to spawn");
                  // check which it is
                  Core.log(StdDraw.colorName(playerColors.get(adj)));
                  // break fixIt;
                }
              }
            }
            break checkIfBroken;
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
      }
    }

    // check if the latest stats show that there hasn't been much change
    public boolean slowProgression(HashMap<Color, ArrayList<Double>> progression) {
      if (progression.get(CoreDraw.RED).size() < 2) {
        return false;
      }
      int lastI = progression.get(CoreDraw.RED).size() - 1;
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
  }

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
    if (_gui == null) {
      _gui = new CoreDraw(w, h, pixls, "NodeGen");
    }

    colors = new ArrayList<Color>(Arrays.asList(new Color[]{
      CoreDraw.RED,
      CoreDraw.BLUE,
      CoreDraw.GREEN,
      CoreDraw.PURPLE
    }));

    // phase 1.1
    ArrayList<Point> ps = generatePoints(w, h, size);
    // phase 1.2
    HashMap<Point, ArrayList<Point>> G = createGraph(ps);
    // phase 2.1
    int players = 4;

    // create population
    int popSize = 100;
    Graph[] pop = new Graph[popSize];
    for (int i = 0; i < popSize; i++) {
      pop[i] = new Graph(G, players);
    }

    // CoreDraw stats = new CoreDraw(56, 20, 20, "Progression");
    // keep track of progression
    // mak colors to arraylist of doubles
    progression = new HashMap<>();
    for (Color c : colors) {
      progression.put(c, new ArrayList<>());
    }
    _gui.enableDoubleBuffering();
    // stats.enableDoubleBuffering();
    int gen = 0;
    // apply genetic algorithm
    genLoop:
    while (true) {
      if (gen > 200) {
        Core.log("failure. Restart.");
        return;
      }
      gen++;
      Core.log("beginning gen number "+gen+"...");
      // for each graph, map to its ratio
      HashMap<Graph, Double> assess = new HashMap<>();
      // Core.freeze("assessing population");
      for (Graph g : pop) {
        assess.put(g, g.cost());
      }
      // now rank based on minimum order
      ArrayList<Graph> ranked = Core.sort(assess, "min");
      // get the top rated
      Graph best = ranked.get(0);
      double bestRatio = assess.get(best);
      Core.log("bestRatio = " + bestRatio);
      // check if ratio is alright
      if (bestRatio < 0.05) {
        // done
        Core.log("victory! : " + bestRatio);
        // draw the successful graph
        best.draw(_gui);
        _gui.show();
        // save data
        _vlis = best.vlis;
        _graph = best.graph;
        _subgraphs = best.subgraphs;
        _playerColors = best.playerColors;
        _spawnEdges = best.spawnEdges;
        break genLoop;
      }
      // we are unsuccessful, we must now enter the mating
      // keep top 10%
      // int top10 = (int)(0.1 * popSize);
      int top10 = 0;
      // for (int i = 0; i < top10; i++) {
      //   pop[i] = ranked.get(i);
      // }
      for (int i = top10; i < popSize; i++) {
        // we must make this a new graph created from the top10
        // for now. lets just make them random
        // pop[i] = new Graph(ranked.get(Core.randInt(top10)));
        // pop[i] = new Graph(pop[Core.randInt(top10)]);
        pop[i] = new Graph(G, players);
        pop[i].lloyd(2);
      }
    }
  }

  public static void drawProgression(CoreDraw stats) {
    // filter down to what can fit on the screen

    // draw up the progression
    int x = progression.get(CoreDraw.RED).size();
    double min = progression.get(CoreDraw.RED).get(0);
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
      stats.clear(CoreDraw.WHITE);
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
