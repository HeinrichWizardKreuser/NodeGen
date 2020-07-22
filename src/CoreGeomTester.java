import java.util.*;
import java.util.function.Function;

public class CoreGeomTester {

  public static void main(String[] args) {
    int size = 10;
    if (args.length != 0) {
      try {
        size = Integer.parseInt(args[0]);
      } catch (Exception e) {
        Core.exit("Please enter a valid size as number");
      }
    }

    int w = 10;
    int h = 10;

    // phase 1.1
    ArrayList<Point> ps = NodeGen.generatePoints(w, h, size);

    // draw the points
    int pixls = 80;
    StdDraw.setCanvasSize(w * pixls, h * pixls);
    StdDraw.setXscale(0, w);
    StdDraw.setYscale(0, h);


    HashMap<Point, ArrayList<Point>> G = NodeGen.createGraph(ps);


    // testDijkstra(G);

    int players = 4;
    HashMap<Point, ArrayList<Point>> graph = NodeGen.generatePlayerSpawns(G, players);

    ArrayList<Point> vlis = new ArrayList<>();
    for (Point p : NodeGen.playerColors.keySet()) {
      vlis.add(p);
    }
    // phase 2.2
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs = CoreGeom.voronoiGraph(vlis, graph);

    // testVoronoiGraph(subgraphs, graph);

    // so voronoi and dijsktra work.
    // Now for centroid
    // absoluteCentre(subgraphs, graph);


    // all is in working condition. New build
    // Core.log("done");
  }


  private static void absoluteCentre(
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs,
    HashMap<Point, ArrayList<Point>> graph) {

    forEachSpawn:
    for (Point spawn : NodeGen.playerColors.keySet()) {

      // get subgraph
      HashMap<Point, ArrayList<Point>> subgraph = subgraphs.get(spawn);

      // pick a random point
      Point rand = Core.randomKey(subgraph);
      // get the node with furthest distance from rand
      HashMap<Point, Double> cost = CoreGeom.dijkstraGraph(rand, subgraph);
      Point v0 = Core.getMaxKey(cost);
      // now get the firthest distance from that node
      cost = CoreGeom.dijkstraGraph(v0, subgraph);
      Point v1 = Core.getMaxKey(cost);
      Core.log("v0.path(v1) = " + cost.get(v1));
      // now get patj from v0 to v1
      // ArrayList<Point> path = dijkstraPath(subgraph, v0, v1);
      ArrayList<Point> path = CoreGeom.aStar(subgraph, v0, v1);

      // draw the full grapg
      StdDraw.clear(StdDraw.WHITE);
      StdDraw.drawGraph(graph, 0.1, StdDraw.LIGHT_GRAY);
      // draw the subgraph
      StdDraw.drawGraph(subgraph, 0.075, StdDraw.BLACK);
      // also the spawn pos
      StdDraw.filledCircle(spawn, 0.05, NodeGen.playerColors.get(spawn));

      // draw the starting node and the ending node
      StdDraw.filledSquare(v0.x + 0.25, v0.y, 0.2, StdDraw.WHITE);
      StdDraw.text(v0.x + 0.25, v0.y, "v0", StdDraw.BLACK);
      StdDraw.filledSquare(v1.x + 0.25, v1.y, 0.2, StdDraw.WHITE);
      StdDraw.text(v1.x + 0.25, v1.y, "v1", StdDraw.BLACK);
      // state the path
      double pathLength = CoreGeom.lineLength(path);
      Core.log("pathLength = " + pathLength);
      Core.freeze();

      // now find the midway point
      // travel along the path until it the remaining cost is less than next road
      double toGo = pathLength / 2;
      ArrayList<Point> line = new ArrayList<>();
      for (Point p : path) {
        line.add(p);
      }
      Point curr = path.remove(0);
      while (toGo > 0) {
        Point next = path.remove(0);
        double nextDist = curr.dist(next);
        // if what we still need to travel is more than next dist, then travel
        if (nextDist < toGo) {
          curr = next;
          toGo -= nextDist;
        } else {
          // then next is will not be reached. from curr to next is the final.
          Point centroid = new Point(curr, next, 1 - toGo / curr.dist(next));
          // return centroid;
          // draw centroid
          StdDraw.filledSquare(centroid, 0.1, NodeGen.playerColors.get(spawn));
          // draw the path
          for (int i = 0; i < line.size() - 1; i++) {
            Point a1 = line.get(i), a2 = line.get(i + 1);
            StdDraw.line(a1, a2, NodeGen.playerColors.get(spawn));
          }
          Core.freeze("Found centroid");
          continue forEachSpawn;
        }
      }

      throw new RuntimeException("Could not path find");
      return null;
    }
  }

  public static ArrayList<Point> aStar(
    HashMap<Point, ArrayList<Point>> graph,
    Point start,
    Point end) {

    ArrayList<Point> curr = new ArrayList<>();
    curr.add(start);
    HashMap<Point, Double> cost = new HashMap<>();

    cost.put(start, start.dist(end));
    // cost.put(start, 0d);

    HashSet<Point> explored = new HashSet<>();
    HashMap<Point, Point> parents = new HashMap<>();
    while (!curr.isEmpty()) {
      // find point with lowest score in cost
      HashMap<Point, Double> map = new HashMap<>();
      for (Point p : curr) {
        map.put(p, cost.get(p));
      }
      Point p = Core.getMinKey(map);
      // check if done
      if (p == end) {
        break;
      }
      curr.remove(p);
      // explore point
      explored.add(p);
      // now go to all adjacent points
      for (Point adj : graph.get(p)) {
        double newDist = cost.get(p) + p.dist(adj) + adj.dist(end);
        // double newDist = cost.get(p) + p.dist(adj);
        if (explored.contains(adj) && newDist > cost.get(adj)) {
          continue;
        }
        curr.add(adj);
        explored.add(adj);
        cost.put(adj, newDist);
        parents.put(adj, p);
      }
    }

    // just get path
    ArrayList<Point> path = new ArrayList<>();
    Point parent = end;
    do {
      path.add(parent);
      parent = parents.get(parent);
    } while (parent != null);
    Collections.reverse(path);
    return path;
  }

  private static void testVoronoiGraph(
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs,
    HashMap<Point, ArrayList<Point>> G) {

    Core.freeze();

    for (Point player : subgraphs.keySet()) {

      StdDraw.clear(StdDraw.WHITE);
      NodeGen.drawGraph(G, StdDraw.LIGHT_GRAY);

      HashMap<Point, ArrayList<Point>> graph = subgraphs.get(player);
      NodeGen.drawGraph(graph, 0.05, StdDraw.BLACK);
      StdDraw.filledSquare(0.2, 0.2, 0.2, NodeGen.playerColors.get(player));

      for (Point p : graph.keySet()) {

        StdDraw.filledCircle(p, 0.05, StdDraw.YELLOW);
        for (Point adj : graph.get(p)) {
          StdDraw.filledCircle(adj, 0.05, NodeGen.playerColors.get(player));
        }
        Core.freeze();

        StdDraw.filledCircle(p, 0.05, StdDraw.BLACK);
        for (Point adj : graph.get(p)) {
          StdDraw.filledCircle(adj, 0.05, StdDraw.BLACK);
        }
      }
    }


    // StdDraw.filledSquare(p.x + 0.5, p.y, 0.2, StdDraw.WHITE);
    // StdDraw.text(p.x + 0.5, p.y, String.format("%.1f", cost.get(p)), StdDraw.BLACK);
    // Core.freeze();
    // StdDraw.filledCircle(p, 0.1, StdDraw.BLACK);
  }

  private static void testDijkstra(HashMap<Point, ArrayList<Point>> graph) {
    Point a = Core.randomKey(graph);

    ArrayList<Point> curr = new ArrayList<>();
    curr.add(a);
    HashMap<Point, Double> cost = new HashMap<>();
    cost.put(a, 0d);
    ArrayList<Point> explored = new ArrayList<>();
    while (!curr.isEmpty()) {
      // find point with lowest score in cost
      HashMap<Point, Double> map = new HashMap<>();
      for (Point p : curr) {
        map.put(p, cost.get(p));
      }
      Point p = Core.getMinKey(map);
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
      }

      // draw new cost next to node p
      StdDraw.filledCircle(p, 0.1, StdDraw.RED);
      StdDraw.filledSquare(p.x + 0.5, p.y, 0.2, StdDraw.WHITE);
      StdDraw.text(p.x + 0.5, p.y, String.format("%.1f", cost.get(p)), StdDraw.BLACK);
      Core.freeze();
      StdDraw.filledCircle(p, 0.1, StdDraw.BLACK);
    }
  }



  /*
  public static Point[] graphCentroid(HashMap<Point, ArrayList<Point>> graph) {
    if (isBroken(graph) != null) {
      Core.exit("broken graph detected");
    }
    // map each node to the value of the longest distance to another node in the graph
    HashMap<Point, Double> furthestCosts = new HashMap<>();
    for (Point p : graph.keySet()) {
      if (graph.get(p).isEmpty()) {
        Core.exit("broke");
      }
      // get the cost matrix
      HashMap<Point, Double> cost = dijkstraGraph(p, graph);
      // get the point with the highest value
      Point furthest = Core.getMaxKey(cost);
      double furthestCost = cost.get(furthest);
      // add to mapping
      furthestCosts.put(p, furthestCost);
    }
    // the node with the minimum max cost will be part of the edge containing the centroid
    Point e0 = Core.getMinKey(furthestCosts);

    HashMap<Point, Double> second = new HashMap<>();
    for (Point adj : graph.get(e0)) {
      second.put(adj, furthestCosts.get(adj));
    }
    Point e1 = Core.getMaxKey(second);
    Core.log("second.size() = " + second.size());
    if (second.size() == 0) {
      StdDraw.circle(e0, 0.5, StdDraw.BLACK);
    }

    // // now search among adjacent nodes for the one with the smallest furthest cost
    // Point e1 = Core.randomKey(graph.get(e0));
    // for (Point adj : graph.get(e0)) {
    //   if (furthestCosts.get(adj) < furthestCosts.get(e1)) {
    //     e1 = adj;
    //   }
    // }

    // now we have the edge. Get weighted ave on that line based on the costs
    double e0Cost = furthestCosts.get(e0), e1Cost = furthestCosts.get(e1);

    double d = e0Cost/(e0Cost + e1Cost);

    Point ave = new Point(e0, e1, d);

    return new Point[]{e0, ave, e1};
  }
  */
}
