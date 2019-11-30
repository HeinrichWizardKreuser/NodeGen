import java.util.*;

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

    testVoronoiGraph(subgraphs, graph);
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
}
