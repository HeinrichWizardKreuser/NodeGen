import java.util.*;
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

    System.out.println("left = " + left);
    System.out.println("right = " + right);
    System.out.println("top = " + top);
    System.out.println("bot = " + bot);

    // find ideal box
    double left_d = x0 + (x1 - x0) * 0.05;
    double right_d = x1 - (x1 - x0) * 0.05;
    double bot_d = y0 + (y1 - y0) * 0.05;
    double top_d = y1 - (y1 - y0) * 0.05;

    System.out.println("left_d = " + left_d);
    System.out.println("right_d = " + right_d);
    System.out.println("top_d = " + top_d);
    System.out.println("bot_d = " + bot_d);


    // find multipliers
    double vertical_multiplier = (top_d - bot_d) / (top - bot);
    double horizontal_multiplier = (right_d - left_d) / (right - left);
    System.out.println("vertical_multiplier = " + vertical_multiplier);
    System.out.println("horizontal_multiplier = " + horizontal_multiplier);


    // find the centre of the current box
    double centre_x = left + (right - left) * 0.5;
    double centre_y = bot + (top - bot) * 0.5;
    System.out.println("centre_x = " + centre_x);
    System.out.println("centre_y = " + centre_y);

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
    System.out.println("centre_x_d = " + centre_x_d);
    System.out.println("centre_y_d = " + centre_y_d);

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

    Core.freeze("applying post fix");

    // post fix
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
            }
          }
        }
      }
    }

    StdDraw.clear(StdDraw.WHITE);
    drawPoints(ps);
    drawEdges(G);


    //ArrayList<Point[]> tris = CoreGeom.getTriangles(del);

    // remove edges connected to vertices of highest degree not equal to 1 or 3




    //HashMap<Point, ArrayList<Point>> map = nodeGen(ps);

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

  private static void drawEdges(HashMap<Point, ArrayList<Point>> del) {
    StdDraw.setPenColor(StdDraw.BLACK);
    for (Point p : del.keySet()) {
      for (Point adj : del.get(p)) {
        StdDraw.line(p, adj);
      }
    }
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
