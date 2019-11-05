import java.util.*;

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

    if (size % 2 != 0) {
      throw new IllegalArgumentException("Size must be even, " + size + " is not");
    }

    // declare space
    int w = 10;
    int h = 10;

    // create a list of random nodes in the space
    ArrayList<Point> ps = CoreGeom.randomPoints(w, h, size, 1);
    // ArrayList<Point> ps = new ArrayList<>(size);
    // Random r = new Random();
    // for (int i = 0; i < size; i++) {
    //   Point p = new Point(r.nextDouble() * w, r.nextDouble() * h);
    //   ps.add(p);
    // }

    // draw the points
    int pixls = 80;
    StdDraw.setCanvasSize(w * pixls, h * pixls);
    StdDraw.setXscale(0, w);
    StdDraw.setYscale(0, h);


    drawPoints(ps);

    HashMap<Point, ArrayList<Point>> del = Delaunay.delaunize(ps);

    //ArrayList<Point[]> tris = CoreGeom.getTriangles(del);

    // remove edges connected to vertices of highest degree not equal to 1 or 3

    loop:
    while (true) {
      // check if we are done
      conditionCheck:
      while (true) {
        for (Point p : del.keySet()) {
          // check if not degree 1 or degree 3
          if (Math.abs(del.get(p).size() - 2) != 1) { // |3-2| = |1-2| = 1
            break conditionCheck;
          }
        }
        // we're done
        break loop;
      }

      // sort verticies based on degree
      HashMap<Point, Double> degrees = new HashMap<>();
      for (Point p : del.keySet()) {
        degrees.put(p, (double)del.get(p).size());
      }
      ArrayList<Point> sorted = Core.sort(degrees, "max");
      Point highestDegree = sorted.get(0);

      // get the highest degree adjacent node
      HashMap<Point, Double> degrees2 = new HashMap<>();
      for (Point adj : del.get(highestDegree)) {
        degrees2.put(adj, degrees.get(adj));
      }
      ArrayList<Point> sorted2 = Core.sort(degrees2, "max");
      Point highestDegree2 = sorted2.get(0);

      // delete edge bro
      del.get(highestDegree).remove(highestDegree2);
      del.get(highestDegree2).remove(highestDegree);


      StdDraw.clear(StdDraw.WHITE);
      drawPoints(ps);
      drawEdges(del);
      Core.freeze();

/*
      stopCondition:
      while (true) {
        for (Point p : del.keySet()) {
          // check if less than 3

          if (Math.abs(del.get(p).size() - 2) != 1) { // |3-2| = |1-2| = 1
            break stopCondition;
          }
        }
        // we're done
        break loop;
      }
*/
      /*
      // based on highest degree, find adjacent with highest degree
      for (Point p : sorted) {
        // get adj point with highest degree
        HashMap<Point, Double> degrees2 = new HashMap<>();
        for (Point adj : del.get(p)) {
          degrees2.put(adj, degrees.get(adj));
        }
        if (degrees2.isEmpty()) {
          Core.exit("degrees2.isEmpty()");
        }
        ArrayList<Point> sorted2 = Core.sort(degrees2, "max");
        // if the highest degree is optimal, this is a failure. How do we solve this?
        // if it is 3, then still remove. only if it is 1, is it a failure
        if (sorted2.get(0) == 1) {
          Core.exit("degrees2: failure");
        }
        // now sorted, get highest degree
        for (adj : sorted2) {
          // so this is the highest,
        }


        int n = degrees.get(adj);

      }
      */
    }



    //HashMap<Point, ArrayList<Point>> map = nodeGen(ps);

  }

  private static void drawPoints(ArrayList<Point> ps) {
    double radius = 0.1;
    StdDraw.setPenColor(StdDraw.BLACK);
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
