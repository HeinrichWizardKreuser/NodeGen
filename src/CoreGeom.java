/*******************************************************************************

Offers methods for computational geomtetry clients

*******************************************************************************/


//http://www.geom.uiuc.edu/~samuelp/del_project.html
//import java.util.ArrayList;
//import java.util.HashMap;
import java.util.*;

//for library & sorting
import java.util.HashMap;
import java.util.ArrayList;
import java.util.HashSet;
//for sorting
import java.util.Collections;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.Map;
import static java.util.stream.Collectors.*;
import static java.util.Map.Entry.*;

public class CoreGeom {

  /***************************************************************************
   *                        VISIBILIY
   ***************************************************************************/
  /**
   * Checks wether the specified corner (one of the corners of the given shape)
   * is visible by the given eye (outside of the shape)
   * @param corner the corner of the shape the eye wants to see
   * @param shape the array of points representing the ordered Corners of some polygon
   * @param eye the Point that we want to check whether the corner is visible to
   */
  public static boolean isVisible(Point corner, Point[] shape, Point eye) {
    Angle eye_corner = new Angle(eye, corner);
    Angle min = null;
    Angle max = null;
    Angle centre = null;
    // get all points that aren't corner on a hashmap mapped to their distances to p
    HashMap<Point, Double> sortMe = new HashMap<>();
    for (Point p : shape) {
      if (p != corner) sortMe.put(p, corner.dist(p));
    }
    sortMe.put(eye, corner.dist(eye));
    ArrayList<Point> sorted = Core.sort(sortMe, "min");
    // iterate over all points in order of closest distance to p
    for (Point p : sorted) {
      if (p == eye) {
        if (centre == null) {
          return true;
        } else {
          return !(min.lessThan(eye_corner) && eye_corner.lessThan(max));
        }
      }
      Angle p_corner = new Angle(p, corner);
      if (min == null) {
        min = p_corner;
      } else if (centre == null) {
        max = p_corner;
        Angle a = min.copy();
        Angle b = max.copy();
        // now find out which one is actually min and which is actually max
        // get centre angle
        centre = new Angle(a, b);
        // create a point in direction of centre
        Point centreP = new Point(centre.angle, 1);
        // if that point is not inside shape, then we need to invert the angle
        if (!centreP.isInside(shape)) {
          centre = new Angle(centre.angle+180);
        }
        if (b.lessThan(centre)) {
          min = b;
          max = a;
        }
        // reset min and max
        // with new angle, shoot in that direction as far as you can
      } else {
        if ((p_corner).lessThan(min)) min = p_corner;
        if ((max).lessThan(p_corner)) max = p_corner;
      }
    }
    return !(min.lessThan(eye_corner) && eye_corner.lessThan(max));
  }

  /***************************************************************************
   *                        EQUATIONS
   ***************************************************************************/
  public static double[] eq(Point a, Point b) {
    return new double[]{
      a.y-b.y,
      b.x-a.x,
      (a.y-b.y)*b.x + (b.x-a.x)*b.y
    };
  }
  public static double[] eq(Point... ab) {
    return eq(ab[0], ab[1]);
  }

  public static Point intersection(double[] ab, double[] cd) {
    double det = ab[0]*cd[1] - cd[0]*ab[1]; //determinant
    if (det == 0) {
      return null; // parallel lines
    }
    return new Point(
      (cd[1]*ab[2] - ab[1]*cd[2])/det,
      (ab[0]*cd[2] - cd[0]*ab[2])/det
    );
  }

  public static Point intersection(Point[] ab, Point[] cd) {
    double[] eqab = eq(ab), eqcd = eq(cd);
    Point intersection = intersection(eqab, eqcd);
    if (intersection == null) {
      return null; // parallel lines
    }
    Point a = ab[0], b = ab[1];
    if (Math.min(a.x, b.x) <= intersection.x && intersection.x <= Math.max(a.x, b.x)
    &&  Math.min(a.y, b.y) <= intersection.y && intersection.y <= Math.max(a.y, b.y)) {
      return intersection; // lies inbetween segment interval
    }
    return null; // equations intersect, but not inside segment interval
  }

  /**
   * Gets the segment bisector. This is the line perpendicular to line ab and
   * passes through the midpoint of ab
   */
  public static double[] bisector(Point a, Point b) {
    Point mid = new Point(a, b);
    double[] eq = eq(a, b);
    return new double[]{
      -eq[1],
      eq[0],
      -eq[1] * mid.x + eq[0] * mid.y
    };
  }

  /**
   * Gets equation of line passing through the corner abc
   * Otherwise known as the angle bisector
   */
  public static double[] bisector(Point a, Point b, Point c) {
    // get the bisector angle af a_b and c_b
    Angle a_b = new Angle(a, b), c_b = new Angle(c, b);
    Angle mid = new Angle(a_b, c_b);
    // get the equation of the line passing through b and one in the direction of mid
    return eq(b, b.directed(mid.angle, 10));
  }

  /**
   * Calculates the circumcentre of a collection of points making up a triangle.
   * @param a,b,c the 3 corners of a triangle
   * @return the coordinate of the circumcentre of the triangle a-b-c
   */
  public static Point cc(Point a, Point b, Point c) {
    double[] abbi = bisector(a, b), bcbi = bisector(b, c);
    return intersection(abbi, bcbi);
  }
  public static Point cc(Point... tri) {
    return cc(tri[0], tri[1], tri[2]);
  }

  // get x from equation
  public static double x(double[] eq, double y) {
    double a = eq[0], b = eq[1], c = eq[2];
    if (a == 0d) {
      return Double.POSITIVE_INFINITY;
    }
    return (c-b*y)/a;
  }

  // get y from equation
  public static double y(double[] eq, double x) {
    double a = eq[0], b = eq[1], c = eq[2];
    if (b == 0d) {
      return Double.POSITIVE_INFINITY;
    }
    return (c-a*x)/b;
  }

  /***************************************************************************
   *                        BIG METHODS
   ***************************************************************************/
  /**
   * Adds points inbetween every point in the corners of the shape so that the
   * distance between all the points in the new shape are less than or equal to
   * the given density
   */
  private static Point[] densify(Point[] shape, double density) {
    ArrayList<Point> densify = new ArrayList<>();
    int n = shape.length;
    int total = 0;
    for (int i = 0; i < n; i++) {
      int j = (i+1 == n) ? 0 : i+1;
      Point start = shape[i];
      Point end = shape[j];
      // get distance between them
      double dist = start.dist(end);
      // check how many we can fit in
      int amnt = (int)Math.ceil(dist / density);
      // get distance between each
      double spacing = dist / (double)amnt;
      // get direction of of the end relative to start
      double degrees = new Angle(end, start).angle;
      // place each point with distance spacing from each other
      for (double space = 0; space < dist; space += spacing) {
        // create a point in direction of degrees
        Point p = start.directed(degrees, space);
        if (Point.epsilon(p.x, (int)p.x)) {
          p.x = (int)p.x;
        }
        if (Point.epsilon(p.y, (int)p.y)) {
          p.y = (int)p.y;
        }
        densify.add(p);
        total++;
      }
    }
    Point[] toReturn = new Point[total];
    for (int i = 0; i < total; i++) {
      toReturn[i] = densify.get(i);
    }
    return toReturn;
  }

  /**
   * Inserts the least amount of points neccesary between s and e in order to
   * connect them and still be legal delaunay triangulation.
   * For when two points are intended to be connected after triangulation but
   * are not, thus we add the minimal amount of points so that they become
   * connected via mutual points via delaunay triangulation
   *
   * @param s is the starting point
   * @param e is the end point
   * @param all is the hashmap hosting all of the points mapped to a list
   *        of all points that they are connected to
   */
  private static void correctConnection(Point s, Point e, HashMap<Point, ArrayList<Point>> all) {
    while (!all.get(s).contains(e)) {
      // find angle left to and right to e relative to s
      Point min = null, max = null;
      double minDiff = 180d, maxDiff = 180d;
      Angle e_s = new Angle(e, s);
      for (Point p : all.get(s)) {
        Angle p_s = new Angle(p, s);
        double absDiff = p_s.absDiff(e_s);
        if ((p_s).lessThan(e_s)) { // min category
          if (absDiff < minDiff) {
            min = p;
            minDiff = absDiff;
          }
        } else { // max category
          if (absDiff < maxDiff) {
            max = p;
            maxDiff = absDiff;
          }
        }
      }
      if (!all.get(min).contains(max)) {
        System.out.println("min must be connected to max");
        System.exit(0);
      }
      // get circumcentre c of min-s-max
      Point c = cc(min, s, max);
      //create point p to insert
      Angle c_s = new Angle(c, s);
      double diff = c_s.absDiff(e_s);
      double d = 0.99*2*c.dist(s)*Angle.cos(diff);
      if (d > s.dist(e)) d = s.dist(e)/2d;
      Point p = s.directed(e_s.angle, d);
      // insert the point
      insert(p, all);
      if (!all.get(s).contains(p)) {
        System.out.println("!s.contains(p)");
        System.exit(0);
      }
      // set the starting node as p
      s = p;
    }
  }

  /**
   * Inserts the given point into the delaunay triangulation and then updates
   * any triangulations neccesary
   *
   * @param toInsert is the point to insert into the triangulation connectAll
   * @param connectAll is the hashmap hosting all of the points mapped to a list
   * of all points that they are connected to
   */
  private static void insert(Point toInsert, HashMap<Point, ArrayList<Point>> connectAll) {
    connectAll.put(toInsert, new ArrayList<Point>(0));
    // determine triangle
    Point[] triangle = getTriangleContaining(toInsert, connectAll);
    // setup connectections from toinsert to all points
    for (Point corner : triangle) {
      connectAll.get(toInsert).add(corner);
      connectAll.get(corner).add(toInsert);
    }
    // get all angles including toInsert
    ArrayList<Point[]> corners = new ArrayList<>(0);
    corners.add(new Point[]{triangle[0], toInsert, triangle[1]});
    corners.add(new Point[]{triangle[1], toInsert, triangle[2]});
    corners.add(new Point[]{triangle[2], toInsert, triangle[0]});
    // update connections as neccesary
    while (!corners.isEmpty()) {
      // create list of next corners we may need to iterate over
      ArrayList<Point[]> nextCorners = new ArrayList<>(0);
      // iterate through all corners
      cornerLoop:
      for (Point[] corner : corners) {
        Point a = corner[0];
        Point b = corner[2];
        // get point c
        Point c = null;
        for (Point p : connectAll.get(a)) {
          if (connectAll.get(p).contains(a) && p != toInsert) {
            Angle a_toInsert = new Angle(a, toInsert);
            Angle b_toInsert = new Angle(b, toInsert);
            Angle min = (a_toInsert.lessThan(b_toInsert)) ? a_toInsert : b_toInsert;
            Angle max = (min == a_toInsert) ? b_toInsert : a_toInsert;
            Angle p_toInsert = new Angle(p, toInsert);
            if (min.lessThan(p_toInsert) && p_toInsert.lessThan(max)) {
              if (c != null) c = c.dist(toInsert) < p.dist(toInsert) ? c : p;
              else c = p;
            }
          }
        }
        if (c == null) continue cornerLoop;
        Point cc = cc(a, toInsert, b);
        // if a-toInsert-b circle contains c
        if (cc.dist(c) < cc.dist(toInsert)) {
          // sever connection a-b
          connectAll.get(a).remove(b);
          connectAll.get(b).remove(a);
          // create connection toInsert-c
          connectAll.get(c).add(toInsert);
          connectAll.get(toInsert).add(c);
          // add new corners created to next list
          nextCorners.add(new Point[]{a, toInsert, c});
          nextCorners.add(new Point[]{c, toInsert, b});
        }
      }
      // set corners list equal to next
      corners = nextCorners;
    }
  }

  /**
   * Returns a length3 Point array representing the triangle which contains the
   * given point
   *
   * @param point is the point we are searching the host triangle for
   * @param all is the hashmap hosting all of the points mapped to a list
   *        of all points that they are connected to
   */
  private static Point[] getTriangleContaining(Point point, HashMap<Point, ArrayList<Point>> all) {
    // get sorted list of all points closest to point
    HashMap<Point, Double> sortMe = new HashMap<>();
    for (Point p : all.keySet()) {
      sortMe.put(p, p.dist(point));
    }
    ArrayList<Point> closest = Core.sort(sortMe, "min");
    if (closest.contains(point)) closest.remove(point);
    int n = closest.size();
    // iterate over all collection of 3 points in order of likely connected
    maxLoop:
    for (int max = 2; max < n; max++) {
      Point c = closest.get(max);
      minLoop:
      for (int min = 0; min < max-1; min++) {
        Point a = closest.get(min);
        if (!all.get(c).contains(a)) continue minLoop;
        jLoop:
        for (int j = min+1; j < max; j++) {
          Point b = closest.get(j);
          if (!all.get(b).contains(a) || !all.get(b).contains(c)) continue jLoop;
          // check if point is contained within triangle:
          double areaSum = Point.area(a,b,point) + Point.area(a,point,c) + Point.area(point,b,c);
          if (Point.epsilon(Point.area(a, b, c), areaSum)) {
            return new Point[]{ a, b, c};
          }
        }
      }
    }
    Core.log("Could not find triangle containing " +point.toString());
    (new int[]{1})[1] = 1;;
    return null;
  }

  public static void createConnection(Point a, Point b, HashMap<Point, ArrayList<Point>> adj) {
    adj.get(a).add(b);
    adj.get(b).add(a);
  }

  public static void severConnection(Point a, Point b, HashMap<Point, ArrayList<Point>> adj) {
    adj.get(a).remove(b);
    adj.get(b).remove(a);
  }

  /**
   * Gets all of the triangles on the triangulated given tiangulated map
   *
   * @param all is the adjacency list of all points mapped to a list of points
   *        that they are connected/adj to
   * @return a list of all triangles(in the form of length3 point arrays)
   */
  public static ArrayList<Point[]> getTriangles(HashMap<Point, ArrayList<Point>> all) {
    ArrayList<Point[]> triangles = new ArrayList<>(0);
    // create copy of all and name it 'adj'
    HashMap<Point, ArrayList<Point>> adj = new HashMap<>();
    for (Point p : all.keySet()) {
      ArrayList<Point> l = new ArrayList<>(all.get(p).size());
      for (Point adjP : all.get(p)) {
        l.add(adjP);
      }
      adj.put(p, l);
    }
    // for each point, create triangle via adj list
    while (!adj.isEmpty()) {
      // get v = point with adj list of shortest length
      Point v = null;
      int vLen = 0;
      // get point with shortest adj list length
      for (Point p : adj.keySet()) {
        int pLen = adj.get(p).size();
        if (v == null || pLen < vLen) {
          v = p;
          vLen = pLen;
        }
      }
      // create triangles with all unique pairs in adjV that are also adj to each other
      ArrayList<Point> adjV = adj.get(v); // pointer-shortcut

      for (Point[] pair : uniquePairs(adjV)) {
        Point s = pair[0];
        Point t = pair[1];
        if (adj.get(s).contains(t)) {
          triangles.add(new Point[]{v, s, t});
        }
      }
      // remove v from all in adjV's adj lists, then remove v from adjmap
      for (Point w : adjV) {
        adj.get(w).remove(v);
      }
      adj.remove(v);
    }
    return triangles;
  }

  public static ArrayList<Point[]> uniquePairs(ArrayList<Point> l) {
    ArrayList<Point[]> uniquePairs = new ArrayList<>(0);
    int n = l.size();
    for (int i = 0; i < n; i++) {
      for (int j = i+1; j < n; j++) {
        uniquePairs.add(new Point[]{l.get(i), l.get(j)});
      }
    }
    return uniquePairs;
  }

  public static ArrayList<Point> randomPoints(double w, double h, int n, double rangeLimit) {
    ArrayList<Point> ps = new ArrayList<Point>(n);
    int count = 0;
    int counter = 0;
    generating:
    for (int i = 0; i < n; i++) {
      if (counter > 100 * n) {
        // then we have been trying this for way too long, just complain
        throw new IllegalArgumentException("Couldn't fit " + n +
          " points into space of " + w + " x " + h);
      }
      Point p = new Point(Math.random() * w, Math.random() * h);
      for (Point o : ps) {
        if (p.dist(o) < rangeLimit) {
          counter++;
          count++;
          i--;
          continue generating;
        }
      }
      ps.add(p);
      counter = 0;
    }
    Core.log("count = " + count);
    return ps;
  }

  public static ArrayList<Point> randomPoints(double w, double h, int n) {
    return randomPoints(w, h, n, 0);
  }

  /**
   * Given a list of nodes on the given graph, it performs voronoi on the graph
   * in such a way that we get all of the nodes mapped to the node in the list
   * that it is closest to
   *
   * @param vlis is the list of the nodes that will represent the sites for the
   *        graph
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adj to in the graph
   * @return a mapping of each point P in vlis to a subgraph. Each subgraph is a
   *         part of the graph, but where all nodes are closer to P than any other
   *         node in vlis. It will also contain new leaf nodes created on edges
   *         where the end points are closer to different nodes in vlis. Thus
   *         these new leaf nodes are created to have more precision
   */
  public static HashMap<Point, HashMap<Point, ArrayList<Point>>> voronoiGraph(
    ArrayList<Point> vlis,
    HashMap<Point, ArrayList<Point>> graph) {

    // voro = map of distance between vlis nodes and their distnaces to other points in the graph
    HashMap<Point, HashMap<Point, Double>> voro = new HashMap<Point, HashMap<Point, Double>>();

    // perform dijkstra for each node in voro
    for (Point a : vlis) {
      voro.put(a, dijkstraGraph(a, graph));
    }

    // closestTo = all points in G mapped to the point in voro they are closest to
    HashMap<Point, Point> closestTo = new HashMap<Point, Point>();
    // edgeConflicts = list of all edges that have endpoints closer to opposing owners
    ArrayList<Point[]> edgeConflicts = new ArrayList<>();
    // edgeFriendly = list of all edges of which both ends belong to the same owners
    ArrayList<Point[]> edgeFriendly = new ArrayList<>();
    // find the owner to which each node is closest to
    for (Point p : graph.keySet()) {
      // find which owner p is closest to
      HashMap<Point, Double> cost = new HashMap<>();
      for (Point a : voro.keySet()) {
        cost.put(a, voro.get(a).get(p));
      }
      closestTo.put(p, Core.getMinKey(cost));

      // check if p is adj to an opposing owner
      for (Point adj : graph.get(p)) {
        // if adj has also been treated
        if (closestTo.containsKey(adj)) {
          // are they closer to different owners?
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
    for (Point[] e : edgeConflicts) {
      double d = e[0].dist(e[1]);
      double len = (voro.get(closestTo.get(e[0])).get(e[0]) +
                    voro.get(closestTo.get(e[1])).get(e[1]) + d)/2 -
                    voro.get(closestTo.get(e[0])).get(e[0]);

      Point border = new Point(e[0], e[1], 1 - len/d);
      borders.put(e, border);
    }

    // create subgraph of each owner
    HashMap<Point, HashMap<Point, ArrayList<Point>>> subgraphs = new HashMap<>();
    for (Point p : voro.keySet()) {
      subgraphs.put(p, new HashMap<>());
    }
    // add each point
    for (Point p : graph.keySet()) {
      // check which owner this belongs to
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

  /**
   * Given a starting node on a given graph, it performs dijkstra's algorithm.
   *
   * @param a the starting node to perform dijkstra from
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adjacent to in the graph
   * @return a cost matrix. It is a mapping of all points in the graph to the
   *         value of their shortest path to a.
   */
  public static HashMap<Point, Double> dijkstraGraph(
    Point a,
    HashMap<Point, ArrayList<Point>> graph) {

    if (a == null) {
      throw new IllegalArgumentException("a is null");
    }
    if (!graph.containsKey(a)) {
      throw new IllegalArgumentException("graph does not contain a");
    }

    ArrayList<Point> curr = new ArrayList<>();
    curr.add(a);
    HashMap<Point, Double> cost = new HashMap<>();
    cost.put(a, 0d);
    HashSet<Point> explored = new HashSet<>();
    while (!curr.isEmpty()) {
      // find point with lowest score in cost
      // HashMap<Point, Double> map = Core.hashMap(curr, (p) -> cost.get(p));
      HashMap<Point, Double> map = new HashMap<>();
      for (Point p : curr) {
        map.put(p, cost.get(p));
      }
      Point p = Core.getMinKey(map);
      curr.remove(p);
      // explore point
      explored.add(p);
      // now go to all adjacent points
      if (p == null) {
        Core.exit("p == null");
      }
      for (Point adj : graph.get(p)) {
        double newDist = cost.get(p) + p.dist(adj);
        if (explored.contains(adj) && newDist > cost.get(adj)) {
          continue;
        }
        curr.add(adj);
        explored.add(adj);
        cost.put(adj, newDist);
      }
    }
    return cost;
  }

  /**
   * Given a starting starting and ending node on a given graph, it performs
   * A* pathfinding. It is an adaption of dijkstra for only finding the path
   * to another node and not necessarily the cost of all nodes from a certain
   * node.
   *
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adjacent to in the graph.
   * @param start the starting node.
   * @param end the node to pathfind to.
   * @return an arraylist that represents the shortest path from start to end on
   *         the graph.
   */
  public static ArrayList<Point> aStar(
    HashMap<Point, ArrayList<Point>> graph,
    Point start,
    Point end) {

    // check that start and end are in the graph
    if (!graph.containsKey(start)) {
      throw new IllegalArgumentException("start node not in graph!");
    }
    if (!graph.containsKey(end)) {
      throw new IllegalArgumentException("end node not in graph!");
    }

    ArrayList<Point> curr = new ArrayList<>();
    curr.add(start);
    HashMap<Point, Double> cost = new HashMap<>();
    cost.put(start, start.dist(end));
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

    // get the path by popping the parents
    ArrayList<Point> path = new ArrayList<>();
    Point parent = end;
    do {
      path.add(parent);
      parent = parents.get(parent);
    } while (parent != null);
    Collections.reverse(path);
    return path;
  }

  /**
   * Finds the absolute center of the graph, assuming that the graph is a tree.
   * It creates a new node to represent this and returns it along with the end
   * points of the edge it exists on.
   *
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adjacent to in the graph
   * @return an ordered size3 array that contains a new node that represents the
   *         centroid at index 1. Index 0 and 2 contain the endpoints of the edge.
   */
  public static Point[] absoluteCentre(
    HashMap<Point, ArrayList<Point>> graph) {

    // pick a random point
    Point rand = Core.randomKey(graph);
    // get the node with furthest distance from rand
    HashMap<Point, Double> cost = dijkstraGraph(rand, graph);
    Point v0 = Core.getMaxKey(cost);
    // now get the firthest distance from that node
    cost = CoreGeom.dijkstraGraph(v0, graph);
    Point v1 = Core.getMaxKey(cost);
    // now get path from v0 to v1
    ArrayList<Point> path = aStar(graph, v0, v1);
    double pathLength = CoreGeom.lineLength(path);
    // now find the midway point
    double toGo = pathLength / 2;
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
        return new Point[]{curr, centroid, next};
      }
    }

    throw new RuntimeException("Could not path find");
  }

  /**
   * Calculates the length of the graph. That is, the sum lengtt of all edges
   * in the graph.
   *
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adjacent to in the graph
   * @return the sum total distance of all edges on the graph.
   */
  public static double graphLength(HashMap<Point, ArrayList<Point>> graph) {
    double graphLength = 0;
    for (Point p : graph.keySet()) {
      for (Point adj : graph.get(p)) {
        graphLength += p.dist(adj);
      }
    }
    return graphLength/2;
  }

  /**
   * Checks that if a point a is connected to b that b is also connect to a
   *
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adjacent to in the graph
   * @return true if the graph is an undirected graph.
   */
  public static boolean isUndirected(HashMap<Point, ArrayList<Point>> graph) {
    // check that all nodes have consensus
    for (Point p : graph.keySet()) {
      for (Point adj : graph.get(p)) {
        if (!graph.get(adj).contains(p)) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Checks if any node can travel to any other node in the graph.
   * Assumes that the graph is undirected.
   *
   * @param graph is all of the points on the graph mapped to all other points
   *        that they are adjacent to in the graph
   * @return true if the graph is a connected graph.
   */
  public static boolean isConnected(HashMap<Point, ArrayList<Point>> graph) {
    if (graph.isEmpty()) {
      throw new IllegalArgumentException("graph is empty");
    }
    Stack<Point> stack = new Stack<>();
    stack.push(Core.randomKey(graph));
    HashSet<Point> explored = new HashSet<>();
    while (stack.isEmpty()) {
      Point p = stack.pop();
      for (Point adj : graph.get(p)) {
        if (!explored.contains(adj)) {
          explored.add(adj);
          stack.push(adj);
        }
      }
    }
    return explored.size() != graph.size();
  }

  public static void sever(Point a, Point b, HashMap<Point, ArrayList<Point>> graph) {
    graph.get(a).remove(b);
    graph.get(b).remove(a);
  }

  public static void connect(Point a, Point b, HashMap<Point, ArrayList<Point>> graph) {
    graph.get(a).add(b);
    graph.get(b).add(a);
  }

  public static void remove(Point toRemove, HashMap<Point, ArrayList<Point>> graph) {
    // sever connections p -> toRemove
    for (Point p : graph.get(toRemove)) {
      graph.get(p).remove(toRemove);
    }
    // remove toRemove (which will remove toRemove as well as its adjacency list)
    graph.remove(toRemove); // removes toRemove -> adjacent points
  }

  public static void addBetween(Point e0, Point e1, HashMap<Point, ArrayList<Point>> graph, Point toAdd) {
    sever(e0, e1, graph);
    graph.put(toAdd, new ArrayList<Point>(0));
    connect(toAdd, e0, graph);
    connect(toAdd, e1, graph);
  }

  public static HashMap<Point, ArrayList<Point>> copy(HashMap<Point, ArrayList<Point>> graph) {
    // make copy of graph
    HashMap<Point, ArrayList<Point>> copy = new HashMap<>();
    for (Point p : graph.keySet()) {
      ArrayList<Point> adj = new ArrayList<>();
      for (Point a : graph.get(p)) {
        adj.add(a);
      }
      copy.put(p, adj);
    }
    return copy;
  }

  public static ArrayList<Point[]> edgeList(HashMap<Point, ArrayList<Point>> graph) {
    // make a copy
    HashMap<Point, ArrayList<Point>> copy = copy(graph);
    // create the edge list
    ArrayList<Point[]> edgeList = new ArrayList<Point[]>();

    for (Point toRemove : graph.keySet()) {
      for (Point adj : copy.get(toRemove)) {
        edgeList.add(new Point[]{toRemove, adj});
        copy.get(adj).remove(toRemove);
      }
    }
    return edgeList;
  }

  public static void fix(HashMap<Point, ArrayList<Point>> graph) {
    // create a copy of the graph
    HashMap<Point, ArrayList<Point>> fixed = new HashMap<>();
    for (Point p : graph.keySet()) {
      fixed.put(p, new ArrayList<Point>());
    }
    // check if p has null
    for (Point p : graph.keySet()) {
      if (graph.get(p) == null) {
        graph.put(p, new ArrayList<Point>());
      }
    }
    // now go through connections and add them
    for (Point p : graph.keySet()) {
      // just add all of the connections both ways
      for (Point adj : graph.get(p)) {
        fixed.get(adj).add(p);
        fixed.get(p).add(p);
      }
    }
    // we have created duplicates but now we will remove them
    for (Point p : graph.keySet()) {
      // check that there are no duplicates
      ArrayList<Point> lis = graph.get(p);
      int N = lis.size();
      outer:
      for (int i = 0; i < N; i++) {
        Point curr = lis.get(i);
        inner:
        for (int j = i + 1; j < N; j++) {
          Point next = lis.get(j);
          // check if they are equal
          if (curr == next) {
            continue outer;
          }
        }
        // then there are no duplicates of curr. add curr
        fixed.get(p).add(curr);
      }
    }
    // check that there are no connections to nodes that don't exist
    for (Point p : fixed.keySet()) {
      // check if this
      for (int i = 0; i < fixed.get(p).size(); i++) {
        // check if adj is in the graph
        Point adj = fixed.get(p).get(i);
        if (!graph.containsKey(adj)) {
          fixed.get(p).remove(adj);
          i--;
        }
      }
    }
    graph = fixed;
  }



  /**
   * Calculates the length of the given line.
   *
   * @param line is a list of points where each point is connected to the points
   *        with index 1 less or 1 more than it.
   * @return the calculated path length of the line.
   */
  public static double lineLength(ArrayList<Point> line) {
    double lineLength = 0;
    int N = line.size();
    for (int i = 0; i < N - 1; i++) {
      Point curr = line.get(i), next = line.get(i + 1);
      lineLength += curr.dist(next);
    }
    return lineLength;
  }
}
