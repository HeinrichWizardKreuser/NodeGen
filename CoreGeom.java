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

}
