/*******************************************************************************
 * This is an implementation of a divide and conquer approach to delauna
 * triangulation of a map. The explenation to the algorithm can be found at
 * http://www.geom.uiuc.edu/~samuelp/del_project.html
 * It is NOT my own original work, this is simply my own original Java
 * implementation of it.
 *
 * DEPENDENCIES:
 *  Point.java
 *  Angle.java
 * Available at: https://github.com/HeinrichWizardKreuser/ComputationalGeometry
 *
 * HOW TO USE:
 *  Create an arraylist of Points(coordinates/locations of x,y values).
 *  Call delaunize() and use your arraylist of points as the parameter ('all')
 *  It will return a HashMap of Points mapped to a list of points that they are
 *  connected to after the triangulation (this is called an adjacency list)
 *
 * @author Heinrich Kreuser
 *
 * Date: 28 May 2019
 *
 *         MIT License
 *
 *         Copyright (c) [2019] [Heinrich Kreuser]
 *
 *         Permission is hereby granted, free of charge, to any person obtaining
 *         a copy of this software and associated documentation files (the
 *         "Software"), to deal in the Software without restriction, including
 *         without limitation the rights to use, copy, modify, merge, publish,
 *         distribute, sublicense, and/or sell copies of the Software, and to
 *         permit persons to whom the Software is furnished to do so, subject to
 *         the following conditions:
 *
 *         The above copyright notice and this permission notice shall be
 *         included in all copies or substantial portions of the Software.
 *
 *         THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 *         EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 *         MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 *         NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 *         BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 *         ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 *         CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 *         SOFTWARE.
 ******************************************************************************/
//for library & sorting
import java.util.HashMap;
import java.util.ArrayList;
//for sorting
import java.util.Collections;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import static java.util.stream.Collectors.*;
import static java.util.Map.Entry.*;

public class Delaunay {

  /**
   * Main operator to create a delaunay-triangulation of all of the given points
   *
   * @param all the list of all Points on some map
   * @return a graph in the form of an adjacency list: A list of all points mapped
   *         to a list of all points the are connected to
   */
  public static HashMap<Point, ArrayList<Point>> delaunize(ArrayList<Point> all) {
    // get sorted list based on x position (if same x, take y)
    int n = all.size();
    Point[] sorted = new Point[all.size()];
    for (int i = 0; i < n; i++) {
      sorted[i] = all.get(i);
    }
    Arrays.sort(sorted);
    //divide into different sets
    return triangulate(sorted);
  }

  /**
   * Takes in a set of points and triangulates them using delaunay triangulation.
   * If the set's size > 3,
   *    it divides it in half and recursively feeds itself both halves.
   * If the size <= 3,
   *    it sets up connections between all points in the set.
   * When it gets two halves back, it will send them through merge() before returning it
   * @param set is the set of points to triangulate, sorted in delaunize()
   * @return the triangulated graph where all points are mapped to a list of points
   *         that they are connected to.
   */
  private static HashMap<Point, ArrayList<Point>> triangulate(Point[] set) {
    //divide if more than three
    int setlen = set.length;
    if (setlen > 3){
      int split = setlen/2;
      Point[] left = Arrays.copyOfRange(set, 0, split);
      Point[] right = Arrays.copyOfRange(set, split, setlen);
  		return merge(triangulate(left), triangulate(right));
    }
    // setup triangulation in small scale
    HashMap<Point, ArrayList<Point>> triangulated = new HashMap<>();
    // check for colinear points if 3
    if (setlen == 3 && Point.orientation_i(set[0], set[1], set[2]) == Point.COLINEAR) {
      // find the point at the center
      Arrays.sort(set); // sorting it will make that the centre point is at index 1 due to colinearcy
      triangulated.put(set[0], new ArrayList<>(Arrays.asList(set[1])));
      triangulated.put(set[1], new ArrayList<>(Arrays.asList(set[0], set[2])));
      triangulated.put(set[2], new ArrayList<>(Arrays.asList(set[1])));
    } else { //regular setup if non colinear points: connect all point to all other points
      int len = set.length;
      for (int i = 0; i < len; i++) { // all indices
        ArrayList<Point> l = new ArrayList<Point>(0);
        for (int j = 0; j < i; j++) { // all before curr index
          l.add(set[j]);
        }
        for (int j = i+1; j < len; j++) { // all after curr index
          l.add(set[j]);
        }
        triangulated.put(set[i], l);
      }
    }
    return triangulated;
  }

  /**
   * Evaluates the current L and R and makes sure that they are the correct
   * choices by using 4 extra theorems
   * @param left,right are the left and right triangulated (so-far) adjacency lists
   * @param l,r are the current L and R points
   * @return the new/confirmed L and R in the form of an Point[]{L, R}
   */
  private static Point[] getLR(
      HashMap<Point, ArrayList<Point>> left,
      HashMap<Point, ArrayList<Point>> right,
      Point l, Point r) {
    ArrayList<Point> usedLs = new ArrayList<>(0);
    ArrayList<Point> usedRs = new ArrayList<>(0);
    //assume l and r are lowest. Test the cases
    loop:
    while (true) {
      //2) If the current l is not visible by the current r, discard edge and look for any new one
      if (!isVisible(l, left, r)) {
        usedLs.add(l);
        // get new L
        l = null;
        for (Point p : left.keySet()) {
          if (!usedLs.contains(p)) {
            if (l == null || p.y < l.y) l = p;
            else if (p.y == l.y && l.x < p.x) l = p;
          }
        }
        continue loop;
      }
      //3) If the current r is not visible by the current l, discard edge and look for any new one
      if (!isVisible(r, right, l)) {
        usedRs.add(r);
        // get new R
        r = null;
        for (Point p : right.keySet()) {
          if (!usedRs.contains(p)) {
            if (r == null || p.y < r.y) r = p;
            else if (p.y == r.y && r.x > p.x) r = p;
          }
        }
        continue loop;
      }
      //4) If there is an r that has an angle less than current r (relative to l), take that and discard previous
      Angle r_l = new Angle(r, l);
      loop4:
      for (Point p : right.keySet()) {
        if (!usedRs.contains(p)) {
          //check if this point has less degrees while being visible.
          Angle p_l = new Angle(p, l);
          if (p_l.lessThan(r_l) && isVisible(p, right, l)) {
            //however if this new one is basically the same angle as previous
            if (r_l.minus(p_l) < 0.25) {
              //but the distance to the new one is further, rather skip
              if (r.dist(l) < p.dist(l)) continue loop4;
            }
            //then we will only take it if it is closer
            usedRs.add(r);
            r = p;
            continue loop;
          }
        }
      }
      //5) If there is an l that has an angle more than current l (relative to r), take that and discard previous
      Angle l_r = new Angle(l, r);
      loop5:
      for (Point p : left.keySet()) {
        if (!usedLs.contains(p)) {
          //check if this point has less degrees while being visible.
          Angle p_r = new Angle(p, r);
          if (l_r.lessThan(p_r) && isVisible(p, left, r)) {
            //however if this new one is basically the same angle as previous
            if (p_r.minus(l_r) < 0.25) {
              //but the distance to the new one is further, rather skip
              if (l.dist(r) < p.dist(r)) continue loop5;
            }
            usedLs.add(l);
            l = p;
            continue loop;
          }
        }
      }
      return new Point[]{l, r};
    }
  }

  /**
   * Checks wether the specified corner (one of the corners of the given shape)
   * is visible by the given eye (outside of the shape)
   * @param corner the corner of the shape the eye wants to see
   * @param shape the array of points representing the ordered Corners of some polygon
   * @param eye the Point that we want to check whether the corner is visible to
   */
  private static boolean isVisible(Point corner, Point[] shape, Point eye) {
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
    ArrayList<Point> sorted = sort(sortMe);
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
  /* Different Parameter version of the above isVisible() */
  private static boolean isVisible(Point corner, HashMap<Point, ArrayList<Point>> shape, Point eye) {
    Point[] arr = new Point[shape.size()];
    int i = 0;
    for (Point p : shape.keySet()) {
      arr[i] = p;
      i++;
    }
    return isVisible(corner, arr, eye);
  }

  /**
   * Takes in two sets that are next to each other (left and right) and merges them. Steps:
   1 Create an edge at the bottom using the lowest(y) of the left and right subsets' points
   2 For each side repsectively, set up a list of candidates and sort the list based on
   * smallest angle difference with the bottom edge
   3 Each side selects a candidate from the sorted list and determines which
   * criterions the candidate meets:
   * C1: Criterion1) Angle less than 180
   * C2: Criterion2) Circumcircle with that candidate via triangle with left,
   * right and candidate points does not contain next candidate on that side
   *
   * Process followed based on critera met:
   * if C1 and C2
   *   Set this candidate as the Final candidate for this side
   * else if !C1
   *   Then no candidates will be chosen from this side (break candidate search loop)
   * else if C1 and !C2
   *   The line drawn from the endpoint (on the side of the candidate) to
   *   the candidate is deleted
   * end if
   *
   * Unless the candidate search loop is broken, this process is repeated until
   * a final candidate for this side is chosen, or all candidates have been exhuasted
   * (In which case, there is no final candidate)
   *
   * After both sides have confirmed their final candidates (if any), the
   * following process is followed:
   - If neither side has submitted a candidate, the merge is complete
   - If only one candidate is submited, we setup a connection between the opposite
   * side's bottom edge endpoint and the final candidate. the new LR lines endpoint
   * on the candidate's side becomes the candidate
   - If both sides have submitted a candidate, the previous if will be followed but
   * with the final candidate being the point not contained by the circumcircle drawn
   * via the bottom Edge's endpoints and the opposite side's final candidate
   *
   * @param left,right are the (so-far) delaunay triangulated sets of the left
   *        and right sides of the map, respectively.
   * @return the delaunized merged version of left and right
   */
  private static HashMap<Point, ArrayList<Point>> merge(
      HashMap<Point, ArrayList<Point>> left,
      HashMap<Point, ArrayList<Point>> right) {
    // create an edge at the bottom using the bottom of the left and right subsets
    HashMap<Point, ArrayList<Point>> merged = new HashMap<>();

    //------------------------------------------------------------------------
    //                            GET LR
    //------------------------------------------------------------------------
    // find bottom of each of left and right
    Point l = null;
    for (Point p : left.keySet()){
      merged.put(p, left.get(p));
      if (l == null || p.y < l.y || (p.y == l.y && p.x > l.x))
        l = p;
    }
    Point r = null;
    for (Point p : right.keySet()){
      merged.put(p, right.get(p));
      if (r == null || p.y < r.y || (p.y == r.y && p.x < p.y))
        r = p;
    }
    if (left.isEmpty() || right.isEmpty())
      return merged;


// check that r and l are in sight
    // ensure that there are no points under the line between r-l from l & r list
    Angle r_l = new Angle(r, l);
    Angle l_r = new Angle(l, r);
    Point[] getLR = getLR(left, right, l, r);
    l = getLR[0];
    r = getLR[1];

    // LOOP
    loop:
    while (true) {
      //------------------------------------------------------------------------
      //                     CANDIDATE SELECTION
      //------------------------------------------------------------------------

      //______________________________LEFT______________________________________
      // create sorted list of all points connected to l based on smallest angle
      HashMap<Point, Double> leftAngles = new HashMap<>();
      r_l = new Angle(r, l);//base line from l's perspective
      for (Point p : merged.get(l)) {
        Angle p_l = new Angle(p, l);//angle of possible candidate relative to l
        double diff = p_l.minus(r_l);
        leftAngles.put(p, diff);
      }
      ArrayList<Point> leftCandidates = sort(leftAngles);
      // search for left candidate
      Point leftFinalCandidate = null;
      leftCandidateSearch:
      for (int i = 0; i < leftCandidates.size(); i++) {
        if (leftFinalCandidate != null) break leftCandidateSearch;
        Point leftCandidate = leftCandidates.get(i);
        if (leftCandidate == r) continue leftCandidateSearch;
        //Now to get the next candidate
        Point nextCandidate = (i+1 < leftCandidates.size()) ? leftCandidates.get(i+1) : null;
        //first criterion: angle < 180
        double angle = leftAngles.get(leftCandidate);
        boolean firstCriterion = angle < 179.99;//180.0;
        //second criterion: circumcircle with that candidate may not contain next candidate on that side
        boolean secondCriterion = true;
        if (nextCandidate != null && angle != 180.0) {
          Point cc = getCc(l, r, leftCandidate);
          if (cc.dist(nextCandidate) < cc.dist(l)) {//if dist to next candidate is outside of radius
            secondCriterion = false;
          }
        }
        //act accordingly to criteria:
        if (firstCriterion && secondCriterion) {
          leftFinalCandidate = leftCandidate;//our final candidate for that side
        } else if (!firstCriterion) {
          break leftCandidateSearch; //Then no candidates will be chosen from that side
        } else if (firstCriterion && !secondCriterion) {
          // the line drawn from the endpoint (on the side of the candidate) to the candidate is deleted
          merged.get(leftCandidate).remove(l);
          merged.get(l).remove(leftCandidate);
        }
      }

      //______________________________RIGHT_____________________________________
      // create sorted list of all points connected to r based on smallest angle
      HashMap<Point, Double> rightAngles = new HashMap<>();
      l_r = new Angle(l, r);//base line from r's perspective
      for (Point p : merged.get(r)) {
        Angle p_r = new Angle(p, r);//angle of possible candidate relative to r
        double diff = l_r.minus(p_r);
        rightAngles.put(p, diff);
      }
      ArrayList<Point> rightCandidates = sort(rightAngles);

      // now we have a sorted list of all possible left and right candidates.
      Point rightFinalCandidate = null;
      rightCandidateSearch:
      for (int i = 0; i < rightCandidates.size(); i++) {
        if (rightFinalCandidate != null) break rightCandidateSearch;
        Point rightCandidate = rightCandidates.get(i);
        if (rightCandidate == l) continue rightCandidateSearch;
        Point nextCandidate = (i+1 < rightCandidates.size()) ? rightCandidates.get(i+1) : null;
        //first criterion: angle < 180
        double angle = rightAngles.get(rightCandidate);
        boolean firstCriterion = angle < 179.99;//180.0;
        //second criterion: circumcircle with that candidate may not contain next candidate on that side
        boolean secondCriterion = true;
        if (nextCandidate != null && angle != 180.0) {
          Point cc = getCc(l, r, rightCandidate);
          if (cc.dist(nextCandidate) < cc.dist(r)) {//if dist to next candidate is outside of radius
            secondCriterion = false;
          }
        }
        //act accordingly to criteria:
        if (firstCriterion && secondCriterion) {
          rightFinalCandidate = rightCandidate;//our final candidate for that side
        } else if (!firstCriterion) {
          break rightCandidateSearch; //Then no candidates will be chosen from that side
        } else if (firstCriterion && !secondCriterion) {
          //the line drawn from the endpoint (on the side of the candidate) to the candidate is deleted
          merged.get(rightCandidate).remove(r);
          merged.get(r).remove(rightCandidate);
        }
      }
      // both sides have finished selecting candidates

      //_______________________CANDIDATE PROCESSING_____________________________
      // If neither side has submitted a candidate
      if (leftFinalCandidate == null && rightFinalCandidate == null) {
        // the merge is complete
        if (!merged.get(l).contains(r)) merged.get(l).add(r);
        if (!merged.get(r).contains(l)) merged.get(r).add(l);
        break loop;
      }

      // only left candidate was chosen
      else if (leftFinalCandidate != null && rightFinalCandidate == null) {
        //check for colinearcy
        double lAngle = leftAngles.get(leftFinalCandidate);
        if (lAngle == 0) {
          //l = lfc
          l = leftFinalCandidate;
        } else {
          //the line that goes from l to final candidate becomes the new LR edge.
          if (!merged.get(l).contains(r)) merged.get(l).add(r);
          if (!merged.get(r).contains(l)) merged.get(r).add(l);
          l = leftFinalCandidate;
        }
      }

      //only right candidate was chosen
      else if (leftFinalCandidate == null && rightFinalCandidate != null) {
        //check for colinearcy
        double rAngle = rightAngles.get(rightFinalCandidate);
        if (rAngle == 0) {
          //r = rfc
          r = rightFinalCandidate;
        } else {
          //the line that goes from l to final candidate becomes the new LR edge.
          if (!merged.get(l).contains(r)) merged.get(l).add(r);
          if (!merged.get(r).contains(l)) merged.get(r).add(l);
          r = rightFinalCandidate;
        }
      }

      //If both sides have submitted a candidate:
      else if (leftFinalCandidate != null && rightFinalCandidate != null) {
        double lAngle = leftAngles.get(leftFinalCandidate);
        double rAngle = rightAngles.get(rightFinalCandidate);

        if (lAngle == 0 && rAngle == 0) {
          // THIS SHOULD NOT BE THE CASE:
          if (leftFinalCandidate == rightFinalCandidate) {
            System.out.println("leftFinalCandidate == rightFinalCandidate");
            System.exit(0);
          }
        }

        if (lAngle == 0) {
          l = leftFinalCandidate;
        } else if (rAngle == 0) {
          r = rightFinalCandidate;
        } else {
          // get the circum circle of both triangles formed with the candidates
          Point lcc = getCc(l, r, leftFinalCandidate);
          Point rcc = getCc(l, r, rightFinalCandidate);//we need only test once
          // the circumcircle that does not contain the other candidate will be the final candidate
          if (lcc.dist(l) > lcc.dist(rightFinalCandidate)){
            //then rightFinalCandidate is the successor
            if (!merged.get(l).contains(r)) merged.get(l).add(r);
            if (!merged.get(r).contains(l)) merged.get(r).add(l);
            r = rightFinalCandidate;
          } else {
            //then leftFinalCandidate is the successor
            if (!merged.get(l).contains(r)) merged.get(l).add(r);
            if (!merged.get(r).contains(l)) merged.get(r).add(l);
            l = leftFinalCandidate;
          }
        }
      }
    }
    return merged;
  }

  /**
   * Calculates the circumcentre of a collection of points making up a triangle.
   * @param a,b,c the 3 corners of a triangle
   * @return the coordinate of the circumcentre of the triangle a-b-c
   */
  private static Point getCc(Point a, Point b, Point c){
    double x1 = a.x, y1 = a.y, x2 = b.x, y2 = b.y, x3 = c.x, y3 = c.y;
    double y2_MINUS_y1 = (y2-y1==0) ? 0.001 : y2-y1;
    double y3_MINUS_y1 = (y3-y1==0) ? 0.001 : y3-y1;
    double x = 0.5*
      ( y3 - y2 + (x3*x3-x1*x1)/(y3_MINUS_y1) + (x1*x1-x2*x2)/(y2_MINUS_y1) )/
      ( (x1-x2)/(y2_MINUS_y1) + (x3-x1)/(y3_MINUS_y1) );
    double y = ((x1-x2)/(y2_MINUS_y1))*x + 0.5*(y1 + y2 + (x2*x2-x1*x1)/(y2_MINUS_y1));
    return new Point(x, y);
  }

  /**
   * Sorts a given mappigng of keys to doubles in order of the mapped values.
   * Pulld from Core.java @https://github.com/HeinrichWizardKreuser/KreuserCore
   *
   * @param toSort HashMap with key of any data type and value as Double
   * @return an arraylist of the keys in order sorted based on doubles they were
   *         originally mapped to in the hashmap
   */
	public static <T> ArrayList<T> sort(HashMap<T, Double> toSort) {
		HashMap<T, Double> sorted = toSort.entrySet().stream().sorted(comparingByValue())
			.collect(toMap(e -> e.getKey(), e -> e.getValue(), (e1, e2) -> e2, LinkedHashMap::new));
		ArrayList<T> toReturn = new ArrayList<>(0);
		for (T t : sorted.keySet()) toReturn.add(t);
		return toReturn;
	}
}
