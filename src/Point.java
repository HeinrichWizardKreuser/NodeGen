/*******************************************************************************
 * This is a coordinate in 2D or 3D space. Offers methods you'd think necessary
 * for classic location data types. Can also act as a vector for some vector
 * functions such as dot and cross product.
 *
 * The methods in the "INTERSECTION" section are not my own work, you can find
 * the methods written in C (originally, I converted them to java) at
 * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
 *
 * @author Heinrich Kreuser
 *
 * Date: 9 June 2019
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
public class Point implements Comparable<Point> {

	/* The x, y and z coordinates of this point */
	public double x = 0d;
	public double y = 0d;
	public double z = 0d;

	public static final Point origin = new Point(0, 0, 0);

	/*****************************************************************************
	 *                           CONSTRUCTORS
	 ****************************************************************************/
	 /* 2d constructor */
	public Point(double x, double y) {
		setTo(x, y);
	}

	/* 3d constructor */
	public Point(double x, double y, double z) {
		setTo(x, y, z);
	}

	/**
	 * Constructs a point from a given string
	 *
	 * @param s the string in the form [x, y] or [x, y, z] where x, y and z are
	 * 				double values describing the coordinates
	 */
	public Point(String s) {
		s.replaceAll("[[]", "");
		s.replaceAll("[]]", "");
		String[] sarr = s.split("[,]");
		int n = sarr.length;
		double[] darr = new double[n];
		for (int i = 0; i < n; i++) {
			darr[i] = Double.parseDouble(sarr[i]);
		}
		if (n == 2) {
			setTo(darr[0], darr[1]);
		} else if (n == 3) {
			setTo(darr[0], darr[1], darr[2]);
		} else {
			throw new IllegalArgumentException("Invalid dimension size " + n + " for "
			 + s + "!!");
		}
	}

	/* Constructor for array */
	public Point(double[] arr) {
	    if (arr.length == 2) {
	      setTo(arr[0], arr[1]);
	    } else if (arr.length == 3) {
	      setTo(arr[0], arr[1], arr[2]);
	    } else {
	      throw new IllegalArgumentException("Invalid array size " + arr.length + "!!");
	    }
	}

	/* Constructor for the average of two points */
	public Point(Point a, Point b) {
		this.x = (a.x + b.x)/2;
		this.y = (a.y + b.y)/2;
		this.z = (a.z + b.z)/2;
	}

	/* Constructer for the weighted average of two points */
	public Point(Point a, Point b, double d) {
		if (d < 0d && !epsilon(d, 0d)) {
			throw new IllegalArgumentException("d("+d+") must be between 0 and 1 inclusive!");
		}
		if (1d < d && !epsilon(d, 1d)) {
			throw new IllegalArgumentException("d("+d+") must be between 0 and 1 inclusive!");
		}
		double i = 1 - d;
		this.x = a.x*d + b.x*i;
		this.y = a.y*d + b.y*i;
		this.z = a.z*d + b.z*i;
	}

	/* Setter methods */
	public void setTo(Point p) {
		this.x = p.x;
		this.y = p.y;
		this.z = p.z;
	}
	public void setTo(double x, double y) {
		this.x = x;
		this.y = y;
	}
	public void setTo(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/*****************************************************************************
	 *                           XYZ MATH
	 ****************************************************************************/
	/* Distance methods*/
	public double dist(Point p) {
		return Math.sqrt( (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y) + (z-p.z)*(z-p.z) );
	}
	public double dist(double x, double y) {
		return dist(new Point(x, y));
	}
	public double dist(double x, double y, double z) {
		return dist(new Point(x, y, z));
	}

	/**
	 * Returns a point with distance d from this point and in direction degrees
	 *
	 * @param degrees the direction of the returned point relative to this point
	 * @param d the distance the returned point must be
	 * @return a point with distance d from this point and in direction degrees
	 */
	public Point directed(double degrees, double d) {
		double x = this.x + d * Math.cos(Math.toRadians(degrees));
		double y = this.y + d * Math.sin(Math.toRadians(degrees));
		return new Point(x, y);
	}

	/* Point arithmatic */
	public Point plus(Point a) {
		return new Point(x + a.x, y + a.y, z + a.z);
	}
	public Point plus(double x, double y) {
		return new Point(this.x + x, this.y + y);
	}

	public Point minus(Point a) {
		return new Point(x - a.x, y - a.y, z - a.z);
	}
	public Point minus(double x, double y) {
		return new Point(this.x - x, this.y - y);
	}

	/* Gets the area of the triangle with corners d, e, f */
	public static double area(Point d, Point e, Point f){
		return Math.abs(d.x*(e.y-f.y)+e.x*(f.y-d.y)+f.x*(d.y-e.y))/2.0;
	}

	/* Gets area of triangle in 3d using dot product method */
	public static double area3d(Point a, Point b, Point c) {
		return a.minus(b).cross(c.minus(b)).len()*0.5;
	}

	/*****************************************************************************
	*                           VECTOR MATH
	****************************************************************************/
 	// determinant
	public double det(Point p) {
		return x * p.y - y * p.x;
	}
	// cross product
	public Point cross(Point p) {
		return new Point(
			y*p.z - z*p.y,
			z*p.x - x*p.z,
			x*p.y - y*p.x
		);
	}

	// dot product
 	public double dot(Point p) {
		return x * p.x + y * p.y + z * p.z;
	}
	// multiplies vector by s (scalar)
	public Point scale(double s) {
		return new Point(s*x, s*y, s*z);
	}

	public double len() {
		return Math.sqrt(x*x + y*y + z*z);
	}
	/**
	 * Gets the projection point if this point on line ab.
	 * NOTE: This is not using point as vector, but a location.
	 */
	public Point proj(Point[] ab) {
		return proj(ab[0], ab[1]);
	}
	public Point proj(Point a, Point b) {
		Point u = b.minus(a);
		Point v = this.minus(a);
		Point proj = u.scale(u.dot(v)/u.dot(u));
		return proj.plus(a);
	}
	/* Gets the shortest distance to line ab */
	public double dist(Point[] ab) {
		return dist(proj(ab[0], ab[1]));
	}
	public double dist(Point a, Point b) {
		return dist(proj(a, b));
	}

	/*****************************************************************************
	 *                           EXTRAS
	 ****************************************************************************/
	@Override
	public String toString() {
		if (z == 0) {
			return "[" + x + ", " + y + "]";
		}
		return "[" + x + ", " + y + ", " + z + "]";
	}

	/**
	 * Converts this point to an array
	 *
	 * @return an array containing this location's x, y and z (if z is nonzero)
	 *				 coordinates in that order
	 */
	public double[] toArray() {
		if (z == 0) {
			return new double[]{x, y};
		}
		return new double[]{x, y, z};
	}

	/* Standard equals methods for equality comparison */
	public boolean equals(Point p) {
		return epsilon(this.x, p.x) && epsilon(this.y, p.y) && epsilon(this.z, p.z);
	}
	public boolean equals(double x, double y) {
		return epsilon(this.x, x) && epsilon(this.y, y);
	}
	public boolean equals(double x, double y, double z) {
		return epsilon(this.x, x) && epsilon(this.y, y) && epsilon(this.z, z);
	}

	/**
	 * Checks whether the given values are within one epsilon difference from
	 * each other. If they are, they are considered to be equal. This is for when
	 * number calculations have a slight error.
	 */
	public static boolean epsilon(double a, double b) {
		return Math.abs(a-b) < EPSILON;
	}
	/*
	 * If two values have < EPSILON diifference from each other, they are
	 * considered to be equal
	 */
	private static final double EPSILON = 1E-6;

 	/**
	 * Retrieves a copy of this point
	 */
	public Point copy() {
		return new Point(this.x, this.y, this.z);
	}

	/**
	 * Override for comparable API. A point is less than another if it's x value
	 * is less. If their x values are equal, then it is smaller if it's y value
	 * is less. Else, they are equal and we return 0.
	 * This idea for how to compare points originates from a divide-and-conquer
	 * approach to delaunay triangulation of a map/collection of points
	 */
	@Override
	public int compareTo(Point p) {
		if (epsilon(this.x, p.x)) {
			if (epsilon(this.y, p.y)) {
				return 0;
			}
			return this.y < p.y ? -1 : +1;
		}
		return this.x < p.x ? -1 : +1;
	}

	/*****************************************************************************
	 *                          UTILITIES
	 ****************************************************************************/
	/**
	 * Creates an array of points with x and y values from the given array
	 * @param coords the list of x and y coordinates to construct an array of
	 * 				points of
	 * @return an array of points of which each point has the x coordinate of every
	 *				 odd index, and y coordinate of every even index
	 */
	public static Point[] arrxy(double... coords) {
		int n = coords.length;
		Point[] arr = new Point[n/2];
		for (int i = 0; i < n; i += 2) {
			arr[i/2] = new Point(coords[i], coords[i+1]);
		}
		return arr;
	}

	/**
	 * Creates an array of points with x, y and z values from the given array
	 * @param coords the list of x, y and z coordinates to construct an array of
	 * 				points of
	 * @return an array of points of which each point has the x, y and z coordinate
	 *				 of every 1st, 2nd and 3rd index respectively
	 * Example: arrxyz(1,2,3,4,5,6) = points (1,2,3) and (4,5,6)
	 */
	public static Point[] arrxyz(double... coords) {
		int n = coords.length;
		Point[] arr = new Point[n/3];
		for (int i = 0; i < n; i += 3) {
			arr[i/3] = new Point(coords[i], coords[i+1], coords[i+2]);
		}
		return arr;
	}

	/*****************************************************************************
	 *                           INTERSECTION
	 * The following are not my own work. You can find out more about the authors
	 * and functions like these on the page I got this from:
	 * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
	 ****************************************************************************/
	/**
	 * Checks whether line segmen a-b and c-d are intersecting
	 * @return true if line segement ab and cd intersects
	 */
	public static boolean intersects(Point a, Point b, Point c, Point d) {
		// Find the four orientations needed for general and
		// special cases
		int o1 = orientation_i(a, b, c);
		int o2 = orientation_i(a, b, d);
		int o3 = orientation_i(c, d, a);
		int o4 = orientation_i(c, d, b);
		// General case
		if (o1 != o2 && o3 != o4) {
			return true;
		}
		// Special Cases
		// a, b and c are colinear and c lies on segment ab
		if (o1 == COLINEAR && onSegment(a, c, b)) {
			return true;
		}
		// a, b and d are colinear and d lies on segment ab
		if (o2 == COLINEAR && onSegment(a, d, b)) {
			return true;
		}
		// c, d and a are colinear and a lies on segment cd
		if (o3 == COLINEAR && onSegment(c, a, d)) {
			return true;
		}
		// c, d and b are colinear and b lies on segment cd
		if (o4 == COLINEAR && onSegment(c, b, d)) {
			return true;
		}
		return false; // Doesn't fall in any of the above cases
	}
	/** Point[] / line parameter version of the above */
	public static boolean intersects(Point[] ab, Point[] cd) {
		return intersects(ab[0], ab[1], cd[0], cd[1]);
	}

	/**
	 * To find orientation of ordered triplet (p, q, r).
	 * See https://www.geeksforgeeks.org/orientation-3-ordered-points/
	 * for details of below formula.
	 *
	 * @param p,q,r are the points to check colinearcy of
	 * @return a string description of the orientation of the points
	 */
	public static String orientation_s(Point p, Point q, Point r) {
		int o = orientation_i(p, q, r);
		switch (o) {
			case COLINEAR:
				return "colinear";
			case CLOCKWISE:
				return "clockwise";
			case COUNTERCLOCKWISE:
				return "counterclockwise";
			default:
				break;
		}
		return null;
	}

	public static int orientation_i(Point p, Point q, Point r) {
		double val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
		if (val == 0) {
			return COLINEAR; // colinear
		}
		return val > 0 ? CLOCKWISE : COUNTERCLOCKWISE;
	}

	public static final int COLINEAR = 0;
	public static final int CLOCKWISE = 1;
	public static final int COUNTERCLOCKWISE = 2;

	/**
	 * Given three colinear points p, q, r, the function checks if
	 * point q lies on line segment 'pr'
	 */
	public static boolean onSegment(Point p, Point q, Point r) {
		return
			q.x <= Math.max(p.x, r.x) && q.x >= Math.min(p.x, r.x) &&
			q.y <= Math.max(p.y, r.y) && q.y >= Math.min(p.y, r.y);
	}

	/**
	 * Checks if this point is inside the given polygon
	 *
	 * @param poly is a collection of points that represent the corners of a polygon
	 * @return true if this point is contained within the given polygon
	 */
	public boolean isInside(Point[] poly) {
		int n = poly.length;
		// There must be at least 3 vertices in polygon[]
		if (n < 3) return false;
		// Create a point for line segment from p to relative infinity
		double extremeX = poly[0].x;
		for (int i = 1; i < n; i++) {
			extremeX = Math.max(extremeX, poly[i].x);
		}
		Point extreme = new Point(extremeX+1, this.y);
		// Count intersections of the above line with sides of polygon
		int count = 0;
		// check with all lines that contain this point
		loop:
		for (int i = 0; i < n; i++) {
			Point a = poly[i];
			if (equals(a)) return true; // if this point is a corner of shape
			int j = i+1 == n ? 0 : i+1;
			Point b = poly[j];
			// atleast one of a or b must be above and under p respectively
			if (a.y < y && b.y < y) continue loop;
			if (a.y > y && b.y > y) continue loop;
			// Check if the line intersects:
			if (intersects(a, b, this, extreme)) count++;
		}
		// Return true if count is odd, false otherwise
		return count%2 == 1;
	}

}
