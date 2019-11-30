import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FileDialog;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.RenderingHints;
import java.awt.Toolkit;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.*;

import java.awt.geom.Arc2D;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Line2D;
import java.awt.geom.Rectangle2D;

import java.awt.image.BufferedImage;
import java.awt.image.DirectColorModel;
import java.awt.image.WritableRaster;

import java.io.File;
import java.io.IOException;

import java.net.MalformedURLException;
import java.net.URL;

import java.util.LinkedList;
import java.util.ArrayList;
import java.util.TreeSet;
import java.util.NoSuchElementException;
import javax.imageio.ImageIO;

import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

public final class StdDraw implements ActionListener, MouseListener, MouseMotionListener, KeyListener, MouseWheelListener {

    public static final Color BLACK = Color.BLACK;
    public static final Color BLUE = Color.BLUE;
    public static final Color CYAN = Color.CYAN;
    public static final Color DARK_GRAY = Color.DARK_GRAY;
    public static final Color GRAY = Color.GRAY;
    public static final Color GREEN  = Color.GREEN;
    public static final Color LIGHT_GRAY = Color.LIGHT_GRAY;
    public static final Color MAGENTA = Color.MAGENTA;
    public static final Color ORANGE = Color.ORANGE;
    public static final Color PINK = Color.PINK;
    public static final Color RED = Color.RED;
    public static final Color WHITE = Color.WHITE;
    public static final Color YELLOW = Color.YELLOW;
    public static final Color BOOK_BLUE = new Color(9, 90, 166);
    public static final Color BOOK_LIGHT_BLUE = new Color(103, 198, 243);
    public static final Color BOOK_RED = new Color(150, 35, 31);
    public static final Color PRINCETON_ORANGE = new Color(245, 128, 37);

//some personal edits
	//title
    public static void setTitle(String newTitle){title = newTitle;}
    public static String title = "Standard Draw";
    //customized colors
    public static final Color CRIMSON = new Color(220,20,60);
    public static final Color MAROON = new Color(128,0,0);
    public static final Color DEEPSKY_BLUE = new Color(0,191,255);
    public static final Color INDIGO = new Color(75,0,130);
    public static final Color LAVENDER = new Color(230,230,250);
    public static final Color AQUAMARINE = new Color(127,255,212);
    public static final Color OLIVE = new Color(128,128,0);
    public static final Color LIGHTSLATE_GRAY = new Color(119,136,153);
    public static final Color GOLD = new Color(255,215,0);
    public static final Color PURPLE = new Color(128,0,128);
    public static final Color SNOW = new Color(255,250,250);
    public static final Color LIGHT_LAVENDER = new Color(204,204,255);
    public static final Color ATOM = new Color(39,44,51);
    public static final Color LIGHT_ATOM = new Color(68,77,89);

    public static String colorName(Color color) {
      if (color.equals(BLACK)) {
        return "BLACK";
      } else if (color.equals(BLUE)) {
        return "BLUE";
      } else if (color.equals(CYAN)) {
        return "CYAN";
      } else if (color.equals(DARK_GRAY)) {
        return "DARK_GRAY";
      } else if (color.equals(GRAY)) {
        return "GRAY";
      } else if (color.equals(GREEN)) {
        return "GREEN";
      } else if (color.equals(LIGHT_GRAY)) {
        return "LIGHT_GRAY";
      } else if (color.equals(MAGENTA)) {
        return "MAGENTA";
      } else if (color.equals(ORANGE)) {
        return "ORANGE";
      } else if (color.equals(PINK)) {
        return "PINK";
      } else if (color.equals(RED)) {
        return "RED";
      } else if (color.equals(WHITE)) {
        return "WHITE";
      } else if (color.equals(YELLOW)) {
        return "YELLOW";
      } else if (color.equals(BOOK_BLUE)) {
        return "BOOK_BLUE";
      } else if (color.equals(BOOK_LIGHT_BLUE)) {
        return "BOOK_LIGHT_BLUE";
      } else if (color.equals(BOOK_RED)) {
        return "BOOK_RED";
      } else if (color.equals(PRINCETON_ORANGE)) {
        return "PRINCETON_ORANGE";
      } else if (color.equals(CRIMSON)) {
        return "CRIMSON";
      } else if (color.equals(MAROON)) {
        return "MAROON";
      } else if (color.equals(DEEPSKY_BLUE)) {
        return "DEEPSKY_BLUE";
      } else if (color.equals(INDIGO)) {
        return "INDIGO";
      } else if (color.equals(LAVENDER)) {
        return "LAVENDER";
      } else if (color.equals(AQUAMARINE)) {
        return "AQUAMARINE";
      } else if (color.equals(OLIVE)) {
        return "OLIVE";
      } else if (color.equals(LIGHTSLATE_GRAY)) {
        return "LIGHTSLATE_GRAY";
      } else if (color.equals(GOLD)) {
        return "GOLD";
      } else if (color.equals(PURPLE)) {
        return "PURPLE";
      } else if (color.equals(SNOW)) {
        return "SNOW";
      } else if (color.equals(LIGHT_LAVENDER)) {
        return "LIGHT_LAVENDER";
      } else if (color.equals(ATOM)) {
        return "ATOM";
      } else if (color.equals(LIGHT_ATOM)) {
        return "LIGHT_ATOM";
      } else {
        // get rgb
        return colorString(color);
      }
    }

    public static String colorString(Color c) {
      return "["+c.getRed()+", "+c.getGreen()+", "+c.getBlue()+"]";
    }

    public static Color randomColor(){
        return new Color( (int)(Math.random()*256), (int)(Math.random()*256), (int)(Math.random()*256) );
    }
    public static Color hue(Color color, double percent){
        double p = percent/100.00;
        return new Color((int)(color.getRed()*p), (int)(color.getGreen()*p), (int)(color.getBlue()*p));
    }

    private static final Color DEFAULT_PEN_COLOR   = BLACK;
    private static final Color DEFAULT_CLEAR_COLOR = WHITE;

    private static Color penColor;// current pen color

    // default canvas size is DEFAULT_SIZE-by-DEFAULT_SIZE
    private static final int DEFAULT_SIZE = 512;
    private static int width  = DEFAULT_SIZE;
    private static int height = DEFAULT_SIZE;

    private static final double DEFAULT_PEN_RADIUS = 0.002;// default pen radius
    private static double penRadius;// current pen radius
    private static boolean defer = false;// show we draw immediately or wait until next show?

    // boundary of drawing canvas, 0% border
    // private static final double BORDER = 0.05;
    private static final double BORDER = 0.00;
    private static final double DEFAULT_XMIN = 0.0;
    private static final double DEFAULT_XMAX = 1.0;
    private static final double DEFAULT_YMIN = 0.0;
    private static final double DEFAULT_YMAX = 1.0;
    private static double xmin, ymin, xmax, ymax;

    // for synchronization
    private static Object mouseLock = new Object();
    private static Object keyLock = new Object();

    // default font
    private static final Font DEFAULT_FONT = new Font("SansSerif", Font.PLAIN, 16);

    // current font
    private static Font font;

    // double buffered graphics
    private static BufferedImage offscreenImage, onscreenImage;
    private static Graphics2D offscreen, onscreen;

    // singleton for callbacks: avoids generation of extra .class files
    private static StdDraw std = new StdDraw();

    // the frame for drawing to the screen
    public static JFrame frame;

    // mouse state
    private static boolean isMousePressed = false;
    private static double mouseX = 0;
    private static double mouseY = 0;

    // queue of typed key characters
    private static LinkedList<Character> keysTyped = new LinkedList<Character>();

    // set of key codes currently pressed down
    private static TreeSet<Integer> keysDown = new TreeSet<Integer>();

    // singleton pattern: client can't instantiate
    private StdDraw() { }

    // static initializer
/*    static {
        init();
    }
*/
    /**
     * Sets the canvas (drawing area) to be 512-by-512 pixels.
     * This also erases the current drawing and resets the coordinate system,
     * pen radius, pen color, and font back to their default values.
     * Ordinarly, this method is called once, at the very beginning
     * of a program.
     */
    public static void setCanvasSize() {
        setCanvasSize(DEFAULT_SIZE, DEFAULT_SIZE);
    }

    /**
     * Sets the canvas (drawing area) to be <em>width</em>-by-<em>height</em> pixels.
     * This also erases the current drawing and resets the coordinate system,
     * pen radius, pen color, and font back to their default values.
     * Ordinarly, this method is called once, at the very beginning
     * of a program.
     *
     * @param  canvasWidth the width as a number of pixels
     * @param  canvasHeight the height as a number of pixels
     * @throws IllegalArgumentException unless both {@code canvasWidth} and
     *         {@code canvasHeight} are positive
     */
    public static void setCanvasSize(int canvasWidth, int canvasHeight) {
        if (canvasWidth <= 0 || canvasHeight <= 0)
            throw new IllegalArgumentException("width and height must be positive");
        width = canvasWidth;
        height = canvasHeight;
        init();
    }

    // init
    private static void init() {
        if (frame != null) frame.setVisible(false);
        frame = new JFrame();
        offscreenImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        onscreenImage  = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        offscreen = offscreenImage.createGraphics();
        onscreen  = onscreenImage.createGraphics();
        setXscale();
        setYscale();
        offscreen.setColor(DEFAULT_CLEAR_COLOR);
        offscreen.fillRect(0, 0, width, height);
        setPenColor();
        setPenRadius();
        setFont();
        clear();

        // add antialiasing
        RenderingHints hints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
                                                  RenderingHints.VALUE_ANTIALIAS_ON);
        hints.put(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY);
        offscreen.addRenderingHints(hints);

        // frame stuff
        ImageIcon icon = new ImageIcon(onscreenImage);
        JLabel draw = new JLabel(icon);

        draw.addMouseListener(std);
        draw.addMouseMotionListener(std);

        draw.addMouseWheelListener(std);//---------------PERSONAL ADD

        frame.setContentPane(draw);
        frame.addKeyListener(std);    // JLabel cannot get keyboard focus
        frame.setResizable(false);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);            // closes all windows
        // frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);      // closes only current window
        frame.setTitle(title);
        frame.setJMenuBar(createMenuBar());
        frame.pack();
        frame.requestFocusInWindow();
        frame.setVisible(true);
    }

    // create the menu bar (changed to private)
    private static JMenuBar createMenuBar() {
        JMenuBar menuBar = new JMenuBar();
        JMenu menu = new JMenu("File");
        menuBar.add(menu);
        JMenuItem menuItem1 = new JMenuItem(" Save...   ");
        menuItem1.addActionListener(std);
        menuItem1.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S,
                                Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()));
        menu.add(menuItem1);
        return menuBar;
    }


   /***************************************************************************
    *  User and screen coordinate systems.
    ***************************************************************************/

    /**
     * Sets the <em>x</em>-scale to be the default (between 0.0 and 1.0).
     */
    public static void setXscale() {
        setXscale(DEFAULT_XMIN, DEFAULT_XMAX);
    }

    /**
     * Sets the <em>y</em>-scale to be the default (between 0.0 and 1.0).
     */
    public static void setYscale() {
        setYscale(DEFAULT_YMIN, DEFAULT_YMAX);
    }

    /**
     * Sets the <em>x</em>-scale and <em>y</em>-scale to be the default
     * (between 0.0 and 1.0).
     */
    public static void setScale() {
        setXscale();
        setYscale();
    }

    /**
     * Sets the <em>x</em>-scale to the specified range.
     *
     * @param  min the minimum value of the <em>x</em>-scale
     * @param  max the maximum value of the <em>x</em>-scale
     * @throws IllegalArgumentException if {@code (max == min)}
     */
    public static void setXscale(double min, double max) {
        double size = max - min;
        if (size == 0.0) throw new IllegalArgumentException("the min and max are the same");
        synchronized (mouseLock) {
            xmin = min - BORDER * size;
            xmax = max + BORDER * size;
        }
    }

    /**
     * Sets the <em>y</em>-scale to the specified range.
     *
     * @param  min the minimum value of the <em>y</em>-scale
     * @param  max the maximum value of the <em>y</em>-scale
     * @throws IllegalArgumentException if {@code (max == min)}
     */
    public static void setYscale(double min, double max) {
        double size = max - min;
        if (size == 0.0) throw new IllegalArgumentException("the min and max are the same");
        synchronized (mouseLock) {
            ymin = min - BORDER * size;
            ymax = max + BORDER * size;
        }
    }

    /**
     * Sets both the <em>x</em>-scale and <em>y</em>-scale to the (same) specified range.
     *
     * @param  min the minimum value of the <em>x</em>- and <em>y</em>-scales
     * @param  max the maximum value of the <em>x</em>- and <em>y</em>-scales
     * @throws IllegalArgumentException if {@code (max == min)}
     */
    public static void setScale(double min, double max) {
        double size = max - min;
        if (size == 0.0) throw new IllegalArgumentException("the min and max are the same");
        synchronized (mouseLock) {
            xmin = min - BORDER * size;
            xmax = max + BORDER * size;
            ymin = min - BORDER * size;
            ymax = max + BORDER * size;
        }
    }

    // helper functions that scale from user coordinates to screen coordinates and back
    public static double  scaleX(double x) { return width  * (x - xmin) / (xmax - xmin); }
    public static double  scaleY(double y) { return height * (ymax - y) / (ymax - ymin); }
    public static double factorX(double w) { return w * width  / Math.abs(xmax - xmin);  }
    public static double factorY(double h) { return h * height / Math.abs(ymax - ymin);  }
    public static double   userX(double x) { return xmin + x * (xmax - xmin) / width;    }
    public static double   userY(double y) { return ymax - y * (ymax - ymin) / height;   }


    /**
     * Clears the screen to the default color (white).
     */
    public static void clear() {
        clear(DEFAULT_CLEAR_COLOR);
    }

    /**
     * Clears the screen to the specified color.
     *
     * @param color the color to make the background
     */
    public static void clear(Color color) {
        offscreen.setColor(color);
        offscreen.fillRect(0, 0, width, height);
        offscreen.setColor(penColor);
        draw();
    }

    /**
     * Returns the current pen radius.
     *
     * @return the current value of the pen radius
     */
    public static double getPenRadius() {
        return penRadius;
    }

    /**
     * Sets the pen size to the default size (0.002).
     * The pen is circular, so that lines have rounded ends, and when you set the
     * pen radius and draw a point, you get a circle of the specified radius.
     * The pen radius is not affected by coordinate scaling.
     */
    public static void setPenRadius() {
        setPenRadius(DEFAULT_PEN_RADIUS);
    }

    /**
     * Sets the radius of the pen to the specified size.
     * The pen is circular, so that lines have rounded ends, and when you set the
     * pen radius and draw a point, you get a circle of the specified radius.
     * The pen radius is not affected by coordinate scaling.
     *
     * @param  radius the radius of the pen
     * @throws IllegalArgumentException if {@code radius} is negative
     */
    public static void setPenRadius(double radius) {
        if (!(radius >= 0)) throw new IllegalArgumentException("pen radius must be nonnegative");
        penRadius = radius;
        float scaledPenRadius = (float) (radius * DEFAULT_SIZE);
        BasicStroke stroke = new BasicStroke(scaledPenRadius, BasicStroke.CAP_ROUND, BasicStroke.JOIN_ROUND);
        // BasicStroke stroke = new BasicStroke(scaledPenRadius);
        offscreen.setStroke(stroke);
    }

    /**
     * Returns the current pen color.
     *
     * @return the current pen color
     */
    public static Color getPenColor() {
        return penColor;
    }

    /**
     * Set the pen color to the default color (black).
     */
    public static void setPenColor() {
        setPenColor(DEFAULT_PEN_COLOR);
    }

    /**
     * Sets the pen color to the specified color.
     * <p>
     * The predefined pen colors are
     * {@code StdDraw.BLACK}, {@code StdDraw.BLUE}, {@code StdDraw.CYAN},
     * {@code StdDraw.DARK_GRAY}, {@code StdDraw.GRAY}, {@code StdDraw.GREEN},
     * {@code StdDraw.LIGHT_GRAY}, {@code StdDraw.MAGENTA}, {@code StdDraw.ORANGE},
     * {@code StdDraw.PINK}, {@code StdDraw.RED}, {@code StdDraw.WHITE}, and
     * {@code StdDraw.YELLOW}.
     *
     * @param color the color to make the pen
     */
    public static void setPenColor(Color color) {
        if (color == null) throw new IllegalArgumentException();
        penColor = color;
        offscreen.setColor(penColor);
    }
    public static void setPencolor(Color color) { setPenColor(color); }

    /**
     * Sets the pen color to the specified RGB color.
     *
     * @param  red the amount of red (between 0 and 255)
     * @param  green the amount of green (between 0 and 255)
     * @param  blue the amount of blue (between 0 and 255)
     * @throws IllegalArgumentException if {@code red}, {@code green},
     *         or {@code blue} is outside its prescribed range
     */
    public static void setPenColor(int red, int green, int blue) {
        if (red   < 0 || red   >= 256) throw new IllegalArgumentException("amount of red must be between 0 and 255");
        if (green < 0 || green >= 256) throw new IllegalArgumentException("amount of green must be between 0 and 255");
        if (blue  < 0 || blue  >= 256) throw new IllegalArgumentException("amount of blue must be between 0 and 255");
        setPenColor(new Color(red, green, blue));
    }

    /**
     * Returns the current font.
     *
     * @return the current font
     */
    public static Font getFont() {
        return font;
    }

    /**
     * Sets the font to the default font (sans serif, 16 point).
     */
    public static void setFont() {
        setFont(DEFAULT_FONT);
    }

    /**
     * Sets the font to the specified value.
     *
     * @param font the font
     */
    public static void setFont(Font font) {
        if (font == null) throw new IllegalArgumentException();
        StdDraw.font = font;
    }


   /***************************************************************************
    *  Drawing geometric shapes.
    ***************************************************************************/

    /**
     * Draws a line segment between (<em>x</em><sub>0</sub>, <em>y</em><sub>0</sub>) and
     * (<em>x</em><sub>1</sub>, <em>y</em><sub>1</sub>).
     *
     * @param  x0 the <em>x</em>-coordinate of one endpoint
     * @param  y0 the <em>y</em>-coordinate of one endpoint
     * @param  x1 the <em>x</em>-coordinate of the other endpoint
     * @param  y1 the <em>y</em>-coordinate of the other endpoint
     */
    public static void line(double x0, double y0, double x1, double y1) {
        offscreen.draw(new Line2D.Double(scaleX(x0), scaleY(y0), scaleX(x1), scaleY(y1)));
        draw();
    }
    public static void line(double x0, double y0, double x1, double y1, Color color){
        setPencolor(color);
        line(x0, y0, x1, y1);
    }
    public static void line(Point p1, Point p2){
        line(p1.x, p1.y, p2.x, p2.y);
    }
    public static void line(Point p1, Point p2, Color color){
        line(p1.x, p1.y, p2.x, p2.y, color);
    }
    /**
     * Draws one pixel at (<em>x</em>, <em>y</em>).
     * This method is private because pixels depend on the display.
     * To achieve the same effect, set the pen radius to 0 and call {@code point()}.
     *
     * @param  x the <em>x</em>-coordinate of the pixel
     * @param  y the <em>y</em>-coordinate of the pixel
     */
    private static void pixel(double x, double y) {
        offscreen.fillRect((int) Math.round(scaleX(x)), (int) Math.round(scaleY(y)), 1, 1);
    }

    /**
     * Draws a point centered at (<em>x</em>, <em>y</em>).
     * The point is a filled circle whose radius is equal to the pen radius.
     * To draw a single-pixel point, first set the pen radius to 0.
     *
     * @param x the <em>x</em>-coordinate of the point
     * @param y the <em>y</em>-coordinate of the point
     */
    public static void point(double x, double y) {
        double xs = scaleX(x);
        double ys = scaleY(y);
        double r = penRadius;
        float scaledPenRadius = (float) (r * DEFAULT_SIZE);

        // double ws = factorX(2*r);
        // double hs = factorY(2*r);
        // if (ws <= 1 && hs <= 1) pixel(x, y);
        if (scaledPenRadius <= 1) pixel(x, y);
        else offscreen.fill(new Ellipse2D.Double(xs - scaledPenRadius/2, ys - scaledPenRadius/2,
                                                 scaledPenRadius, scaledPenRadius));
        draw();
    }
    public static void point(Point p, Color color){
        setPencolor(color);
        point(p.x, p.y);
    }
    public static void point(Point p){
        point(p.x, p.y);
    }

    /**
     * Draws a circle of the specified radius, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the circle
     * @param  y the <em>y</em>-coordinate of the center of the circle
     * @param  radius the radius of the circle
     * @throws IllegalArgumentException if {@code radius} is negative
     */
    public static void circle(double x, double y, double radius) {
        if (!(radius >= 0)) throw new IllegalArgumentException("radius must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*radius);
        double hs = factorY(2*radius);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.draw(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }
    public static void circle(Point p, double radius) {
        circle(p.x, p.y, radius);
    }
    public static void circle(Point p, double radius, Color color) {
        setPenColor(color);
        circle(p.x, p.y, radius);
    }

    public static void circle(double x, double y, double radius, double penRadiusMultiplier) {
        double OGpenRadius = getPenRadius();
        setPenRadius(DEFAULT_PEN_RADIUS*penRadiusMultiplier);
        circle(x, y, radius);
        setPenRadius(OGpenRadius);
    }

    /**
     * Draws a filled circle of the specified radius, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the circle
     * @param  y the <em>y</em>-coordinate of the center of the circle
     * @param  radius the radius of the circle
     * @throws IllegalArgumentException if {@code radius} is negative
     */
    public static void filledCircle(double x, double y, double radius) {
        if (!(radius >= 0)) throw new IllegalArgumentException("radius must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*radius);
        double hs = factorY(2*radius);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.fill(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }

    public static void filledCircle(Point p, double radius) {
      filledCircle(p.x, p.y, radius);
    }

    public static void filledCircle(double x, double y, double radius, Color color) {
      setPenColor(color);
      filledCircle(x, y, radius);
    }

    public static void filledCircle(Point p, double radius, Color color) {
      setPenColor(color);
      filledCircle(p.x, p.y, radius);
    }


    /**
     * Draws an ellipse with the specified semimajor and semiminor axes,
     * centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the ellipse
     * @param  y the <em>y</em>-coordinate of the center of the ellipse
     * @param  semiMajorAxis is the semimajor axis of the ellipse
     * @param  semiMinorAxis is the semiminor axis of the ellipse
     * @throws IllegalArgumentException if either {@code semiMajorAxis}
     *         or {@code semiMinorAxis} is negative
     */
    public static void ellipse(double x, double y, double semiMajorAxis, double semiMinorAxis) {
        if (!(semiMajorAxis >= 0)) throw new IllegalArgumentException("ellipse semimajor axis must be nonnegative");
        if (!(semiMinorAxis >= 0)) throw new IllegalArgumentException("ellipse semiminor axis must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*semiMajorAxis);
        double hs = factorY(2*semiMinorAxis);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.draw(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }

    /**
     * Draws an ellipse with the specified semimajor and semiminor axes,
     * centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the ellipse
     * @param  y the <em>y</em>-coordinate of the center of the ellipse
     * @param  semiMajorAxis is the semimajor axis of the ellipse
     * @param  semiMinorAxis is the semiminor axis of the ellipse
     * @throws IllegalArgumentException if either {@code semiMajorAxis}
     *         or {@code semiMinorAxis} is negative
     */
    public static void filledEllipse(double x, double y, double semiMajorAxis, double semiMinorAxis) {
        if (!(semiMajorAxis >= 0)) throw new IllegalArgumentException("ellipse semimajor axis must be nonnegative");
        if (!(semiMinorAxis >= 0)) throw new IllegalArgumentException("ellipse semiminor axis must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*semiMajorAxis);
        double hs = factorY(2*semiMinorAxis);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.fill(new Ellipse2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }


    /**
     * Draws a circular arc of the specified radius,
     * centered at (<em>x</em>, <em>y</em>), from angle1 to angle2 (in degrees).
     *
     * @param  x the <em>x</em>-coordinate of the center of the circle
     * @param  y the <em>y</em>-coordinate of the center of the circle
     * @param  radius the radius of the circle
     * @param  angle1 the starting angle. 0 would mean an arc beginning at 3 o'clock.
     * @param  angle2 the angle at the end of the arc. For example, if
     *         you want a 90 degree arc, then angle2 should be angle1 + 90.
     * @throws IllegalArgumentException if {@code radius} is negative
     */
    public static void arc(double x, double y, double radius, double angle1, double angle2) {
        if (radius < 0) throw new IllegalArgumentException("arc radius must be nonnegative");
        while (angle2 < angle1) angle2 += 360;
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*radius);
        double hs = factorY(2*radius);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.draw(new Arc2D.Double(xs - ws/2, ys - hs/2, ws, hs, angle1, angle2 - angle1, Arc2D.OPEN));
        draw();
    }

    /**
     * Draws a square of side length 2r, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the square
     * @param  y the <em>y</em>-coordinate of the center of the square
     * @param  halfLength one half the length of any side of the square
     * @throws IllegalArgumentException if {@code halfLength} is negative
     */
    public static void square(double x, double y, double halfLength) {
        if (!(halfLength >= 0)) throw new IllegalArgumentException("half length must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*halfLength);
        double hs = factorY(2*halfLength);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.draw(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }
    public static void square(double x, double y, double halfLength, Color color){
        setPencolor(color);
        square(x, y, halfLength);
    }

    /**
     * Draws a filled square of the specified size, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the square
     * @param  y the <em>y</em>-coordinate of the center of the square
     * @param  halfLength one half the length of any side of the square
     * @throws IllegalArgumentException if {@code halfLength} is negative
     */
    public static void filledSquare(double x, double y, double halfLength) {
        if (!(halfLength >= 0)) throw new IllegalArgumentException("half length must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*halfLength);
        double hs = factorY(2*halfLength);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.fill(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }
    public static void filledSquare(double x, double y, double halfLength, Color color){
        setPencolor(color);
        filledSquare(x, y, halfLength);
    }
    public static void filledSquare(Point p, double halfLength){
        filledSquare(p.x, p.y, halfLength);
    }
    public static void filledSquare(Point p, double halfLength, Color color){
        setPencolor(color);
        filledSquare(p.x, p.y, halfLength);
    }

    /**
     * Draws a rectangle of the specified size, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the rectangle
     * @param  y the <em>y</em>-coordinate of the center of the rectangle
     * @param  halfWidth one half the width of the rectangle
     * @param  halfHeight one half the height of the rectangle
     * @throws IllegalArgumentException if either {@code halfWidth} or {@code halfHeight} is negative
     */
    public static void rectangle(double x, double y, double halfWidth, double halfHeight) {
        if (!(halfWidth  >= 0)) throw new IllegalArgumentException("half width must be nonnegative");
        if (!(halfHeight >= 0)) throw new IllegalArgumentException("half height must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*halfWidth);
        double hs = factorY(2*halfHeight);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.draw(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }
    public static void rectangle(double x, double y, double halfWidth, double halfHeight, Color color){
        setPencolor(color);
        rectangle(x, y, halfWidth, halfHeight);
    }

    /**
     * Draws a filled rectangle of the specified size, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the center of the rectangle
     * @param  y the <em>y</em>-coordinate of the center of the rectangle
     * @param  halfWidth one half the width of the rectangle
     * @param  halfHeight one half the height of the rectangle
     * @throws IllegalArgumentException if either {@code halfWidth} or {@code halfHeight} is negative
     */
    public static void filledRectangle(double x, double y, double halfWidth, double halfHeight) {
        if (!(halfWidth  >= 0)) throw new IllegalArgumentException("half width must be nonnegative");
        if (!(halfHeight >= 0)) throw new IllegalArgumentException("half height must be nonnegative");
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(2*halfWidth);
        double hs = factorY(2*halfHeight);
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else offscreen.fill(new Rectangle2D.Double(xs - ws/2, ys - hs/2, ws, hs));
        draw();
    }
    public static void filledRectangle(double x, double y, double halfWidth, double halfHeight, Color color){
        setPencolor(color);
        filledRectangle(x, y, halfWidth, halfHeight);
    }


    /**
     * Draws a polygon with the vertices
     * (<em>x</em><sub>0</sub>, <em>y</em><sub>0</sub>),
     * (<em>x</em><sub>1</sub>, <em>y</em><sub>1</sub>), ...,
     * (<em>x</em><sub><em>n</em>–1</sub>, <em>y</em><sub><em>n</em>–1</sub>).
     *
     * @param  x an array of all the <em>x</em>-coordinates of the polygon
     * @param  y an array of all the <em>y</em>-coordinates of the polygon
     * @throws IllegalArgumentException unless {@code x[]} and {@code y[]}
     *         are of the same length
     */
    public static void polygon(double[] x, double[] y) {
        if (x == null) throw new IllegalArgumentException("x-coordinate array is null");
        if (y == null) throw new IllegalArgumentException("y-coordinate array is null");
        int n1 = x.length;
        int n2 = y.length;
        if (n1 != n2) throw new IllegalArgumentException("arrays must be of the same length");
        int n = n1;
        if (n == 0) return;

        GeneralPath path = new GeneralPath();
        path.moveTo((float) scaleX(x[0]), (float) scaleY(y[0]));
        for (int i = 0; i < n; i++)
            path.lineTo((float) scaleX(x[i]), (float) scaleY(y[i]));
        path.closePath();
        offscreen.draw(path);
        draw();
    }

    public static void polygon(double[] x, double[] y, Color newColor){
        setPencolor(newColor);
        polygon(x, y);
    }

    public static void polygon(Point[] poly){
    //convert to collection of xs and ys and pass as polygon to be drawn
        double xs[] = new double[poly.length];
        double ys[] = new double[poly.length];
        for (int i = 0; i < poly.length; i++){
            xs[i] = poly[i].x;
            ys[i] = poly[i].y;
        }
        polygon(xs, ys);
    }

    public static void polygon(Point[] poly, Color newColor){
        setPencolor(newColor);
        polygon(poly);
    }

    public static void polygon(ArrayList<Point> poly, Color newColor){
        setPencolor(newColor);
    //convert to collection of xs and ys and pass as polygon to be drawn
        double xs[] = new double[poly.size()];
        double ys[] = new double[poly.size()];
        for (int i = 0; i < poly.size(); i++){
            xs[i] = poly.get(i).x;
            ys[i] = poly.get(i).y;
        }
        polygon(xs, ys);
    }


    /**
     * Draws a polygon with the vertices
     * (<em>x</em><sub>0</sub>, <em>y</em><sub>0</sub>),
     * (<em>x</em><sub>1</sub>, <em>y</em><sub>1</sub>), ...,
     * (<em>x</em><sub><em>n</em>–1</sub>, <em>y</em><sub><em>n</em>–1</sub>).
     *
     * @param  x an array of all the <em>x</em>-coordinates of the polygon
     * @param  y an array of all the <em>y</em>-coordinates of the polygon
     * @throws IllegalArgumentException unless {@code x[]} and {@code y[]}
     *         are of the same length
     */
    public static void filledPolygon(double[] x, double[] y) {
        if (x == null) throw new IllegalArgumentException("x-coordinate array is null");
        if (y == null) throw new IllegalArgumentException("y-coordinate array is null");
        int n1 = x.length;
        int n2 = y.length;
        if (n1 != n2) throw new IllegalArgumentException("arrays must be of the same length");
        int n = n1;
        if (n == 0) return;

        GeneralPath path = new GeneralPath();
        path.moveTo((float) scaleX(x[0]), (float) scaleY(y[0]));
        for (int i = 0; i < n; i++)
            path.lineTo((float) scaleX(x[i]), (float) scaleY(y[i]));
        path.closePath();
        offscreen.fill(path);
        draw();
    }

    public static void filledPolygon(double[] x, double[] y, Color newColor){
        setPencolor(newColor);
        filledPolygon(x, y);
    }

    public static void filledPolygon(Point[] poly, Color newColor){
        setPencolor(newColor);
    //convert to collection of xs and ys and pass as polygon to be drawn
        double xs[] = new double[poly.length];
        double ys[] = new double[poly.length];
        for (int i = 0; i < poly.length; i++){
            xs[i] = poly[i].x;
            ys[i] = poly[i].y;
        }
        filledPolygon(xs, ys);
    }

    public static void filledPolygon(ArrayList<Point> poly, Color newColor){
        setPencolor(newColor);
    //convert to collection of xs and ys and pass as polygon to be drawn
        double xs[] = new double[poly.size()];
        double ys[] = new double[poly.size()];
        for (int i = 0; i < poly.size(); i++){
            xs[i] = poly.get(i).x;
            ys[i] = poly.get(i).y;
        }
        filledPolygon(xs, ys);
    }

   /***************************************************************************
    *  Drawing images.
    ***************************************************************************/
    // get an image from the given filename
    private static Image getImage(String filename) {
        if (filename == null) throw new IllegalArgumentException();

        // to read from file
        ImageIcon icon = new ImageIcon(filename);

        // try to read from URL
        if ((icon == null) || (icon.getImageLoadStatus() != MediaTracker.COMPLETE)) {
            try {
                URL url = new URL(filename);
                icon = new ImageIcon(url);
            }
            catch (MalformedURLException e) {
                /* not a url */
            }
        }

        // in case file is inside a .jar (classpath relative to StdDraw)
        if ((icon == null) || (icon.getImageLoadStatus() != MediaTracker.COMPLETE)) {
            URL url = StdDraw.class.getResource(filename);
            if (url != null)
                icon = new ImageIcon(url);
        }

        // in case file is inside a .jar (classpath relative to root of jar)
        if ((icon == null) || (icon.getImageLoadStatus() != MediaTracker.COMPLETE)) {
            URL url = StdDraw.class.getResource("/" + filename);
            if (url == null) throw new IllegalArgumentException("image " + filename + " not found");
            icon = new ImageIcon(url);
        }

        return icon.getImage();
    }

   /***************************************************************************
    * [Summer 2016] Should we update to use ImageIO instead of ImageIcon()?
    *               Seems to have some issues loading images on some systems
    *               and slows things down on other systems.
    *               especially if you don't call ImageIO.setUseCache(false)
    *               One advantage is that it returns a BufferedImage.
    ***************************************************************************/
/*
    private static BufferedImage getImage(String filename) {
        if (filename == null) throw new IllegalArgumentException();

        // from a file or URL
        try {
            URL url = new URL(filename);
            BufferedImage image = ImageIO.read(url);
            return image;
        }
        catch (IOException e) {
            // ignore
        }

        // in case file is inside a .jar (classpath relative to StdDraw)
        try {
            URL url = StdDraw.class.getResource(filename);
            BufferedImage image = ImageIO.read(url);
            return image;
        }
        catch (IOException e) {
            // ignore
        }

        // in case file is inside a .jar (classpath relative to root of jar)
        try {
            URL url = StdDraw.class.getResource("/" + filename);
            BufferedImage image = ImageIO.read(url);
            return image;
        }
        catch (IOException e) {
            // ignore
        }
        throw new IllegalArgumentException("image " + filename + " not found");
    }
*/
    /**
     * Draws the specified image centered at (<em>x</em>, <em>y</em>).
     * The supported image formats are JPEG, PNG, and GIF.
     * As an optimization, the picture is cached, so there is no performance
     * penalty for redrawing the same image multiple times (e.g., in an animation).
     * However, if you change the picture file after drawing it, subsequent
     * calls will draw the original picture.
     *
     * @param  x the center <em>x</em>-coordinate of the image
     * @param  y the center <em>y</em>-coordinate of the image
     * @param  filename the name of the image/picture, e.g., "ball.gif"
     * @throws IllegalArgumentException if the image filename is invalid
     */
    public static void picture(double x, double y, String filename) {
        // BufferedImage image = getImage(filename);
        Image image = getImage(filename);
        double xs = scaleX(x);
        double ys = scaleY(y);
        // int ws = image.getWidth();    // can call only if image is a BufferedImage
        // int hs = image.getHeight();
        int ws = image.getWidth(null);
        int hs = image.getHeight(null);
        if (ws < 0 || hs < 0) throw new IllegalArgumentException("image " + filename + " is corrupt");

        offscreen.drawImage(image, (int) Math.round(xs - ws/2.0), (int) Math.round(ys - hs/2.0), null);
        draw();
    }

    /**
     * Draws the specified image centered at (<em>x</em>, <em>y</em>),
     * rotated given number of degrees.
     * The supported image formats are JPEG, PNG, and GIF.
     *
     * @param  x the center <em>x</em>-coordinate of the image
     * @param  y the center <em>y</em>-coordinate of the image
     * @param  filename the name of the image/picture, e.g., "ball.gif"
     * @param  degrees is the number of degrees to rotate counterclockwise
     * @throws IllegalArgumentException if the image filename is invalid
     */
    public static void picture(double x, double y, String filename, double degrees) {
        // BufferedImage image = getImage(filename);
        Image image = getImage(filename);
        double xs = scaleX(x);
        double ys = scaleY(y);
        // int ws = image.getWidth();    // can call only if image is a BufferedImage
        // int hs = image.getHeight();
        int ws = image.getWidth(null);
        int hs = image.getHeight(null);
        if (ws < 0 || hs < 0) throw new IllegalArgumentException("image " + filename + " is corrupt");

        offscreen.rotate(Math.toRadians(-degrees), xs, ys);
        offscreen.drawImage(image, (int) Math.round(xs - ws/2.0), (int) Math.round(ys - hs/2.0), null);
        offscreen.rotate(Math.toRadians(+degrees), xs, ys);

        draw();
    }

    /**
     * Draws the specified image centered at (<em>x</em>, <em>y</em>),
     * rescaled to the specified bounding box.
     * The supported image formats are JPEG, PNG, and GIF.
     *
     * @param  x the center <em>x</em>-coordinate of the image
     * @param  y the center <em>y</em>-coordinate of the image
     * @param  filename the name of the image/picture, e.g., "ball.gif"
     * @param  scaledWidth the width of the scaled image (in screen coordinates)
     * @param  scaledHeight the height of the scaled image (in screen coordinates)
     * @throws IllegalArgumentException if either {@code scaledWidth}
     *         or {@code scaledHeight} is negative
     * @throws IllegalArgumentException if the image filename is invalid
     */
    public static void picture(double x, double y, String filename, double scaledWidth, double scaledHeight) {
        Image image = getImage(filename);
        if (scaledWidth  < 0) throw new IllegalArgumentException("width  is negative: " + scaledWidth);
        if (scaledHeight < 0) throw new IllegalArgumentException("height is negative: " + scaledHeight);
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(scaledWidth);
        double hs = factorY(scaledHeight);
        if (ws < 0 || hs < 0) throw new IllegalArgumentException("image " + filename + " is corrupt");
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else {
            offscreen.drawImage(image, (int) Math.round(xs - ws/2.0),
                                       (int) Math.round(ys - hs/2.0),
                                       (int) Math.round(ws),
                                       (int) Math.round(hs), null);
        }
        draw();
    }

    public static void picture(double x, double y, Image image, double scaledWidth, double scaledHeight) {
        if (scaledWidth  < 0) throw new IllegalArgumentException("width  is negative: " + scaledWidth);
        if (scaledHeight < 0) throw new IllegalArgumentException("height is negative: " + scaledHeight);
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(scaledWidth);
        double hs = factorY(scaledHeight);
        if (ws < 0 || hs < 0) throw new IllegalArgumentException("image is corrupt");
        if (ws <= 1 && hs <= 1) pixel(x, y);
        else {
            offscreen.drawImage(image, (int) Math.round(xs - ws/2.0),
                                       (int) Math.round(ys - hs/2.0),
                                       (int) Math.round(ws),
                                       (int) Math.round(hs), null);
        }
        draw();
    }


    /**
     * Draws the specified image centered at (<em>x</em>, <em>y</em>), rotated
     * given number of degrees, and rescaled to the specified bounding box.
     * The supported image formats are JPEG, PNG, and GIF.
     *
     * @param  x the center <em>x</em>-coordinate of the image
     * @param  y the center <em>y</em>-coordinate of the image
     * @param  filename the name of the image/picture, e.g., "ball.gif"
     * @param  scaledWidth the width of the scaled image (in screen coordinates)
     * @param  scaledHeight the height of the scaled image (in screen coordinates)
     * @param  degrees is the number of degrees to rotate counterclockwise
     * @throws IllegalArgumentException if either {@code scaledWidth}
     *         or {@code scaledHeight} is negative
     * @throws IllegalArgumentException if the image filename is invalid
     */
    public static void picture(double x, double y, String filename, double scaledWidth, double scaledHeight, double degrees) {
        if (scaledWidth < 0) throw new IllegalArgumentException("width is negative: " + scaledWidth);
        if (scaledHeight < 0) throw new IllegalArgumentException("height is negative: " + scaledHeight);
        Image image = getImage(filename);
        double xs = scaleX(x);
        double ys = scaleY(y);
        double ws = factorX(scaledWidth);
        double hs = factorY(scaledHeight);
        if (ws < 0 || hs < 0) throw new IllegalArgumentException("image " + filename + " is corrupt");
        if (ws <= 1 && hs <= 1) pixel(x, y);

        offscreen.rotate(Math.toRadians(-degrees), xs, ys);
        offscreen.drawImage(image, (int) Math.round(xs - ws/2.0),
                                   (int) Math.round(ys - hs/2.0),
                                   (int) Math.round(ws),
                                   (int) Math.round(hs), null);
        offscreen.rotate(Math.toRadians(+degrees), xs, ys);

        draw();
    }

   /***************************************************************************
    *  Drawing text.
    ***************************************************************************/

    /**
     * Write the given text string in the current font, centered at (<em>x</em>, <em>y</em>).
     *
     * @param  x the center <em>x</em>-coordinate of the text
     * @param  y the center <em>y</em>-coordinate of the text
     * @param  text the text to write
     */
    public static void text(double x, double y, String text) {
        if (text == null) throw new IllegalArgumentException();
        if (text.equals("")) return;
        offscreen.setFont(font);
        FontMetrics metrics = offscreen.getFontMetrics();
        double xs = scaleX(x);
        double ys = scaleY(y);
        int ws = metrics.stringWidth(text);
        int hs = metrics.getDescent();
        offscreen.drawString(text, (float) (xs - ws/2.0), (float) (ys + hs));
        draw();
    }
    public static void text(double x, double y, String text, Color color) {
        setPencolor(color);
        text(x, y, text);
    }
    public static void text(Point p, String text) {
      text(p.x, p.y, text);
    }
    public static void text(Point p, String text, Color color) {
      setPencolor(color);
      text(p.x, p.y, text);
    }

    //returns the string width as we would input it into StdDraw functions
    public static double getTextWidth(String text){
        offscreen.setFont(font);
        return offscreen.getFontMetrics().stringWidth(text)*0.01;

    }

    public static double getTextHeight(String text){
        offscreen.setFont(font);
        return offscreen.getFontMetrics().getHeight()*0.01;

    }

    /**
     * Write the given text string in the current font, centered at (<em>x</em>, <em>y</em>) and
     * rotated by the specified number of degrees.
     * @param  x the center <em>x</em>-coordinate of the text
     * @param  y the center <em>y</em>-coordinate of the text
     * @param  text the text to write
     * @param  degrees is the number of degrees to rotate counterclockwise
     */
    public static void text(double x, double y, String text, double degrees) {
        if (text == null) throw new IllegalArgumentException();
        if (text.equals("")) return;
        double xs = scaleX(x);
        double ys = scaleY(y);
        offscreen.rotate(Math.toRadians(-degrees), xs, ys);
        text(x, y, text);
        offscreen.rotate(Math.toRadians(+degrees), xs, ys);
    }


    /**
     * Write the given text string in the current font, left-aligned at (<em>x</em>, <em>y</em>).
     * @param  x the <em>x</em>-coordinate of the text
     * @param  y the <em>y</em>-coordinate of the text
     * @param  text the text
     */
    public static void textLeft(double x, double y, String text) {
        if (text == null) throw new IllegalArgumentException();
        if (text.equals("")) return;
        offscreen.setFont(font);
        FontMetrics metrics = offscreen.getFontMetrics();
        double xs = scaleX(x);
        double ys = scaleY(y);
        int hs = metrics.getDescent();
        offscreen.drawString(text, (float) xs, (float) (ys + hs));
        draw();
    }
    public static void textLeft(double x, double y, String text, Color color){
        setPencolor(color);
        textLeft(x, y, text);
    }

    /**
     * Write the given text string in the current font, right-aligned at (<em>x</em>, <em>y</em>).
     *
     * @param  x the <em>x</em>-coordinate of the text
     * @param  y the <em>y</em>-coordinate of the text
     * @param  text the text to write
     */
    public static void textRight(double x, double y, String text) {
        if (text == null) throw new IllegalArgumentException();
        if (text.equals("")) return;
        offscreen.setFont(font);
        FontMetrics metrics = offscreen.getFontMetrics();
        double xs = scaleX(x);
        double ys = scaleY(y);
        int ws = metrics.stringWidth(text);
        int hs = metrics.getDescent();
        offscreen.drawString(text, (float) (xs - ws), (float) (ys + hs));
        draw();
    }
    public static void textRight(double x, double y, String text, Color color){
        setPencolor(color);
        textRight(x, y, text);
    }



    /**
     * Copies the offscreen buffer to the onscreen buffer, pauses for t milliseconds
     * and enables double buffering.
     * @param t number of milliseconds
     * @deprecated replaced by {@link #enableDoubleBuffering()}, {@link #show()}, and {@link #pause(int t)}
     */
    @Deprecated
    public static void show(int t) {
        show();
        pause(t);
        enableDoubleBuffering();
    }

    /**
     * Pause for t milliseconds. This method is intended to support computer animations.
     * @param t number of milliseconds
     */
    public static void pause(int t) {
        try {
            Thread.sleep(t);
        }
        catch (InterruptedException e) {
            System.out.println("Error sleeping");
        }
    }

    /**
     * Copies offscreen buffer to onscreen buffer. There is no reason to call
     * this method unless double buffering is enabled.
     */
    public static void show() {
        onscreen.drawImage(offscreenImage, 0, 0, null);
        frame.repaint();
    }

    // draw onscreen if defer is false
    private static void draw() {
        if (!defer) show();
    }

    /**
     * Enable double buffering. All subsequent calls to
     * drawing methods such as {@code line()}, {@code circle()},
     * and {@code square()} will be deffered until the next call
     * to show(). Useful for animations.
     */
    public static void enableDoubleBuffering() {
        defer = true;
    }

    /**
     * Disable double buffering. All subsequent calls to
     * drawing methods such as {@code line()}, {@code circle()},
     * and {@code square()} will be displayed on screen when called.
     * This is the default.
     */
    public static void disableDoubleBuffering() {
        defer = false;
    }


   /***************************************************************************
    *  Save drawing to a file.
    ***************************************************************************/

    /**
     * Saves the drawing to using the specified filename.
     * The supported image formats are JPEG and PNG;
     * the filename suffix must be {@code .jpg} or {@code .png}.
     *
     * @param  filename the name of the file with one of the required suffixes
     */
    public static void save(String filename) {
        Core.println("save("+filename+")");
        if (filename == null) throw new IllegalArgumentException();
        File file = new File(filename);
        String suffix = filename.substring(filename.lastIndexOf('.') + 1);

        // png files
        if ("png".equalsIgnoreCase(suffix)) {
            try {
                ImageIO.write(onscreenImage, suffix, file);
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }

        // need to change from ARGB to RGB for JPEG
        // reference: http://archives.java.sun.com/cgi-bin/wa?A2=ind0404&L=java2d-interest&D=0&P=2727
        else if ("jpg".equalsIgnoreCase(suffix)) {
            WritableRaster raster = onscreenImage.getRaster();
            WritableRaster newRaster;
            newRaster = raster.createWritableChild(0, 0, width, height, 0, 0, new int[] {0, 1, 2});
            DirectColorModel cm = (DirectColorModel) onscreenImage.getColorModel();
            DirectColorModel newCM = new DirectColorModel(cm.getPixelSize(),
                                                          cm.getRedMask(),
                                                          cm.getGreenMask(),
                                                          cm.getBlueMask());
            BufferedImage rgbBuffer = new BufferedImage(newCM, newRaster, false,  null);
            try {
                ImageIO.write(rgbBuffer, suffix, file);
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }

        else {
            System.out.println("Invalid image file type: " + suffix);
        }
    }


    /**
     * This method cannot be called directly.

     Personal further research @KreuserCorp
        This method is called when the File/Save... Button is pressed. It is called due to an action listener
     */
    @Override
    public void actionPerformed(ActionEvent e) {
        Core.println("actionPerformed()");
        FileDialog chooser = new FileDialog(StdDraw.frame, "Use a .png or .jpg extension", FileDialog.SAVE);
        chooser.setVisible(true);
        String filename = chooser.getFile();
        if (filename != null) {
            StdDraw.save(chooser.getDirectory() + File.separator + chooser.getFile());
        }
    }


   /***************************************************************************
    *  Mouse interactions.
    ***************************************************************************/

    /**
     * Returns true if the mouse is being pressed.
     *
     * @return {@code true} if the mouse is being pressed; {@code false} otherwise
     */
    public static boolean isMousePressed() {
        synchronized (mouseLock) {
            return isMousePressed;
        }
    }
    //@Deprecated
    public static boolean mousePressed() {
        synchronized (mouseLock) {
            return isMousePressed;
        }
    }
    public static double mouseX() {
        synchronized (mouseLock) {
            return mouseX;
        }
    }
    public static double mouseY() {
        synchronized (mouseLock) {
            return mouseY;
        }
    }


//--------------------PERSONAL MOUSE SCROLL ADDIN----------------------

//Check if the mousewheel has been moved
    public static boolean mouseWheelMoved(){
        if (mouseScrolls != 0) return true;
        return false;
    }
    private static boolean mouseScrolled = false;

//Return the amount of units by which the mouse has been scrolle
    @Override
    public void mouseWheelMoved(MouseWheelEvent e){
        mouseScrolls = -e.getWheelRotation();
    }

    private static int mouseScrolls;
    public static int mouseScrolls(){
        int temp = mouseScrolls;
        mouseScrolls = 0;
        return temp;
    }

    public static Point mouseLoc(){
        return new Point(mouseX(), mouseY());
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mouseClicked(MouseEvent e) {
        // this body is intentionally left empty
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mouseEntered(MouseEvent e) {
        // this body is intentionally left empty
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mouseExited(MouseEvent e) {
        // this body is intentionally left empty
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mousePressed(MouseEvent e) {
        synchronized (mouseLock) {
            mouseX = StdDraw.userX(e.getX());
            mouseY = StdDraw.userY(e.getY());
            isMousePressed = true;
        }
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mouseReleased(MouseEvent e) {
        synchronized (mouseLock) {
            isMousePressed = false;
        }
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mouseDragged(MouseEvent e)  {
        synchronized (mouseLock) {
            mouseX = StdDraw.userX(e.getX());
            mouseY = StdDraw.userY(e.getY());
        }
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void mouseMoved(MouseEvent e) {
        synchronized (mouseLock) {
            mouseX = StdDraw.userX(e.getX());
            mouseY = StdDraw.userY(e.getY());
        }
    }


   /***************************************************************************
    *  Keyboard interactions.
    ***************************************************************************/

    /**
     * Returns true if the user has typed a key (that has not yet been processed).
     *
     * @return {@code true} if the user has typed a key (that has not yet been processed
     *         by {@link #nextKeyTyped()}; {@code false} otherwise
     */
    public static boolean hasNextKeyTyped() {
        synchronized (keyLock) {
            return !keysTyped.isEmpty();
        }
    }

    /**
     * Returns the next key that was typed by the user (that your program has not already processed).
     * This method should be preceded by a call to {@link #hasNextKeyTyped()} to ensure
     * that there is a next key to process.
     * This method returns a Unicode character corresponding to the key
     * typed (such as {@code 'a'} or {@code 'A'}).
     * It cannot identify action keys (such as F1 and arrow keys)
     * or modifier keys (such as control).
     *
     * @return the next key typed by the user (that your program has not already processed).
     * @throws NoSuchElementException if there is no remaining key
     */
    public static char nextKeyTyped() {
        synchronized (keyLock) {
            if (keysTyped.isEmpty()) {
                throw new NoSuchElementException("your program has already processed all keystrokes");
            }
            return keysTyped.remove(keysTyped.size() - 1);
            // return keysTyped.removeLast();
        }
    }

    /**
     * Returns true if the given key is being pressed.
     * <p>
     * This method takes the keycode (corresponding to a physical key)
    *  as an argument. It can handle action keys
     * (such as F1 and arrow keys) and modifier keys (such as shift and control).
     * See {@link KeyEvent} for a description of key codes.
     *
     * @param  keycode the key to check if it is being pressed
     * @return {@code true} if {@code keycode} is currently being pressed;
     *         {@code false} otherwise
     */
    public static boolean isKeyPressed(int keycode) {
        synchronized (keyLock) {
            return keysDown.contains(keycode);
        }
    }


    /**
     * This method cannot be called directly.
     */
    @Override
    public void keyTyped(KeyEvent e) {
        synchronized (keyLock) {
            keysTyped.addFirst(e.getKeyChar());
        }
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void keyPressed(KeyEvent e) {
        synchronized (keyLock) {
            keysDown.add(e.getKeyCode());
        }
    }

    /**
     * This method cannot be called directly.
     */
    @Override
    public void keyReleased(KeyEvent e) {
        synchronized (keyLock) {
            keysDown.remove(e.getKeyCode());
        }
    }




    /**
     * Test client.
     *
     * @param args the command-line arguments
     */
    public static void main(String[] args) {
        StdDraw.square(0.2, 0.8, 0.1);
        StdDraw.filledSquare(0.8, 0.8, 0.2);
        StdDraw.circle(0.8, 0.2, 0.2);

        StdDraw.setPenColor(StdDraw.BOOK_RED);
        StdDraw.setPenRadius(0.02);
        StdDraw.arc(0.8, 0.2, 0.1, 200, 45);

        // draw a blue diamond
        StdDraw.setPenRadius();
        StdDraw.setPenColor(StdDraw.BOOK_BLUE);
        double[] x = { 0.1, 0.2, 0.3, 0.2 };
        double[] y = { 0.2, 0.3, 0.2, 0.1 };
        StdDraw.filledPolygon(x, y);

        // text
        StdDraw.setPenColor(StdDraw.BLACK);
        StdDraw.text(0.2, 0.5, "black text");
        StdDraw.setPenColor(StdDraw.WHITE);
        StdDraw.text(0.8, 0.8, "white text");
    }

}
