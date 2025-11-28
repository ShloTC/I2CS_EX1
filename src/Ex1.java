/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
	/** Epsilon value for numerical computation, it serves as a "close enough" threshold. */
	public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
	/** The zero polynomial function is represented as an array with a single (0) entry. */
	public static final double[] ZERO = {0};
	/**
	 * Computes the f(x) value of the polynomial function at x.
	 * @param poly - polynomial function
	 * @param x
	 * @return f(x) - the polynomial function value at x.
	 */
	public static double f(double[] poly, double x) {
		double ans = 0;
		for(int i=0;i<poly.length;i++) {
			double c = Math.pow(x, i);
			ans += c*poly[i];
		}
		return ans;
	}
	/** Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
	 * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps, 
	 * assuming p(x1)*p(x2) <= 0.
	 * This function should be implemented recursively.
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
	 */
	public static double root_rec(double[] p, double x1, double x2, double eps) {
		double f1 = f(p,x1);
		double x12 = (x1+x2)/2;
		double f12 = f(p,x12);
		if (Math.abs(f12)<eps) {return x12;}
		if(f12*f1<=0) {return root_rec(p, x1, x12, eps);}
		else {return root_rec(p, x12, x2, eps);}
	}
	/**
	 * This function computes a polynomial representation from a set of 2D points on the polynom.
	 * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
	 * Note: this function only works for a set of points containing up to 3 points, else returns null.
	 * @param xx - array containing every x coordinate for each 2d point
	 * @param yy - array containing every y coordinate for each 2d point
	 * @return an array of doubles representing the coefficients of the polynom.
	 */
	public static double[] PolynomFromPoints(double[] xx, double[] yy) {
		double [] ans = null;
		int lx = xx.length;
		int ly = yy.length;
		if(xx!=null && yy!=null && lx==ly && lx>1 && lx<4) {
		/* add you code below */
            if (lx == 2) {
                //f(x) = a*x + b
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];

                double a = (y2 - y1) / (x2 - x1);
                double b = y1 - a*x1;

                ans = new double[]{b, a};
            }
            else {
                // f(x) = a*x^2 + b*x + c
                double x1 = xx[0], y1 = yy[0];
                double x2 = xx[1], y2 = yy[1];
                double x3 = xx[2], y3 = yy[2];

                double denom = (x1-x2)*(x1-x3)*(x2-x3);

                double a = ( y1*(x2-x3) + y2*(x3-x1) + y3*(x1-x2) ) / denom;

                double b = ( y1*(x3*x3 - x2*x2)
                        + y2*(x1*x1 - x3*x3)
                        + y3*(x2*x2 - x1*x1) ) / denom;

                double c = ( y1*(x2*x3*(x2-x3))
                        + y2*(x3*x1*(x3-x1))
                        + y3*(x1*x2*(x1-x2)) ) / denom;

                ans = new double[]{c, b, a};
            }
		/* /////////////////// */

		}
		return ans;
	}
	/** Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
	 * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
	 * @param p1 first polynomial function
	 * @param p2 second polynomial function
	 * @return true iff p1 represents the same polynomial function as p2 (within epsilon)
	 */
	public static boolean equals(double[] p1, double[] p2) {
		boolean ans = true;
        /* add you code below */
        if (p1 == null || p2 == null) return false;
        // highest exponent
        int n = Math.max(p1.length, p2.length) - 1;
        // test at n+1 points: x = 0,...,n (two polynoms of degree <=n that agree on n+1 points are identical)
        for (int x = 0; x <= n; x++) {
            double v1 = f(p1, x);
            double v2 = f(p2, x);
            if (Math.abs(v1 - v2) > EPS) {
                return false; //if current exponent has a difference of more than epsilon between both polynoms
            }
        }
         /////////////////// */
		return ans;
	}

	/** 
	 * Computes a String representing the polynomial function.
	 * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
	 * @param poly the polynomial function represented as an array of doubles
	 * @return String representing the polynomial function:
	 */
	public static String poly(double[] poly) {
		String ans = "";
		if(poly.length==0) {ans="0";}
		else {
            /* add you code below*/
            StringBuilder sb = new StringBuilder(); //build alterable string
            boolean first = true;

            // loop from highest power to the lowest
            for (int i = poly.length - 1; i >= 0; i--) {
                double c = poly[i];

                if (Math.abs(c) < EPS) continue; // skip zero terms

                // Handle sign
                if (first) {
                    // first term: write as-is
                    sb.append(c);
                } else {
                    if (c >= 0) sb.append(" +").append(c);
                    else sb.append(" ").append(c); // c is already negative
                }

                // Handle power
                if (i == 1) sb.append("x");
                else if (i >= 2) sb.append("x^").append(i);

                first = false;
            }

            // If all coefficients were zero
            if (first) return "0";

            return sb.toString();

             /////////////////// */
		}
		return ans;
	}
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
	 * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0. (find x within range where both functions have the same y)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
	 * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps. (x within range where both functions have the same y)
	 */
	public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        /* add you code below */
        double g1 = f(p1, x1) - f(p2, x1);
        double mid = (x1 + x2) / 2.0;
        double gmid = f(p1, mid) - f(p2, mid);

        if (Math.abs(gmid) < eps) return mid;

        if (g1 * gmid <= 0) {
            return sameValue(p1, p2, x1, mid, eps);
        } else {
            return sameValue(p1, p2, mid, x2, eps);
        }

    }
         /////////////////// */

	/**
	 * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
	 * This function computes an approximation of the length of the function between f(x1) and f(x2) 
	 * using n inner sample points and computing the segment-path between them.
	 * assuming x1 < x2. 
	 * This function should be implemented iteratively (none recursive).
	 * @param p - the polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfSegments - (A positive integer value (1,2,...).
	 * @return the length approximation of the function between f(x1) and f(x2).
	 */
	public static double length(double[] p, double x1, double x2, int numberOfSegments) {
		double ans = x1;
        /** add you code below

         /////////////////// */
		return ans;
	}
	
	/**
	 * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
	 * This function computes an approximation of the area between the polynomial functions within the x-range.
	 * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
	 * @param p1 - first polynomial function
	 * @param p2 - second polynomial function
	 * @param x1 - minimal value of the range
	 * @param x2 - maximal value of the range
	 * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
	 * @return the approximated area between the two polynomial functions within the [x1,x2] range.
	 */
	public static double area(double[] p1,double[]p2, double x1, double x2, int numberOfTrapezoid) {
		double ans = 0;
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the array representation of a polynomial function from a String
	 * representation. Note:given a polynomial function represented as a double array,
	 * getPolynomFromString(poly(p)) should return an array equals to p.
	 * 
	 * @param p - a String representing polynomial function.
	 * @return
	 */
	public static double[] getPolynomFromString(String p) {
		double [] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /* add you code below */

        //replace " -" with "+-" so we can split by "+"
        p = p.replace(" -", "+-");
        String[] parts = p.split("\\+");

        // First pass: find max power
        int maxPower = 0;
        for (String term : parts) {
            term = term.trim();
            if (term.isEmpty()) continue;

            if (term.contains("x^")) {
                int pow = Integer.parseInt(term.substring(term.indexOf("^") + 1));
                if (pow > maxPower) maxPower = pow;
            }
            else if (term.contains("x")) {
                if (maxPower < 1) maxPower = 1;
            }
        }

        ans = new double[maxPower + 1];

        // Second pass: fill coefficients
        for (String term : parts) {
            term = term.trim();
            if (term.isEmpty()) continue;

            if (term.contains("x^")) {
                // general term: ax^n
                int pow = Integer.parseInt(term.substring(term.indexOf("^") + 1));
                double coef = Double.parseDouble(term.substring(0, term.indexOf("x")));
                ans[pow] = coef;
            }
            else if (term.contains("x")) {
                // linear term: ax
                double coef = Double.parseDouble(term.substring(0, term.indexOf("x")));
                ans[1] = coef;
            }
            else {
                // constant term: c
                double coef = Double.parseDouble(term);
                ans[0] = coef;
            }
        }

        return ans;
    }
         /////////////////// */
	/**
	 * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] add(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
	 * @param p1
	 * @param p2
	 * @return
	 */
	public static double[] mul(double[] p1, double[] p2) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
		return ans;
	}
	/**
	 * This function computes the derivative of the p0 polynomial function.
	 * @param po
	 * @return
	 */
	public static double[] derivative (double[] po) {
		double [] ans = ZERO;//
        /** add you code below

         /////////////////// */
		return ans;
	}
}
