
#include "utils.h"

#include <vector>

using namespace std;


//------------------------------------------------
// check if a point is within ellipse
// x and y are the coordinates of the query point. f1 and f2 are the coordinates
// of the two foci of the ellipse and e is the eccentricity (0 = circle, 1 =
// line).
bool point_intersects_ellipse(double x, double y,
                              double f1x, double f1y, double f2x, double f2y, double e) {
  
  // deal with limiting values of eccentricity
  if (e < 0.0 | e > 1.0) {
    Rcpp::stop("error in point_intersects_ellipse(): e outside range [0,1]");
  }
  if (e == 0.0) {  // infinitely large circle
    return true;
  }
  if (e == 1.0) { // line between foci
    return false;
  }
  
  // get linear eccentricity and semi-major axis of ellipse
  double c = 0.5*dist_euclid_2d(f1x, f1y, f2x, f2y);
  double a = c/e;
  
  // get summed distance from point to each focus
  double d = dist_euclid_2d(x, y, f1x, f1y) + dist_euclid_2d(x, y, f2x, f2y);
  
  // return whether point is inside ellipse
  return (d <= 2*a);
}

//------------------------------------------------
// check if a line of the form y = mx + k intersects an ellipse
// m and k are the gradient and intercept of the line. f1 and f2 are the
// coordinates of the two foci of the ellipse and e is the eccentricity (0 =
// circle, 1 = line). If the line does intersect then in addition to returning
// true, the x- and y- coordinates of the two intersection points are stored in
// [h1x, h1y] and [h2x, h2y].
bool abline_intersects_ellipse(double m, double k,
                               double f1x, double f1y, double f2x, double f2y, double e,
                               double &h1x, double &h1y, double &h2x, double &h2y) {
  
  // deal with limiting values of eccentricity
  if (e < 0.0 | e > 1.0) {
    Rcpp::stop("error in abline_intersects_ellipse(): e outside range [0,1]");
  }
  if (e == 0.0) {  // infinitely large circle
    return true;
  }
  if (e == 1.0) { // line between foci
    return false;
  }
  
  // get further useful properties of the ellipse
  double c = 0.5*dist_euclid_2d(f1x, f1y, f2x, f2y);  // linear eccentricity
  double a = c/e;  // semi-major axis
  double b = sqrt(sq(a) - sq(c)); // semi-minor axis
  double qx = 0.5*(f1x + f2x);  // coordinates of centre
  double qy = 0.5*(f1y + f2y);
  double theta = atan2(f2y - f1y, f2x - f1x);  // angle of ellipse, in radians measured counter-clockwise from the positive x-axis
  
  // define some terms for convenience
  double a2 = sq(a);
  double b2 = sq(b);
  double qx2 = sq(qx);
  double qy2 = sq(qy);
  double sint = sin(theta);
  double cost = cos(theta);
  double sin2t = sq(sint);
  double cos2t = sq(cost);
  
  // the general equation of an ellipse centred at [qx,qy] and rotated by theta is:
  // ((x-qx)cos(theta)+(y-qy)sin(theta))^2/a^2 + ((x-qx)sin(theta)-(y-qy)cos(theta))^2/b^2 = 1
  // this can be expressed in a polynomial of the form:
  // Ax^2 + Bx + Cy^2 + Dy + Exy + F = 0
  // The values of the constants are:
  double A = b2*cos2t + a2*sin2t;
  double B = -2*b2*qx*cos2t - 2*a2*qx*sin2t -2*qy*cost*sint*(b2-a2);
  double C = b2*sin2t + a2*cos2t;
  double D = -2*b2*qy*sin2t - 2*a2*qy*cos2t -2*qx*cost*sint*(b2-a2);
  double E = 2*sint*cost*(b2-a2);
  double F = cos2t*(b2*qx2 + a2*qy2) + sin2t*(b2*qy2 + a2*qx2) + 2*qx*qy*cost*sint*(b2-a2) - a2*b2;
  
  // substituting y = mx + k into the polynomial we obtain a quadratic:
  // Gx^2 + Hx + I = 0
  // The values of the constants are
  double G = A + C*sq(m) + E*m;
  double H = B + 2*C*m*k + D*m + E*k;
  double I = C*sq(k) + D*k + F;
  
  // if the descrimant of the quadratic is negative then the line does not
  // intersect, hence return false
  double desc = sq(H) - 4*G*I;
  if (desc < 0) {
    return false;
  }
  
  // if the descriminant is positive then calculate the two solutions,
  // corresponding to the two points of intersection
  double root_desc = sqrt(desc);
  h1x = (-H + root_desc)/(2*G);
  h2x = (-H - root_desc)/(2*G);
  h1y = m*h1x + k;
  h2y = m*h2x + k;
  
  return true;
}

//------------------------------------------------
// check if a line connecting two coordinates intersects an ellipse
// l1 and l2 are the start- and end-points of a line. f1 and f2 are the
// coordinates of the two foci of the ellipse and e is the eccentricity (0 =
// circle, 1 = line).
bool line_intersects_ellipse(double l1x, double l1y, double l2x, double l2y,
                             double f1x, double f1y, double f2x, double f2y, double e) {
  
  // ensure coordinates are increasing in x. Swap if not
  if (l2x < l1x) {
    double tmp = l1x;
    l1x = l2x;
    l2x = tmp;
    tmp = l1y;
    l1y = l2y;
    l2y = tmp;
  }
  
  // check if start or end of line intersects the ellipse. If so then return
  // true
  if (point_intersects_ellipse(l1x, l1y, f1x, f1y, f2x, f2y, e)) {
    return true;
  }
  if (point_intersects_ellipse(l2x, l2y, f1x, f1y, f2x, f2y, e)) {
    return true;
  }
  
  // the input coordinates are a subset of the infinite line y = mx + k.
  // Determine the values of the intercept and gradient of this line
  double m = (l2y - l1y)/(l2x - l1x);
  double k = l1y - m*l1x;
  
  // determine whether the infinite line intersects the ellipse. If not then
  // return
  double h1x, h1y, h2x, h2y;
  bool ab_intersect = abline_intersects_ellipse(m, k, f1x, f1y, f2x, f2y, e, h1x, h1y, h2x, h2y);
  if (!ab_intersect) {
    return false;
  }
  
  // determine whether either of the points of intersection lie between the
  // original coordinates
  if (h1x > l1x & h1x < l2x) {
    return true;
  }
  if (h2x > l1x & h2x < l2x) {
    return true;
  }
  
  return false;
}

//------------------------------------------------
// check if a hex intersects an ellipse
// [hx,hy] are the coordinates of the centroid of a hex (with points at top and
// bottom) and hs is the width of this hex. f1 and f2 are the coordinates of the
// two foci of the ellipse and e is the eccentricity (0 = circle, 1 = line).
bool hex_intersects_ellipse(double hx, double hy, double hs,
                            double f1x, double f1y, double f2x, double f2y, double e) {
  
  // manually define 7 points: the 6 points of the hex (starting at the top)
  // plus the first point repeated to complete the polygon
  double c = 1/sqrt(3);
  vector<double> px = {hx, hx+hs/2, hx+hs/2, hx, hx-hs/2, hx-hs/2, hx};
  vector<double> py = {hy+hs*c, hy+hs/2*c, hy-hs/2*c, hy-hs*c, hy-hs/2*c, hy+hs/2*c, hy+hs*c};
  
  // if any line of the hex polygon intersects the ellipse then return true
  for (unsigned int i = 0; i < 6; ++i) {
    if (line_intersects_ellipse(px[i], py[1], px[i+1], py[i+1], f1x, f1y, f2x, f2y, e)) {
      return true;
    }
  }
  
  return false;
}