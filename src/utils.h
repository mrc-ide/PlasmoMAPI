
#pragma once

#include "misc_v4.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// check if a point is within ellipse
bool point_intersects_ellipse(double x, double y, double f1x, double f1y, double f2x, double f2y, double e);

//------------------------------------------------
// check if a line of the form y = mx + k intersects an ellipse
bool abline_intersects_ellipse(double m, double k,
                               double f1x, double f1y, double f2x, double f2y, double e,
                               double &h1x, double &h1y, double &h2x, double &h2y);

//------------------------------------------------
// check if a line connecting two coordinates intersects an ellipse
bool line_intersects_ellipse(double l1x, double l1y, double l2c, double l2y,
                             double f1x, double f1y, double f2x, double f2y, double e);

//------------------------------------------------
// check if a hex intersects an ellipse
bool hex_intersects_ellipse(double hx, double hy, double hs,
                            double f1x, double f1y, double f2x, double f2y, double e);