
#pragma once

#include "misc_v9.h"

#ifdef RCPP_ACTIVE
#include <Rcpp.h>
#endif

//------------------------------------------------
// test for collision between a line of the form y = mx + k and a finite line
// between two points
bool collision_test_abline_line(double m, double k,
                                double l1x, double l1y, double l2x, double l2y,
                                double &hx, double &hy);

//------------------------------------------------
// test for collision between a rectangle and a circle
bool collision_test_circle_rect(double rx, double ry, double w, double h, double theta,
                                double cx, double cy, double r);

//------------------------------------------------
// test for collision between a point and an ellipse
bool collision_test_point_ellipse(double x, double y, double f1x, double f1y, double f2x, double f2y, double e);

//------------------------------------------------
// test for collision between a line of the form y = mx + k and an ellipse
bool collision_test_abline_ellipse(double m, double k,
                                   double f1x, double f1y, double f2x, double f2y, double e,
                                   double &h1x, double &h1y, double &h2x, double &h2y);

//------------------------------------------------
// test for collision between a line between two points and an ellipse
bool collision_test_line_ellipse(double l1x, double l1y, double l2c, double l2y,
                                 double f1x, double f1y, double f2x, double f2y, double e);

//------------------------------------------------
// test for collision between a hexagon and an ellipse
bool collision_test_hex_ellipse(double hx, double hy, double w,
                                double f1x, double f1y, double f2x, double f2y, double e);
