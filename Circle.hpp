/**
 * @file Circle.hpp
 *
 * This file defines Circle Class.
 *
 * Class member & explaination:
 * - Point<double> center: it stores the center of Circle
 * - double& a, &b: legacy purpose, it references to center.x and center.y
 * - double r: radius of Circle
 * - double s: the estimate of sigma (standard deviation)
 * - double g: the norm of the gradient of the objective function
 * - int i: number of used outer iteration for circle fitting
 * - int k: number of used inner iteration for circle fitting
 */
#pragma once
#include "Point.hpp"
#include "helper_funcs.hpp"

class Circle {
public:
  // The fields of a Circle
  double r, s, g, Gx, Gy;
  double &a, &b;
  Point<double> center;
  int i, j;

  // constructors
  Circle() : center(), a(center.x), b(center.y) {
    a = 0.;
    b = 0.;
    r = 1.;
    s = 0.;
    i = 0;
    j = 0;
  }
  Circle(double cx, double cy, double radius)
      : center(cx, cy), a(center.x), b(center.y) {
    r = radius;
  }

  Circle &operator=(const Circle &c) {
    center = c.center;
    r = c.r;
    s = c.s;
    i = c.i;
    j = c.j;

    return *this;
  }

  // routines
  void print(void) {
    cout << endl;
    cout << setprecision(10) << "center (" << a << "," << b << ")  radius " << r
         << "  sigma " << s << "  gradient " << g << "  iter " << i
         << "  inner " << j << endl;
  }
};