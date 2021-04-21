#pragma once
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace std;

const double REAL_MAX = numeric_limits<double>::max();
const double REAL_MIN = numeric_limits<double>::min();
const double REAL_EPSILON = numeric_limits<double>::epsilon();

template <typename T> static inline T SQR(T t) { return t * t; };

static double pythag(double a, double b) {
  a = abs(a), b = abs(b);

  if (a > b)
    return a * sqrt(1 + SQR(b / a));
  else
    return (b == 0.0 ? 0.0 : b * sqrt(1 + SQR(a / b)));
}