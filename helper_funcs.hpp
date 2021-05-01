#pragma once
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

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

static void split(const string &s, vector<string> &tokens,
                  const string &delimiters = " ") {
  string::size_type lastPos = s.find_first_not_of(delimiters, 0);
  string::size_type pos = s.find_first_of(delimiters, lastPos);
  while (string::npos != pos || string::npos != lastPos) {
    tokens.push_back(
        s.substr(lastPos, pos - lastPos)); // use emplace_back after C++11
    lastPos = s.find_first_not_of(delimiters, pos);
    pos = s.find_first_of(delimiters, lastPos);
  }
}