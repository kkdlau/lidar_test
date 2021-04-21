#pragma once
#include "Point.hpp"
#include <cmath>
#include <string>

class PolarVector {
public:
  double rad;
  double mag;

  PolarVector() {
    rad = .0;
    mag = .0;
  }

  PolarVector(const double r, const double m) {
    rad = r;
    mag = m;
  }

  static PolarVector fromCartesian(const Point<double> &p) {
    return PolarVector(atan(p.y / p.x), sqrt(p.x * p.x + p.y * p.y));
  }

  static PolarVector inDegree(const double degree, const double mag) {
    return PolarVector(degree / 180 * M_PI, mag);
  }

  Point<double> toCartesian() const {
    return Point<double>(cos(rad) * mag, sin(rad) * mag);
  }

  std::string toString() {
    return "angle: " +
           std::to_string(rad / (2 * 3.1415926) * 360) +
          ", mag: " + std::to_string(mag);
  }

    std::string toCSVString() {
    return std::to_string(rad / (2 * M_PI) * 360) +
          "," + std::to_string(mag);
  }
};

template <typename T> Point<T> Point<T>::fromPolar(const PolarVector &p) {
  return Point<T>(cos(p.rad) * p.mag, sin(p.rad) * p.mag);
}