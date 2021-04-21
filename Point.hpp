#pragma once
#include <string>

class PolarVector;

template <typename T> class Point {
public:
  T x;
  T y;

  Point() {
    x = 0;
    y = 0;
  }

  Point(T x, T y) {
    this->x = x;
    this->y = y;
  }

  Point<T> &operator=(const Point<T> &p) {
    x = p.x;
    y = p.y;
    return *this;
  }
  Point<T> operator+(Point<T> &p) const {
    return Point<T>(this->x + p.x, this->y + p.y);
  }
  Point<T> operator/(double d) const {
    return Point<T>(this->x / d, this->y / d);
  }

  static Point<T> fromPolar(const PolarVector &p);

  std::string toString() const {
    return std::string("x: ") + std::to_string(x) + ", " + std::string("y: ") +
           std::to_string(y);
  }

  std::string toCSVString() const {
    return std::to_string(x) + ", " + std::to_string(y);
  }

  float slope(const Point &p) const {
    return (this->y - p.y) / (this->x - p.x);
  }
};