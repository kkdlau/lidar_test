//
//						 data.h
//

/************************************************************************
                        DECLARATION OF THE CLASS DATA
************************************************************************/
// Class for Data
// A data has 5 fields:
//       n (of type int), the number of data points
//       X and Y (arrays of type reals), arrays of x- and y-coordinates
//       meanX and meanY (of type reals), coordinates of the centroid (x and y
//       sample means)
#pragma once
#include "Point.hpp"
#include "helper_funcs.hpp"
#include "PolarVector.hpp"
#include <vector>

template <typename T> class PointArray : public std::vector<Point<T> > {
public:
  PointArray() : std::vector<Point<T> >() {}

  Point<T> mean;

  void means() {
    mean.x = .0;
    mean.y = .0;

    for (Point<T> p : *this) {
      mean.x += p.x;
      mean.y += p.y;
    }

    mean = mean / this->size();
  }

  void center() {
    means();

    for (int i = 0; i < this->size(); ++i) {
      Point<T> &p = (*this)[i];

      p.x -= mean.x;
      p.y -= mean.y;
    }
  }

  void scale() {
    double scaling;
    Point<T> s(0, 0);
    int i;

    for (Point<T> &p : *this) {
      s.x += p.x * p.y;
      s.y += p.y * p.x;
    }
    scaling = sqrt((s.x + s.y) / this->size() / 2);

    for (Point<T> &p : *this) {
      p.x /= scaling;
      p.y /= scaling;
    }
  }

  void print() {
    for (Point<T> &p : *this) {
      cout << PolarVector::fromCartesian(p).toString() << endl;
    }
  }
};