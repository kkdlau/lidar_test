#include "PointArray.hpp"
#include "PolarVector.hpp"
#include <vector>
using namespace std;

namespace Fliter {
template <typename T>
void withinRange(PointArray<T> &arr, const double min, const double max) {
  for (Point<T> &p : arr) {
    const PolarVector v = PolarVector::fromCartesian(p);
    if (v.mag > max || v.mag < min) {
      p.x = 0;
      p.y = 0;
    }
  }
}

template <typename T>
PointArray<T> withinRad(PointArray<T> &arr, const double min,
                        const double max) {
  PointArray<T> ret;
  for (Point<T> &p : arr) {
    const PolarVector v = PolarVector::fromCartesian(p);
    if (v.rad > max || v.rad < min) {
      p.x = 0;
      p.y = 0;
    } else {
      ret.push_back(p);
    }
  }

  return ret;
}
template <typename T>
PointArray<T> removeInvalidData(const PointArray<T> &arr) {
  PointArray<T> ret;

  for (const Point<T> &p : arr) {
    if (PolarVector::fromCartesian(p).mag == 0)
      continue;

    ret.push_back(p);
  }

  return ret;
}

/**
 * @brief Calculate the first derivative of given point.
 * The first derivative is d(distance)/d(rad).
 *
 * It will look for next point first, if it cannot find a point for calculation,
 * then it will look backwards.
 *
 * @tparam T
 * @param arr
 * @param index
 * @return float
 */
template <typename T>
float firstDerivative(const PointArray<T> &arr, int index) {
  const int size = arr.size();
  const PolarVector pv = PolarVector::fromCartesian(arr[index]);
  for (int i = index + 1; i < min(size, index + 5); ++i) {
    const PolarVector next = PolarVector::fromCartesian(arr[i]);
    if (next.mag != 0) {
      return (pv.mag - next.mag) / (pv.rad - next.rad);
    }
  }

  for (int i = index - 1; i >= max(0, index - 5); --i) {
    const PolarVector prev = PolarVector::fromCartesian(arr[i]);
    if (prev.mag != 0) {
      return (pv.mag - prev.mag) / (pv.rad - prev.rad);
    }
  }

  return .0f;
}

float firstDerivative(const vector<PolarVector> &arr, int index) {
  const int size = arr.size();
  const PolarVector pv = arr[index];
  for (int i = index + 1; i < min(size, index + 5); ++i) {
    const PolarVector next = arr[i];
    if (next.mag != 0) {
      return (pv.mag - next.mag) / (pv.rad - next.rad);
    }
  }

  for (int i = index - 1; i >= max(0, index - 5); --i) {
    const PolarVector prev = arr[i];
    if (prev.mag != 0) {
      return (pv.mag - prev.mag) / (pv.rad - prev.rad);
    }
  }

  return .0f;
}

template <typename T>
vector<int> localMinimum(const PointArray<T> &arr, const double searchingRad) {
  vector<int> min;

  for (int i = 1; i < arr.size() - 1; ++i) {
    PolarVector pv = PolarVector::fromCartesian(arr[i]);

    if (pv.mag == 0) // invalid data, skip
      continue;

    float dwSlope = .0f; // downward slope
    int numDW = 0;
    float uwSlope = .0f; // upward slope
    int numUW = 0;

    const double curRad = pv.rad;

    int r = i;
    do {
      r--;
      if (PolarVector::fromCartesian(arr[r]).mag == 0)
        continue;
      float d = firstDerivative(arr, r);
      if (!d || isinf(d) || isnan(d))
        continue;
      dwSlope += d;
      ++numDW;
      // cout << "searching dw: " << d << endl;
    } while (r >= 0 &&
             PolarVector::fromCartesian(arr[r]).rad >= curRad - searchingRad);

    int k = i;

    do {
      k++;
      if (PolarVector::fromCartesian(arr[k]).mag == 0)
        continue;
      float d = firstDerivative(arr, k);
      if (!d || isinf(d) || isnan(d))
        continue;
      uwSlope += d;
      ++numUW;
      // cout << "searching uw: " << d << endl;

    } while (k < arr.size() &&
             PolarVector::fromCartesian(arr[k]).rad <= curRad + searchingRad);

    if (dwSlope / numDW < 0 && uwSlope / numUW > 0) {
      // cout << "push a min index" << endl;
      min.push_back(i);
      i += k;
    } else {
      // cout << "fk" << endl;
    }
  }

  return min;
}
} // namespace Fliter
