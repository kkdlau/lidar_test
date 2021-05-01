#include "PointArray.hpp"
#include "PolarVector.hpp"
#include <functional>
#include <vector>
using namespace std;

enum class Trend { DOWNWARD, NONE, UPWARD };
inline Trend slopeTrend(const double slope) {
  if (isinf(slope) || isnan(slope))
    return Trend::NONE;
  else if (slope < 0)
    return Trend::DOWNWARD;
  else
    return Trend::UPWARD;
}
namespace Filter {
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
double avgSlope(const PointArray<T> &arr, const int START, const int dir,
                function<bool(double)> shouldStop) {
  double v = 0;
  unsigned count = 0;

  int i = START;
  PolarVector pv;

  do {
    double d = firstDerivative(arr, i);
    if (d != 0) {
      v += d;
      count++;
    }
    i += 1 * dir;
    pv = PolarVector::fromCartesian(arr[i]);
  } while (!shouldStop(pv.rad) && i >= 0 && i < arr.size());

  return v / count;
}

template <typename T>
vector<int> localMinimum(const PointArray<T> &arr, double searchingRad) {
  const double POT_RADIUS = 160.0;
  const double MIN_SEARCHING_PLANE_SIZE = 40.0;
  vector<int> min;

  for (int i = 1; i < arr.size() - 1; ++i) {
    PolarVector pv = PolarVector::fromCartesian(arr[i]);

    if (pv.mag == 0) // invalid data, skip
      continue;

    const double curRad = pv.rad;
    for (int plane = 150; plane >= MIN_SEARCHING_PLANE_SIZE; plane -= 10) {
      searchingRad = atan2(plane, pv.mag + POT_RADIUS);

      double lslope = avgSlope(arr, i, -1, [curRad, searchingRad](double rad) {
        return rad < curRad - searchingRad;
      });

      double rslope = avgSlope(arr, i, 1, [curRad, searchingRad](double rad) {
        return rad > curRad + searchingRad;
      });
#define to_deg(r) ((r) / M_PI * 180)
      if (slopeTrend(lslope) == Trend::DOWNWARD &&
          slopeTrend(rslope) == Trend::UPWARD) {
        cout << "min on:" << to_deg(pv.rad) << ", searching from "
             << to_deg(pv.rad - searchingRad) << " to "
             << to_deg(pv.rad + searchingRad) << endl;
        min.push_back(i);
        while (PolarVector::fromCartesian(arr[i]).rad < curRad + searchingRad)
          i++;

        break;
      }
    }
  }

  return min;
}
} // namespace Filter
