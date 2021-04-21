#include "Circle.hpp"
#include "PointArray.hpp"
#include "Statistics.hpp"
#include "helper_funcs.hpp"

namespace CircleFittingModel {
Circle kasaFitting(PointArray<double> &data)
/*
      Circle fit to a given set of data points (in 2D)

      This is an algebraic fit, disovered and rediscovered by many people.
      One of the earliest publications is due to Kasa:

      I. Kasa, "A curve fitting procedure and its error analysis",
      IEEE Trans. Inst. Meas., Vol. 25, pages 8-14, (1976)

      Input:  data     - the class of data (contains the given points):

              data.n   - the number of data points
              data.X[] - the array of X-coordinates
              data.Y[] - the array of Y-coordinates

     Output:
               circle - parameters of the fitting circle:

               circle.a - the X-coordinate of the center of the fitting circle
               circle.b - the Y-coordinate of the center of the fitting circle
               circle.r - the radius of the fitting circle
               circle.s - the root mean square error (the estimate of sigma)
               circle.j - the total number of iterations

     The method is based on the minimization of the function

                 F = sum [(x-a)^2 + (y-b)^2 - R^2]^2

     This is perhaps the simplest and fastest circle fit.

     It works well when data points are sampled along an entire circle
     or a large part of it (at least half circle).

     It does not work well when data points are sampled along a small arc
     of a circle. In that case the method is heavily biased, it returns
     circles that are too often too small.

     It is the oldest algebraic circle fit (first published in 1972?).
     For 20-30 years it has been the most popular circle fit, at least
     until the more robust Pratt fit (1987) and Taubin fit (1991) were invented.

       Nikolai Chernov  (September 2012)
*/
{
  int i;

  double Xi, Yi, Zi;
  double Mxy, Mxx, Myy, Mxz, Myz;
  double B, C, G11, G12, G22, D1, D2;

  Circle circle;

  data.means(); // Compute x- and y- sample means (via a function in the class
                // "data")

  //     computing moments

  Mxx = Myy = Mxy = Mxz = Myz = 0.;

  for (i = 0; i < data.size(); i++) {
    Xi = data[i].x - data.mean.x; //  centered x-coordinates
    Yi = data[i].y - data.mean.y; //  centered y-coordinates
    Zi = Xi * Xi + Yi * Yi;

    Mxx += Xi * Xi;
    Myy += Yi * Yi;
    Mxy += Xi * Yi;
    Mxz += Xi * Zi;
    Myz += Yi * Zi;
  }
  Mxx /= data.size();
  Myy /= data.size();
  Mxy /= data.size();
  Mxz /= data.size();
  Myz /= data.size();

  //    solving system of equations by Cholesky factorization

  G11 = sqrt(Mxx);
  G12 = Mxy / G11;
  G22 = sqrt(Myy - G12 * G12);

  D1 = Mxz / G11;
  D2 = (Myz - D1 * G12) / G22;

  //    computing paramters of the fitting circle

  C = D2 / G22 / 2;
  B = (D1 - G12 * C) / G11 / 2;

  //       assembling the output

  circle.a = B + data.mean.x;
  circle.b = C + data.mean.y;
  circle.r = sqrt(B * B + C * C + Mxx + Myy);
  circle.s = Sigma(data, circle);
  circle.i = 0;
  circle.j = 0;

  return circle;
}
} // namespace CircleFittingModel
