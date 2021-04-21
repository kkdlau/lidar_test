#include "Circle.hpp"
#include "PointArray.hpp"
#include "Statistics.hpp"
#include "helper_funcs.hpp"

namespace CircleFittingModel {

/**
 * @brief Levenberg Marquardt Full Circle Fitting Model.
 *
 * @param data array of points.
 * @param circleIni Initial circle, it will be used as initial guess.
 * @param LambdaIni The initial value of the control parameter "lambda" for the
 * Levenberg-Marquardt procedure. It is suggested to use a smaller number (e.g.
 * 0.001)
 * @param circle Output parameter. The guessed circle will be returned via this
 * parameter.
 * @return int 0 = successfully found a circle,
 * 1 = exceeding the limit of outer iteration,
 * 2 = exceeding the limit of inner iteration,
 * 3 = the points are divergence, too difficult to find a circle
 */
int levenbergMarquardtFullFitting(PointArray<double> &data, Circle &circleIni,
                                  double LambdaIni, Circle &circle)
/*
       Geometric circle fit to a given set of data points (in 2D)

       Output:
               integer function value is a code:
                          0:  normal termination, the best fitting circle is
                              successfully found
                          1:  the number of outer iterations exceeds the limit
   (99) (indicator of a possible divergence) 2:  the number of inner iterations
   exceeds the limit (99) (another indicator of a possible divergence) 3:  the
   coordinates of the center are too large (a strong indicator of divergence)

               circle - parameters of the fitting circle ("best fit")

               circle.a - the X-coordinate of the center of the fitting circle
               circle.b - the Y-coordinate of the center of the fitting circle
               circle.r - the radius of the fitting circle
               circle.s - the root mean square error (the estimate of sigma)
               circle.i - the total number of outer iterations (updating the
   parameters) circle.j - the total number of inner iterations (adjusting
   lambda)

       Algorithm:  Levenberg-Marquardt running over the full parameter space
   (a,b,r)

       See a detailed description in Section 4.5 of the book by Nikolai Chernov:
       "Circular and linear regression: Fitting circles and lines by least
   squares" Chapman & Hall/CRC, Monographs on Statistics and Applied
   Probability, volume 117, 2010.

                Nikolai Chernov,  February 2014
*/
{
  int code, i, iter, inner, IterMAX = 99;

  double factorUp = 10., factorDown = 0.04, lambda, ParLimit = 1.e+6;
  double dx, dy, ri, u, v;
  double Mu, Mv, Muu, Mvv, Muv, Mr, UUl, VVl, Nl, F1, F2, F3, dX, dY, dR;
  double epsilon = 3.e-8;
  double G11, G22, G33, G12, G13, G23, D1, D2, D3;

  Circle Old, New;

  //       starting with the given initial circle (initial guess)

  New = circleIni;

  //       compute the root-mean-square error via function Sigma; see
  //       Utilities.cpp

  New.s = Sigma(data, New);

  //       initializing lambda, iteration counters, and the exit code

  lambda = LambdaIni;
  iter = inner = code = 0;

NextIteration:

  Old = New;
  if (++iter > IterMAX) {
    code = 1;
    goto enough;
  }

  //       computing moments

  Mu = Mv = Muu = Mvv = Muv = Mr = 0.;

  for (i = 0; i < data.size(); i++) {
    dx = data[i].x - Old.a;
    dy = data[i].y - Old.b;
    ri = sqrt(dx * dx + dy * dy);
    u = dx / ri;
    v = dy / ri;
    Mu += u;
    Mv += v;
    Muu += u * u;
    Mvv += v * v;
    Muv += u * v;
    Mr += ri;
  }
  Mu /= data.size();
  Mv /= data.size();
  Muu /= data.size();
  Mvv /= data.size();
  Muv /= data.size();
  Mr /= data.size();

  //       computing matrices

  F1 = Old.a + Old.r * Mu - data.mean.x;
  F2 = Old.b + Old.r * Mv - data.mean.y;
  F3 = Old.r - Mr;

  Old.g = New.g = sqrt(F1 * F1 + F2 * F2 + F3 * F3);

try_again:

  UUl = Muu + lambda;
  VVl = Mvv + lambda;
  Nl = 1 + lambda;

  //         Cholesly decomposition

  G11 = sqrt(UUl);
  G12 = Muv / G11;
  G13 = Mu / G11;
  G22 = sqrt(VVl - G12 * G12);
  G23 = (Mv - G12 * G13) / G22;
  G33 = sqrt(Nl - G13 * G13 - G23 * G23);

  D1 = F1 / G11;
  D2 = (F2 - G12 * D1) / G22;
  D3 = (F3 - G13 * D1 - G23 * D2) / G33;

  dR = D3 / G33;
  dY = (D2 - G23 * dR) / G22;
  dX = (D1 - G12 * dY - G13 * dR) / G11;

  if ((abs(dR) + abs(dX) + abs(dY)) / (1 + Old.r) < epsilon)
    goto enough;

  //       updating the parameters

  New.a = Old.a - dX;
  New.b = Old.b - dY;

  if (abs(New.a) > ParLimit || abs(New.b) > ParLimit) {
    code = 3;
    goto enough;
  }

  New.r = Old.r - dR;

  if (New.r <= 0.) {
    lambda *= factorUp;
    if (++inner > IterMAX) {
      code = 2;
      goto enough;
    }
    goto try_again;
  }

  //       compute the root-mean-square error via function Sigma; see
  //       Utilities.cpp

  New.s = Sigma(data, New);

  //       check if improvement is gained

  if (New.s < Old.s) //   yes, improvement
  {
    lambda *= factorDown;
    goto NextIteration;
  } else //   no improvement
  {
    if (++inner > IterMAX) {
      code = 2;
      goto enough;
    }
    lambda *= factorUp;
    goto try_again;
  }

  //       exit

enough:

  Old.i = iter;  // total number of outer iterations (updating the parameters)
  Old.j = inner; // total number of inner iterations (adjusting lambda)

  circle = Old;

  return code;
}
} // namespace CircleFittingModel
