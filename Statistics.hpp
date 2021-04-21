#pragma once

#include "Circle.hpp"
#include "PointArray.hpp"

#include "helper_funcs.hpp"

static void RandomNormalPair(double &x, double &y)
/*
    Generator of pseudo-random numbers
    with standard normal distribution
    based on Box-Muller transformation

    "double" can be replaced by "float",
    "double", or "long double"; or it
    can be predefined as 1 of these types

    Input:  n1
    Output:  2 real values, x and y,
    that are random, independent, and
    have standard normal distribution
    (with mean 0 and variance 1)

    Call:

        RandomNormalPair(x,y);

    Uses standard C++ random generator rand()

    To reseed the generator rand(), call

        srand ( (unsigned)time(NULL) );

    before you start calling this function

       Nikolai Chernov, November 2011
*/
{
  double rand1, rand2, wrand;
  /*
  //       version 1, by direct calculation (slower)

      double pi=3.141592653589793;
      rand1 = (double)rand()/RAND_MAX;
      rand2 = (double)rand()/RAND_MAX;
      x = sqrt(-2*log(rand1))*cos(2*pi*rand2);
      y = sqrt(-2*log(rand1))*sin(2*pi*rand2);
  */
  //       version 2, in polar form
  //         (faster and more stable numerically)

  do {
    rand1 = 2 * rand() / RAND_MAX - 1;
    rand2 = 2 * rand() / RAND_MAX - 1;
    wrand = rand1 * rand1 + rand2 * rand2;
  } while (wrand >= 1);
  wrand = sqrt((-2 * log(wrand)) / wrand);
  x = rand1 * wrand;
  y = rand2 * wrand;
}

//*********************** SimulateArc ******************************

static void SimulateArc(PointArray<double> &data, double a, double b, double R,
                        double theta1, double theta2, double sigma)
/*
          Simulate data points equally spaced along a circular arc with Gaussian
  noise

  input:
          a,b         the coordinates of the center of the circle
          R           the radius of circle
          theta1      first  endpoint of the arc (in radians)
          theta2      second endpoint of the arc (in radians)
          sigma       noise level (standard deviation of residuals)
*/
{
  int N = data.size();
  double theta, dx, dy;

  for (int i = 0; i < N; i++) {
    theta = theta1 + (theta2 - theta1) * i / (N - 1);

    //			isotropic Gaussian noise

    RandomNormalPair(dx, dy);
    data[i].x = a + R * cos(theta) + sigma * dx;
    data[i].y = b + R * sin(theta) + sigma * dy;
  }
}

//********************* SimulateRandom ****************************

static void SimulateRandom(PointArray<double> &data, double Window)
/*
          Simulate data points with uniform distribution
          in the square |x|<Window, |y|<Window

  input:
          nPoints  the number of data points
*/
{
  // PointArray<double> data(nPoints);

  for (int i = 0; i < data.size(); i++) {
    data[i].x = Window * (2 * rand() / RAND_MAX - 1);
    data[i].y = Window * (2 * rand() / RAND_MAX - 1);
  }
}

//****************** Sigma ************************************
//
//   estimate of Sigma = square root of RSS divided by N
//   gives the root-mean-square error of the geometric circle fit

static double Sigma(PointArray<double> &data, Circle &circle) {
  double sum = 0., dx, dy;

  for (int i = 0; i < data.size(); i++) {
    dx = data[i].x - circle.a;
    dy = data[i].y - circle.b;
    sum += SQR(sqrt(dx * dx + dy * dy) - circle.r);
  }
  return sqrt(sum / data.size());
}

//****************** SigmaReduced ************************************
//
//   estimate of Sigma = square root of RSS divided by N
//   gives the root-mean-square error of the geometric circle fit
//
//   uses only the center of the circle (a,b), not the radius
//   the function computes the optimal radius here

static double SigmaReduced(PointArray<double> &data, Circle &circle) {
  int i, n = data.size();
  double sum = 0., dx, dy, r, D[n];

  for (i = 0; i < n; i++) {
    dx = data[i].x - circle.a;
    dy = data[i].y - circle.b;
    D[i] = sqrt(dx * dx + dy * dy);
    sum += D[i];
  }
  r = sum / n;

  for (sum = 0., i = 0; i < n; i++)
    sum += SQR(D[i] - r);

  return sqrt(sum / n);
}

//****************** SigmaReducedNearLinearCase ****************
//
//   estimate of Sigma = square root of RSS divided by N

static double SigmaReducedNearLinearCase(PointArray<double> &data,
                                         Circle &circle) {
  int i, n = data.size();
  double a0, b0, del, s, c, x, y, z, p, t, g, W, Z;

  a0 = circle.a - data.mean.x;
  b0 = circle.b - data.mean.y;
  del = 1 / sqrt(a0 * a0 + b0 * b0);
  s = b0 * del;
  c = a0 * del;

  for (W = Z = 0., i = 0; i < n; i++) {
    x = data[i].x - data.mean.x;
    y = data[i].y - data.mean.y;
    z = x * x + y * y;
    p = x * c + y * s;
    t = del * z - 2 * p;
    g = t / (1 + sqrt(1 + del * t));
    W += (z + p * g) / (2 + del * g);
    Z += z;
  }
  W /= n;
  Z /= n;

  return sqrt(Z - W * (2 + del * del * W));
}

//****************** SigmaReducedForCenteredScaled ****************
//
//   estimate of Sigma = square root of RSS divided by N

static double SigmaReducedForCenteredScaled(PointArray<double> &data,
                                            Circle &circle) {
  int i, n = data.size();
  double sum = 0., dx, dy, r;

  for (i = 0; i < n; i++) {
    dx = data[i].x - circle.a;
    dy = data[i].y - circle.b;
    sum += sqrt(dx * dx + dy * dy);
  }
  r = sum / n;

  return sqrt(SQR(circle.a) + SQR(circle.b) - r * r + 2);
}

//****************** OptimalRadius ******************************
//
//     compute the optimal radius of a circle, given its center (a,b)

static double OptimalRadius(PointArray<double> &data, Circle &circle) {
  double Mr = 0., dx, dy;

  for (int i = 0; i < data.size(); i++) {
    dx = data[i].x - circle.a;
    dy = data[i].y - circle.b;
    Mr += sqrt(dx * dx + dy * dy);
  }
  return Mr / data.size();
}
