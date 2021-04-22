/*
 *  RPLIDAR
 *  Ultra Simple Data Grabber Demo App
 *
 *  Copyright (c) 2009 - 2014 RoboPeak Team
 *  http://www.robopeak.com
 *  Copyright (c) 2014 - 2019 Shanghai Slamtec Co., Ltd.
 *  http://www.slamtec.com
 *
 */
/*
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "Circle.hpp"
#include "CircleFitByKasa.hpp"
#include "CircleFitByLevenbergMarquardtFull.hpp"
#include "CircleFitByTaubin.hpp"
#include "PointArray.hpp"
#include "filter.hpp"

#include <fstream>
#include <signal.h>
#include <string>
#include <vector>
using namespace std;

template <typename T>
PointArray<T> partition(PointArray<T> arr, int index,
                        unsigned simpleSize = 20) {
  PointArray<T> ret;
  const unsigned half = simpleSize << 1;

  for (int i = max(0u, (unsigned)(index - half));
       i < min((unsigned)arr.size(), (unsigned)(index + half)); ++i) {
    if (PolarVector::fromCartesian(arr[i]).mag == 0)
      continue;
    ret.push_back(arr[i]);
  }

  return ret;
}

void split(const string &s, vector<string> &tokens,
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

int main(int argc, const char *argv[]) {
  ifstream f;
  f.open("data.csv", ios::in);
  if (!f.is_open())
    return 1;
  vector<PolarVector> pvs;
  PointArray<double> pts;

  string line;

  getline(f, line); // reading the header

  while (getline(f, line)) {
    vector<string> tokens;
    split(line, tokens, ",");
    float dist = stof(tokens[1]);
    if (dist < 500 || dist > 1700)
      dist = 0;
    PolarVector pv = PolarVector::inDegree(stof(tokens[0]), dist);
    cout << pv.toCSVString() << endl;

    // cout << pv.toString() << endl;
    pvs.push_back(pv);
    pts.push_back(pv.toCartesian());
  }

  vector<int> min = Filter::localMinimum(pts, 0.03);

  cout << "found local minimum:" << endl;
  for (int m : min) {
    cout << m << ": " << pvs[m].toString() << endl;

    PointArray<double> circleData = partition(pts, m);

    for (const auto &p : circleData) {
      cout << PolarVector::fromCartesian(p).toString() << endl;
    }

    cout << "circle approximation:" << endl;

    cout << "\napproximation by taubinFitting:" << endl;
    Circle c = CircleFittingModel::taubinFitting(circleData);

    c.print();

    cout << "\napproximation by levenbergMarquardtFullFitting:" << endl;
    CircleFittingModel::levenbergMarquardtFullFitting(circleData, c, 0.0001, c);

    c.print();
  }

  f.close();
  return 0;
}
