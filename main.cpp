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
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//
#include "Circle.hpp"
#include "CircleFitByKasa.hpp"
#include "CircleFitByLevenbergMarquardtFull.hpp"
#include "CircleFitByTaubin.hpp"
#include "MainBoardCMD.hpp"
#include "PointArray.hpp"
#include "filter.hpp"
//
#include "serialib.cpp"
#include "serialib.hpp"
//
#include "sdk/include/rplidar.h"
#include <signal.h>
#include <string>
#include <vector>
using namespace std;

// read this:
// https://lucidar.me/en/serialib/cross-plateform-rs232-serial-library/

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

bool get_mb_cmd(serialib &serial, BytewiseRequestLidarCMD &cmd) {
  char tmp;
  char success = 0;
  do {
    serial.readChar(&tmp);
  } while (tmp != COMMR_DEFAULT_START_BYTE);

  success = serial.readChar(&tmp, COMMR_DEFAULT_TIMEOUT);

  if (!success) {
    cout << "failed to read num_byte" << endl;
    return false;
  }

  short num_btye = tmp;
  cout << "start receiving bytes from main board - num btye: " << num_btye
       << endl;

  for (int i = 0; i < num_btye; ++i) {
    success = serial.readChar(&cmd.bytes[i], COMMR_DEFAULT_TIMEOUT);
    if (!success) {
      cout << "failed to read data segment." << endl;
      return false;
    }
  }

  return true;
}

void transmit_circle(serialib &serial, const vector<BytewiseEstCricle> cs) {
  serial.writeChar(COMMR_DEFAULT_START_BYTE);
  if (cs.size() > 8) {
    // special case: msg is too large
  } else {
    serial.writeChar(sizeof(BytewiseEstCricle) * cs.size()); // num_bytes
    for (const BytewiseEstCricle &c : cs) {
      serial.writeBytes(c.bytes, sizeof(BytewiseEstCricle));
      usleep(COMMR_DEFAULT_TIMEOUT / 2);
    }
  }
}

int main(int argc, const char *argv[]) {

  serialib serial;
  char success = serial.openDevice("/dev/tty.usbserial-1420", 115200);

  if (!success) {
    cout << "failed to connect TTL" << endl;
    return 1;
  }

  while (true) {
    // BytewiseRequestLidarCMD cmd;
    // bool ok = get_mb_cmd(serial, cmd); // get read lidar cmd from main board

    // if (!ok)
    //   continue;

    // cout << cmd << endl;

    vector<BytewiseEstCricle> cs = {};

    char tmp;
    cin >> tmp;

    Circle c(300, 700, 1000);
    Circle c2(1000, 2000, 3000);
    Circle c3(1000, 2000, 3000);
    Circle c4(1000, 2000, 3000);
    Circle c5(1000, 2000, 3000);
    cs.push_back(c.toSendableCircle());
    cs.push_back(c2.toSendableCircle());
    cs.push_back(c3.toSendableCircle());
    cs.push_back(c4.toSendableCircle());
    cs.push_back(c5.toSendableCircle());

    transmit_circle(serial, cs);
  }

  ifstream f;
  f.open("Sample1/data.csv", ios::in);
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
    if (dist < 3000 || dist > 6000)
      dist = 0;
    PolarVector pv = PolarVector::inDegree(stof(tokens[0]), dist);
    cout << pv.toCSVString() << endl;

    // cout << pv.toString() << endl;
    pvs.push_back(pv);
    pts.push_back(pv.toCartesian());
  }

  vector<int> min = Filter::localMinimum(pts, 0.004);

  cout << "found local minimum:" << endl;
  for (int m : min) {
    cout << m << ": " << pvs[m].toString() << endl;

    Circle best;

    best.r = 10000000;

    for (int s = 3; s < 50; s += 1) {
      PointArray<double> circleData = partition(pts, m, s);
      Circle c = CircleFittingModel::taubinFitting(circleData);
      CircleFittingModel::levenbergMarquardtFullFitting(circleData, c, 0.0001,
                                                        c);
      if (abs(160 - c.r) < abs(160 - best.r)) {
        best = c;
      }
    }

    best.print();
  }

  f.close();
  return 0;
}
