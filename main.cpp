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
using namespace rp::standalone::rplidar;

#ifndef _countof
#define _countof(_Array) (int)(sizeof(_Array) / sizeof(_Array[0]))
#endif

// #define USE_MOCK_DATA 1

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
    serial.writeChar(sizeof(BytewiseEstCricle) * 8);

    for (int i = 0; i < 8; ++i) {
      const BytewiseEstCricle &c = cs[i];

      cout << "send a circle" << endl;
      serial.writeBytes(c.bytes, sizeof(BytewiseEstCricle));
      usleep(COMMR_DEFAULT_TIMEOUT / 2);
    }
  } else {
    serial.writeChar(sizeof(BytewiseEstCricle) * cs.size()); // num_bytes
    for (const BytewiseEstCricle &c : cs) {
      cout << "send a circle" << endl;
      serial.writeBytes(c.bytes, sizeof(BytewiseEstCricle));
      usleep(COMMR_DEFAULT_TIMEOUT / 2);
    }
  }
}

const char *opt_com_path = NULL;
_u32 baudrateArray[2] = {115200, 256000};
_u32 opt_com_baudrate = 0;
u_result op_result;
bool useArgcBaudrate = false;
RPlidarDriver *drv = nullptr;

void my_handler(int s) {
#ifndef USE_MOCK_DATA
  drv->stop();
  drv->stopMotor();
  RPlidarDriver::DisposeDriver(drv);
  drv = NULL;
  cout << "Ctrl + C Event received, stop LiDAR." << endl;
#endif
  exit(1);
}

void connect_lidar(RPlidarDriver *&drv) {
  if (!drv) {
    fprintf(stderr, "insufficent memory, exit\n");
    exit(-2);
  }
  rplidar_response_device_info_t devinfo;
  bool connectSuccess = false;

  // make connection...
  if (useArgcBaudrate) {
    if (!drv)
      drv = RPlidarDriver::CreateDriver(DRIVER_TYPE_SERIALPORT);
    if (IS_OK(drv->connect(opt_com_path, opt_com_baudrate))) {
      op_result = drv->getDeviceInfo(devinfo);

      if (IS_OK(op_result)) {
        connectSuccess = true;
      } else {
        delete drv;
        drv = NULL;
      }
    }
  } else {
    size_t baudRateArraySize =
        (sizeof(baudrateArray)) / (sizeof(baudrateArray[0]));
    for (size_t i = 0; i < baudRateArraySize; ++i) {
      if (!drv)
        drv = RPlidarDriver::CreateDriver(DRIVER_TYPE_SERIALPORT);
      if (IS_OK(drv->connect(opt_com_path, baudrateArray[i]))) {
        op_result = drv->getDeviceInfo(devinfo);

        if (IS_OK(op_result)) {
          connectSuccess = true;
          break;
        } else {
          delete drv;
          drv = NULL;
        }
      }
    }
  }
  if (!connectSuccess) {

    fprintf(stderr, "Error, cannot bind to the specified serial port %s.\n",
            opt_com_path);
    exit(1);
  }

  // print out the device serial number, firmware and hardware version number..
  printf("RPLIDAR S/N: ");
  for (int pos = 0; pos < 16; ++pos) {
    printf("%02X", devinfo.serialnum[pos]);
  }

  printf("\n"
         "Firmware Ver: %d.%02d\n"
         "Hardware Rev: %d\n",
         devinfo.firmware_version >> 8, devinfo.firmware_version & 0xFF,
         (int)devinfo.hardware_version);
}

bool checkRPLIDARHealth(RPlidarDriver *drv) {
  u_result op_result;
  rplidar_response_device_health_t healthinfo;

  op_result = drv->getHealth(healthinfo);
  if (IS_OK(op_result)) { // the macro IS_OK is the preperred way to judge
                          // whether the operation is succeed.
    printf("RPLidar health status : %d\n", healthinfo.status);
    if (healthinfo.status == RPLIDAR_STATUS_ERROR) {
      fprintf(stderr, "Error, rplidar internal error detected. Please reboot "
                      "the device to retry.\n");
      // enable the following code if you want rplidar to be reboot by software
      // drv->reset();
      return false;
    } else {
      return true;
    }

  } else {
    fprintf(stderr, "Error, cannot retrieve the lidar health code: %x\n",
            op_result);
    return false;
  }
}

int main(int argc, const char *argv[]) {

  signal(SIGINT, my_handler);
  serialib serial;
  char success = serial.openDevice("/dev/tty.usbserial-1430", 115200);

  if (!success) {
    cout << "failed to connect TTL" << endl;
    return 1;
  }

  drv = RPlidarDriver::CreateDriver(DRIVER_TYPE_SERIALPORT);

  if (!opt_com_path) {
#ifdef _WIN32
    // use default com port
    opt_com_path = "\\\\.\\com57";
#elif __APPLE__
    opt_com_path = "/dev/tty.usbserial-0001";
#else
    opt_com_path = "/dev/ttyUSB0";
#endif
  }
  vector<PolarVector> pvs;
  PointArray<double> pts;

#ifndef USE_MOCK_DATA

  connect_lidar(drv);

  if (!checkRPLIDARHealth(drv)) {
    cout << "Lidar health check failed" << endl;
    exit(1);
  }

  drv->startMotor();
  drv->setMotorPWM(800);
  // start scan...
  drv->startScan(0, 1);

#else
  ifstream f;
  f.open("Sample3/data.csv", ios::in);
  if (!f.is_open())
    return 1;

  string line;

  getline(f, line); // reading the header

  vector<string> tokens;
  while (getline(f, line)) {
    tokens.push_back(line);
  }

  f.close();

#endif

  while (true) {
    BytewiseRequestLidarCMD cmd = {.cmd = {.search_max_dist = 7000,
                                           .search_min_dist = 2000,
                                           .serach_min_range = 0,
                                           .serach_max_range = 180}};
    // bool ok = get_mb_cmd(serial, cmd); // get read lidar cmd from main board

    // if (!ok)
    //   continue;

    cout << cmd << endl;
#define within(d, _min, _max) (d >= _min && d <= _max)
#ifndef USE_MOCK_DATA
    rplidar_response_measurement_node_hq_t nodes[8192];
    size_t count = _countof(nodes);

    op_result = drv->grabScanDataHq(nodes, count);

    if (IS_OK(op_result)) {
      drv->ascendScanData(nodes, count);
      ofstream plotFile;
      plotFile.open("plot.csv");

      plotFile << "angle,distance" << endl;

      for (int pos = 0; pos < (int)count; ++pos) {

        double dist = nodes[pos].dist_mm_q2 / 4.0f;
        double angle = nodes[pos].angle_z_q14 * 90.f / (1 << 14);
        double r = angle / 360 * 2 * M_PI;
        if (angle > 180)
          break;

        if (!within(dist, cmd.cmd.search_min_dist, cmd.cmd.search_max_dist) ||
            !within(angle, cmd.cmd.serach_min_range,
                    cmd.cmd.serach_max_range)) {
          dist = 0;
        }

        PolarVector pv = PolarVector::inDegree(angle, dist);

        plotFile << pv.toCSVString() << endl;

        pvs.push_back(pv);
        pts.push_back(pv.toCartesian());
      }
      plotFile.close();
    }
#else
    for (string &s : tokens) {
      vector<string> splitted;

      split(s, splitted, ",");
      float angle = stof(splitted[0]);
      float dist = stof(splitted[1]);
      if (!within(dist, cmd.cmd.search_min_dist, cmd.cmd.search_max_dist) ||
          !within(angle, cmd.cmd.serach_min_range, cmd.cmd.serach_max_range)) {
        dist = 0;
      }

      PolarVector pv = PolarVector::inDegree(angle, dist);
      pvs.push_back(pv);
      pts.push_back(pv.toCartesian());
    }
#endif

    vector<int> min = Filter::localMinimum(pts, 0.004);
    vector<BytewiseEstCricle> cs = {};

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
      cs.push_back(best.toSendableCircle());
    }

    for (BytewiseEstCricle &c : cs) {
      c.c.center.x *= -1; // flip the sign
    }
    transmit_circle(serial, cs);

    pvs.clear();
    pts.clear();
  }
  return 0;
}
