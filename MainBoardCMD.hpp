#pragma once
#include <iostream>
using namespace std;

#define COMMR_DEFAULT_START_BYTE 100u
#define COMMR_DEFAULT_TIMEOUT 30

typedef struct {
  float search_min_dist; // unit in mm
  float search_max_dist;
  float serach_min_range; // unit in degree
  float serach_max_range;
} RequestLidarCMD;

typedef union {
  char bytes[sizeof(RequestLidarCMD)];
  RequestLidarCMD cmd;
} BytewiseRequestLidarCMD;

static ostream &operator<<(ostream &cout, const BytewiseRequestLidarCMD &cmd) {
  cout << "Lidar CMD:" << endl;
  cout << "\tmin dist: " << cmd.cmd.search_min_dist << endl;
  cout << "\tmax dist: " << cmd.cmd.search_max_dist << endl;
  cout << "\tmin range: " << cmd.cmd.serach_min_range << endl;
  cout << "\tmax rnage: " << cmd.cmd.serach_max_range << endl;

  return cout;
}
