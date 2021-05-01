#pragma once
#include "serialib.hpp"
#include <string>

using namespace std;

void receiver_routine(void *args) {
  string device_path = *(string *)args;

  while (true) {
    /* code */
  }
}