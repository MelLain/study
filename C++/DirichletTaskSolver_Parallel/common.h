#pragma once

#include <cmath>

#include <exception>
#include <iostream>
#include <iomanip>
#include <string>

namespace DTS {

struct ProcBounds {
  bool is_up;
  bool is_low;
  bool is_left;
  bool is_right;
};

const double EPS = 1e-4;
const double RAND_CONST = 0.05;

const std::string OUT_VALUE_FILE = "OUT_VALUE_FILE.txt";
const std::string OUT_TRUE_FILE = "OUT_TRUE_FILE.txt";

enum FlagType {
  START_ITER,
  TERMINATE,
};

}  // namespace DTS
