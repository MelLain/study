#pragma once

#include <cmath>

#include <exception>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

/*
Uncomment and set num threads to use OpenMP
#include "omp.h"
*/

namespace DTS {

struct ProcBounds {
  bool is_up;
  bool is_low;
  bool is_left;
  bool is_right;
};

const double EPS = 1e-4;
const double RAND_CONST = 0.05;

enum FlagType {
  START_ITER,
  TERMINATE,
};

}  // namespace DTS
