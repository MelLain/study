#pragma once

#include <cmath>

#include <exception>
#include <iostream>
#include <string>

namespace DTS {

enum ProcType {
  GLOBAL_PROC,
  UPPER_PROC,
  CENTER_PROC,
  LOWER_PROC,
};

const double EPS = 1e-5;
const double RAND_CONST = 0.05;

const std::string OUT_VALUE_FILE = "OUT_VALUE_FILE.txt";
const std::string OUT_TRUE_FILE = "OUT_TRUE_FILE.txt";
const std::string OUT_ERROR_FILE = "OUT_ERROR_FILE.txt";

enum FlagType {
  START_ITER,
  TERMINATE,
};

}  // namespace DTS
