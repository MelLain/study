#pragma once

#include <cmath>

#include <exception>
#include <iostream>
#include <string>

namespace DTS {

const int ARG_PARSE_ERROR = -1;

const double EPS = 1e-20;
const double RAND_CONST = 0.05;

const std::string OUT_VALUE_FILE = "OUT_VALUE_FILE.txt";
const std::string OUT_TRUE_FILE = "OUT_TRUE_FILE.txt";
const std::string OUT_ERROR_FILE = "OUT_ERROR_FILE.txt";

}  // namespace DTS