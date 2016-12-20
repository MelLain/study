#pragma once

#include <cmath>

#include <stdexcept>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <mpi.h>
#include <cuda_runtime.h>

namespace DTS {

#define CheckCuda(code) {\
  if (code != cudaSuccess) {\
    std::stringstream ss;\
    ss << __FILE__ << ", line " << __LINE__ << " : " << cudaGetErrorString(code);\
    throw std::runtime_error(ss.str().c_str());\
  }\
}\

__device__ double main_func_k(double r, double c);
__device__ double bound_func_k(double r, double c);
__device__ double true_func_k(double r, double c);
__host__ __device__ double step_func_k(double value, double q);

__device__ double matrix_get(double* data, size_t row_index, size_t col_index, size_t num_cols);
__device__ void matrix_set(double* data, size_t row_index, size_t col_index, double value, size_t num_cols);

__device__ double grid_r_step(double* data_x, double* data_y, size_t row, size_t num_cols);
__device__ double grid_c_step(double* data_x, double* data_y, size_t col, size_t num_cols);

struct GridData {
  double r_lower_bound;
  double r_upper_bound;
  size_t r_num_points;

  double c_lower_bound;
  double c_upper_bound;
  size_t c_num_points;

  double q;
};

struct ProcBounds {
  ProcBounds(bool q, bool w, bool e, bool r): is_up(q), is_low(w), is_left(e), is_right(r) { }
  
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

} // namespace DTS
