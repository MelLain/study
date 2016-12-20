#include "common.h"

namespace DTS {

__device__ double main_func_k(double r, double c) {
  return 8 - 12 * pow(r, 2) - 12 * pow(c, 2);
}

__device__ double bound_func_k(double r, double c) {
  return pow((1 - pow(r, 2)), 2) + pow((1 - pow(c, 2)), 2);
}

__device__ double true_func_k(double r, double c) {
  return pow((1 - pow(r, 2)), 2) + pow((1 - pow(c, 2)), 2);
}

__host__ __device__ double step_func_k(double value, double q) {
  return (pow(1 + value, q) - 1) / (pow(2.0, q) - 1);
}


__device__ double matrix_get(double* data, size_t row_index, size_t col_index, size_t num_cols) {
  return data[num_cols * row_index + col_index];
}

__device__ void matrix_set(double* data, size_t row_index, size_t col_index, double value, size_t num_cols) {
  data[num_cols * row_index + col_index] = value;
}

__device__ double grid_r_step(double* data_x, double* data_y, size_t row, size_t num_cols) {
  return row > 0 ? matrix_get(data_x, row, 0, num_cols) - matrix_get(data_x, row - 1, 0, num_cols) : 0.0;
}

__device__ double grid_c_step(double* data_x, double* data_y, size_t col, size_t num_cols) {
  return col > 0 ? matrix_get(data_y, 0, col, num_cols) - matrix_get(data_y, 0, col - 1, num_cols) : 0.0;
}

} // namespace DTS
