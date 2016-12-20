#pragma once

#include <utility>

#include "common.h"
#include "cuda_operations.h"

namespace DTS {

class GradientDescent {
 public:
  GradientDescent(const GridData& grid_data, const ProcBounds& proc_bounds,
                  size_t num_procs, size_t proc_rank, size_t num_row_procs, size_t num_points,
                  size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx);

  ~GradientDescent() { 
    Clear_CUDA(grid_x_cuda_);
    Clear_CUDA(grid_y_cuda_);
    Clear_CUDA(values_cuda_);
    Clear_CUDA(values_laplass_cuda_);
    Clear_CUDA(old_values_cuda_);
    Clear_CUDA(residuals_cuda_);
    Clear_CUDA(residuals_laplass_cuda_);
    Clear_CUDA(gradients_cuda_);
    Clear_CUDA(gradients_laplass_cuda_);
    Clear_CUDA(temp_matrix_cuda_);
  }

  std::pair<size_t, double> FitModel();

 private:
  //save_results_file();

  ProcBounds proc_bounds_;
  size_t num_procs_;
  size_t proc_rank_;
  size_t num_points_;

  std::pair<bool, bool> first_send_;
  std::pair<int, int> left_right_proc_;

  size_t num_rows_;
  size_t num_cols_;

  // cuda pointers to all matrices
  double* grid_x_cuda_;
  double* grid_y_cuda_;
  double* values_cuda_;
  double* values_laplass_cuda_;
  double* old_values_cuda_;
  double* residuals_cuda_;
  double* residuals_laplass_cuda_;
  double* gradients_cuda_;
  double* gradients_laplass_cuda_;
  double* temp_matrix_cuda_;
};

}  // namespace DTS
