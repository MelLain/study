#include <cmath>
#include <string>
#include <utility>

#include "mpi_helpers.h"
#include "gradient_descent.h"

namespace DTS {

GradientDescent::GradientDescent(const GridData& grid_data, const ProcBounds& proc_bounds,
                                 size_t num_procs, size_t proc_rank, size_t num_row_procs, size_t num_points,
                                 size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx)
  : proc_bounds_(proc_bounds)
  , num_procs_(num_procs)
  , proc_rank_(proc_rank)
  , num_points_(num_points)
  , first_send_(std::make_pair((proc_rank % num_row_procs) % 2 == 1, (proc_rank / num_row_procs) % 2 == 0))
  , left_right_proc_(std::make_pair(proc_bounds.is_left ? -1 : proc_rank - num_row_procs,
                                    proc_bounds.is_right ? -1 : proc_rank + num_row_procs))
  , num_rows_(0)
  , num_cols_(0)
{
  CheckCuda(cudaSetDevice(0));

  InitGrid_CUDA(&grid_x_cuda_, &grid_y_cuda_, &num_rows_, &num_cols_, grid_data, proc_bounds,
                start_row_idx, end_row_idx, start_col_idx, end_col_idx);
  InitValues_CUDA(&values_cuda_, grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_, proc_bounds);

  InitMatrix_CUDA(&values_laplass_cuda_, num_rows_, num_cols_);
  InitMatrix_CUDA(&old_values_cuda_, num_rows_, num_cols_);
  InitMatrix_CUDA(&residuals_cuda_, num_rows_, num_cols_);
  InitMatrix_CUDA(&residuals_laplass_cuda_, num_rows_, num_cols_);
  InitMatrix_CUDA(&gradients_cuda_, num_rows_, num_cols_);
  InitMatrix_CUDA(&gradients_laplass_cuda_, num_rows_, num_cols_);
  InitMatrix_CUDA(&temp_matrix_cuda_, num_rows_, num_cols_);
}

std::pair<size_t, double> GradientDescent::FitModel() {
  size_t processed_iter = 0;
  double error = 1e+10;
  while (true) {
    FlagType flag;
    if (proc_rank_ == 0) {
      send_flag_to_all(num_procs_, START_ITER);
    } else {
      receive_flag(&flag, 0, proc_rank_);
      if (flag == TERMINATE) {
        double pre_err = CountPreError_CUDA(values_cuda_, temp_matrix_cuda_,
                                            grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
        send_value(pre_err, 0, proc_rank_);
        break;
      }
    }

    ExchangeRowsCols_CUDA(values_cuda_, num_rows_, num_cols_, proc_bounds_, first_send_, left_right_proc_, proc_rank_);

    // step 1: count residuals
    FivePointsLaplass_CUDA(values_laplass_cuda_, values_cuda_, grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
    ExchangeRowsCols_CUDA(values_laplass_cuda_, num_rows_, num_cols_, proc_bounds_, first_send_, left_right_proc_, proc_rank_);

    CountResiduals_CUDA(residuals_cuda_, values_laplass_cuda_, grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_, proc_bounds_);

    FivePointsLaplass_CUDA(residuals_laplass_cuda_, residuals_cuda_, grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
    ExchangeRowsCols_CUDA(residuals_laplass_cuda_, num_rows_, num_cols_, proc_bounds_, first_send_, left_right_proc_, proc_rank_);

    // step 2: count alpha
    double alpha_den = ProductByPointAndSum_CUDA(gradients_laplass_cuda_, gradients_cuda_,
                                                 grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
    double alpha_nom = ProductByPointAndSum_CUDA(residuals_laplass_cuda_, gradients_cuda_,
                                                 grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);

    if (proc_rank_ == 0) {
      alpha_den = collect_value_from_all(num_procs_) + alpha_den;
      send_value_to_all(num_procs_, alpha_den);
      alpha_nom = collect_value_from_all(num_procs_) + alpha_nom;
      send_value_to_all(num_procs_, alpha_nom);
    } else {
      send_value(alpha_den, 0, proc_rank_);
      receive_value(&alpha_den, 0, proc_rank_);
      send_value(alpha_nom, 0, proc_rank_);
      receive_value(&alpha_nom, 0, proc_rank_);
    }
    double alpha = alpha_den > 0.0 ? alpha_nom / alpha_den : 0.0;

    // step 3: count new gradients
    CountGradients_CUDA(gradients_cuda_, residuals_cuda_, temp_matrix_cuda_,
                        alpha, num_rows_, num_cols_, processed_iter == 0); 
    processed_iter++;

    // step 4: count new gradients laplass
    FivePointsLaplass_CUDA(gradients_laplass_cuda_, gradients_cuda_,
                           grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);

    // step 5: count new tau
    double tau_den = ProductByPointAndSum_CUDA(gradients_laplass_cuda_, gradients_cuda_,
                                               grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
    double tau_nom = ProductByPointAndSum_CUDA(residuals_cuda_, gradients_cuda_,
                                               grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
    if (proc_rank_ == 0) {
      tau_den = collect_value_from_all(num_procs_) + tau_den;
      send_value_to_all(num_procs_, tau_den);
      tau_nom = collect_value_from_all(num_procs_) + tau_nom;
      send_value_to_all(num_procs_, tau_nom);
    } else {
      send_value(tau_den, 0, proc_rank_);
      receive_value(&tau_den, 0, proc_rank_);
      send_value(tau_nom, 0, proc_rank_);
      receive_value(&tau_nom, 0, proc_rank_);
    }
    double tau = tau_den > 0.0 ? tau_nom / tau_den : 0.0;

    // step 6: count new values and process error
    CountNewValues_CUDA(values_cuda_, old_values_cuda_, gradients_cuda_, temp_matrix_cuda_, tau, num_rows_, num_cols_);
    double val_diff = CountValuesDifference_CUDA(values_cuda_, old_values_cuda_, temp_matrix_cuda_,
                                            grid_x_cuda_, grid_y_cuda_, num_rows_, num_cols_);
    if (proc_rank_ == 0) {
      double difference = sqrt(collect_value_from_all(num_procs_) + val_diff);
      if (difference < EPS) {
        send_flag_to_all(num_procs_, TERMINATE);
        double pre_error = CountPreError_CUDA(values_cuda_, temp_matrix_cuda_, grid_x_cuda_,
                                               grid_y_cuda_, num_rows_, num_cols_);
        error = sqrt(collect_value_from_all(num_procs_) + pre_error);
        break;
      }      
    } else {
      send_value(val_diff, 0, proc_rank_);
    }
  }
  return std::make_pair(processed_iter, error);
}

/*
void GradientDescent::save_results_file() {
  std::ofstream out_value_file;
  std::ofstream out_true_file;

  size_t start_row_shift = proc_bounds_.is_up ? 0 : 1;
  size_t end_row_shift = proc_bounds_.is_low ? 0 : 1;
  size_t start_col_shift = proc_bounds_.is_left ? 0 : 1;
  size_t end_col_shift = proc_bounds_.is_right ? 0 : 1;

  out_value_file.open("VALUE_PART_POINTS_" + std::to_string(num_points_) + "_PROC_" + std::to_string(proc_rank_));
  out_true_file.open("TRUE_PART_POINTS_" + std::to_string(num_points_) + "_PROC_" + std::to_string(proc_rank_));

  for (size_t i = start_row_shift; i < values_->num_rows() - end_row_shift; ++i) {
    for (size_t j = start_col_shift; j < values_->num_cols() - end_col_shift; ++j) {
      out_value_file << (*values_)(i, j) << ", ";
      out_true_file << true_func((*grid_)(i, j)) << ", ";                                                                                       
    }
    out_value_file << std::endl;
    out_true_file << std::endl;
  }

  out_value_file.close();
  out_true_file.close();
}*/

}  // namespace DTS
