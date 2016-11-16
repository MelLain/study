#include <cmath>
#include <string>
#include <utility>

#include "matrix_operations.h"
#include "mpi_helpers.h"
#include "gradient_descent.h"

namespace DTS {

GradientDescent::GradientDescent(const GridData& grid_data, const Functions& functions, const ProcBounds& proc_bounds,
                                 size_t proc_rank, size_t num_row_procs, size_t num_points,
                                 size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx)
  : functions_(functions)
  , proc_bounds_(proc_bounds)
  , proc_rank_(proc_rank)
  , first_send_(std::make_pair((proc_rank % num_row_procs) % 2 == 1, (proc_rank / num_row_procs) % 2 == 0))
  , left_right_proc_(std::make_pair(proc_bounds.is_left ? 0 : proc_rank - num_row_procs, proc_bounds.is_right ? 0 : proc_rank + num_row_procs))
  , num_points_(num_points)
{
  clear();
  init_grid(grid_data, start_row_idx, end_row_idx, start_col_idx, end_col_idx);
  init_values();

  size_t w = grid_->num_rows();
  size_t h = grid_->num_cols();

  old_values_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  gradients_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  gradients_laplass_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
}

void GradientDescent::FitModel() {
  bool first_iter = true;
  while (true) {
    FlagType flag;
    receive_flag(&flag, 0, proc_rank_);

    if (flag == TERMINATE) {
      send_value(count_pre_error(), 0, proc_rank_);
      save_results_file();
      break;
    }

    exchange_mirror_rows(values_);

    // step 1: count residuals
    auto residuals = count_residuals();
    auto residuals_lap = FivePointsLaplass(*residuals, *grid_);
    exchange_mirror_rows(residuals_lap);

    // step 2: count alpha
    double alpha_den_part = ProductByPointAndSum(*gradients_laplass_, *gradients_, *grid_);
    double alpha_nom_part =  ProductByPointAndSum(*residuals_lap, *gradients_, *grid_);

    send_value(alpha_den_part, 0, proc_rank_);
    double alpha_den;
    receive_value(&alpha_den, 0, proc_rank_);

    send_value(alpha_nom_part, 0, proc_rank_);
    double alpha_nom;
    receive_value(&alpha_nom, 0, proc_rank_);

    double alpha = alpha_den > 0.0 ? alpha_nom / alpha_den : 0.0;

    // step 3: count new gradients
    if (first_iter) {
      *gradients_ = *residuals;
      first_iter = false;
    } else {
      gradients_ = *residuals - *gradients_ * alpha;
    }

    // step 4: count new gradients laplass
    gradients_laplass_ = FivePointsLaplass(*gradients_, *grid_);
    exchange_mirror_rows(gradients_laplass_);

    // step 5: count tau
    double tau_den_part = ProductByPointAndSum(*gradients_laplass_, *gradients_, *grid_);
    double tau_nom_part = ProductByPointAndSum(*residuals, *gradients_, *grid_);

    send_value(tau_den_part, 0, proc_rank_);
    double tau_den;
    receive_value(&tau_den, 0, proc_rank_);

    send_value(tau_nom_part, 0, proc_rank_);
    double tau_nom;
    receive_value(&tau_nom, 0, proc_rank_);

    double tau = tau_den > 0.0 ? tau_nom / tau_den : 0.0;

    residuals->clear();  // free memory
    residuals_lap->clear();

    // step 6: count new values
    count_new_values(tau);

    // step 7: send difference to master and finish iter 
    send_value(values_difference(), 0, proc_rank_);
  }
}

std::shared_ptr<DM> GradientDescent::count_residuals() {
  auto residuals = std::shared_ptr<DM>(new DM(values_->num_rows(), values_->num_cols(), 0.0));

  auto value_lap = FivePointsLaplass(*values_, *grid_);
  exchange_mirror_rows(value_lap);

  for (size_t i = 0; i < residuals->num_rows(); ++i) {
    for (size_t j = 0; j < residuals->num_cols(); ++j) {
      (*residuals)(i, j) = grid_->is_bound_point({ static_cast<double>(i), static_cast<double>(j) }, proc_bounds_) ? 0.0 :
        (*value_lap)(i, j) - functions_.main_func((*grid_)(i, j));
    }
  }

  return residuals;
}

void GradientDescent::count_new_values(double tau) {
  *old_values_ = *values_;
  values_ = *values_ - *gradients_ * tau;
}

double GradientDescent::count_pre_error() const {
  auto psi = *values_;

  for (size_t i = 0; i < psi.num_rows(); ++i) {
    for (size_t j = 0; j < psi.num_cols(); ++j) {
      psi(i, j) = functions_.true_func((*grid_)(i, j)) - (*values_)(i, j);
    }
  }

  return ProductByPointAndSum(psi, psi, *grid_);
}

double GradientDescent::values_difference() const {
  auto difference = *values_ - *old_values_;
  return ProductByPointAndSum(*difference, *difference, *grid_);
}

void GradientDescent::init_grid(const GridData& grid_data, size_t start_row_idx, size_t end_row_idx,
                                size_t start_col_idx, size_t end_col_idx) {
  start_row_idx -= !proc_bounds_.is_up ? 1 : 0;
  end_row_idx += !proc_bounds_.is_low ? 1 : 0;
  start_col_idx -= !proc_bounds_.is_left ? 1 : 0;
  end_col_idx += !proc_bounds_.is_right ? 1 : 0;

  grid_ = std::shared_ptr<Grid>(new Grid(end_row_idx - start_row_idx, end_col_idx - start_col_idx));

  for (size_t i = 0; i < grid_data.r_num_points; ++i) {
    if (i >= start_row_idx && i < end_row_idx) {
      for (size_t j = 0; j < grid_data.c_num_points; ++j) {
        if (j >= start_col_idx && j < end_col_idx) {
          double cur_r_value = grid_data.r_upper_bound * functions_.step_func(static_cast<double>(i) / (grid_data.r_num_points - 1), grid_data.q) +
            grid_data.r_lower_bound * (1 - functions_.step_func(static_cast<double>(i) / (grid_data.r_num_points - 1), grid_data.q));
          double cur_c_value = grid_data.c_upper_bound * functions_.step_func(static_cast<double>(j) / (grid_data.c_num_points - 1), grid_data.q) +
            grid_data.c_lower_bound * (1 - functions_.step_func(static_cast<double>(j) / (grid_data.c_num_points - 1), grid_data.q));
          (*grid_)(i - start_row_idx, j - start_col_idx) = { cur_r_value, cur_c_value };
        }
      }
    }
  }
}

void GradientDescent::init_values() {
  values_ = std::shared_ptr<DM>(new DM(grid_->num_rows(), grid_->num_cols(), RAND_CONST));

  if (proc_bounds_.is_up) {
    for (size_t i = 0; i < grid_->num_cols(); ++i) {
      (*values_)(0, i) = functions_.bound_func((*grid_)(0, i));
    }
  }

  if (proc_bounds_.is_low) {
    for(size_t i = 0; i < grid_->num_cols(); ++i) {
      (*values_)(grid_->num_rows() - 1, i) = functions_.bound_func((*grid_)(grid_->num_rows() - 1, i));
    }
  }

  if (proc_bounds_.is_left) {
    for(size_t i = 0; i < grid_->num_rows(); ++i) {
      (*values_)(i, 0) = functions_.bound_func((*grid_)(i, 0));
    }
  }

  if (proc_bounds_.is_right) {
    for(size_t i = 0; i < grid_->num_rows(); ++i) {
      (*values_)(i, grid_->num_cols() - 1) = functions_.bound_func((*grid_)(i, grid_->num_cols() - 1));
    }
  }
}

void GradientDescent::exchange_mirror_rows(std::shared_ptr<DM> values) {
  if (proc_bounds_.is_up && proc_bounds_.is_low && proc_bounds_.is_left && proc_bounds_.is_right) {
    return;
  }

  // exchange rows
  if (first_send_.first) {
    if (!proc_bounds_.is_low) {
      send_vector(values->get_row(values->num_rows() - 2), proc_rank_ + 1, proc_rank_);
    }
    if (!proc_bounds_.is_up) { 
      send_vector(values->get_row(1), proc_rank_ - 1, proc_rank_);
    }
 
    if (!proc_bounds_.is_low) {
      receive_vector(&values->get_row_non_const(values->num_rows() - 1), proc_rank_ + 1, proc_rank_ + 1);
    }
    if (!proc_bounds_.is_up) {
      receive_vector(&values->get_row_non_const(0), proc_rank_ - 1, proc_rank_ - 1);
    }
  } else {
    if (!proc_bounds_.is_low) { 
      receive_vector(&values->get_row_non_const(values->num_rows() - 1), proc_rank_ + 1, proc_rank_ + 1);
    }
    if (!proc_bounds_.is_up) {
      receive_vector(&values->get_row_non_const(0), proc_rank_ - 1, proc_rank_ - 1);
    }

    if (!proc_bounds_.is_low) {
      send_vector(values->get_row(values->num_rows() - 2), proc_rank_ + 1, proc_rank_);
    }
    if (!proc_bounds_.is_up) {
      send_vector(values->get_row(1), proc_rank_ - 1, proc_rank_);
    }
  }

  // exchange cols (maybe need to remove copy-paste)
  std::vector<double> data(values->num_rows(), 0.0);
  if (first_send_.second){
    if (!proc_bounds_.is_left) {
      for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (*values)(i, 1);
      }
      send_vector(data, left_right_proc_.first, proc_rank_);
    }
    if (!proc_bounds_.is_right) {
      for (size_t i = 0; i < data.size(); ++i) {
	data[i]= (*values)(i, values->num_cols() - 2);
      }
      send_vector(data, left_right_proc_.second, proc_rank_);
    }

    if (!proc_bounds_.is_left) {
      receive_vector(&data, left_right_proc_.first, left_right_proc_.first);
      for (size_t i = 0; i < data.size(); ++i) {
        (*values)(i, 0) = data[i];
      }
    }
    if (!proc_bounds_.is_right) {
      receive_vector(&data, left_right_proc_.second, left_right_proc_.second);
      for (size_t i = 0; i < data.size(); ++i) {
	(*values)(i, values->num_cols() - 1) = data[i];
      }
    }
  } else {
    if (!proc_bounds_.is_left) {
      receive_vector(&data, left_right_proc_.first, left_right_proc_.first);
      for (size_t i = 0; i < data.size(); ++i) {
	(*values)(i, 0) = data[i];
      } 
    }
    if (!proc_bounds_.is_right) {
      receive_vector(&data, left_right_proc_.second, left_right_proc_.second);
      for (size_t i = 0; i < data.size(); ++i) {
        (*values)(i, values->num_cols() - 1) = data[i];
      }
    }

    if (!proc_bounds_.is_left) {
      for (size_t i = 0; i < data.size(); ++i) {
        data[i] = (*values)(i, 1);
      }
      send_vector(data, left_right_proc_.first, proc_rank_);
    }
    if (!proc_bounds_.is_right) {
      for (size_t i = 0; i < data.size(); ++i) {
        data[i]= (*values)(i, values->num_cols() - 2);
      }
      send_vector(data, left_right_proc_.second, proc_rank_);
    }
  }
}

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
      out_true_file << functions_.true_func((*grid_)(i, j)) << ", ";                                                                                       
    }
    out_value_file << std::endl;
    out_true_file << std::endl;
  }

  out_value_file.close();
  out_true_file.close();
}

}  // namespace DTS
