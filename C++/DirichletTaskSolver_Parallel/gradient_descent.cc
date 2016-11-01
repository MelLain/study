#include <cmath>

#include "matrix_operations.h"
#include "mpi_helpers.h"
#include "gradient_descent.h"

namespace DTS {

GradientDescent::GradientDescent(const GridData& grid_data, const Functions& functions,
  ProcType proc_type, size_t proc_rank, size_t start_idx, size_t end_idx)
  : functions_(functions)
  , proc_type_(proc_type)
  , proc_rank_(proc_rank)
{
  clear();
  init_grid(grid_data, start_idx, end_idx);
  init_values();

  size_t h = grid_->num_cols();
  size_t w = grid_->num_rows();

  old_values_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  gradients_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  gradients_laplass_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  old_gradients_laplass_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
}

void GradientDescent::FitModel() {
  bool first_iter = true;
  while (true) {
    FlagType flag;
    receive_flag(&flag, 0, proc_rank_);

    if (flag == TERMINATE) {
      send_value(count_pre_error(), 0, proc_rank_);
      send_matrix(*values_, 0, proc_rank_);
      send_grid(*grid_, 0, proc_rank_);
      break;
    }

    exchange_mirror_rows(values_);
    *old_gradients_laplass_ = *gradients_laplass_;

    // step 1: count residuals
    auto residuals = count_residuals();
    auto residuals_lap = FivePointsLaplass(*residuals, *grid_);
    exchange_mirror_rows(residuals_lap);

    // step 2: count alpha
    double alpha_den_part = ProductByPointAndSum(*old_gradients_laplass_, *gradients_, *grid_);
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
      (*residuals)(i, j) = grid_->is_bound_point({ static_cast<double>(i), static_cast<double>(j) }, proc_type_) ? 0.0 :
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

void GradientDescent::init_grid(const GridData& grid_data, size_t start_index, size_t end_index) {
  if (proc_type_ != GLOBAL_PROC) {
    start_index -= proc_type_ == UPPER_PROC ? 0 : 1;
    end_index += proc_type_ == LOWER_PROC ? 0 : 1;
  }
  grid_ = std::shared_ptr<Grid>(new Grid(end_index - start_index, grid_data.num_points));

  double denominator = 0.0;
  for (size_t i = 0; i < grid_data.num_points - 1; ++i) {
    denominator += pow(grid_data.multiplier, i);
  }
  double step_0 = denominator > 0.0 ? fabs(grid_data.upper_bound - grid_data.lower_bound) / denominator : 0.0;

  double cur_r_value = grid_data.lower_bound;
  for (size_t i = 0; i < grid_data.num_points; ++i) {
    if (i >= start_index && i < end_index) {
      double cur_c_value = grid_data.lower_bound;
      for (size_t j = 0; j < grid_->num_cols(); ++j) {
        (*grid_)(i - start_index, j) = { cur_r_value, cur_c_value };
        cur_c_value += step_0 * pow(grid_data.multiplier, j);
      }
    }
    cur_r_value += step_0 * pow(grid_data.multiplier, i);
  }
}

void GradientDescent::init_values() {
  values_ = std::shared_ptr<DM>(new DM(grid_->num_rows(), grid_->num_cols(), RAND_CONST));

  if (proc_type_ != CENTER_PROC) {
    for (size_t i = 0; i < grid_->num_cols(); ++i) {
      if (proc_type_ != LOWER_PROC) {
        (*values_)(0, i) = functions_.bound_func((*grid_)(0, i));
      }

      if (proc_type_ != UPPER_PROC) {
        (*values_)(grid_->num_rows() - 1, i) = functions_.bound_func((*grid_)(grid_->num_rows() - 1 , i));
      }
    }
  }

  for (size_t i = 0; i < grid_->num_rows(); ++i) {
    (*values_)(i, 0) = functions_.bound_func((*grid_)(i, 0));
    (*values_)(i, grid_->num_cols() - 1) = functions_.bound_func((*grid_)(i, grid_->num_cols() - 1));
  }
}

void GradientDescent::exchange_mirror_rows(std::shared_ptr<DM> values) {
  // proc_rank_ mod 2 == 1 -> first send, then receive
  //                  == 0 -> first receive, then send
  if (proc_type_ == GLOBAL_PROC) {
    return;
  }

  if (proc_rank_ % 2 == 1) {
    if (proc_type_ != LOWER_PROC) {
      send_vector(values->get_row(values->num_rows() - 2), proc_rank_ + 1, proc_rank_);
    }
    if (proc_type_ != UPPER_PROC) {
      send_vector(values->get_row(1), proc_rank_ - 1, proc_rank_);
    }
    if (proc_type_ != LOWER_PROC) {
      receive_vector(&values->get_row_non_const(values->num_rows() - 1), proc_rank_ + 1, proc_rank_ + 1);
    }
    if (proc_type_ != UPPER_PROC) {
      receive_vector(&values->get_row_non_const(0), proc_rank_ - 1, proc_rank_ - 1);
    }
  } else {
    if (proc_type_ != LOWER_PROC) {
      receive_vector(&values->get_row_non_const(values->num_rows() - 1), proc_rank_ + 1, proc_rank_ + 1);
    }
    if (proc_type_ != UPPER_PROC) {
      receive_vector(&values->get_row_non_const(0), proc_rank_ - 1, proc_rank_ - 1);
    }
    if (proc_type_ != LOWER_PROC) {
      send_vector(values->get_row(values->num_rows() - 2), proc_rank_ + 1, proc_rank_);
    }
    if (proc_type_ != UPPER_PROC) {
      send_vector(values->get_row(1), proc_rank_ - 1, proc_rank_);
    }
  }
}

}  // namespace DTS
