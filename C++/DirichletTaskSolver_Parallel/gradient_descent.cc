#include "mpi_helpers.h"

#include "gradient_descent.h"

namespace DTS {

GradientDescent::GradientDescent(const Grid& grid, const Functions& functions,
				 ProcType proc_type, size_t proc_rank, size_t start_idx, size_t end_idx)
  : grid_(grid)
  , functions_(functions)
  , proc_type_(proc_type)
  , proc_rank_(proc_rank)
  , start_index_(start_idx)
  , end_index_(end_idx)
{
  clear();
  values_ = init_values();

  size_t h = grid.num_height_points();
  size_t w = end_index_ - start_index_;

  old_values_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  gradients_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
  old_gradients_ = std::shared_ptr<DM>(new DM(w, h, 0.0));
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
      break;
    }

    // step 1: count residuals
    auto residuals = count_residuals();
    auto residuals_lap = DM::FivePointsLaplass(*residuals, grid_, start_index_);
    exchange_mirror_rows(residuals_lap);

    // step 2: count alpha
    double alpha_den_part = DM::ProductByPointAndSum(*old_gradients_laplass_, *old_gradients_, grid_, proc_type_, start_index_);
    double alpha_nom_part =  DM::ProductByPointAndSum(*residuals_lap, *old_gradients_, grid_, proc_type_, start_index_);

    send_value(alpha_den_part, 0, proc_rank_);
    double alpha_den;
    receive_value(&alpha_den, 0, proc_rank_);

    send_value(alpha_nom_part, 0, proc_rank_);
    double alpha_nom;
    receive_value(&alpha_nom, 0, proc_rank_);

    double alpha = alpha_den > 0.0 ? alpha_nom / alpha_den : 0.0;

    // step 3: count new gradients
    old_gradients_.swap(gradients_);
    if (first_iter) {
      *gradients_ = *residuals;
      first_iter = false;
    } else {
      gradients_ = *residuals - *residuals_lap * alpha;
    }

    // step 4: count new gradients laplass
    *old_gradients_laplass_ = *gradients_laplass_;
    gradients_laplass_ = DM::FivePointsLaplass(*gradients_, grid_, start_index_);
    exchange_mirror_rows(gradients_laplass_);
    
    // step 5: count tau
    double tau_den_part = DM::ProductByPointAndSum(*gradients_laplass_, *gradients_, grid_, proc_type_, start_index_);
    double tau_nom_part = DM::ProductByPointAndSum(*residuals, *gradients_, grid_, proc_type_, start_index_);

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

  auto value_lap = DM::FivePointsLaplass(*values_, grid_, start_index_);
  exchange_mirror_rows(value_lap);

  for (size_t i = 0; i < residuals->num_rows(); ++i) {
    for (size_t j = 0; j < residuals->num_cols(); ++j) {
      (*residuals)(i, j) = grid_.isBoundPoint({ static_cast<double>(i + start_index_), static_cast<double>(j) }) ? 0.0 :
        (*value_lap)(i, j) - functions_.main_func(grid_(i + start_index_, j));
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
      psi(i, j) = functions_.true_func(grid_(i + start_index_, j)) - (*values_)(i, j);
    }
  }

  return DM::ProductByPointAndSum(psi, psi, grid_, proc_type_, start_index_);
}

double GradientDescent::values_difference() const {
  auto difference = *values_ - *old_values_;
  return DM::ProductByPointAndSum(*difference, *difference, grid_, proc_type_, start_index_);
}

std::shared_ptr<DM> GradientDescent::init_values() {
  if (proc_type_ != GLOBAL_PROC) {
    start_index_ -= proc_type_ == UPPER_PROC ? 0 : 1;
    end_index_ += proc_type_ == LOWER_PROC ? 0 : 1;
  }
  size_t width_size = end_index_ - start_index_;
  auto retval = std::shared_ptr<DM>(new DM(width_size, grid_.num_height_points(), RAND_CONST));

  if (proc_type_ != CENTER_PROC) {
    for (size_t i = 0; i < grid_.num_height_points(); ++i) {
      if (proc_type_ != LOWER_PROC) {
        (*retval)(0, i) = functions_.bound_func(grid_(start_index_, i));
      }

      if (proc_type_ != UPPER_PROC) {
        (*retval)(width_size - 1, i) = functions_.bound_func(grid_(end_index_ - 1 , i));
      }
    }
  }

  for (size_t i = 0; i < width_size; ++i) {
    (*retval)(i, 0) = functions_.bound_func(grid_(i + start_index_, 0));
    (*retval)(i, grid_.num_height_points() - 1) = functions_.bound_func(grid_(i + start_index_, grid_.num_height_points() - 1));
  }

  return retval;
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
