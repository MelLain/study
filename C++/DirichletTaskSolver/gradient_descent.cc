#include "gradient_descent.h"

namespace DTS {

  namespace {
    std::shared_ptr<DM> CreateInitMatrix(const Grid& grid, Func bound_func) {
      auto retval = std::shared_ptr<DM>(new DM(grid.num_width_points(), grid.num_height_points(), RAND_CONST));

      for (int i = 0; i < grid.num_height_points(); ++i) {
        (*retval)(0, i) = bound_func(grid(0, i));
        (*retval)(grid.num_width_points() - 1, i) = bound_func(grid(grid.num_width_points() - 1, i));
      }

      for (int i = 0; i < grid.num_width_points(); ++i) {
        (*retval)(i, 0) = bound_func(grid(i, 0));
        (*retval)(i, grid.num_height_points() - 1) = bound_func(grid(i, grid.num_height_points() - 1));
      }

      return retval;
    }
  }

GradientDescent::GradientDescent(const Grid& grid, const Functions& functions)
  : grid_(grid)
  , functions_(functions)
{
  clear();
  values_ = CreateInitMatrix(grid, functions.bound_func);

  int h = grid.num_height_points();
  int w = grid.num_width_points();
  old_values_ = std::shared_ptr<DM>(new DM(h, w, 0.0));
  gradients_ = std::shared_ptr<DM>(new DM(h, w, 0.0));
  old_gradients_ = std::shared_ptr<DM>(new DM(h, w, 0.0));
  gradients_laplass_ = std::shared_ptr<DM>(new DM(h, w, 0.0));
  old_gradients_laplass_ = std::shared_ptr<DM>(new DM(h, w, 0.0));
}

void GradientDescent::FitModel() {
  bool first_iter = true;
  while (true) {
    auto residuals = count_residuals();
    auto residuals_lap = DM::FivePointsLaplass(*residuals, grid_);

    count_new_gradients(*residuals, *residuals_lap, first_iter);
    first_iter = false;

    *old_gradients_laplass_ = *gradients_laplass_;
    gradients_laplass_ = DM::FivePointsLaplass(*gradients_, grid_);
    
    double tau = count_tau(*residuals);
    residuals->clear();  // free memory
    count_new_values(tau);

    error_by_iter_.push_back(count_error());

    if (finished())
      break;
  }
}

double GradientDescent::count_tau(const DM& res) const {
  double denominator = DM::ProductByPointAndSum(*gradients_laplass_, *gradients_, grid_);
  return denominator > 0.0 ? DM::ProductByPointAndSum(res, *gradients_, grid_) / denominator : 0.0;
}

double GradientDescent::count_alpha(const DM& res_lap) const {
  double denominator = DM::ProductByPointAndSum(*old_gradients_laplass_, *old_gradients_, grid_);
  return denominator > 0.0 ? DM::ProductByPointAndSum(res_lap, *old_gradients_, grid_) / denominator : 0.0;
}

std::shared_ptr<DM> GradientDescent::count_residuals() const {
  auto value_lap = DM::FivePointsLaplass(*values_, grid_);
  auto residuals = std::shared_ptr<DM>(new DM(grid_.num_height_points(), grid_.num_width_points(), 0.0));

  for (int i = 0; i < residuals->num_rows(); ++i) {
    for (int j = 0; j < residuals->num_cols(); ++j) {
      (*residuals)(i, j) = grid_.isBoundPoint({ static_cast<double>(i), static_cast<double>(j) }) ? 0.0 :
        (*value_lap)(i, j) - functions_.main_func(grid_(i, j));
    }
  }

  return residuals;
}

void GradientDescent::count_new_gradients(const DM& res, const DM& res_lap, bool first_iter) {
  old_gradients_.swap(gradients_);
  
  if (first_iter) {
    *gradients_ = res;
    return;
  }

  gradients_ = res - res_lap * count_alpha(res_lap);
}

void GradientDescent::count_new_values(double tau) {
  *old_values_ = *values_;
  values_ = *values_ - *gradients_ * tau;
}

double GradientDescent::count_error() const {
  auto psi = *values_;

  for (int i = 0; i < psi.num_rows(); ++i) {
    for (int j = 0; j < psi.num_cols(); ++j) {
      psi(i, j) = functions_.true_func(grid_(i, j)) - psi(i, j);
    }
  }

  return sqrt(DM::ProductByPointAndSum(psi, psi, grid_));
}

bool GradientDescent::finished() const {
  auto difference = *values_ - *old_values_;
  return sqrt(DM::ProductByPointAndSum(*difference, *difference, grid_)) < EPS;
}


}  // namespace DTS