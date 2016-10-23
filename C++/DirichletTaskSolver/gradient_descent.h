#pragma once

#include <memory>
#include <vector>

#include "common.h"
#include "matrix.h"
#include "grid.h"

namespace DTS {

typedef double(*Func)(const Point& point);
typedef Matrix<double> DM;

struct Functions {
  Func main_func;
  Func bound_func;
  Func true_func;
};

class GradientDescent {
 public:
  GradientDescent(const Grid& grid, const Functions& functions);

  void FitModel();
  const DM& GetCurrentMatrix() const { return *values_; }
  const std::vector<double>& GetErrors() const { return error_by_iter_; }

  void clear() {
    values_ = nullptr;
    old_values_ = nullptr;
    gradients_ = nullptr;
    gradients_laplass_ = nullptr;
    old_gradients_laplass_ = nullptr;
    error_by_iter_.clear();
  }

 private:
   double count_tau(const DM& res) const;
   double count_alpha(const DM& res_lap) const;
   std::shared_ptr<DM> count_residuals() const;
   void count_new_gradients(const DM& res, const DM& res_lap, bool first_iter);
   void count_new_values(double tau);

   double count_error() const;
   bool finished() const;
   size_t num_processed_iters() const { return error_by_iter_.size(); }

  Grid grid_;
  Functions functions_;

  std::shared_ptr<DM> values_;
  std::shared_ptr<DM> old_values_;
  std::shared_ptr<DM> gradients_;
  std::shared_ptr<DM> old_gradients_;
  std::shared_ptr<DM> gradients_laplass_;
  std::shared_ptr<DM> old_gradients_laplass_;

  std::vector<double> error_by_iter_;
};

}  // namespace DTS
