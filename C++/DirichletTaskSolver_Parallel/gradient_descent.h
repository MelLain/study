#pragma once

#include <memory>
#include <vector>

#include "common.h"
#include "matrix.h"
#include "grid.h"

namespace DTS {

typedef double(*Func)(const Point& point);

struct Functions {
  Func main_func;
  Func bound_func;
  Func true_func;
};

class GradientDescent {
 public:
  GradientDescent(const Grid& grid, const Functions& functions,
		  ProcType proc_type, size_t proc_rank, size_t start_idx, size_t end_idx);

  void FitModel();
  const DM& GetCurrentMatrix() const { return *values_; }

  void clear() {
    values_.reset();
    old_values_.reset();
    gradients_.reset();
    gradients_laplass_.reset();
    old_gradients_laplass_.reset();
  }

 private:
  std::shared_ptr<DM> count_residuals();
  void count_new_values(double tau);

  double count_pre_error() const;
  double values_difference() const;

  std::shared_ptr<DM> init_values();
  void exchange_mirror_rows(std::shared_ptr<DM> values);

  Grid grid_;
  Functions functions_;
  ProcType proc_type_;
  size_t proc_rank_;
  size_t start_index_;
  size_t end_index_;

  std::shared_ptr<DM> values_;
  std::shared_ptr<DM> old_values_;
  std::shared_ptr<DM> gradients_;
  std::shared_ptr<DM> gradients_laplass_;
  std::shared_ptr<DM> old_gradients_laplass_;
};

}  // namespace DTS
