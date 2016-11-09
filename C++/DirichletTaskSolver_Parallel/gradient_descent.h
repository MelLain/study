#pragma once

#include <memory>
#include <vector>
#include <utility>

#include "common.h"
#include "matrix.h"
#include "grid.h"

namespace DTS {

typedef double(*Func)(const Point& point);
typedef double(*StepFunc)(double value, double q);

struct Functions {
  Func main_func;
  Func bound_func;
  Func true_func;
  StepFunc step_func;
};

class GradientDescent {
 public:
  GradientDescent(const GridData& grid_data, const Functions& functions, const ProcBounds& proc_bounds,
                  size_t proc_rank, std::pair<bool, bool> first_send, std::pair<size_t, size_t> left_right_proc,
                  size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx);

  void FitModel();
  const DM& GetCurrentMatrix() const { return *values_; }

  void clear() {
    grid_.reset();
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

  void init_grid(const GridData& grid_data, const Functions& functions,
                 size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx);
  void init_values();
  void exchange_mirror_rows(std::shared_ptr<DM> values);

  Functions functions_;
  ProcBounds proc_bounds_;
  size_t proc_rank_;

  std::pair<bool, bool> first_send_;
  std::pair<size_t, size_t> left_right_proc_;

  std::shared_ptr<Grid> grid_;

  std::shared_ptr<DM> values_;
  std::shared_ptr<DM> old_values_;
  std::shared_ptr<DM> gradients_;
  std::shared_ptr<DM> gradients_laplass_;
  std::shared_ptr<DM> old_gradients_laplass_;
};

}  // namespace DTS
