#pragma once

#include <memory>

#include "common.h"
#include "matrix.h"

namespace DTS {

struct Point {
  double r_value;
  double c_value;

  friend std::ostream& operator<<(std::ostream& out, const Point& p) {
    return out << "(" << p.r_value << ", " << p.c_value << ")";
  }
};

struct GridData {
  double lower_bound;
  double upper_bound;
  size_t num_points;
  double multiplier;
};

class Grid {
 public:
  Grid(size_t num_rows, size_t num_cols);

  size_t num_rows() const { return data_.num_rows(); }
  size_t num_cols() const { return data_.num_cols(); }

  double r_step(size_t row) const;
  double c_step(size_t col) const;

  const Point& operator()(size_t row, size_t col) const;
  Point& operator()(size_t row, size_t col);

  bool is_bound_point(const Point& p, ProcType proc_type) const;

  void debug_print() const;

 private:
  Matrix<Point> data_;
};

}  // namespace DTS
