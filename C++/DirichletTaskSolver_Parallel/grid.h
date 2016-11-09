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
  double r_lower_bound;
  double r_upper_bound;
  size_t r_num_points;

  double c_lower_bound;
  double c_upper_bound;
  size_t c_num_points;

  double q;
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

  bool is_bound_point(const Point& p, const ProcBounds& proc_bounds) const;

  void debug_print() const;

 private:
  Matrix<Point> data_;
};

}  // namespace DTS
