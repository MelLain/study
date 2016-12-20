#include "grid.h"

namespace DTS {

Grid::Grid(size_t num_rows, size_t num_cols)
  : data_(Matrix<Point>(num_rows, num_cols, { 0.0, 0.0 })) { }

double Grid::r_step(size_t row) const {
  return row > 0 ? data_(row, 0).r_value - data_(row - 1, 0).r_value : 0.0;
}

double Grid::c_step(size_t col) const {
  return col > 0 ? data_(0, col).c_value - data_(0, col - 1).c_value : 0.0;
}

const Point& Grid::operator()(size_t row, size_t col) const {
  return data_(row, col);
}

Point& Grid::operator()(size_t row, size_t col) {
  return data_(row, col);
}

bool Grid::is_bound_point(const Point& p, const ProcBounds& proc_bounds) const {
  // This method is correct if number of rows and cols per processor is >= 2
  bool left_bound = (p.c_value == 0 && proc_bounds.is_left);
  bool right_bound = (p.c_value == data_.num_cols() - 1 && proc_bounds.is_right);
  bool up_bound = (p.r_value == 0 && proc_bounds.is_up);
  bool low_bound = (p.r_value == data_.num_rows() - 1 && proc_bounds.is_low);

  return left_bound || right_bound || up_bound || low_bound;
}

void Grid::debug_print() const {
  for (size_t i = 0; i < data_.num_rows(); ++i) {
    for (size_t j = 0; j < data_.num_cols(); ++j) {
      std::cout << (*this)(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

}  // namespace DTS
