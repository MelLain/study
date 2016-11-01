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

bool Grid::is_bound_point(const Point& p, ProcType proc_type) const {
  // This method is correct if number of rows per processor is >= 2

  bool up = (proc_type == GLOBAL_PROC || proc_type == UPPER_PROC);
  bool lw = (proc_type == GLOBAL_PROC || proc_type == LOWER_PROC);

  bool row_bound = (p.r_value == 0 && up) || (p.r_value == data_.num_rows() - 1 && lw);
  bool col_bound = p.c_value == 0 || p.c_value == data_.num_cols() - 1;
  return row_bound || col_bound;
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
