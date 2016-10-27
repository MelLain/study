#include <cmath>

#include "grid.h"

namespace DTS {

Grid::Grid(const Bounds& bounds, size_t num_width_points, size_t num_height_points, double multiplier)
  : bounds_(bounds)
  , num_width_points_(num_width_points + 1)
  , num_height_points_(num_height_points + 1)
  , multiplier_(multiplier)
  , step_width_0_(0)
  , step_height_0_(0)
{
  double denominator = 0.0;
  for (size_t i = 0; i < num_width_points; ++i) {
    denominator += pow(multiplier, i);
  }
  step_width_0_ = denominator > 0 ? (bounds.width_upper_bound - bounds.width_lower_bound) / denominator : 0.0;

  denominator = 0.0;
  for (size_t i = 0; i < num_height_points; ++i) {
    denominator+= pow(multiplier, i);
  }
  step_height_0_ = denominator > 0 ? (bounds.height_upper_bound - bounds.height_lower_bound) / denominator : 0.0;
}

double Grid::step_height(size_t col) const {
  return col > 0 ? step_height_0_ * pow(multiplier_, col - 1) : 0.0;
}

double Grid::step_width(size_t row) const {
  return row > 0 ? step_width_0_ * pow(multiplier_, row - 1) : 0.0;
}

Point Grid::operator()(size_t row, size_t col) const {
  double row_value = 0.0;
  for (size_t i = 0; i < row; ++i) {
    row_value += step_width_0_ * pow(multiplier_, i);
  }

  double col_value = 0.0;
  for (size_t i = 0; i < col; ++i) {
    col_value += step_height_0_ * pow(multiplier_, i);
  }

  return{ (width_lower_bound() + row_value), (height_lower_bound() + col_value) };
}

bool Grid::isBoundPoint(const Point& p) const {
  bool height_bound = p.height == 0 || p.height == num_height_points_ - 1;
  bool width_bound = p.width == 0 || p.width == num_width_points_ - 1;
  return height_bound || width_bound;
}

void Grid::debug_print() const {
  for (size_t i = 0; i < num_width_points(); ++i) {
    for (size_t j = 0; j < num_height_points(); ++j) {
      std::cout << (*this)(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

}  // namespace DTS
