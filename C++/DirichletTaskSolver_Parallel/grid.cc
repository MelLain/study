#include "grid.h"

namespace DTS {

Grid::Grid(const Bounds& bounds, size_t num_height_points, size_t num_width_points)
  : bounds_(bounds)
  , num_height_points_(num_height_points + 1)
  , num_width_points_(num_width_points + 1)
  , step_height_((bounds_.height_upper_bound - bounds_.height_lower_bound) / (num_height_points))
  , step_width_((bounds_.width_upper_bound - bounds_.width_lower_bound) / (num_height_points)) { }

Point Grid::operator()(size_t row, size_t col) const {
  return{ (height_lower_bound() + step_height_ * row), (width_lower_bound() + step_width_ * col) };
}

bool Grid::isBoundPoint(const Point& p) const {
  bool height_bound = p.height == 0 || p.height == num_height_points_ - 1;
  bool width_bound = p.width == 0 || p.width == num_width_points_ - 1;
  return height_bound || width_bound;
}

void Grid::debug_print() const {
  for (size_t i = 0; i < num_height_points(); ++i) {
    for (size_t j = 0; j < num_width_points(); ++j) {
      std::cout << (*this)(i, j) << " ";
    }
    std::cout << std::endl;
  }
}

}  // namespace DTS
