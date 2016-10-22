#include "grid.h"

Grid::Grid(const Bounds& bounds, int num_height_points, int num_width_points)
  : bounds_(bounds)
  , num_height_points_(num_height_points + 1)
  , num_width_points_(num_width_points + 1)
  , step_height_((bounds_.height_upper_bound - bounds_.height_lower_bound) / (num_height_points))
  , step_width_((bounds_.width_upper_bound - bounds_.width_lower_bound) / (num_height_points)) { }

Point Grid::operator()(int row, int col) const {
  return{ (height_lower_bound() + step_height_ * row), (width_lower_bound() + step_width_ * col) };
}

bool Grid::isBoundPoint(const Point& p) const {
  bool height_bound = p.height == 0 || p.height == num_height_points_ - 1;
  bool width_bound = p.width == 0 || p.width == num_width_points_ - 1;
  return height_bound || width_bound;
}

void Grid::debug_print() const {
  for (int i = 0; i < num_height_points(); ++i) {
    for (int j = 0; j < num_width_points(); ++j) {
      std::cout << (*this)(i, j) << " ";
    }
    std::cout << std::endl;
  }
}
