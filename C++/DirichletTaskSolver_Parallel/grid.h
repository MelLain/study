#pragma once

#include "common.h"

namespace DTS {

struct Point {
  double width;
  double height;

  friend std::ostream& operator<<(std::ostream& out, const Point& p) {
    return out << "(" << p.width << ", " << p.height << ")";
  }
};

struct Bounds {
  double width_lower_bound;
  double width_upper_bound;
  double height_lower_bound;
  double height_upper_bound;
};

class Grid {
 public:
  Grid(const Bounds& bounds, size_t num_width_points, size_t num_height_points, double multiplier);

  double height_lower_bound() const { return bounds_.height_lower_bound; };
  double height_upper_bound() const { return bounds_.height_upper_bound; };
  double width_lower_bound() const { return bounds_.width_lower_bound; };
  double width_upper_bound() const { return bounds_.width_upper_bound; };

  size_t num_width_points() const { return num_width_points_; }
  size_t num_height_points() const { return num_height_points_; }

  double step_height(size_t col) const;
  double step_width(size_t row) const;

  Point operator()(size_t row, size_t col) const;

  bool isBoundPoint(const Point& p) const;

  void debug_print() const;

 private:
  Bounds bounds_;
  size_t num_width_points_;
  size_t num_height_points_;
  double multiplier_;
  double step_width_0_;
  double step_height_0_;
};

}  // namespace DTS
