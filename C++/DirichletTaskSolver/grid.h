#pragma once

#include "common.h"

namespace DTS {

struct Point {
  double height;
  double width;

  friend std::ostream& operator<<(std::ostream& out, const Point& p) {
    return out << "(" << p.height << ", " << p.width << ")";
  }
};

struct Bounds {
  double height_lower_bound;
  double height_upper_bound;
  double width_lower_bound;
  double width_upper_bound;
};

class Grid {
 public:
   Grid(const Bounds& bounds, size_t num_height_points, size_t num_width_points);

  double height_lower_bound() const { return bounds_.height_lower_bound; };
  double height_upper_bound() const { return bounds_.height_upper_bound; };
  double width_lower_bound() const { return bounds_.width_lower_bound; };
  double width_upper_bound() const { return bounds_.width_upper_bound; };

  size_t num_width_points() const { return num_width_points_; }
  size_t num_height_points() const { return num_height_points_; }

  double step_height() const { return step_height_; }
  double step_width() const { return step_width_; }

  Point operator()(size_t row, size_t col) const;

  bool isBoundPoint(const Point& p) const;

  void debug_print() const;

 private:
  Bounds bounds_;
  size_t num_height_points_;
  size_t num_width_points_;
  double step_height_;
  double step_width_;
};

}  // namespace DTS
