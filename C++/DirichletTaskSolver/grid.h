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
   Grid(const Bounds& bounds, int num_height_points, int num_width_points);

  double height_lower_bound() const { return bounds_.height_lower_bound; };
  double height_upper_bound() const { return bounds_.height_upper_bound; };
  double width_lower_bound() const { return bounds_.width_lower_bound; };
  double width_upper_bound() const { return bounds_.width_upper_bound; };

  int num_width_points() const { return num_width_points_; }
  int num_height_points() const { return num_height_points_; }

  double step_height() const { return step_height_; }
  double step_width() const { return step_width_; }

  Point operator()(int row, int col) const;

  bool isBoundPoint(const Point& p) const;

  void debug_print() const;

 private:
  Bounds bounds_;
  int num_height_points_;
  int num_width_points_;
  double step_height_;
  double step_width_;
};

}  // namespace DTS