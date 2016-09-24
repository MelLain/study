#pragma once

#include <cmath>

#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <iterator>

namespace NLabAlgOne {

const double INF = 1e+20;
const double EPS = 1e-20;
const int MAX_POINTS_TO_BRUTE_FORCE = 3;

typedef std::vector<std::pair<double, double> > TPoints;
typedef TPoints::iterator TPsIter;

struct TResult {
  TResult() { reset(); }

  void operator=(const TResult& rhs) {
    points_1.clear();
    points_2.clear();
    for (int i = 0; i < rhs.points_1.size(); ++i) {
      points_1.push_back(rhs.points_1[i]);
      points_2.push_back(rhs.points_2[i]);
    }
    dist = rhs.dist;
  }

  void reset() {
    points_1.clear();
    points_2.clear();
    dist = INF;
  }

  TPoints points_1;
  TPoints points_2;
  double dist;
};

TPsIter MergeSort(TPsIter in_iter, TPsIter temp_iter, int left_bound, int right_bound);
void PrintResult(const TResult& result, const std::string& data_path);
double EuclidDist(const std::pair<double, double>& point_1, const std::pair<double, double>& point_2);
TResult BruteForce(TPsIter x_y_begin, TPsIter x_y_end);
TResult FindMinDist(TPsIter x_y_begin, const TPoints& y_x_points, int left_bound, int right_bound);

}  // namespace NLabAlgOne