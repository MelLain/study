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

const double INF = 1e+40;
const int MAX_POINTS_TO_BRUTE_FORCE = 3;

typedef std::vector<std::pair<float, float> > TPoints;
typedef TPoints::iterator TPsIter;

struct TResult {
  std::pair<float, float> point_1;
  std::pair<float, float> point_2;
  float dist;
};

TPsIter MergeSort(TPsIter in_iter, TPsIter temp_iter, int left_bound, int right_bound);
void PrintResult(const TResult& result, const std::string& data_path);
float EuclidDist(const std::pair<float, float>& point_1, const std::pair<float, float>& point_2);
TResult BruteForce(TPsIter x_y_begin, int size);
TResult FindMinDist(TPsIter x_y_begin, TPsIter y_x_begin, TPsIter y_x_end, int left_bound, int right_bound);

}  // namespace NLabAlgOne