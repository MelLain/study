#include "helper.h"

namespace NLabAlgOne {

TPsIter MergeSort(TPsIter in_iter, TPsIter temp_iter, int left_bound, int right_bound) {

  if (left_bound == right_bound) {
    *(temp_iter + left_bound) = *(in_iter + left_bound);
    return temp_iter;
  }

  int middle = static_cast<int>((left_bound + right_bound) / 2);
  auto left_array_iter = MergeSort(in_iter, temp_iter, left_bound, middle);
  auto right_array_iter = MergeSort(in_iter, temp_iter, middle + 1, right_bound);

  auto out_array_iter = (&(*left_array_iter) == &(*in_iter)) ? temp_iter : in_iter;

  int cur_left_index = left_bound;
  int cur_right_index = middle + 1;

  for (int index = left_bound; index <= right_bound; ++index) {
    if (cur_left_index <= middle && cur_right_index <= right_bound) {
     if (*(left_array_iter + cur_left_index) < *(right_array_iter + cur_right_index)) {
        *(out_array_iter + index) = *(left_array_iter + cur_left_index++);
      } else {
        *(out_array_iter + index) = *(right_array_iter + cur_right_index++);
      }
    } else {
      if (cur_left_index <= middle) {
        *(out_array_iter + index) = *(left_array_iter + cur_left_index++);
      } else {
        *(out_array_iter + index) = *(right_array_iter + cur_right_index++);
      }
    }
  }

  return out_array_iter;
}

void PrintResult(const TResult& result, const std::string& data_path) {
  std::cout << "File = " << data_path << std::endl;
  for (int i = 0; i < result.points_1.size(); ++i) {
    std::cout << "[ (" << result.points_1[i].first << ", " << result.points_1[i].second << "), ("
              << result.points_2[i].first << ", " << result.points_2[i].second << ")]" << std::endl;
  }
  std::cout << "Value = " << result.dist << std::endl;
}

double EuclidDist(const std::pair<double, double>& point_1, const std::pair<double, double>& point_2) {
  return std::sqrt(std::pow(point_1.first - point_2.first, 2) + std::pow(point_1.second - point_2.second, 2));
}

TResult BruteForce(TPsIter x_y_begin, TPsIter x_y_end) {
  TResult result;
  for (auto iter = x_y_begin; iter != x_y_end; ++iter) {
    auto local_iter = iter + 1;
    for (; local_iter != x_y_end; ++local_iter) {
      double dist = EuclidDist(*iter, *local_iter);
      if (fabs(dist - result.dist) < EPS) {
        result.points_1.push_back(*iter);
        result.points_2.push_back(*local_iter);
      } else if (dist < result.dist) {
        result.reset();
        result.points_1.push_back(*iter);
        result.points_2.push_back(*local_iter);
        result.dist = dist;
      }
    }
  }
  return result;
}

TResult FindMinDist(TPsIter x_y_begin, const TPoints& y_x_points, int left_bound, int right_bound) {
  int size = right_bound - left_bound;
  if (size <= MAX_POINTS_TO_BRUTE_FORCE) {
    return BruteForce(x_y_begin + left_bound, x_y_begin + right_bound);
  }
  int middle = left_bound + static_cast<int>(size / 2);

  TResult res_left = FindMinDist(x_y_begin, y_x_points, left_bound, middle);
  TResult res_right = FindMinDist(x_y_begin, y_x_points, middle, right_bound);
  
  TResult result;
  if (fabs(res_left.dist - res_right.dist) < EPS) {
    
    result = res_right;
    for (int i = 0; i < res_left.points_1.size(); ++i) {
      result.points_1.push_back(res_left.points_1[i]);
      result.points_2.push_back(res_left.points_2[i]);
    }
  } else {
    result = res_left.dist > res_right.dist ? res_right : res_left;
  }

  for (int index = 0; index < y_x_points.size() - 1; ++index) {
    if (y_x_points[index].second <= (x_y_begin + right_bound - 1)->first &&
        y_x_points[index].second >= (x_y_begin + left_bound)->first) {
      for (auto local_index = index + 1; local_index != (index + 7); ++local_index) {
        if (local_index >= y_x_points.size())
          break;

        if ((y_x_points[local_index].second <= (x_y_begin + middle - 1)->first &&
             y_x_points[index].second <= (x_y_begin + middle - 1)->first) ||
            (y_x_points[local_index].second >= (x_y_begin + middle)->first &&
             y_x_points[index].second >= (x_y_begin + middle)->first)) {
          continue;  // don't combine pairs from one part of array
        }

        double cur_dist = EuclidDist(y_x_points[index], y_x_points[local_index]);
        if (fabs(cur_dist - result.dist) < EPS) {
          result.points_1.push_back(std::make_pair(y_x_points[index].second, y_x_points[index].first));
          result.points_2.push_back(std::make_pair(y_x_points[local_index].second, y_x_points[local_index].first));
        } 
        else if (cur_dist < result.dist) {
          result.reset();
          result.points_1.push_back(std::make_pair(y_x_points[index].second, y_x_points[index].first));
          result.points_2.push_back(std::make_pair(y_x_points[local_index].second, y_x_points[local_index].first));
          result.dist = cur_dist;
        }
      }
    }
  }
  return result;
}

}  // namespace NLabAlgOne