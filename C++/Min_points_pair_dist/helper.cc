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

  auto out_array_iter = (left_array_iter == in_iter) ? temp_iter : in_iter;

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
  std::cout << "File = " << data_path << ", ";
  std::cout << "(" << result.point_1.first << ", " << result.point_1.second << "), ("
                   << result.point_2.first << ", " << result.point_2.second << "), ";
  std::cout << "Value = " << result.dist << std::endl;
}

float EuclidDist(const std::pair<float, float>& point_1, const std::pair<float, float>& point_2) {
  return std::sqrt(std::pow(point_1.first - point_2.first, 2) + std::pow(point_1.second - point_2.second, 2));
}

TResult BruteForce(TPsIter x_y_begin, int size) {
  float min_dist = INF;
  std::pair<float, float> min_point_1 = std::make_pair(-1.0f, -1.0f);
  std::pair<float, float> min_point_2 = std::make_pair(-1.0f, -1.0f);

  for (auto iter = x_y_begin; iter != (x_y_begin + size); ++iter) {
    auto local_iter = iter + 1;
    for (; local_iter != (x_y_begin + size); ++local_iter) {
      float dist = EuclidDist(*iter, *local_iter);
      if (dist < min_dist) {
        min_point_1 = *iter;
        min_point_2 = *local_iter;
        min_dist = dist;
      }
    }
  }
  return {min_point_1, min_point_2, min_dist};
}

TResult FindMinDist(TPsIter x_y_begin, TPsIter y_x_begin, TPsIter y_x_end, int left_bound, int right_bound) {
  int size = right_bound - left_bound;
  if (size <= MAX_POINTS_TO_BRUTE_FORCE) {
    return BruteForce(x_y_begin, size);
  }
  int middle = left_bound + static_cast<int>(size / 2);

  TResult res_left = FindMinDist(x_y_begin, y_x_begin, y_x_end, left_bound, middle);
  TResult res_right = FindMinDist(x_y_begin, y_x_begin, y_x_end, middle, right_bound);
  TResult cur_best_res = res_left.dist > res_right.dist ? res_right : res_left;

  for (auto iter = y_x_begin; iter != y_x_end; ++iter) {
    if (iter->second < (x_y_begin + right_bound)->first) {
      for (auto inner_iter = (iter + 1); inner_iter != (iter + 7) && inner_iter != y_x_end; ++inner_iter) {
        float cur_dist = EuclidDist(*iter, *inner_iter);
        if (cur_dist < cur_best_res.dist) {
          cur_best_res.point_1 = std::make_pair(iter->second, iter->first);
          cur_best_res.point_2 = std::make_pair(inner_iter->second, inner_iter->first);
          cur_best_res.dist = cur_dist;
        }
      }
    }
  }
  return cur_best_res;
}

}  // namespace NLabAlgOne