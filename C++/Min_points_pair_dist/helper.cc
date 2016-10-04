#include "helper.h"

namespace NLabAlgOne {

void ReadPoints(TPoints* points_x_y, std::string data_path) {
  std::ifstream file_stream;
  file_stream.open(data_path);

  int num_points = 0;
  while (!file_stream.eof()) {
    std::string str;
    std::getline(file_stream, str);

    if (str.empty())
      continue;

    std::istringstream iss(str);
    std::vector<std::string> coordinates{ (std::istream_iterator<std::string>(iss)),
      (std::istream_iterator<std::string>()) };
    if (coordinates.size() == 1) {
      num_points = std::stoi(coordinates[0]);
    }
    else if (coordinates.size() == 2){
      points_x_y->push_back(std::make_pair(std::make_pair(std::stof(coordinates[0]), std::stof(coordinates[1])), 0));
    }
    else {
      num_points = std::stoi(coordinates[0]);

      for (int i = 1; i < coordinates.size(); i += 2) {
        points_x_y->push_back(std::make_pair(std::make_pair(std::stof(coordinates[i]), std::stof(coordinates[i + 1])), 0));
      }
      break;
    }
  }
}

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
    std::cout << "[ (" << result.points_1[i].first.first << ", " << result.points_1[i].first.second << "), ("
              << result.points_2[i].first.first << ", " << result.points_2[i].first.second << ")]" << std::endl;
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
      double dist = EuclidDist(iter->first, local_iter->first);
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

  TPoints current_delta_strip;
  for (int index = 0; index < y_x_points.size(); ++index) {
    if (std::fabs((x_y_begin + middle)->first.first - y_x_points[index].first.second) <= result.dist) {
      current_delta_strip.push_back(y_x_points[index]);
    }
  }

  for (int index = 0; index < current_delta_strip.size(); ++index) {
    if (current_delta_strip[index].second < right_bound && current_delta_strip[index].second >= left_bound) {
      for (auto local_index = index + 1; local_index != (index + 7); ++local_index) {
        if (local_index >= current_delta_strip.size())
          break;

        if ((current_delta_strip[local_index].second < middle && current_delta_strip[index].second < middle) ||
          (current_delta_strip[local_index].second >= middle && current_delta_strip[index].second >= middle)) {
          continue;  // don't combine pairs from one part of array
        }

        double cur_dist = EuclidDist(current_delta_strip[index].first, current_delta_strip[local_index].first);
        if (fabs(cur_dist - result.dist) < EPS) {
          result.points_1.push_back(std::make_pair(std::make_pair(
            current_delta_strip[index].first.second, current_delta_strip[index].first.first), current_delta_strip[index].second));
          result.points_2.push_back(std::make_pair(std::make_pair(
            current_delta_strip[local_index].first.second, current_delta_strip[local_index].first.first),
            current_delta_strip[local_index].second));
        } 
        else if (cur_dist < result.dist) {
          result.reset();
          result.points_1.push_back(std::make_pair(std::make_pair(
            current_delta_strip[index].first.second, current_delta_strip[index].first.first), current_delta_strip[index].second));
          result.points_2.push_back(std::make_pair(std::make_pair(
            current_delta_strip[local_index].first.second, current_delta_strip[local_index].first.first),
            current_delta_strip[local_index].second));
          result.dist = cur_dist;
        }
      }
    }
  }
  return result;
}

}  // namespace NLabAlgOne