#include "helper.h"

using namespace NLabAlgOne;

int main(int argc, char* argv[]) {
  try {
    if (argc != 2) {
      throw std::runtime_error("Incorrect parameters. Try ./srcmain <data_filename>");
    }

    TPoints points_x_y;
    ReadPoints(&points_x_y, argv[1]);
    int num_points = points_x_y.size();

    // if there're too few points find the solution using brute force
    if (num_points <= MAX_POINTS_TO_BRUTE_FORCE) {
      PrintResult(BruteForce(points_x_y.begin(), points_x_y.end()), argv[1]);
      return 0;
    }

    // sort of the coordinates
    TPoints temp_x(num_points, std::make_pair(std::make_pair(0.0f, 0.0f), 0));
    auto sort_result_x_y = MergeSort(points_x_y.begin(), temp_x.begin(), 0, num_points - 1);
    
    TPoints points_y_x;
    for (int i = 0; i < num_points; ++i) {
      sort_result_x_y[i].second = i;
      points_y_x.push_back(std::make_pair(std::make_pair(
        sort_result_x_y[i].first.second, sort_result_x_y[i].first.first), sort_result_x_y[i].second));
    }
    
    TPoints temp_y(num_points, std::make_pair(std::make_pair(0.0f, 0.0f), 0));
    auto sort_result_y_x = MergeSort(points_y_x.begin(), temp_y.begin(), 0, num_points - 1);

    points_x_y.assign(sort_result_x_y, sort_result_x_y + num_points);
    points_y_x.assign(sort_result_y_x, sort_result_y_x + num_points);

    // find and print the common solution
    PrintResult(FindMinDist(points_x_y.begin(), points_y_x, 0, num_points), argv[1]);

  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
return 0;
}