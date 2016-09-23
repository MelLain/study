#include "helper.h"

using namespace NLabAlgOne;

int main(int argc, char* argv[]) {
  try {
    if (argc != 2) {
      throw std::runtime_error("Incorrect parameters. Try ./srcmain <data_filename>");
    }

    // read points
    std::ifstream myfile;
    myfile.open(argv[1]);

    std::string str;
    std::getline(myfile, str);
    std::istringstream iss(str);

    std::vector<std::string> coordinates{(std::istream_iterator<std::string>(iss)),
                                         (std::istream_iterator<std::string>())};

    int num_points = std::stoi(coordinates[0]);
    TPoints points_x_y(num_points, std::make_pair(0.0f, 0.0f));
    TPoints points_y_x(num_points, std::make_pair(0.0f, 0.0f));

    int point_id = 0;
    for (int i = 1; i < coordinates.size(); i += 2, ++point_id) {
      points_x_y[point_id] = std::make_pair(std::stof(coordinates[i]), std::stof(coordinates[i + 1]));
      points_y_x[point_id] = std::make_pair(std::stof(coordinates[i + 1]), std::stof(coordinates[i]));
    }

    // if there're too few points find the solution using brute force
    if (num_points <= MAX_POINTS_TO_BRUTE_FORCE) {
      PrintResult(BruteForce(points_x_y.begin(), num_points), argv[1]);
      return 0;
    }

    // sort of the coordinates
    TPoints temp_x(num_points, std::make_pair(0.0f, 0.0f));
    auto sort_result_y_x = MergeSort(points_y_x.begin(), temp_x.begin(), 0, num_points - 1);
    TPoints temp_y(num_points, std::make_pair(0.0f, 0.0f));
    auto sort_result_x_y = MergeSort(points_x_y.begin(), temp_y.begin(), 0, num_points - 1);

    points_x_y.assign(sort_result_x_y, sort_result_x_y + num_points);
    points_y_x.assign(sort_result_y_x, sort_result_y_x + num_points);

    // find and print the common solution
    PrintResult(FindMinDist(points_x_y.begin(), points_y_x.begin(), points_y_x.end(), 0, num_points), argv[1]);

  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
return 0;
}