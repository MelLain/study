#include "helper.h"

using namespace NLabAlgOne;

int main(int argc, char* argv[]) {
  try {
    if (argc != 2) {
      throw std::runtime_error("Incorrect parameters. Try ./srcmain <data_filename>");
    }

    // read points
    std::ifstream file_stream;
    file_stream.open(argv[1]);

    int num_points = 0;
    TPoints points_x_y;
    TPoints points_y_x;
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
        points_x_y.push_back(std::make_pair(std::stof(coordinates[0]), std::stof(coordinates[1])));
        points_y_x.push_back(std::make_pair(std::stof(coordinates[1]), std::stof(coordinates[0])));
      }
      else {
        num_points = std::stoi(coordinates[0]);

        for (int i = 1; i < coordinates.size(); i += 2) {
          points_x_y.push_back(std::make_pair(std::stof(coordinates[i]), std::stof(coordinates[i + 1])));
          points_y_x.push_back(std::make_pair(std::stof(coordinates[i + 1]), std::stof(coordinates[i])));
        }
        break;
      }
    }
    // if there're too few points find the solution using brute force
    if (num_points <= MAX_POINTS_TO_BRUTE_FORCE) {
      PrintResult(BruteForce(points_x_y.begin(), points_x_y.end()), argv[1]);
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
    PrintResult(FindMinDist(points_x_y.begin(), points_y_x, 0, num_points), argv[1]);

  } catch (std::exception& e) {
    std::cout << e.what() << std::endl;
  }
return 0;
}