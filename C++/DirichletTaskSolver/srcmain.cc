#include <ctime>

#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include <string>

#include "common.h"
#include "matrix.h"
#include "grid.h"
#include "gradient_descent.h"

using namespace DTS;

void PrintResults(const DM& values, const Grid& grid, Functions functions, const std::vector<double>& errors) {
  std::ofstream out_file;
  out_file.open(OUT_VALUE_FILE);
  for (size_t i = 0; i < values.num_rows(); ++i) {
    for (size_t j = 0; j < values.num_cols(); ++j) {
      out_file << std::setw(10) << values(i, j);
    }
    out_file << std::endl;
  }
  out_file.close();

  out_file.open(OUT_TRUE_FILE);
  for (size_t i = 0; i < grid.num_height_points(); ++i) {
    for (size_t j = 0; j < grid.num_width_points(); ++j) {
      out_file << std::setw(10) << functions.true_func(grid(i, j));
    }
    out_file << std::endl;
  }
  out_file.close();

  out_file.open(OUT_ERROR_FILE);
  for (double e : errors) {
    out_file << e << std::endl;
  }
  out_file.close();
}

const std::string HELP_STR = "usage: ./srcmain <grid_size>";

int main(int argc, char* argv[]) {
  const clock_t begin_time = clock();
  try {
    int grid_size = 0;
    if (argc != 2) {
      std::cout << HELP_STR << std::endl;
      return ARG_PARSE_ERROR;
    } else {
      grid_size = std::stoi(argv[1]);
    }

    Grid grid = Grid({ 0, 2, 0, 2 }, grid_size, grid_size);
    Functions functions = { [](const Point& p){ return (p.width * p.width +
                                                        p.height * p.height) * sin(p.height * p.width); },
                            [](const Point& p){ return 1.0 + sin(p.height * p.width); },
                            [](const Point& p){ return 1.0 + sin(p.height * p.width); } };

    auto model = GradientDescent(grid, functions);
    model.FitModel();
    PrintResults(model.GetCurrentMatrix(), grid, functions, model.GetErrors());

  } catch (const std::exception& e) {
    std::cout << e.what();
  }

  std::cout << "Finished! Elapsed time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << " sec." << std::endl;
  std::cout << "Press any key to continue..." << std::endl;
  getchar();
  return 0;
}
