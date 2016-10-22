#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>

#include "common.h"
#include "matrix.h"
#include "grid.h"
#include "gradient_descent.h"

void PrintResults(const DM& values, const Grid& grid, Functions functions, const std::vector<double>& errors) {
  std::ofstream out_file;
  out_file.open(OUT_VALUE_FILE);
  for (int i = 0; i < values.num_rows(); ++i) {
    for (int j = 0; j < values.num_cols(); ++j) {
      out_file << std::setw(10) << values(i, j);
    }
    out_file << std::endl;
  }
  out_file.close();

  out_file.open(OUT_TRUE_FILE);
  for (int i = 0; i < grid.num_height_points(); ++i) {
    for (int j = 0; j < grid.num_width_points(); ++j) {
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

int main(int argc, char* argv[]) {
  try {
    Grid grid = Grid({ 0, 2, 0, 2 }, 20, 20);
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

  std::cout << "Finished! Press any key to continue..." << std::endl;
  getchar();
  return 0;
}
