#include <cmath>

#include <memory>
#include <vector>
#include <string>

#include "mpi.h"

#include "common.h"
#include "grid.h"
#include "gradient_descent.h"
#include "mpi_helpers.h"

using namespace DTS;

namespace {
  int is_power_of_two(int number) {
    if (number <= 0) {
      return -1;
    }

    int retval = 0;
    while (number % 2 == 0) {
      ++retval;
      number >>= 1;
    }

    if (number >> 1) {
      return -1;
    } else {
      return retval;
    } 
  }

  int num_row_processors(int num_row_points, int num_col_points, int num_processors) {
    float row_var = static_cast<float>(num_row_points);
    float col_var = static_cast<float>(num_col_points);
    size_t num_r = 0;
    size_t num_c = 0;
    
    for (size_t i = 0; i < num_processors; ++i) {
      if (row_var >= col_var) {
        row_var *= 0.5;
        ++num_r;
      } else {
        col_var *= 0.5;
        ++num_c;
      }
      if (num_r * num_c == num_processors)
        break;
    }

    return num_r;
  }

  void save_info_file(size_t num_procs, double error, size_t num_processed_iter, size_t num_row_procs, size_t num_points, double e_time) {
  std::ofstream out_info_file;
  out_info_file.open("INFO_MODEL_" + std::to_string(num_points) + "_POINTS_" +  std::to_string(num_procs) + "_PROCS");
  out_info_file << "Num slave processors:\n" << num_procs << std::endl;
  out_info_file << "Final error:\n" << error << std::endl;
  out_info_file << "Num processed iters:\n" << num_processed_iter << std::endl;
  out_info_file << "Elapsed time:\n" << e_time << std::endl;
  out_info_file << "Num row processors:\n" << num_row_procs << std::endl;
  out_info_file << "Num points:\n" << num_points << std::endl;
  out_info_file.close();
}

} // namespace

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  double begin_time = MPI_Wtime();

  int num_processors = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

  try {
    int power_of_two = is_power_of_two(num_processors);

    if (num_processors < 1 || power_of_two == -1) {
      throw std::runtime_error("Invalid number of procs (not 2^q, q >= 0)");
    }

    size_t grid_size = 0;
    if (argc != 2) {
      throw std::runtime_error("Incorrect usage! Use <run_command> srcmain <grid_size>");
    } else {
      grid_size = std::stoi(argv[1]);
    }

    if (num_processors > static_cast<int>(grid_size * 0.5)) {
      throw std::runtime_error("Invalid number of procs (too many)");
    }

    size_t num_row_procs = num_row_processors(grid_size, grid_size,
      static_cast<int>((num_processors) * ((power_of_two % 2 == 0) ? 1 : 0.5)));
    num_row_procs *= power_of_two % 2 == 0 ? 1 : 2;

    size_t r_grid_size = grid_size;
    size_t c_grid_size = grid_size;
    GridData grid_data = { 0.0, 1.0, r_grid_size, 0.0, 1.0, c_grid_size, 1.5 };

    Functions functions = { [](const Point& p) { return 8 - 12 * pow(p.r_value, 2) - 12 * pow(p.c_value, 2); },
                            [](const Point& p) { return pow((1 - pow(p.r_value, 2)), 2) + pow((1 - pow(p.c_value, 2)), 2); },
                            [](const Point& p) { return pow((1 - pow(p.r_value, 2)), 2) + pow((1 - pow(p.c_value, 2)), 2); },
                            [](double value, double q) { return (pow(1 + value, q) - 1) / (pow(2.0, q) - 1); } };

    size_t sub_rank = (rank % num_row_procs) + 1;
    size_t int_num_rows = r_grid_size / num_row_procs;
    size_t mod_num_rows = r_grid_size % num_row_procs;
    size_t num_rows_to_process = int_num_rows + (sub_rank <= mod_num_rows ? 1 : 0);

    size_t r_start_idx = 0;
    if (sub_rank <= mod_num_rows) {
      r_start_idx = (sub_rank - 1) * (int_num_rows + 1);
    } else {
      r_start_idx = (int_num_rows + 1) * mod_num_rows + (int_num_rows * (sub_rank - mod_num_rows - 1));
    }

    size_t num_col_procs = num_processors / num_row_procs;
    sub_rank = rank / num_row_procs + 1;
    size_t int_num_cols = c_grid_size / num_col_procs;
    size_t mod_num_cols = c_grid_size % num_col_procs;
    size_t num_cols_to_process = int_num_cols + (sub_rank <= mod_num_cols ? 1 : 0);

    size_t c_start_idx = 0;
    if (sub_rank <= mod_num_cols) {
      c_start_idx = (sub_rank - 1) * (int_num_cols + 1);
    } else {
      c_start_idx = (int_num_cols + 1) * mod_num_cols + (int_num_cols * (sub_rank - mod_num_cols - 1));
    }

    ProcBounds proc_bounds = { false, false, false, false };
    if (num_processors > 2) {
      if (rank % num_row_procs == 0) {
        proc_bounds.is_up = true;
      }
      if ((rank + 1) % num_row_procs == 0) {
        proc_bounds.is_low = true;
      }
      if (rank < num_row_procs) {
        proc_bounds.is_left = true;
      }
      if (rank >= num_row_procs * (num_col_procs - 1)) {
        proc_bounds.is_right = true;
      }
    } else {
      proc_bounds = { true, true, true, true };
    }

    // each submatrix has both sides rows and cols mirroring
    auto model = GradientDescent(grid_data,
                                 functions,
                                 proc_bounds,
                                 num_processors,
                                 rank,
                                 num_row_procs,
                                 grid_size,
                                 r_start_idx,
                                 r_start_idx + num_rows_to_process,
                                 c_start_idx,
                                 c_start_idx + num_cols_to_process);
    model.FitModel();

    if (rank == 0) {
      std::cout << rank  << ": Finished! Elapsed time: " << MPI_Wtime() - begin_time << " sec." << std::endl
                << "Final error: " << model.Error() << std::endl << "Num iters processed: " << model.NumProcessedIter() << std::endl;

      save_info_file(num_processors, model.Error(), model.NumProcessedIter(), num_row_procs, grid_size, MPI_Wtime() - begin_time);
    }
  } catch (const std::exception& e) {
    std::cout << rank << " : " << e.what() << std::endl;
  }

  MPI_Finalize();
  return 0;
}
