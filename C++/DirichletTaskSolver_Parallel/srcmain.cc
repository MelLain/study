#include <cmath>
#include <ctime>

#include <fstream>
#include <memory>
#include <vector>
#include <string>
#include <utility>

#include "mpi.h"

#include "common.h"
#include "matrix.h"
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
    int num_r = 0;
    int num_c = 0;
    
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

} // namespace

void print_results(const std::vector<std::shared_ptr<DM> >& values, const std::vector<std::shared_ptr<Grid> >& grids, Functions functions) {
  /*
  std::ofstream out_value_file;
  std::ofstream out_true_file;
  out_value_file.open(OUT_VALUE_FILE);
  out_true_file.open(OUT_TRUE_FILE);

  for (size_t i = 0; i < values.size(); ++i) {
    size_t start_shift = i > 0 ? 1 : 0;
    size_t end_shift = i < values.size() - 1 ? 1 : 0;
    
    for (size_t r = start_shift; r < values[i]->num_rows() - end_shift; ++r) {
      for (size_t c = 0; c < values[i]->num_cols(); ++c) {
        out_value_file << std::setw(15) << (*values[i])(r, c);
        out_true_file << std::setw(15) << functions.true_func((*grids[i])(r, c));
        //out_value_file << (*values[i])(r, c) << ", ";
        //out_true_file << functions.true_func(grid(r, c)) << ", ";
      }
      out_value_file << std::endl;
      out_true_file << std::endl;
    }
  }
  out_value_file.close();
  out_true_file.close();
  */
}

int main(int argc, char* argv[]) {
  double begin_time = MPI_Wtime();
  int rank = 0;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int num_processors = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

  try {
    int power_of_two = is_power_of_two(num_processors - 1);
    if (num_processors < 2 || power_of_two == -1) {
      throw std::runtime_error("Invalid number of procs (< 2 or not 2^q + 1)");
    }

    size_t grid_size = 0;
    if (argc != 2) {
      throw std::runtime_error("Incorrect usage! Use ./srcmain <grid_size>");
    } else {
      grid_size = std::stoi(argv[1]);
    }

    if (num_processors - 1 > static_cast<int>(grid_size * 0.5)) {
      throw std::runtime_error("Invalid number of procs (too many)");
    }

    size_t r_grid_size = grid_size;
    size_t c_grid_size = grid_size;
    GridData grid_data = { 0.0, 1.0, r_grid_size, 0.0, 1.0, c_grid_size, 1.0 };

    Functions functions = { [](const Point& p) { return 8 - 12 * pow(p.r_value, 2) - 12 * pow(p.c_value, 2); },
                            [](const Point& p) { return pow((1 - pow(p.r_value, 2)), 2) + pow((1 - pow(p.c_value, 2)), 2); },
                            [](const Point& p) { return pow((1 - pow(p.r_value, 2)), 2) + pow((1 - pow(p.c_value, 2)), 2); },
                            [](double value, double q) { return (pow(1 + value, q) - 1) / (pow(2.0, q) - 1); } };

    if (rank == 0) {
      size_t num_processed_iter = 0;
      double error = 1;
      while (true) {
        ++num_processed_iter;
        // start new iter
        send_flag_to_all(num_processors, START_ITER);

        // wait for parts of alpha, sum them and send back
        double alpha_den = collect_value_from_all(num_processors);
        send_value_to_all(num_processors, alpha_den);

        double alpha_nom = collect_value_from_all(num_processors);
        send_value_to_all(num_processors, alpha_nom);

	// wait for parts of tau, sum them and send back
        double tau_den = collect_value_from_all(num_processors);
        send_value_to_all(num_processors, tau_den);

        double tau_nom = collect_value_from_all(num_processors);
        send_value_to_all(num_processors, tau_nom);

        // wait for the end of the iteration and collect errors
        double difference = sqrt(collect_value_from_all(num_processors));

        if (difference < EPS) {
	  send_flag_to_all(num_processors, TERMINATE);

          // receive results
          error = sqrt(collect_value_from_all(num_processors));

          std::vector<std::shared_ptr<DM> > all_values;
          for (int i = 1; i < num_processors; ++i) {
            all_values.push_back(receive_matrix(i, i));
          }

          std::vector<std::shared_ptr<Grid> > all_grids;
          for (int i = 1; i < num_processors; ++i) {
            all_grids.push_back(receive_grid(i, i));
          }

          print_results(all_values, all_grids, functions);
	  break;
        }
      }

      std::cout << rank  << ": Finished! Elapsed time: " << MPI_Wtime() - begin_time << " sec." << std::endl
                << "Final error: " << error << std::endl << "Num iters processed: " << num_processed_iter << std::endl;
    } else {
      size_t num_row_procs = num_row_processors(grid_size, grid_size, static_cast<int>((num_processors - 1) * ((power_of_two % 2 == 0) ? 1 : 0.5)));
      num_row_procs *= power_of_two % 2 == 0 ? 1 : 2;
      size_t sub_rank = (rank - 1) % num_row_procs + 1;
      size_t int_num_rows = r_grid_size / num_row_procs;
      size_t mod_num_rows = r_grid_size % num_row_procs;
      size_t num_rows_to_process = int_num_rows + (sub_rank <= mod_num_rows ? 1 : 0);

      size_t r_start_idx = 0;
      if (sub_rank <= mod_num_rows) {
        r_start_idx = (sub_rank - 1) * (int_num_rows + 1);
      } else {
        r_start_idx = (int_num_rows + 1) * mod_num_rows + (int_num_rows * (sub_rank - mod_num_rows - 1));
      }

      size_t num_col_procs = (num_processors - 1) / num_row_procs;
      sub_rank = (rank - 1) / num_row_procs + 1;
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
        if (rank % num_row_procs == 1) {
          proc_bounds.is_up = true;
        }
        if (rank % num_row_procs == 0) {
          proc_bounds.is_low = true;
        }
        if (rank <= num_row_procs) {
          proc_bounds.is_left = true;
        }
        if (rank > num_row_procs * (num_col_procs - 1)) {
	  proc_bounds.is_right = true;
        }
      } else {
        proc_bounds = { true, true, true, true };
      }

      //std::cout << "R: " << rank << " * " << r_start_idx << " - " << r_start_idx + num_rows_to_process << " | " << c_start_idx << " - " << c_start_idx + num_cols_to_process << std::endl;

      // determines the order of sending-receiving mirror (rows, cols)
      auto first_send = std::make_pair((rank % num_row_procs) % 2 == 1, (rank / num_row_procs) % 2 == 0);

      // determines left and right neighbours of the current processor
      auto left_right_procs = std::make_pair(proc_bounds.is_left ? 0 : rank - num_row_procs,
                                             proc_bounds.is_right ? 0 : rank + num_row_procs);

      // each submatrix has both sides rows and cols mirroring
      auto model = GradientDescent(grid_data,
                                   functions,
                                   proc_bounds,
                                   rank,
                                   first_send,
                                   left_right_procs,
                                   r_start_idx,
                                   r_start_idx + num_rows_to_process,
                                   c_start_idx,
                                   c_start_idx + num_cols_to_process);
      model.FitModel();
    }
  } catch (const std::exception& e) {
    std::cout << rank << " : " << e.what() << std::endl;
  }

  MPI_Finalize();
  return 0;
}
