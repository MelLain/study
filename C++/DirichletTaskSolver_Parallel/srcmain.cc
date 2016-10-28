#include <cmath>
#include <ctime>

#include <iomanip>
#include <fstream>
#include <memory>
#include <vector>
#include <string>

#include "mpi.h"

#include "common.h"
#include "matrix.h"
#include "grid.h"
#include "gradient_descent.h"
#include "mpi_helpers.h"

using namespace DTS;

void print_results(const std::vector<std::shared_ptr<DM> >& values, const Grid& grid, Functions functions) {
  std::ofstream out_file;
  out_file.open(OUT_VALUE_FILE);

  for (size_t i = 0; i < values.size(); ++i) {
    size_t start_shift = i > 0 ? 1 : 0;
    size_t end_shift = i < values.size() - 1 ? 1 : 0;
    
    for (size_t r = start_shift; r < values[i]->num_rows() - end_shift; ++r) {
      for (size_t c = 0; c < values[i]->num_cols(); ++c) {
        out_file << std::setw(15) << (*values[i])(r, c);
      }
      out_file << std::endl;
    }
  }
  out_file.close();

  out_file.open(OUT_TRUE_FILE);
  for (size_t i = 0; i < grid.num_height_points(); ++i) {
    for (size_t j = 0; j < grid.num_width_points(); ++j) {
      out_file << std::setw(15) << functions.true_func(grid(i, j));
    }
    out_file << std::endl;
  }
  out_file.close();
}

int main(int argc, char* argv[]) {
  const clock_t begin_time = clock();
  int rank = 0;

  try {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_processors = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);

    if (num_processors < 2) {
      throw std::runtime_error("Invalid number of procs (< 2)");
    }

    int grid_size = 0;
    if (argc != 2) {
      throw std::runtime_error("Incorrect usage! Use ./srcmain <grid_size>");
    } else {
      grid_size = std::stoi(argv[1]);
    }

    if (num_processors > (grid_size + 2)) {
      throw std::runtime_error("Invalid number of procs (too many)");
    }

    Grid grid = Grid({ 0, 1, 0, 1 }, grid_size, grid_size, 1.5);
    Functions functions = { [](const Point& p){ return 8 - 12 * pow(p.width, 2) - 12 * pow(p.height, 2); },
                            [](const Point& p){ return pow((1 - pow(p.width, 2)), 2) + pow((1 - pow(p.height, 2)), 2); },
                            [](const Point& p){ return pow((1 - pow(p.width, 2)), 2) + pow((1 - pow(p.height, 2)), 2); } };

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
          print_results(all_values, grid, functions);
          break;
        }
      }

      std::cout << rank  << ": Finished! Elapsed time: " << float(clock() - begin_time) / CLOCKS_PER_SEC << " sec." << std::endl
                << "Final error: " << error << std::endl << "Num iters processed: " << num_processed_iter << std::endl;
    } else {
      ProcType proc_type = GLOBAL_PROC;
      if (num_processors > 2) {
        proc_type = rank == 1 ? UPPER_PROC : (rank == num_processors - 1 ? LOWER_PROC : CENTER_PROC);
      }

      size_t int_num_rows = (grid_size + 1) / (num_processors - 1);
      size_t mod_num_rows = (grid_size + 1) % (num_processors - 1);
      size_t num_rows_to_process = int_num_rows + (rank <= mod_num_rows ? 1 : 0);

      size_t start_idx = 0;
      if (rank <= mod_num_rows) {
        start_idx = (rank - 1) * (int_num_rows + 1);
      } else {
        start_idx = (int_num_rows + 1) * mod_num_rows + (int_num_rows * (rank - mod_num_rows - 1));
      }

      // each submatrix has both sides rows mirroring 
      auto model = GradientDescent(grid, functions, proc_type, rank, start_idx, start_idx + num_rows_to_process);
      model.FitModel();
    }
  } catch (const std::exception& e) {
    std::cout << rank << " : " << e.what() << std::endl;
  }

  return 0;
}
