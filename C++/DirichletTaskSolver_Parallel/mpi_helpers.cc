#include "mpi_helpers.h"

namespace DTS {

void send_matrix(const DM& values, int receiver_rank, int tag) {
  int num_rows = values.num_rows();
  int num_cols = values.num_cols();
  MPI_Send(&num_rows, 1, MPI_INT, receiver_rank, tag, MPI_COMM_WORLD);
  MPI_Send(&num_cols, 1, MPI_INT, receiver_rank, tag, MPI_COMM_WORLD);
    for (int i = 0; i < num_rows; ++i) {
      send_vector(values.get_row(i), receiver_rank, tag);
    }
}

std::shared_ptr<DM> receive_matrix(int sender_rank, int tag) {
  MPI_Status status;
  int num_rows;
  int num_cols;
  MPI_Recv(&num_rows, 1, MPI_INT, sender_rank, tag, MPI_COMM_WORLD, &status);
  MPI_Recv(&num_cols, 1, MPI_INT, sender_rank, tag, MPI_COMM_WORLD, &status);

  auto values = std::shared_ptr<DM>(new DM(num_rows, num_cols, 0.0));
  for (int i = 0; i < num_rows; ++i) {
    receive_vector(&(values->get_row_non_const(i)), sender_rank, tag);
  }
  return values;
}

void send_grid(const Grid& values, int receiver_rank, int tag) {
  int num_rows = values.num_rows();
  int num_cols = values.num_cols();
  MPI_Send(&num_rows, 1, MPI_INT, receiver_rank, tag, MPI_COMM_WORLD);
  MPI_Send(&num_cols, 1, MPI_INT, receiver_rank, tag, MPI_COMM_WORLD);
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      send_value(values(i, j).r_value, receiver_rank, tag);
      send_value(values(i, j).c_value, receiver_rank, tag);
    }
  }
}

std::shared_ptr<Grid> receive_grid(int sender_rank, int tag) {
  MPI_Status status;
  int num_rows;
  int num_cols;
  MPI_Recv(&num_rows, 1, MPI_INT, sender_rank, tag, MPI_COMM_WORLD, &status);
  MPI_Recv(&num_cols, 1, MPI_INT, sender_rank, tag, MPI_COMM_WORLD, &status);

  auto values = std::shared_ptr<Grid>(new Grid(num_rows, num_cols));
  double row;
  double col;
  for (int i = 0; i < num_rows; ++i) {
    for (int j = 0; j < num_cols; ++j) {
      receive_value(&row, sender_rank, tag);
      receive_value(&col, sender_rank, tag);
      (*values)(i, j) = { row, col };
    }
  }
  return values;
}

void send_vector(const std::vector<double>& values, int receiver_rank, int tag) {
  int size = values.size();
  MPI_Send(&size, 1, MPI_INT, receiver_rank, tag, MPI_COMM_WORLD);
  MPI_Send(&values[0], size, MPI_DOUBLE, receiver_rank, tag, MPI_COMM_WORLD);
}

void receive_vector(std::vector<double>* values, int sender_rank, int tag) {
  MPI_Status status;
  int size;
  MPI_Recv(&size, 1, MPI_INT, sender_rank, tag, MPI_COMM_WORLD, &status);
  values->resize(size);
  MPI_Recv(&(*values)[0], size, MPI_DOUBLE, sender_rank, tag, MPI_COMM_WORLD, &status);
}

void send_value(double value, int receiver_rank, int tag) {
  MPI_Send(&value, 1, MPI_DOUBLE, receiver_rank, tag, MPI_COMM_WORLD);
}

void receive_value(double* value, int sender_rank, int tag) {
  MPI_Status status;
  MPI_Recv(value, 1, MPI_DOUBLE, sender_rank, tag, MPI_COMM_WORLD, &status);
}

void send_flag(FlagType flag, int receiver_rank, int tag) {
  MPI_Send(&flag, 1, MPI_INT, receiver_rank, tag, MPI_COMM_WORLD);
}

void receive_flag(FlagType* flag, int sender_rank, int tag) {
  MPI_Status status;
  MPI_Recv(flag, 1, MPI_DOUBLE, sender_rank, tag, MPI_COMM_WORLD, &status);
}

double collect_value_from_all(int num_procs) {
  double value = 0;
  for (size_t i = 1; i < num_procs; ++i) {
    double value_part;
    receive_value(&value_part, i, i);
    value += value_part;
  }
  return value;
}

void send_value_to_all(int num_procs, double value) {
  for (size_t i = 1; i < num_procs; ++i) {
    send_value(value, i, i);
  }
}

void send_flag_to_all(int num_procs, FlagType flag) {
  for (size_t i = 1; i < num_procs; ++i) {
    send_flag(flag, i, i);
  }
}

} // namespace DTS
