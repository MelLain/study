#include "mpi_helpers.h"

namespace DTS {

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

void send_receive_vector(const std::vector<double>& to_send, std::vector<double>* to_recv,
                         int send_rank, int recv_rank, int send_tag, int recv_tag) {
  MPI_Status status;
  int size = to_send.size();
  MPI_Sendrecv(&to_send[0], size, MPI_DOUBLE, recv_rank, send_tag,
               &(*to_recv)[0], size, MPI_DOUBLE, send_rank, recv_tag, MPI_COMM_WORLD, &status);
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
