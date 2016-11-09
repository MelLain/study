#pragma once

#include <memory>
#include <vector>

#include <mpi.h>

#include "common.h"
#include "matrix.h"
#include "grid.h"

namespace DTS {

void send_matrix(const DM& values, int receiver_rank, int tag);
std::shared_ptr<DM> receive_matrix(int sender_rank, int tag);

void send_grid(const Grid& values, int receiver_rank, int tag);
std::shared_ptr<Grid> receive_grid(int sender_rank, int tag);

void send_vector(const std::vector<double>& values, int receiver_rank, int tag);
void receive_vector(std::vector<double>* values, int sender_rank, int tag);

void send_value(double value, int receiver_rank, int tag);
void receive_value(double* value, int sender_rank, int tag);

void send_flag(FlagType flag, int receiver_rank, int tag);
void receive_flag(FlagType* flag, int sender_rank, int tag);

double collect_value_from_all(int num_procs);
void send_value_to_all(int num_procs, double value);
void send_flag_to_all(int num_procs, FlagType flag);

}  // namespace DTS
