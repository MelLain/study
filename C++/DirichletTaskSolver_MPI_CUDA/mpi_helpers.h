#pragma once

#include <memory>
#include <vector>

#include <mpi.h>

#include "common.h"

namespace DTS {

void send_vector(const std::vector<double>& values, int receiver_rank, int tag);
void receive_vector(std::vector<double>* values, int sender_rank, int tag);

 void send_receive_vector(const std::vector<double>& to_send, std::vector<double>* to_recv,
			  int send_rank, int recv_rank, int send_tag, int recv_tag);

void send_value(double value, int receiver_rank, int tag);
void receive_value(double* value, int sender_rank, int tag);

void send_flag(FlagType flag, int receiver_rank, int tag);
void receive_flag(FlagType* flag, int sender_rank, int tag);

double collect_value_from_all(int num_procs);
void send_value_to_all(int num_procs, double value);
void send_flag_to_all(int num_procs, FlagType flag);

}  // namespace DTS
