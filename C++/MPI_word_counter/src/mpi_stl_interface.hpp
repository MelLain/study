#ifndef MPI_STL_INTERFACE_HPP_
#define MPI_STL_INTERFACE_HPP_

#include <vector>
#include <map>
#include <string>

#include "common.hpp"

void SendVectorViaMpi(std::vector<std::string>& buffer, int recv_rank, int tag);
void RecvVectorViaMpi(std::vector<std::string>* buffer, int send_rank, int tag);

void SendMapViaMpi(WordsCountMap& buffer, int recv_rank, int tag);
void RecvMapViaMpi(WordsCountMap* buffer, int send_rank, int tag);

#endif // MPI_STL_INTERFACE_HPP_
