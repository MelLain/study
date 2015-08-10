#include "mpi.h"

#include "file_processor.hpp"

void SendVectorViaMpi(std::vector<std::string>& buffer, int recv_rank, int tag) {
  MPI_Status status;
  MPI_Request request;
  int size = buffer.size();
  MPI_Isend(&size, 1, MPI_INT, recv_rank,
      tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);

  for (auto& elem : buffer) {
    MPI_Isend(const_cast<char*>(elem.c_str()), elem.size(), MPI_CHAR, recv_rank,
        tag, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);
  }
}

void RecvVectorViaMpi(std::vector<std::string>* buffer, int send_rank, int tag) {
  MPI_Status status;
  MPI_Request request;
  int size;
  MPI_Irecv(&size, 1, MPI_INT, send_rank,
      tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  
  buffer->clear();
  for (int i = 0; i < size; ++i) {
    MPI_Probe(send_rank, tag, MPI_COMM_WORLD, &status);
    int len;
    MPI_Get_count(&status, MPI_CHAR, &len);
    char* buf = new char[len];
	
    MPI_Irecv(buf, len, MPI_CHAR, send_rank,
        tag, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);    
    
    std::string elem(buf, len);
    delete [] buf;
    buffer->push_back(elem);
  }
}

void SendMapViaMpi(WordsCountMap& buffer, int recv_rank, int tag) {
  MPI_Status status;
  MPI_Request request;
  int size = buffer.size();
  MPI_Isend(&size, 1, MPI_INT, recv_rank,
      tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);

  for (auto& elem : buffer) {
    MPI_Isend(const_cast<char*>(elem.first.c_str()), elem.first.size(), MPI_CHAR, recv_rank,
        tag, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);

    MPI_Isend(&(elem.second), 1, MPI_INT, recv_rank,
        tag, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);
  }	
}

void RecvMapViaMpi(WordsCountMap* buffer, int send_rank, int tag) {
  MPI_Status status;
  MPI_Request request;
  int size;
  MPI_Irecv(&size, 1, MPI_INT, send_rank,
      tag, MPI_COMM_WORLD, &request);
  MPI_Wait(&request, &status);
  
  buffer->clear();
  for (int i = 0; i < size; ++i) {
    MPI_Probe(send_rank, tag, MPI_COMM_WORLD, &status);
    int len;
    MPI_Get_count(&status, MPI_CHAR, &len);
    char* buf = new char[len];
	
    MPI_Irecv(buf, len, MPI_CHAR, send_rank,
        tag, MPI_COMM_WORLD, &request);
    MPI_Wait(&request, &status);    
    
    std::pair<std::string, int> elem(std::string(buf, len), -1);
    delete [] buf;
    
    int value;
    MPI_Irecv(&value, 1, MPI_INT, send_rank,
        tag, MPI_COMM_WORLD, &request);    
    MPI_Wait(&request, &status); 
        
    elem.second = value;
    buffer->insert(elem);
  }
}
