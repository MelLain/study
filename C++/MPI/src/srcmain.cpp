#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "mpi.h"

#include "boost/date_time/posix_time/posix_time.hpp"

#include "common.hpp"
#include "file_processor.hpp"
#include "filename_reader.hpp"
#include "mpi_stl_interface.hpp"

namespace bfs = boost::filesystem;

void MoveData(std::vector<WordsCountMap>& sender, WordsCountMap& reciever);
void PrintData(WordsCountMap& source);
void SaveData(WordsCountMap& source);

int main(int argc, char* argv[]) {
  auto time_start = boost::posix_time::microsec_clock::local_time();

  MPI_Init(&argc, &argv);
  int rank = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int no_processors = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &no_processors);

  std::string data_path;
  int document_part_size;
  int no_top_words;

  if (no_processors < 2) {
    throw InvalidNumberOfProcessors(std::string("Number of ") +
        "processors should be less or equal 2\n");
  }

  switch (argc) {
   case 1 :
    data_path = DefaultDataPath;
    document_part_size = DefaultDocumentPartSize;
    no_top_words = DefaultNoTopWords;
    break;
   case 4:
    data_path = argv[1];
    document_part_size = std::atoi(argv[2]);
    no_top_words = std::atoi(argv[3]);
    break;
   default:
    throw InvalidNumberOfArguments(std::string("Usage: mpirun -n ") +
        std::string("<number of processes> ./srcmain <data pathname> ") +
        "<size of document batch> <number of top words>");
  }

  if (rank == 0) {
  // MASTER'S CODE
    if (!bfs::exists(data_path) || !bfs::is_directory(data_path)) {
      throw IncorrectPath(
          "Given path's is not exists, or it's name is incorrect\n");
    }
	  
    FilenameReader filename_reader(data_path);
    // will contain info about all words counts in collection
    WordsCountMap global_text_info;
    // will contain info about top-words counts in collection
    WordsCountMap words_occurrence;

    for (;;) {
      bool no_more_files = false;
      if (!filename_reader.has_more_documents()) {
        no_more_files = true;
      }

      MPI_Status status[no_processors - 1];
      MPI_Request requests[no_processors - 1];
      int continue_work = !no_more_files;

      // this messages will tell workers, if the process of word 
      // counting is over or no
      for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
        MPI_Isend(&continue_work, 1, MPI_INT, rank_idx,
            ContTag, MPI_COMM_WORLD, &(requests[rank_idx - 1]));
      }
      MPI_Waitall(no_processors - 1, requests, status);

      if (!no_more_files) {
        for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
          // this loop reads next part of filenames to process and send
          // them to next worker
          Filenames cur_filenames;
          if (filename_reader.has_more_documents()) {
            cur_filenames =
                filename_reader.GetNextFilenames(document_part_size);
          }
          SendVectorViaMpi(cur_filenames, rank_idx, SendTag);
        }

        // this vector will contain all counts from every worker, 
        // gained on this itertion
        std::vector<WordsCountMap> local_text_info;
        for (int i = 1; i < no_processors; ++i) {
          local_text_info.push_back(WordsCountMap());
        }

        // receive info from every worker
        for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
          RecvMapViaMpi(&(local_text_info[rank_idx - 1]), rank_idx, RecvTag);
        }

        // reduce all data into global storage
        MoveData(local_text_info, global_text_info);
      } else {
        // if all files are passed, go to second step
        break;
      }
    }

    // PrintData(global_text_info);
    std::cout << "Number of unique words: " <<
        global_text_info.size() << std::endl;
    auto time_finish = boost::posix_time::microsec_clock::local_time();
    auto time_for_first_step = (time_finish - time_start).total_milliseconds();
    std::cout << "Elapsed time on step #1: " <<  time_for_first_step
        << " msec.\n";

    int no_unique_words = global_text_info.size();
    if (no_top_words > no_unique_words) {
      std::cout << "Given number of top words is too large, all (" <<
          no_unique_words << ") unique words in collection "
          << "will be used\n";
      no_top_words = no_unique_words;
    }

    // this vectors will contain all words and their counts to find 
    // most frequent ones
    TopWords words;
    std::vector<int> words_counts;

    // this vector will contain most frequent words
    TopWords top_words;

    if (global_text_info.size() >= 1) {
      for (auto& elem : global_text_info) {
        words.push_back(elem.first);
        words_counts.push_back(elem.second);
      }

      // having vector with words and vector with corresponding counts,
      // function will sort second by values, and first --- by indices of
      // of sorting of second
      SortTwoVectorsByFirst<int, std::string>(
          words_counts, words, 0, words.size() - 1);

      for (int i = 0; i < no_top_words; ++i) {
        top_words.push_back(words[words.size() - 1 - i]);
      }
    }

    // top-words are sent to workers
    for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
      SendVectorViaMpi(top_words, rank_idx, SendTag); 	
    }

    filename_reader.reset();

    for (;;) {
      bool no_more_files = false;
      if (!filename_reader.has_more_documents()) {
        no_more_files = true;
      }

      MPI_Status status[no_processors - 1];
      MPI_Request requests[no_processors - 1];
      int continue_work = !no_more_files;

      // this messages will tell workers, if the process of word 
      // counting is over or no
      for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
        MPI_Isend(&continue_work, 1, MPI_INT, rank_idx, 
            ContTag, MPI_COMM_WORLD, &(requests[rank_idx - 1]));
      }
      MPI_Waitall(no_processors - 1, requests, status);

      if (!no_more_files) {
          // this loop reads next part of filenames to process and send
          // them to next worker
        for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
          Filenames cur_filenames;
          if (filename_reader.has_more_documents()) {
            cur_filenames =
              filename_reader.GetNextFilenames(document_part_size);
          }
          SendVectorViaMpi(cur_filenames, rank_idx, SendTag);
        }

        // this vector will contain all counts from every worker, 
        // gained on this itertion
        std::vector<WordsCountMap> local_text_info;
        for (int i = 1; i < no_processors; ++i) {
          local_text_info.push_back(WordsCountMap());
        }

        // receive info from every worker
        for (int rank_idx = 1; rank_idx < no_processors; ++rank_idx) {
          RecvMapViaMpi(&(local_text_info[rank_idx - 1]), rank_idx, RecvTag);
        }

        // reduce all data into global storage
        MoveData(local_text_info, words_occurrence);
      } else {
        // if all files are passed, finish processing
        break;
      }
    }
    // PrintData(words_occurrence);
    SaveData(words_occurrence);
    time_finish = boost::posix_time::microsec_clock::local_time();
    std::cout << "Elapsed time on step #2: " <<
        (time_finish - time_start).total_milliseconds() << " msec.\n";
    std::cout << "Total time: " << time_for_first_step +
        (time_finish - time_start).total_milliseconds() << " msec.\n";
  } else {
    // WORKER'S CODE
    FileProcessor file_processor(data_path);
    
    for (;;) {		
      MPI_Request request;
      MPI_Status status;
      int continue_work;
      Filenames cur_filenames;

      // receive message, which tells, if this step is over, or no
      MPI_Irecv(&continue_work, 1, MPI_INT, 0, ContTag, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);

      if (!continue_work) {
        // if all files are passed, go to second step
        break;
      }

      // recieve names of files to process
      RecvVectorViaMpi(&cur_filenames, 0, SendTag);
      // call processor and get the result
      auto& cur_text_info = file_processor.CountWords(cur_filenames);
      // send the result to master
      SendMapViaMpi(cur_text_info, 0, RecvTag);
    }

    TopWords top_words; 
    
    // receive top-words to process second step
    RecvVectorViaMpi(&top_words, 0, SendTag);

    for (;;) {
      MPI_Request request;
      MPI_Status status;
      int continue_work;
      Filenames cur_filenames;

      // receive message, which tells, if this step is over, or no
      MPI_Irecv(&continue_work, 1, MPI_INT, 0, ContTag, MPI_COMM_WORLD, &request);
      MPI_Wait(&request, &status);

      if (!continue_work) {
        // if all files are passed, finish processing
        break;
      }

      // recieve names of files to process
      RecvVectorViaMpi(&cur_filenames, 0, SendTag);

      // call processor and get the final result --- couples of 
      // top-words with co-occurrences
      auto& cur_text_info = file_processor.CountWordCouplesOccurrence(
          cur_filenames, top_words);
      // send the result to master
      SendMapViaMpi(cur_text_info, 0, RecvTag);
    }
  }
  
  MPI_Finalize(); 
  return 0;
}

// this function aggregates data from several maps into one
void MoveData(std::vector<WordsCountMap>& sender, WordsCountMap& reciever) {
  for (WordsCountMap& elem : sender) {
    for (auto iter = elem.begin(); iter != elem.end(); ++iter) {
      auto find_iter = reciever.find(iter->first);
      if (find_iter != reciever.end()) {
        find_iter->second += iter->second;
      } else {
        reciever.insert(std::pair<std::string, int>(
          iter->first, iter->second));
      }
    }
  }
}

void PrintData(WordsCountMap& source) {
  for (auto& elem : source) {
    std::cout << "Key '" << elem.first << "', value "
      << elem.second << std::endl;
  }
}

void SaveData(WordsCountMap& source) {
  std::ofstream res_file;
  res_file.open ("results.txt");
  for (auto& elem : source) {
    res_file << elem.first << ' ' << elem.second << std::endl;
  }
  res_file.close();
}
