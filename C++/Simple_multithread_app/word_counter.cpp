
#include <memory>
#include <fstream>
#include <iostream>

#include "boost/bind.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem.hpp"

#include "word_counter.h"

using namespace word_counter_ns;

WordCounter::WordCounter(int processors_count, const std::string& documents_path_name)
    : processors_count_(processors_count),
      documents_path_name_(documents_path_name) ,
      words_frequency_(),
      no_words_(std::make_shared<int>(0)),
      no_unique_words_(0),
      documents_queue_(),
      processors_() {}

void WordCounter::Reconfigure(int processors_count, const std::string& documents_path_name) {
  processors_count_ = processors_count;
  documents_path_name_ = documents_path_name;
}

ResultType* WordCounter::ProceedWordCount() {
  if (processors_count_ <= 0) {
    throw IncorrectProcessorsCount("The number of processors must be positive int!\n");
  }

  CreateDocumentsQueue();

  while (static_cast<int>(processors_.size()) < processors_count_) {
    processors_.push_back(boost::thread(boost::bind(&ThreadFunction, this)));
  }

  for (auto iter = processors_.begin(); iter != processors_.end(); ++iter) {
    iter->join();
  }
  no_unique_words_ = words_frequency_.keys().size();

  std::cout << "Word counting process had finished!\n";
  return &words_frequency_;
}
  
void WordCounter::CreateDocumentsQueue() {
  boost::filesystem::path documents_path_(documents_path_name_);
  boost::filesystem::directory_iterator end_iter;

  if (!boost::filesystem::exists(documents_path_) ||
      !boost::filesystem::is_directory(documents_path_)) {
    throw IncorrectPath("Given path name is incorrect or not exists!\n");
  }

  bool has_any_document = false;
  for (boost::filesystem::directory_iterator dir_iter(documents_path_); 
      dir_iter != end_iter; ++dir_iter) {
    if (boost::filesystem::is_regular_file(*dir_iter) && 
        dir_iter->path().extension() == ".txt") {
      has_any_document = true;
      documents_queue_.push(dir_iter->path().filename().string());    
    }
  }
  if (!has_any_document) {
    throw HasNoDocuments("Given path doesn't contain any .txt file!");
  }  
}

void WordCounter::ThreadFunction(WordCounter* word_counter) {
  try {
    while (word_counter->documents_queue_.size() > 0) {
      std::string document_name;
      if (word_counter->documents_queue_.try_pop(&document_name)) {
        std::ifstream file_stream;
        file_stream.open((word_counter->documents_path_name_ + "\\" + document_name).c_str());

        if (file_stream.fail()) {
          throw FailedOpenFile("Can't open text file " + document_name + " \n");
        }

        std::string word;
        while (file_stream >> word) {
          bool has_key = word_counter->words_frequency_.has_key(word);
          if (has_key) {
            (*(word_counter->words_frequency_.get(word)))++;
          } else {
            word_counter->words_frequency_.set(word, std::make_shared<int>(1));
          }
          
          (*(word_counter->no_words_.get()))++;
        }
        if (!file_stream.eof()) {
          throw FileReadingError("Error while reading " + document_name + " had occurred!\n");
        }
      }
      boost::this_thread::sleep(boost::posix_time::millisec(1000));
    }
  } catch (FailedOpenFile &obj) {
    std::cout << obj.what();
  } catch (FileReadingError &obj) {
    std::cout << obj.what();
  } catch (...) {
    std::cout << "Some exception in ThreadFunction()!\n";
  }
  std::cout << "Thread had finished it's work!\n";
}
