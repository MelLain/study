#pragma once

#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <iostream>
#include <unordered_map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <fstream>
#include <thread>

#include "boost/lexical_cast.hpp"
#include <boost/filesystem.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/program_options.hpp>
#include "boost/thread.hpp"
#include "boost/thread/mutex.hpp"

const int CHECK_UPDATES_PERIOD = 1;
const std::string TEMP_TOKEN = "__temp__@#%^&__";

const int DEFAULT_NUM_TOPICS = 10;
const int DEFAULT_NUM_O_ITER = 10;
const int DEFAULT_NUM_I_ITER = 1;
const int DEFAULT_NUM_TOP_TOKENS = 10;
const int DEFAULT_NUM_PROCESSORS = 1;
const std::string DEFAULT_DATA_PATH = "";

const std::string DICTIONARY_NAME = "DICTIONARY.dict";
const std::string OUTPUT_FILENAME = "TOP_TOKENS.tokens";

typedef std::string Token;

class TFilenameQueue {
public:
  TFilenameQueue(const std::vector<std::string>& filenames) : data_(filenames) { }

  TFilenameQueue(const TFilenameQueue&) = delete;

  std::string pop_back() {
    boost::lock_guard<boost::mutex> guard(locker_);
    if (data_.size() < 1) {
      throw std::runtime_error("TFilenameQueue::pop_back: queue is empty");
    }
    std::string value = data_[data_.size() - 1];
    data_.pop_back();
    return value;
  }

  std::string push_back(const std::string& filename) {
    boost::lock_guard<boost::mutex> guard(locker_);
    data_.push_back(filename);
  }

  int size() {
    boost::lock_guard<boost::mutex> guard(locker_);
    return data_.size();
  }

  bool is_empty() {
    boost::lock_guard<boost::mutex> guard(locker_);
    return data_.empty();
  }

private:
  std::vector<std::string> data_;
  boost::mutex locker_;
};