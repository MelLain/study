#include <fstream>
#include <vector>
#include <string>

#include "file_processor.hpp"
    
WordsCountMap& FileProcessor::CountWords(Filenames& filenames) {
  text_info_->clear();
  for (std::string& name : filenames) {
    std::ifstream file_stream;
    // open next file to process
    file_stream.open((path_ + "/" + name).c_str());
        
    if (file_stream.fail()) {
      throw FailedOpenFile("Can't open text file " + name + " \n");
    }
        
    std::string cur_word;       
    while (file_stream >> cur_word) {
      // read next word and increas it's count or create it, if 
      // that doesn't exist
      auto elem_iter = text_info_->find(cur_word);
      if (elem_iter != text_info_->end()) {
        elem_iter->second++;
      } else {
        text_info_->insert(std::pair<std::string, int>(cur_word, 1));
      }
    }        
    if (!file_stream.eof()) {
      throw FileReadingError("Error reading file" + name + "\n");
    }
  }
  return *text_info_;
}

WordsCountMap& FileProcessor::CountWordCouplesOccurrence(
    Filenames& filenames, TopWords top_words) {
  words_occurrence_->clear();
  int no_words = top_words.size();
  
  for (std::string& name : filenames) {
    // this map will ne used to get info about top-words from this file,
    // every position is initialized by one of top-words and zero
    WordsCountMap top_words_map;
    for (auto& elem : top_words) {
      top_words_map.insert(std::pair<std::string, int>(elem, 0));
    }
    // open next file to process
    std::ifstream file_stream;
    file_stream.open((path_ + "/" + name).c_str());
        
    if (file_stream.fail()) {
      throw FailedOpenFile("Can't open text file " + name + " \n");
    }

    // read words and fill top-words map, it is better for search by 
    // key than vector
    std::string cur_word;     
    while (file_stream >> cur_word) {
      auto iter = top_words_map.find(cur_word);
      if (iter != top_words_map.end()) {
        iter->second++;
      }
    }

    // map can't be used for sorting by values, so info about top-words
    // moves to vector
    std::vector<int> top_words_counts;
    for (auto& elem : top_words) {
      auto iter = top_words_map.find(elem);
      int value = 0;
      if (iter != top_words_map.end()) {
        value = iter->second;
	  }
	  top_words_counts.push_back(value);
    }
    // having vector with words and vector with corresponding counts,
    // function will sort second by values, and first --- by indices of
    // of sorting of second
    SortTwoVectorsByFirst<int, std::string>(
        top_words_counts, top_words, 0, no_words - 1);

    // such lookup allow to find all co-occurrences relatively easy
    for (int i = 0; i < no_words - 1; ++i) {
      for (int j = i + 1; j < no_words; ++j) {
        std::string key;
        // store names in lexicographical order to prevent duplication
        if (top_words[i] > top_words[j]) {
          key = top_words[i] + " " + top_words[j];
        } else {
          key = top_words[j] + " " + top_words[i];
        }
        int value = std::min(top_words_counts[i], top_words_counts[j]);
        
        // fill counts from this file into general result 
        auto iter = words_occurrence_->find(key);
        if (iter == words_occurrence_->end()) {
          words_occurrence_->insert(StrIntPair(key, value));
        } else {
          iter->second += value;
        }
      }
    } 
          
    if (!file_stream.eof()) {
      throw FileReadingError("Error reading file" + name + "\n");
    }
  }
  return *words_occurrence_;
}
