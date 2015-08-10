#ifndef SRC_FILE_PROCESSOR_HPP_
#define SRC_FILE_PROCESSOR_HPP_

#include <memory>
#include <vector>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/filesystem/path.hpp"

#include "common.hpp"

class FileProcessor {
 public:
  FileProcessor(std::string path) 
    : text_info_(new WordsCountMap()), 
      words_occurrence_(new WordsCountMap()),
      path_(path) { }
          
  ~FileProcessor() { }
    
  WordsCountMap& CountWords(Filenames& filenames);
  WordsCountMap& CountWordCouplesOccurrence(
      Filenames& filenames, TopWords top_words);

 private:
  // contains info about collection's words and their counts
  std::shared_ptr<WordsCountMap> text_info_;
  // contains pair of top-words with counts of their co-occurrences
  std::shared_ptr<WordsCountMap> words_occurrence_;
  std::string path_;
};

#endif // SRC_FILE_PROCESSOR_HPP_
