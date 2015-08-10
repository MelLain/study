#ifndef SRC_FILENAME_READER_HPP_
#define SRC_FILENAME_READER_HPP_

#include <vector>
#include <string>

#include "boost/filesystem.hpp"
#include "boost/filesystem/path.hpp"

namespace bfs = boost::filesystem;

class FilenameReader {
 public:
  FilenameReader(std::string path);
  ~FilenameReader() { }
    
  Filenames& GetNextFilenames(int documents_part_size);
  bool has_more_documents() { return has_more_documents_; }
  void reset();

 private:
  bfs::path path_;
  bfs::directory_iterator cur_dir_iter_;
  bfs::directory_iterator end_dir_iter_;
  Filenames current_filenames_;  
  bool has_more_documents_;
};

#endif // SRC_FILENAME_READER_HPP_
