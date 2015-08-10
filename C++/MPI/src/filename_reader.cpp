#include <vector>
#include <string>

#include "common.hpp"
#include "filename_reader.hpp"

namespace bfs = boost::filesystem;

FilenameReader::FilenameReader(std::string path) : path_(path) {
  end_dir_iter_ = bfs::directory_iterator();
  cur_dir_iter_ = bfs::directory_iterator(path_);
  has_more_documents_ = (end_dir_iter_ != cur_dir_iter_);
}

void FilenameReader::reset() {
  end_dir_iter_ = bfs::directory_iterator();
  cur_dir_iter_ = bfs::directory_iterator(path_);
  has_more_documents_ = (end_dir_iter_ != cur_dir_iter_);
}

Filenames& FilenameReader::GetNextFilenames(
    int documents_part_size) {
  current_filenames_.clear();
  bool has_any_document = false;
  for (int counter = 0; 
       cur_dir_iter_ != end_dir_iter_ && counter < documents_part_size; 
       ++cur_dir_iter_) {
            
    if (bfs::is_regular_file(*cur_dir_iter_) && 
        cur_dir_iter_->path().extension() == ".txt") {
      has_any_document = true;
      current_filenames_.push_back(
        cur_dir_iter_->path().filename().string());
      ++counter;    
    }
  }

  if (!has_any_document) {
      has_more_documents_ = false;
  }
  return current_filenames_;
}
