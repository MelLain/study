
#ifndef WC_WORD_COUNTER_H_
#define WC_WORD_COUNTER_H_

#include <string>

#include "boost/thread.hpp"

#include "common.h"
#include "thread_safe_holder.h"

namespace word_counter_ns {

class WordCounter : boost::noncopyable {
public:
  WordCounter(int processors_count, const std::string& documents_path_name);
  ~WordCounter() {}
  void Reconfigure(int processors_count, const std::string& documents_path_name);

  ResultType* ProceedWordCount();
  int GetNoWords() const { return *no_words_.get().get(); }
  int GetNoUniqueWords() const { return no_unique_words_; }

 private:
  int processors_count_;
  std::string documents_path_name_;
  ResultType words_frequency_;
  word_counter_ns::ThreadSafeHolder<int> no_words_;
  int no_unique_words_;
  word_counter_ns::ThreadSafeQueue<std::string> documents_queue_;
  std::vector<boost::thread> processors_;

  void CreateDocumentsQueue();
  static void WordCounter::ThreadFunction(WordCounter* word_counter);
};

} // namespace word_counter_ns

#endif // WC_WORD_COUNTER_H_