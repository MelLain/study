#include "word_counter.h"

int main(int argc, char* argv[]) {
  try {
    auto time_start = boost::posix_time::microsec_clock::local_time();

    word_counter_ns::WordCounter word_counter(1, argv[0]);

    std::cout << "Start!\n";
    auto retval = word_counter.ProceedWordCount();
    auto time_finish = boost::posix_time::microsec_clock::local_time();
    std::cout << "============================================\n";
    std::cout << "Total number of words  = " << word_counter.GetNoWords() <<  std::endl;
    std::cout << "Number of unique words = " << word_counter.GetNoUniqueWords() << std::endl;
    std::cout << "Elapsed time           = " << 
      (time_finish - time_start).total_milliseconds() << " milliseconds.\n";
    std::cout << "============================================\n";

    //auto keys = retval->keys();
    //for (auto iter = keys.begin(); iter != keys.end(); ++iter) {
    //  std::cout << *iter << ' ' << *(retval->get(*iter)) << std::endl;
    //}
  } catch (IncorrectProcessorsCount &obj) {
    std::cout << obj.what();
  } catch (IncorrectPath &obj) {
    std::cout << obj.what();
  } catch (HasNoDocuments &obj) {
    std::cout << obj.what();
  } catch (...) {
    std::cout << "Some exception in main()!\n";
  }

  char c;
  std::cin >> c;
  return 0;
}