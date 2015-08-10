#ifndef SRC_COMMON_HPP_
#define SRC_COMMON_HPP_

#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

#include "boost/bind.hpp"
#include "boost/function.hpp"

#define DEFINE_EXCEPTION_TYPE(Type, BaseType)                  \
class Type : public BaseType {                                 \
 public:                                                       \
  explicit Type() : BaseType("") { }                           \
  explicit Type(std::string message) : BaseType(message) { }   \
};

DEFINE_EXCEPTION_TYPE(IncorrectPath, std::runtime_error);
DEFINE_EXCEPTION_TYPE(FailedOpenFile, std::runtime_error);
DEFINE_EXCEPTION_TYPE(FileReadingError, std::runtime_error);
DEFINE_EXCEPTION_TYPE(InvalidNumberOfArguments, std::runtime_error);
DEFINE_EXCEPTION_TYPE(InvalidNumberOfProcessors, std::runtime_error);

typedef std::pair<std::string, int> StrIntPair;
typedef std::map<std::string, int> WordsCountMap;
typedef std::vector<std::string> Filenames;
typedef std::vector<std::string> TopWords;

const std::string DefaultDataPath = "../data/";
const int DefaultDocumentPartSize = 3;
const int DefaultNoTopWords = 1;
const int ContTag = 0;
const int SendTag = 1;
const int RecvTag = 2;

// uses to sort the first vector by values and the second --- by the indices
// of sorting of first one (based on quick sort algorithm)
template<typename T_1, typename T_2>
  void SortTwoVectorsByFirst(std::vector<T_1>& src_1, 
      std::vector<T_2>& src_2, int first, int last) {
    
    int i = first;
    int j = last;
    int x = src_1[(first + last) / 2];
    
    do {
      while (src_1[i] < x) i++;
      while (src_1[j] > x) j--;
      
      if (i <= j) {
        if (src_1[i] > src_1[j]) {
          auto temp_1 = src_1[i];
          src_1[i] = src_1[j];
          src_1[j] = temp_1;
          
          auto temp_2 = src_2[i];
          src_2[i] = src_2[j];
          src_2[j] = temp_2;          
        }
        i++;
        j--;
      }
    } while (i <= j);

  if (i < last)
    SortTwoVectorsByFirst(src_1, src_2, i, last);
  if (first < j)
    SortTwoVectorsByFirst(src_1, src_2, first, j);
}

#endif // SRC_COMMON_HPP_
