
#ifndef WC_COMMON_H_
#define WC_COMMON_H_

#define DEFINE_EXCEPTION_TYPE(Type, BaseType)                  \
class Type : public BaseType {                                 \
 public:                                                       \
  explicit Type() : BaseType("") {}                            \
  explicit Type(std::string message) : BaseType(message) {}    \
};

DEFINE_EXCEPTION_TYPE(IncorrectPath, std::runtime_error);
DEFINE_EXCEPTION_TYPE(IncorrectProcessorsCount, std::runtime_error);
DEFINE_EXCEPTION_TYPE(HasNoDocuments, std::runtime_error);
DEFINE_EXCEPTION_TYPE(FailedOpenFile, std::runtime_error);
DEFINE_EXCEPTION_TYPE(FileReadingError, std::runtime_error);

#define ResultType word_counter_ns::ThreadSafeCollectionHolder<std::string, int>

#endif // WC_COMMON_H_