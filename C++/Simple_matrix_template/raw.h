#ifndef MATRIX_PROD_RAW_H
#define MATRIX_PROD_RAW_H

namespace matrix_prod {
typedef unsigned int uint;
const int INCREASE_RESERVE_STEP = 2;

template<typename T>
class Raw {
 public:
  Raw() : no_elem_(0), no_reserved_(0), data_(nullptr) { };

  Raw(const Raw<T>& rhs)
    : no_elem_(rhs.no_elem()),
      no_reserved_(rhs.no_reserved_()) {
    data_ = new T[rhs.no_reserved_()];
    for (int i = 0; i < rhs.no_elem(); ++i) {
      data_[i] = new T();
      *data[i] = rhs[i]
    }
  }

  Raw(uint no_elem, T init_value)
    : no_elem_(no_elem),
      no_reserved_(no_elem) {
    data_ = new T*[no_elem];
    for (int i = 0; i < no_elem; ++i) {
      data_[i] = new T();      
      *data_[i] = init_value;
    }
  }
  
  const T& operator[](uint index) const { return *data_[index]; }
  T& operator[](uint index) { return *data_[index]; }
  
  Raw<T>& operator=(const Raw<T>& rhs) {
    no_elem_ = rhs.no_elem();
    no_reserved_ = rhs.no_reserved();
    data_ = new T[no_reserved_];
    for (int i = 0; i < no_elem_; ++i) {
      data_[i] = new T();
      *data_[i] = rhs[i];
    }
    return *this;
  }

  ~Raw() {
    for (int i = 0; i < no_elem_; ++i) delete data_[i];
    delete[] data_;
  }

  inline uint no_elem() const { return no_elem_; }
  inline uint no_reserved() const { return no_reserved_; }
  
  void add(const T& elem) {
    if (no_elem_ == no_reserved_) {
      no_reserved_ += INCREASE_RESERVE_STEP;
      T** buffer = new T*[no_reserved_];
      for (int i = 0; i < no_elem_; ++i) {
        buffer[i] = data_[i];
      }
      if (no_elem_ != 0) delete[] data_;
      data_ = buffer;
    }
    data_[no_elem_++] = new T();
    *data_[no_elem_ - 1] = elem;
  }

  void clear() {
    for (int i = 0; i < no_elem_; ++i) delete data_[i];
    delete[] data_;
    no_elem_ = 0;
    no_reserved_ = 0;
    data_ = nullptr;
  }

 private:
   uint no_elem_;
   uint no_reserved_;
   T** data_;
};
}  // matrix_prod

#endif // MATRIX_PROD_RAW_H