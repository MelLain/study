#ifndef MATRIX_PROD_MATRIX_H
#define MATRIX_PROD_MATRIX_H

#include "raw.h"

namespace matrix_prod {

template<typename T>
class Matrix {
 public:
  Matrix() : no_raws_(0), no_cols_(0), no_reserved_(0), data_(nullptr) { };

  Matrix(const Matrix<T>& rhs)
    : no_raws_(rhs.no_raws()),
      no_cols_(rhs.no_cols()) {
    data_ = new Raw<T>*[rhs.no_raws()];
    for (int i = 0; i < rhs.no_raws(); ++i) {
      data_[i] = new Raw<T>(rhs.no_cols(), rhs[i][0]);
      for (int j = 0; j < rhs.no_cols(); ++j) {
        (*data_[i])[j] = rhs[i][j];
      }
    }
  }

  Matrix(uint no_raws, uint no_cols, T init_value)
    : no_raws_(no_raws),
      no_cols_(no_cols) {
    data_ = new Raw<T>*[no_raws];
    for (int i = 0; i < no_raws; ++i) {
      data_[i] = new Raw<T>(no_cols, init_value);
    }
  }

  const Raw<T>& operator[](uint index) const { return *data_[index]; }
  Raw<T>& operator[](uint index) { return *data_[index]; }

  Matrix<T>& operator=(const Matrix<T>& rhs) {
    no_raws_ = rhs.no_raws();
    no_cols_ = rhs.no_cols();
    data_ = new Raw<T>*[no_raws_];
    for (int i = 0; i < no_raws_; ++i) {
      data_[i] = new Raw<T>(no_cols_, rhs[i][0]);
      for (int j = 0; j < no_cols_; ++j) {
        (*data_[i])[j] = rhs[i][j];
      }
    }
    return *this;
  }

  ~Matrix() {
    for (int i = 0; i < no_raws_; ++i) delete data_[i];
    delete[] data_;
  }

  inline uint no_raws() const { return no_raws_; }
  inline uint no_cols() const { return no_cols_; }
  inline uint no_reserved() const { return no_reserved_; }

  void add_raw(const Raw<T>& raw) {
    if (no_cols_ == 0) {
      no_cols_ = raw.no_elem();
    } else {
      if (no_cols_ != raw.no_elem())
        throw std::exception("Matrix.add(): inconsistent sizes");
    }
  
    if (no_raws_ == no_reserved_) {
      no_reserved_ += INCREASE_RESERVE_STEP;
      Raw<T>** buffer = new Raw<T>*[no_reserved_];
      for (int i = 0; i < no_raws_; ++i) {
        buffer[i] = data_[i];
      }
      if (no_raws_ != 0) delete[] data_;
      data_ = buffer;
    }
    data_[no_raws_++] = new Raw<T>(raw.no_elem(), raw[0]);
    for (int i = 0 ; i < raw.no_elem(); ++i) {
      (*data_[no_raws_ - 1])[i] = raw[i];
    }
  }

 private:
   uint no_raws_;
   uint no_cols_;
   uint no_reserved_;
   Raw<T>** data_;
};
}  // matrix_prod

#endif // MATRIX_PROD_MATRIX_H