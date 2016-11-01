#pragma once

#include <memory>
#include <vector>

#include "common.h"

namespace DTS {

template <typename T>
class Matrix {
 public:
  Matrix() { clear(); }

  Matrix(size_t num_rows, size_t num_cols, const T& default_value = T())
    : num_rows_(num_rows)
    , num_cols_(num_cols)
  {
    std::vector<T> temp(num_cols_, default_value);

    for (size_t i = 0; i < num_rows; ++i) {
      data_.push_back(temp);
    }
  }

  Matrix(const Matrix<T>& rhs) { assign(rhs); }
  void operator=(const Matrix<T>& rhs) { assign(rhs); }

  void clear() {
    num_rows_ = 0;
    num_cols_ = 0;
    data_.clear();
  }

  size_t num_rows() const { return num_rows_; }
  size_t num_cols() const { return num_cols_; }

  T& operator()(size_t row_index, size_t col_index) { return data_[row_index][col_index]; }
  const T& operator()(size_t row_index, size_t col_index) const { return data_[row_index][col_index]; }

  const std::vector<T>& get_row(size_t index) const { return data_[index]; }
  std::vector<T>& get_row_non_const(size_t index) { return data_[index]; }

  void debug_print() const {
    for (size_t i = 0; i < num_rows(); ++i) {
      for (size_t j = 0; j < num_cols(); ++j) {
        std::cout << data_[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }

  std::shared_ptr<Matrix<T> > operator+(std::shared_ptr<Matrix<T> > rhs) const {
    return operator+(*rhs);
  }

  std::shared_ptr<Matrix<T> > operator+(const Matrix<T>& rhs) const {
    if (num_cols() != rhs.num_cols() || num_rows() != rhs.num_rows()) {
      throw std::runtime_error("Matrix::operator+: inconsistent sizes of matrices");
    }
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (size_t i = 0; i < num_rows(); ++i) {
      for (size_t j = 0; j < num_cols(); ++j) {
        (*retval)(i, j) = data_[i][j] + rhs(i, j);
      }
    }
    return retval;
  }

  std::shared_ptr<Matrix<T> > operator-(std::shared_ptr<Matrix<T> > rhs) const {
    return operator-(*rhs);
  }

  std::shared_ptr<Matrix<T> > operator-(const Matrix<T>& rhs) const {
    if (num_cols() != rhs.num_cols() || num_rows() != rhs.num_rows()) {
      throw std::runtime_error("Matrix::operator-: inconsistent sizes of matrices");
    }
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (size_t i = 0; i < num_rows(); ++i) {
      for (size_t j = 0; j < num_cols(); ++j) {
        (*retval)(i, j) = data_[i][j] - rhs(i, j);
      }
    }
    return retval;
  }

  std::shared_ptr<Matrix<T> > operator*(std::shared_ptr<Matrix<T> > rhs) const {
    return operator*(*rhs);
  }

  std::shared_ptr<Matrix<T> > operator*(const Matrix<T>& rhs) const {
    if (num_cols() != rhs.num_cols() || num_rows() != rhs.num_rows()) {
      throw std::runtime_error("Matrix::operator*: inconsistent sizes of matrices");
    }
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (size_t i = 0; i < num_rows(); ++i) {
      for (size_t j = 0; j < num_cols(); ++j) {
        (*retval)(i, j) = data_[i][j] * rhs(i, j);
      }
    }
    return retval;
  }

  std::shared_ptr<Matrix<T> > operator*(const T& val) const {
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (size_t i = 0; i < num_rows(); ++i) {
      for (size_t j = 0; j < num_cols(); ++j) {
        (*retval)(i, j) = data_[i][j] * val;
      }
    }
    return retval;
  }

 private:
  void assign(const Matrix<T>& rhs) {
    clear();

    num_rows_ = rhs.num_rows();
    num_cols_ = rhs.num_cols();
    data_ = rhs.data_;
  }

  size_t num_rows_;
  size_t num_cols_;
  std::vector<std::vector<T> > data_;
};

typedef Matrix<double> DM;

}  // namespace DTS
