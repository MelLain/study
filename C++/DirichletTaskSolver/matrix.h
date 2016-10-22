#pragma once

#include <memory>
#include <vector>

#include "common.h"
#include "grid.h"

template <typename T>
class Matrix {
 public:
  Matrix() { clear(); }

  Matrix(int num_rows, int num_cols, const T& default_value = T())
    : num_rows_(num_rows)
    , num_cols_(num_cols)
  {
    std::vector<T> temp(num_cols_, default_value);

    for (int i = 0; i < num_rows; ++i) {
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

  int num_rows() const { return num_rows_; }
  int num_cols() const { return num_cols_; }

  T& operator()(int row_index, int col_index) { return data_[row_index][col_index]; }
  const T& operator()(int row_index, int col_index) const { return data_[row_index][col_index]; }

  void debug_print() const {
    for (int i = 0; i < num_rows(); ++i) {
      for (int j = 0; j < num_cols(); ++j) {
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
      throw std::exception("Matrix::operator+: inconsistent sizes of matrices");
    }
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (int i = 0; i < num_rows(); ++i) {
      for (int j = 0; j < num_cols(); ++j) {
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
      throw std::exception("Matrix::operator-: inconsistent sizes of matrices");
    }
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (int i = 0; i < num_rows(); ++i) {
      for (int j = 0; j < num_cols(); ++j) {
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
      throw std::exception("Matrix::operator*: inconsistent sizes of matrices");
    }
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (int i = 0; i < num_rows(); ++i) {
      for (int j = 0; j < num_cols(); ++j) {
        (*retval)(i, j) = data_[i][j] * rhs(i, j);
      }
    }
    return retval;
  }

  std::shared_ptr<Matrix<T> > operator*(const T& val) const {
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix(num_rows(), num_cols(), 0.0));

    for (int i = 0; i < num_rows(); ++i) {
      for (int j = 0; j < num_cols(); ++j) {
        (*retval)(i, j) = data_[i][j] * val;
      }
    }
    return retval;
  }

  static std::shared_ptr<Matrix<T> > FivePointsLaplass(const Matrix<T>& src, const Grid& grid) {
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix<T>(src.num_rows(), src.num_cols(), 0.0));

    for (int i = 1; i < retval->num_rows() - 1; ++i) {
      for (int j = 1; j < retval->num_cols() - 1; ++j) {
        double part_1 = (src(i, j) - src(i - 1, j)) / grid.step_height() -
                        (src(i + 1, j) - src(i, j)) / grid.step_height();
        double part_2 = (src(i, j) - src(i, j - 1)) / grid.step_width() -
                        (src(i, j + 1) - src(i, j)) / grid.step_width();
        (*retval)(i, j) = part_1 / grid.step_height() + part_2 / grid.step_width();
      }
    }

    return retval;
  }

  static double ProductByPointAndSum(const Matrix<double>& src_1, const Matrix<double>& src_2, const Grid& grid) {
    double retval = 0.0;

    if (src_1.num_cols() != src_2.num_cols() || src_1.num_rows() != src_2.num_rows()) {
      throw std::exception("Matrix::ProductByPointAndSum: inconsistent sizes of matrices");
    }

    for (int i = 0; i < src_1.num_cols(); ++i) {
      for (int j = 0; j < src_2.num_rows(); ++j) {
        retval += grid.step_height() * grid.step_width() * src_1(i, j) * src_2(i, j);
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

  int num_rows_;
  int num_cols_;
  std::vector<std::vector<T> > data_;
};
