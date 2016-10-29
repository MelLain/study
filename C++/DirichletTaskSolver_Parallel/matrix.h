#pragma once

#include <memory>
#include <vector>

#include "common.h"
#include "grid.h"

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

  static std::shared_ptr<Matrix<T> > FivePointsLaplass(const Matrix<T>& src, const Grid& grid, size_t start_index) {
    auto retval = std::shared_ptr<Matrix<T> >(new Matrix<T>(src.num_rows(), src.num_cols(), 0.0));

    for (size_t i = 1; i < retval->num_rows() - 1; ++i) {
      for (size_t j = 1; j < retval->num_cols() - 1; ++j) {
        double part_1 = (src(i, j) - src(i - 1, j)) / grid.step_width(i + start_index) -
                        (src(i + 1, j) - src(i, j)) / grid.step_width(i + start_index + 1);
        double part_2 = (src(i, j) - src(i, j - 1)) / grid.step_height(j) -
                        (src(i, j + 1) - src(i, j)) / grid.step_height(j + 1);

        (*retval)(i, j) = 2 * part_1 / (grid.step_width(i + start_index) + grid.step_width(i + start_index + 1)) +
                          2 * part_2 / (grid.step_height(j) + grid.step_height(j + 1));
      }
    }

    return retval;
  }

  static double ProductByPointAndSum(const Matrix<double>& src_1, const Matrix<double>& src_2,
				     const Grid& grid, ProcType proc_type = GLOBAL_PROC, size_t start_index = 0) {
    double retval = 0.0;

    if (src_1.num_cols() != src_2.num_cols() || src_1.num_rows() != src_2.num_rows()) {
      throw std::runtime_error("Matrix::ProductByPointAndSum: inconsistent sizes of matrices");
    }

    size_t start_shift = (proc_type == UPPER_PROC || proc_type == GLOBAL_PROC) ? 0 : 1;
    size_t end_shift = (proc_type == LOWER_PROC || proc_type == GLOBAL_PROC) ? 0 : 1;
    for (size_t i = start_shift; i < src_1.num_rows() - end_shift; ++i) {
      for (size_t j = 0; j < src_1.num_cols(); ++j) {
        size_t i_real = i + start_index + 1 < grid.num_width_points() - 1 ? i + start_index + 1 : i + start_index;
        size_t j_real = j + 1 < grid.num_height_points() - 1 ? j + 1 : j;
        retval += 0.25 * (grid.step_width(i + start_index) + grid.step_width(i_real)) * (grid.step_height(j) + grid.step_height(j_real)) * src_1(i, j) * src_2(i, j);
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
