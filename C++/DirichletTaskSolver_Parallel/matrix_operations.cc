#include "matrix_operations.h"

namespace DTS {

std::shared_ptr<DM> FivePointsLaplass(const DM& src, const Grid& grid) {
  auto retval = std::shared_ptr<DM>(new DM(src.num_rows(), src.num_cols(), 0.0));

  for (size_t i = 1; i < retval->num_rows() - 1; ++i) {
    for (size_t j = 1; j < retval->num_cols() - 1; ++j) {
      double part_1 = (src(i, j) - src(i - 1, j)) / grid.r_step(i) - (src(i + 1, j) - src(i, j)) / grid.r_step(i + 1);
      double part_2 = (src(i, j) - src(i, j - 1)) / grid.c_step(j) - (src(i, j + 1) - src(i, j)) / grid.c_step(j + 1);

      (*retval)(i, j) = 2 * part_1 / (grid.r_step(i) + grid.r_step(i + 1)) + 2 * part_2 / (grid.c_step(j) + grid.c_step(j + 1));
    }
  }

  return retval;
}

double ProductByPointAndSum(const DM& src_1, const DM& src_2, const Grid& grid) {
  double retval = 0.0;

  if (src_1.num_cols() != src_2.num_cols() || src_1.num_rows() != src_2.num_rows()) {
    throw std::runtime_error("ProductByPointAndSum: inconsistent sizes of matrices");
  }

  for (size_t i = 1; i < src_1.num_rows() - 1; ++i) {
    for (size_t j = 1; j < src_1.num_cols() - 1; ++j) {
      size_t i_real = i + 1 < grid.num_rows() - 1 ? i + 1 : i;
      size_t j_real = j + 1 < grid.num_cols() - 1 ? j + 1 : j;

      retval += 0.25 * (grid.r_step(i) + grid.r_step(i_real)) * (grid.c_step(j) + grid.c_step(j_real)) * src_1(i, j) * src_2(i, j);
    }
  }

  return retval;
}

}  // namespace DTS
