#include <utility>

#include "common.h"

namespace DTS {

  void Clear_CUDA(double* src);

  void PrintMatrix_CUDA(double* src, size_t num_rows, size_t num_cols);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void CopyMatrix_CUDA(double* dst, double* src, size_t num_rows, size_t num_cols);
  __global__ void CopyMatrixImpl(double* dst, double* src, size_t num_elems);

  void SumMatrices_CUDA(double* dst, double* src_1, double* src_2, size_t num_rows, size_t num_cols);
  __global__ void SumMatricesImpl(double* dst, double* src_1, double* src_2, size_t num_elems);

  void DiffMatrices_CUDA(double* dst, double* src_1, double* src_2, size_t num_rows, size_t num_cols);
  __global__ void DiffMatricesImpl(double* dst, double* src_1, double* src_2, size_t num_elems);

  void ProdMatrixByScalar_CUDA(double* dst, double* src, double alpha, size_t num_rows, size_t num_cols);
  __global__ void ProdMatrixByScalarImpl(double* dst, double* src, double alpha, size_t num_elems);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void FivePointsLaplass_CUDA(double* dst, double* src, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);
  __global__ void FivePointsLaplassImpl(double* dst, double* src, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);

  double ProductByPointAndSum_CUDA(double* src_1, double* src_2, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);
  
  __global__ void ProductByPointAndSumImpl(double* pre_retval, double* src_1, double* src_2,
                                           double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void InitValues_CUDA(double** values, double* grid_x, double* grid_y,
                       size_t num_rows, size_t num_cols, const ProcBounds& proc_bounds);
  
  __global__ void InitValuesImpl(double* dst, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols,
                                 bool is_up, bool is_low, bool is_left, bool is_right);

  void InitGrid_CUDA(double** grid_x, double** grid_y, size_t* num_rows, size_t* num_cols,
                     const GridData& grid_data, const ProcBounds& proc_bounds,
                     size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx);
  
  void InitMatrix_CUDA(double** matrix, size_t num_rows, size_t num_cols);
  __global__ void InitMatrixImpl(double* dst, size_t num_elems);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void ExchangeRowsCols_CUDA(double* matrix, size_t num_rows, size_t num_cols, const ProcBounds& proc_bounds,
                             std::pair<bool, bool> first_send, std::pair<int, int> left_right_proc, size_t proc_rank);

  __global__ void ExtractCol_CUDA(double* dst, double* src, size_t col, size_t num_rows, size_t num_cols);
  __global__ void InsertCol_CUDA(double* dst, double* src, size_t col, size_t num_rows, size_t num_cols);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void CountResiduals_CUDA(double* residuals, double* values_laplass, double* grid_x, double* grid_y,
                           size_t num_rows, size_t num_cols, const ProcBounds& proc_bounds);

  __global__ void CountResidualsImpl(double* residuals, double* values_laplass, double* grid_x, double* grid_y,
                                     size_t num_rows, size_t num_cols, bool is_up, bool is_low, bool is_left, bool is_right);

  double CountPreError_CUDA(double* values, double* temp_matrix, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);

  __global__ void CountPreErrorImpl(double* psi, double* values, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void CountNewValues_CUDA(double* values, double* old_values, double* gradients,
                           double* temp_matrix, double tau, size_t num_rows, size_t num_cols);

  double CountValuesDifference_CUDA(double* values, double* old_values, double* temp_matrix,
                                    double* grid_x, double* grid_y, size_t num_rows, size_t num_cols);

  void CountGradients_CUDA(double* gradients, double* residuals, double* temp_matrix,
                           double alpha, size_t num_rows, size_t num_cols, bool first_iter);

} // namespace DTS
