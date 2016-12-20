#include <vector>

#include "mpi_helpers.h"
#include "cuda_operations.h"

namespace DTS {

void Clear_CUDA(double* src) {
  CheckCuda(cudaFree(src));
}

void PrintMatrix_CUDA(double* src, size_t num_rows, size_t num_cols) {
  size_t num_elems = num_rows * num_cols;
  double* dst = new double[num_elems];
  CheckCuda(cudaMemcpy(dst, src, num_elems * sizeof(double), cudaMemcpyDeviceToHost));
  for (size_t i = 0; i < num_rows; ++i) {
    for (size_t j = 0; j < num_cols; ++j) {
      std::cout << dst[i * num_cols + j] << "  ";
    }
    std::cout << std::endl;
  }
  delete[] dst;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CopyMatrix_CUDA(double* dst, double* src, size_t num_rows, size_t num_cols) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int num_elems = num_rows * num_cols;
  int temp = (int)(devProp.maxThreadsPerBlock / 2);
  int num_threads = (temp > num_elems) ? num_elems : temp;
  int num_blocks = (num_elems - 1) / num_threads + 1;
  CopyMatrixImpl<<<num_blocks, num_threads>>>(dst, src, num_elems);
}

__global__ void CopyMatrixImpl(double* dst, double* src, size_t num_elems) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (i < num_elems) {
    dst[i] = src[i];
  }
}

void SumMatrices_CUDA(double* dst, double* src_1, double* src_2, size_t num_rows, size_t num_cols) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int num_elems = num_rows * num_cols;
  int temp = (int)(devProp.maxThreadsPerBlock / 2);
  int num_threads = (temp > num_elems) ? num_elems : temp;
  int num_blocks = (num_elems - 1) / num_threads + 1;
  SumMatricesImpl<<<num_blocks, num_threads>>>(dst, src_1, src_2, num_elems);
}

__global__ void SumMatricesImpl(double* dst, double* src_1, double* src_2, size_t num_elems) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (i < num_elems) {
    dst[i] = src_1[i] + src_2[i];
  }
}

void DiffMatrices_CUDA(double* dst, double* src_1, double* src_2, size_t num_rows, size_t num_cols) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int num_elems = num_rows * num_cols;
  int temp = (int)(devProp.maxThreadsPerBlock / 2);
  int num_threads = (temp > num_elems) ? num_elems : temp;
  int num_blocks = (num_elems - 1) / num_threads + 1;
  DiffMatricesImpl<<<num_blocks, num_threads>>>(dst, src_1, src_2, num_elems);
}

__global__ void DiffMatricesImpl(double* dst, double* src_1, double* src_2, size_t num_elems) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (i < num_elems) {
    dst[i] = src_1[i] - src_2[i];
  }
}

void ProdMatrixByScalar_CUDA(double* dst, double* src, double alpha, size_t num_rows, size_t num_cols) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int num_elems = num_rows * num_cols;
  int temp = (int)(devProp.maxThreadsPerBlock / 2);
  int num_threads = (temp > num_elems) ? num_elems : temp;
  int num_blocks = (num_elems - 1) / num_threads + 1;
  ProdMatrixByScalarImpl<<<num_blocks, num_threads>>>(dst, src, alpha, num_elems);
}

__global__ void ProdMatrixByScalarImpl(double* dst, double* src, double alpha, size_t num_elems) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (i < num_elems) {
    dst[i] = src[i] * alpha;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void FivePointsLaplass_CUDA(double* dst, double* src, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int temp = (int)(sqrt(devProp.maxThreadsPerBlock) / 2);
  int num_threads_x = (temp > num_rows) ? num_rows : temp;
  int num_threads_y = (temp > num_cols) ? num_cols : temp;
  dim3 num_threads(num_threads_x, num_threads_y);
  dim3 num_blocks((num_rows - 1) / num_threads_x + 1, (num_cols - 1) / num_threads_y + 1);

  FivePointsLaplassImpl<<<num_blocks, num_threads>>>(dst, src, grid_x, grid_y, num_rows, num_cols);
}

__global__ void FivePointsLaplassImpl(double* dst, double* src, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if (i < (num_rows - 1) && j < (num_cols - 1) && i > 0 && j > 0) {
    double part_1 = (matrix_get(src, i, j, num_cols) - matrix_get(src, i - 1, j, num_cols)) / grid_r_step(grid_x, grid_y, i, num_cols) -
      (matrix_get(src, i + 1, j, num_cols) - matrix_get(src, i, j, num_cols)) / grid_r_step(grid_x, grid_y, i + 1, num_cols);
    double part_2 = (matrix_get(src, i, j, num_cols) - matrix_get(src, i, j - 1, num_cols)) / grid_c_step(grid_x, grid_y, j, num_cols) -
      (matrix_get(src, i, j + 1, num_cols) - matrix_get(src, i, j, num_cols)) / grid_c_step(grid_x, grid_y, j + 1, num_cols);
    matrix_set(dst, i, j, 2 * part_1 / (grid_r_step(grid_x, grid_y, i, num_cols) + grid_r_step(grid_x, grid_y, i + 1, num_cols)) +
                          2 * part_2 / (grid_c_step(grid_x, grid_y, j, num_cols) + grid_c_step(grid_x, grid_y, j + 1, num_cols)), num_cols);
  }

  if (((i == 0 || i == (num_rows - 1)) && (j < num_cols)) || ((j == 0 || j == (num_cols - 1)) && (i < num_rows))) {
    matrix_set(dst, i, j, 0.0, num_cols);
  }
}

double ProductByPointAndSum_CUDA(double* src_1, double* src_2, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  // http://mc.stanford.edu/cgi-bin/images/5/55/Darve_cme343_cuda_4.pdf
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  double* pre_retval;
  double* pre_retval_cuda;

  int num_elems = num_cols * num_rows;
  int temp = (int)(devProp.maxThreadsPerBlock / 2);
  int num_threads = (temp > num_elems) ? num_elems : temp;
  int num_blocks = (num_elems - 1) / num_threads + 1;

  pre_retval = new double[num_blocks + 1];
  CheckCuda(cudaMalloc((void**)(&pre_retval_cuda), (num_blocks + 1) * sizeof(double)));

  ProductByPointAndSumImpl<<<num_blocks, num_threads, num_threads * sizeof(double)>>>(pre_retval_cuda, src_1, src_2,
                                                                                      grid_x, grid_y, num_rows, num_cols);

  CheckCuda(cudaMemcpy(pre_retval, pre_retval_cuda, num_blocks * sizeof(double), cudaMemcpyDeviceToHost));

  double retval = 0.0;
  for (size_t i = 0; i < num_blocks; ++i) {
    retval += pre_retval[i];
  }
  return retval;
}

__global__ void ProductByPointAndSumImpl(double* pre_retval, double* src_1, double* src_2,
                                         double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  extern __shared__ double shared_mem[];
  int tid = threadIdx.x;
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  int num_elems = num_rows * num_cols;
  int col = (i % num_cols);
  int row = (i / num_cols);
  shared_mem[tid] = ((i > num_cols) && (i < (num_elems - num_cols)) && (col > 0) && (col < (num_cols - 1))) ?
                    src_1[i] * src_2[i] * 0.25 *
    (grid_r_step(grid_x, grid_y, row, num_cols) + grid_r_step(grid_x, grid_y, row + 1, num_cols)) *
    (grid_c_step(grid_x, grid_y, col, num_cols) + grid_c_step(grid_x, grid_y, col + 1, num_cols)) : 0.0;
  __syncthreads();

  for (int k = blockDim.x / 2; k > 0; k /= 2) {
    if (tid < k) {
      shared_mem[tid] += shared_mem[tid + k];
    }
    __syncthreads();
  }

  if (tid == 0) {
    pre_retval[blockIdx.x] = shared_mem[0];
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void InitValues_CUDA(double** values, double* grid_x, double* grid_y,
                     size_t num_rows, size_t num_cols, const ProcBounds& proc_bounds) {
  int num_elems = num_rows * num_cols;
  CheckCuda(cudaMalloc(values, num_elems * sizeof(double)));

  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int temp = (int)(sqrt(devProp.maxThreadsPerBlock) / 2);
  int num_threads_x = (temp > num_rows) ? num_rows : temp;
  int num_threads_y = (temp > num_cols) ? num_cols : temp;
  dim3 num_threads(num_threads_x, num_threads_y);
  dim3 num_blocks((num_rows - 1) / num_threads_x + 1, (num_cols - 1) / num_threads_y + 1);

  InitValuesImpl<<<num_blocks, num_threads>>>(*values, grid_x, grid_y, num_rows, num_cols,
                                              proc_bounds.is_up, proc_bounds.is_low, proc_bounds.is_left, proc_bounds.is_right);
}

__global__ void InitValuesImpl(double* dst, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols,
                               bool is_up, bool is_low, bool is_left, bool is_right) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if (i < num_rows && j < num_cols) {
    if ((is_up && i == 0) || (is_low && i == (num_rows - 1)) || (is_left && j == 0) || (is_right && j == num_cols - 1)) {
      matrix_set(dst, i, j, bound_func_k(matrix_get(grid_x, i, j, num_cols), matrix_get(grid_y, i, j, num_cols)), num_cols);
    } else {
      matrix_set(dst, i, j, RAND_CONST, num_cols);
    }
  }
}

void InitGrid_CUDA(double** grid_x, double** grid_y, size_t* num_rows, size_t* num_cols,
                   const GridData& grid_data, const ProcBounds& proc_bounds,
                   size_t start_row_idx, size_t end_row_idx, size_t start_col_idx, size_t end_col_idx) {
  start_row_idx -= !proc_bounds.is_up ? 1 : 0;
  end_row_idx += !proc_bounds.is_low ? 1 : 0;
  start_col_idx -= !proc_bounds.is_left ? 1 : 0;
  end_col_idx += !proc_bounds.is_right ? 1 : 0;

  *num_rows = end_row_idx - start_row_idx;
  *num_cols = end_col_idx - start_col_idx;

  int num_elems = (*num_rows) * (*num_cols);
  CheckCuda(cudaMalloc(grid_x, num_elems * sizeof(double)));
  CheckCuda(cudaMalloc(grid_y, num_elems * sizeof(double)));

  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  double* grid_x_cuda = new double[num_elems];
  double* grid_y_cuda = new double[num_elems];

  for (size_t i = 0; i < grid_data.r_num_points; ++i) {
    if (i >= start_row_idx && i < end_row_idx) {
      for (size_t j = 0; j < grid_data.c_num_points; ++j) {
        if (j >= start_col_idx && j < end_col_idx) {
          double cur_r_value = grid_data.r_upper_bound * step_func_k(static_cast<double>(i) / (grid_data.r_num_points - 1), grid_data.q) +
            grid_data.r_lower_bound * (1 - step_func_k(static_cast<double>(i) / (grid_data.r_num_points - 1), grid_data.q));
          double cur_c_value = grid_data.c_upper_bound * step_func_k(static_cast<double>(j) / (grid_data.c_num_points - 1), grid_data.q) +
            grid_data.c_lower_bound * (1 - step_func_k(static_cast<double>(j) / (grid_data.c_num_points - 1), grid_data.q));
          grid_x_cuda[(i - start_row_idx) * (*num_cols) + (j - start_col_idx)] = cur_r_value;
          grid_y_cuda[(i - start_row_idx) * (*num_cols) + (j - start_col_idx)] = cur_c_value;
        }
      }
    }
  }

  CheckCuda(cudaMemcpy(*grid_x, grid_x_cuda, num_elems * sizeof(double), cudaMemcpyHostToDevice));
  CheckCuda(cudaMemcpy(*grid_y, grid_y_cuda, num_elems * sizeof(double), cudaMemcpyHostToDevice));
  delete[] grid_x_cuda;
  delete[] grid_y_cuda;
}

void InitMatrix_CUDA(double** matrix, size_t num_rows, size_t num_cols) {
  int num_elems = num_rows * num_cols;
  CheckCuda(cudaMalloc(matrix, num_elems * sizeof(double)));

  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int temp = (int)(devProp.maxThreadsPerBlock / 2);
  int num_threads = (temp > num_elems) ? num_elems : temp;
  int num_blocks = (num_elems - 1) / num_threads + 1;

  InitMatrixImpl<<<num_blocks, num_threads>>>(*matrix, num_elems);
}

__global__ void InitMatrixImpl(double* dst, size_t num_elems) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (i < num_elems) {
    dst[i] = 0.0;
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ExchangeRowsCols_CUDA(double* matrix, size_t num_rows, size_t num_cols, const ProcBounds& proc_bounds,
                           std::pair<bool, bool> first_send, std::pair<int, int> left_right_proc, size_t proc_rank) {
  if (proc_bounds.is_up && proc_bounds.is_low && proc_bounds.is_left && proc_bounds.is_right) {
    return;
  }

  // exchange rows
  std::vector<double> data_send(num_cols, 0.0);
  std::vector<double> data_recv(num_cols, 0.0);

  if (!proc_bounds.is_low) {    
    CheckCuda(cudaMemcpy(&(data_send[0]), &(matrix[num_cols * (num_rows - 2)]), num_cols * sizeof(double), cudaMemcpyDeviceToHost));
    send_receive_vector(data_send,
                        &data_recv,
                        proc_rank + 1,
                        proc_rank + 1,
                        proc_rank,
                        proc_rank + 1);
    CheckCuda(cudaMemcpy(&(matrix[num_cols * (num_rows - 1)]), &(data_recv[0]), num_cols * sizeof(double), cudaMemcpyHostToDevice));
  }

  if (!proc_bounds.is_up) {
    CheckCuda(cudaMemcpy(&(data_send[0]), &(matrix[num_cols]), num_cols * sizeof(double), cudaMemcpyDeviceToHost));
    send_receive_vector(data_send,
                        &data_recv,
                        proc_rank - 1,
                        proc_rank - 1,
                        proc_rank,
                        proc_rank - 1);
    CheckCuda(cudaMemcpy(&(matrix[0]), &(data_recv[0]), num_cols * sizeof(double), cudaMemcpyHostToDevice)); 
  }

  // exchange cols
  data_send.assign(num_rows, 0.0);
  data_recv.assign(num_rows, 0.0);

  if (!proc_bounds.is_left) {
    cudaDeviceProp devProp;
    CheckCuda(cudaGetDeviceProperties(&devProp, 0));
    int num_threads = (int)(devProp.maxThreadsPerBlock / 2);
    int num_blocks = (num_rows * num_cols - 1) / num_threads + 1;

    double* temp;
    CheckCuda(cudaMalloc((void**)(&temp), num_rows * sizeof(double)));
    ExtractCol_CUDA<<<num_blocks, num_threads>>>(temp, matrix, 1, num_rows, num_cols);
    CheckCuda(cudaMemcpy(&(data_send[0]), temp, num_rows * sizeof(double), cudaMemcpyDeviceToHost));

    send_receive_vector(data_send,
                        &data_recv,
                        left_right_proc.first,
                        left_right_proc.first,
                        proc_rank,
                        left_right_proc.first);

    CheckCuda(cudaMemcpy(temp, &(data_recv[0]), num_rows * sizeof(double), cudaMemcpyHostToDevice));
    InsertCol_CUDA<<<num_blocks, num_threads>>>(matrix, temp, 0, num_rows, num_cols);
    CheckCuda(cudaFree(temp));
  }

  if (!proc_bounds.is_right) {
    cudaDeviceProp devProp;
    CheckCuda(cudaGetDeviceProperties(&devProp, 0));
    int num_threads = (int)(devProp.maxThreadsPerBlock / 2);
    int num_blocks = (num_rows * num_cols - 1) / num_threads + 1;

    double* temp;
    CheckCuda(cudaMalloc((void**)(&temp), num_rows * sizeof(double)));
    ExtractCol_CUDA<<<num_blocks, num_threads>>>(temp, matrix, num_cols - 2, num_rows, num_cols);
    CheckCuda(cudaMemcpy(&(data_send[0]), temp, num_rows * sizeof(double), cudaMemcpyDeviceToHost));

    send_receive_vector(data_send,
                        &data_recv,
                        left_right_proc.second,
                        left_right_proc.second,
                        proc_rank,
                        left_right_proc.second);

    CheckCuda(cudaMemcpy(temp, &(data_recv[0]), num_rows * sizeof(double), cudaMemcpyHostToDevice));
    InsertCol_CUDA<<<num_blocks, num_threads>>>(matrix, temp, num_cols - 1, num_rows, num_cols);
    CheckCuda(cudaFree(temp));
  }
}

__global__ void ExtractCol_CUDA(double* dst, double* src, size_t col, size_t num_rows, size_t num_cols) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  
  if (i < num_rows * num_cols) {
    if (i % num_cols == col) {
      dst[i / num_cols] = src[i];
    }
  }
}

__global__ void InsertCol_CUDA(double* dst, double* src, size_t col, size_t num_rows, size_t num_cols) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;

  if (i < num_rows * num_cols) {
    if (i % num_cols == col) {
      dst[i] = src[i / num_cols];
    }
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CountResiduals_CUDA(double* residuals, double* values_laplass, double* grid_x, double* grid_y,
                         size_t num_rows, size_t num_cols, const ProcBounds& proc_bounds) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int temp = (int)(sqrt(devProp.maxThreadsPerBlock) / 2);
  int num_threads_x = (temp > num_rows) ? num_rows : temp;
  int num_threads_y = (temp > num_cols) ? num_cols : temp;
  dim3 num_threads(num_threads_x, num_threads_y);
  dim3 num_blocks((num_rows - 1) / num_threads_x + 1, (num_cols - 1) / num_threads_y + 1);

  CountResidualsImpl<<<num_blocks, num_threads>>>(residuals, values_laplass, grid_x, grid_y, num_rows, num_cols,
                                                  proc_bounds.is_up, proc_bounds.is_low, proc_bounds.is_left, proc_bounds.is_right);
}

__global__ void CountResidualsImpl(double* residuals, double* values_laplass, double* grid_x, double* grid_y,
                                   size_t num_rows, size_t num_cols, bool is_up, bool is_low, bool is_left, bool is_right) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if (i < num_rows && j < num_cols) {
    if ((j == 0 && is_left) || (j == (num_cols - 1) && is_right) || (i == 0 && is_up) || (i == (num_rows - 1) && is_low)) {
      matrix_set(residuals, i, j, 0.0, num_cols);
    } else {
      double r = matrix_get(grid_x, i, j, num_cols);
      double c = matrix_get(grid_y, i, j, num_cols);
      matrix_set(residuals, i, j, matrix_get(values_laplass, i, j, num_cols) - main_func_k(r, c), num_cols);
    }
  }
}

double CountPreError_CUDA(double* values, double* temp_matrix, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  cudaDeviceProp devProp;
  CheckCuda(cudaGetDeviceProperties(&devProp, 0));

  int temp = (int)(sqrt(devProp.maxThreadsPerBlock) / 2);
  int num_threads_x = (temp > num_rows) ? num_rows : temp;
  int num_threads_y = (temp > num_cols) ? num_cols : temp;
  dim3 num_threads(num_threads_x, num_threads_y);
  dim3 num_blocks((num_rows - 1) / num_threads_x + 1, (num_cols - 1) / num_threads_y + 1);

  CountPreErrorImpl<<<num_blocks, num_threads>>>(temp_matrix, values, grid_x, grid_y, num_rows, num_cols);

  return ProductByPointAndSum_CUDA(temp_matrix, temp_matrix, grid_x, grid_y, num_rows, num_cols);
}

__global__ void CountPreErrorImpl(double* psi, double* values, double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = (blockIdx.y * blockDim.y) + threadIdx.y;

  if (i < num_rows && j < num_cols) {
    double r = matrix_get(grid_x, i, j, num_cols);
    double c = matrix_get(grid_y, i, j, num_cols);
    matrix_set(psi, i, j, matrix_get(values, i, j, num_cols) - true_func_k(r, c), num_cols);
  }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void CountNewValues_CUDA(double* values, double* old_values, double* gradients,
                         double* temp_matrix, double tau, size_t num_rows, size_t num_cols) {
  CopyMatrix_CUDA(old_values, values, num_rows, num_cols);
  ProdMatrixByScalar_CUDA(temp_matrix, gradients, tau, num_rows, num_cols);
  DiffMatrices_CUDA(values, old_values, temp_matrix, num_rows, num_cols);
}

double CountValuesDifference_CUDA(double* values, double* old_values, double* temp_matrix,
                                  double* grid_x, double* grid_y, size_t num_rows, size_t num_cols) {
  DiffMatrices_CUDA(temp_matrix, values, old_values, num_rows, num_cols);
  return ProductByPointAndSum_CUDA(temp_matrix, temp_matrix, grid_x, grid_y, num_rows, num_cols);
}

void CountGradients_CUDA(double* gradients, double* residuals, double* temp_matrix,
                         double alpha, size_t num_rows, size_t num_cols, bool first_iter) {
  if (first_iter) {
    CopyMatrix_CUDA(gradients, residuals, num_rows, num_cols);
  } else {
    ProdMatrixByScalar_CUDA(temp_matrix, gradients, alpha, num_rows, num_cols);
    DiffMatrices_CUDA(gradients, residuals, temp_matrix, num_rows, num_cols);
  }
}

} // namespace DTS
