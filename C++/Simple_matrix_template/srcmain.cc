#include <exception>
#include <iostream>
#include <sstream>
#include <string>

#include "raw.h"
#include "matrix.h"

namespace mp = matrix_prod;

mp::Matrix<float> prod(const mp::Matrix<float>& mat_fst, const mp::Matrix<float>& mat_snd) {
  mp::Matrix<float> retval(mat_fst.no_raws(), mat_snd.no_cols(), 0.0f);

  for (int i = 0; i < mat_fst.no_raws(); ++i) {
    for (int j = 0; j < mat_snd.no_cols(); ++j) {
      for (int k = 0; k < mat_fst.no_cols(); ++k)
        retval[i][j] += mat_fst[i][k] * mat_snd[k][j];
    }
  }
  return retval;
}

void read_matrix(mp::Matrix<float>* matrix) {
  bool first_line = true;
  mp::Raw<float> current_raw;
  
  for (std::string line; std::getline(std::cin, line) && !line.empty();) {
    mp::Raw<std::string> buffer;
    std::stringstream ss(line);
    std::string elem;
    if (first_line) {
      ss >> elem;
      ss >> elem;
      first_line = false;
    }
    while (ss >> elem) {
      current_raw.add(std::stof(elem));
    }

    matrix->add_raw(current_raw);
    current_raw.clear();
  }
}

int main(int argc, char* argv[]) {
  try {
    mp::Matrix<float> A;
    mp::Matrix<float> B;
  
    read_matrix(&A);
    read_matrix(&B);

    if (A.no_cols() < 1) throw std::exception("Incorrect number of columns in first matrix\n");
    if (A.no_raws() < 1) throw std::exception("Incorrect number of raws in first matrix\n");
    if (B.no_cols() < 1) throw std::exception("Incorrect number of columns in second matrix\n");
    if (B.no_raws() < 1) throw std::exception("Incorrect number of raws in second matrix\n");
  
    if (A.no_cols() != B.no_raws())
      throw std::exception((std::string("Number of columns in first matrix should be equal to") +
                            std::string("number of raws in second matrix\n")).c_str());

    mp::Matrix<float> C = prod(A, B);
    for (int i = 0; i < C.no_raws(); ++i) {
      for (int j = 0; j < C.no_cols(); ++j) {
        std::cout << C[i][j] << ' ';
      }
      std::cout << std::endl;
    }

} catch (const std::exception& err) {
  std::cout << err.what() << std::endl;
}

  std::cout << "Press any key to continue...\n";
  char c = getchar();
  return 0;
}