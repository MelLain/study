#include <cfloat>
#include <cmath>
#include <ctime>

#include <iostream>
#include <fstream>
#include <sstream> 
#include <string>

#include "lagrange.h"

const int INCORRECT_INPUT = -1;
const int SUCCESS = 0;
const std::string error_str = std::string("Usage: interpolation /fp:strict [<test_case>]") +
                              std::string(",\n <test_case>:\n 0 --- basic test\n 1 --- ") +
                              std::string("test #1\n 2 --- test #2\n 3 --- test #3");
void test_0();
void test_1();
void test_2();
void test_3();

inline int factorial(int n) {
  if (n == 1 || n == 0) return 1;
  return n * factorial(n - 1);
}

int main(int argc, char* argv[])
{
  int test_case = 0;
  if (argc > 2) {
    std::cout << error_str << std::endl;
    return INCORRECT_INPUT;
  }

  if (argc == 2) {
    test_case = strtol(argv[1], nullptr, 10);
    if (test_case != 0 &&test_case != 1 && test_case != 2 && test_case != 3) {
      std::cout << error_str << std::endl;
      return INCORRECT_INPUT;
    }
  }

  switch(test_case) {
   case 0: test_0(); break;
   case 1: test_1(); break;
   case 2: test_2(); break;
   case 3: test_3(); break;
  }
  return SUCCESS;
}

void test_0() {
  std::ofstream output_data_file("output_test_0.txt");

  for (int n = 2; n <= 10; n += 1) {
    FunctionArgs args;
    FunctionVals vals;
    double left = 0.0;
    double right = 0.9;

    double delta = (right - left) / n;

    double x = 0.0;
    for (int i = 0; i < n; ++i, x += delta) {
      args.push_back(x);
      vals.push_back(exp(x));
    }
    Lagrange lagrange_polynomial;
    auto coefficients = lagrange_polynomial.FindLagrangePolynomial(args, vals);

    FunctionArgs nodes = args;
    args.clear();
    vals.clear();

    int no_steps = 1000;
    double delta_step = (right - left) / no_steps;
    double current_step = left;
    for (int i = 0; i < no_steps; ++i, current_step += delta_step) {
      args.push_back(current_step);
      vals.push_back(exp(current_step));
    }
    auto found_vals = *lagrange_polynomial.FindInterpolatedValues(args);
    
    Errors final_count_error = lagrange_polynomial.GetValuesCountError();
    double mean_analytic_error = 0.0;
    double mean_count_error = 0.0;
    double mean_error = 0.0;
    double max_error = 0.0;
    double mean_limit_error = 0.0;

    for (int i = 0; i < no_steps; ++i) {
      double cur_real_error = fabs(vals[i] - found_vals[i]);
      double cur_analytic_error = exp(right);
      for (int j = 0; j < n; ++j) {
        cur_analytic_error *= args[i] - nodes[j];
      }
      cur_analytic_error /= factorial(n);
      mean_error += cur_real_error;
      mean_analytic_error += cur_analytic_error;
      mean_count_error += final_count_error[i];
      if (cur_real_error > max_error) {
        max_error = cur_real_error;
      }
      mean_limit_error += final_count_error[i] + cur_analytic_error;
    }

    mean_error /= no_steps;
    mean_limit_error /= no_steps;
    mean_analytic_error /= no_steps;
    mean_count_error /= no_steps;

    output_data_file << mean_error << ' ' << mean_limit_error << ' ' <<
        mean_analytic_error << ' ' << mean_count_error << std::endl;

    std::cout << "n: " << n << std::endl;
    std::cout << "Mean error: " << mean_error << std::endl;
    std::cout << "Mean analitical limit error: " << mean_analytic_error << std::endl;
    std::cout << "Mean count limit error: " << mean_count_error << std::endl;
    std::cout << "Mean limit error: " << mean_limit_error << std::endl;
    std::cout << "Maximal error: " << max_error << std::endl;
  }
  output_data_file.close();
  std::cout << "Press any key to continue..." << std::endl;
  char c = getchar();
}

void test_1() {
  std::ofstream output_data_file("output_test_1.txt");

  for (int n = 2; n <= 10; n += 1) {
    FunctionArgs args;
    FunctionVals vals;
    double left = -0.8;
    double right = 0.9;

    double delta = (right - left) / n;

    double x = 0.0;
    for (int i = 0; i < n; ++i, x += delta) {
      args.push_back(x);
      vals.push_back(5 * cos(x));
    }
    Lagrange lagrange_polynomial;
    auto coefficients = lagrange_polynomial.FindLagrangePolynomial(args, vals);

    FunctionArgs nodes = args;
    args.clear();
    vals.clear();

    int no_steps = 1000;
    double delta_step = (right - left) / no_steps;
    double current_step = left;
    for (int i = 0; i < no_steps; ++i, current_step += delta_step) {
      args.push_back(current_step);
      vals.push_back(5 * cos(current_step));
    }
    auto found_vals = *lagrange_polynomial.FindInterpolatedValues(args);
    
    Errors final_count_error = lagrange_polynomial.GetValuesCountError();
    double mean_analytic_error = 0.0;
    double mean_count_error = 0.0;
    double mean_error = 0.0;
    double max_error = 0.0;
    double mean_limit_error = 0.0;

    for (int i = 0; i < no_steps; ++i) {
      double cur_real_error = fabs(vals[i] - found_vals[i]);

      double cur_analytic_error = 5.0;
      if (n % 2 == 1) 
        cur_analytic_error *= sin(right);

      for (int j = 0; j < n; ++j) {
        cur_analytic_error *= args[i] - nodes[j];
      }
      cur_analytic_error /= factorial(n);
      mean_error += cur_real_error;
      mean_analytic_error += cur_analytic_error;
      mean_count_error += final_count_error[i];
      if (cur_real_error > max_error) {
        max_error = cur_real_error;
      }
      mean_limit_error += final_count_error[i] + cur_analytic_error;
    }

    mean_error /= no_steps;
    mean_limit_error /= no_steps;
    mean_analytic_error /= no_steps;
    mean_count_error /= no_steps;

    output_data_file << mean_error << ' ' << mean_limit_error << ' ' <<
        mean_analytic_error << ' ' << mean_count_error << std::endl;

    std::cout << "n: " << n << std::endl;
    std::cout << "Mean error: " << mean_error << std::endl;
    std::cout << "Mean analitical limit error: " << mean_analytic_error << std::endl;
    std::cout << "Mean count limit error: " << mean_count_error << std::endl;
    std::cout << "Mean limit error: " << mean_limit_error << std::endl;
    std::cout << "Maximal error: " << max_error << std::endl;
  }
  output_data_file.close();
  std::cout << "Press any key to continue..." << std::endl;
  char c = getchar();
}

void test_2() {
  std::ofstream output_data_file("output_test_2.txt");

  for (int n = 2; n <= 10; n += 1) {
    FunctionArgs args;
    FunctionVals vals;
    double left = 0.0;
    double right = 1.0;

    double delta = (right - left) / n;

    double x = 0.0;
    for (int i = 0; i < n; ++i, x += delta) {
      args.push_back(x);
      vals.push_back(4 * x * x + exp(2 * x));
    }
    Lagrange lagrange_polynomial;
    auto coefficients = lagrange_polynomial.FindLagrangePolynomial(args, vals);

    FunctionArgs nodes = args;
    args.clear();
    vals.clear();

    int no_steps = 1000;
    double delta_step = (right - left) / no_steps;
    double current_step = left;
    for (int i = 0; i < no_steps; ++i, current_step += delta_step) {
      args.push_back(current_step);
      vals.push_back(4 * current_step * current_step + exp(2 * current_step));
    }
    auto found_vals = *lagrange_polynomial.FindInterpolatedValues(args);
    
    Errors final_count_error = lagrange_polynomial.GetValuesCountError();
    double mean_analytic_error = 0.0;
    double mean_count_error = 0.0;
    double mean_error = 0.0;
    double max_error = 0.0;
    double mean_limit_error = 0.0;

    for (int i = 0; i < no_steps; ++i) {
      double cur_real_error = fabs(vals[i] - found_vals[i]) / fabs(vals[i]);
      
      double cur_analytic_error = 8.0 + 2.0 * exp(2 * right);
      if (n > 2) {
        double k = 2.0;
        for (int j = 0; j < n - 2; ++j) k *= 2.0;
        cur_analytic_error = k * exp(2 * right);
      }

      for (int j = 0; j < n; ++j) {
        cur_analytic_error *= args[i] - nodes[j];
      }
      cur_analytic_error /= factorial(n);
      mean_error += cur_real_error;
      mean_analytic_error += cur_analytic_error;
      mean_count_error += final_count_error[i];
      if (cur_real_error > max_error) {
        max_error = cur_real_error;
      }
      mean_limit_error += final_count_error[i] + cur_analytic_error;
    }

    mean_error /= no_steps;
    mean_limit_error /= no_steps;
    mean_analytic_error /= no_steps;
    mean_count_error /= no_steps;

    output_data_file << mean_error << ' ' << mean_limit_error << ' ' <<
        mean_analytic_error << ' ' << mean_count_error << std::endl;

    std::cout << "n: " << n << std::endl;
    std::cout << "Mean error: " << mean_error << std::endl;
    std::cout << "Mean analitical limit error: " << mean_analytic_error << std::endl;
    std::cout << "Mean count limit error: " << mean_count_error << std::endl;
    std::cout << "Mean limit error: " << mean_limit_error << std::endl;
    std::cout << "Maximal error: " << max_error << std::endl;
  }
  output_data_file.close();
  std::cout << "Press any key to continue..." << std::endl;
  char c = getchar();
}

void test_3() {
  int n = 4;
  int no_experiments = 20000;
  srand(time(NULL));
  std::ofstream output_data_file("output_test_3.txt");
  for (int experiment = 0; experiment < no_experiments; ++experiment) {
    FunctionArgs args;
    FunctionVals vals;

    double a_0 = rand() / (RAND_MAX + 1.0);
    double a_1 = rand() / (RAND_MAX + 1.0);
    double b = rand() % 10;
    double c = rand() % 5;
    double d = rand() / (RAND_MAX + 1.0) * 3;

    double left = 0.0;
    double right = 0.9;
    double delta = (right - left) / n;
    
    double x = 0.0;
    for (int i = 0; i < n; ++i, x += delta) {
      args.push_back(x);
      vals.push_back(4 * x * x + exp(2 * x));
    }
    Lagrange lagrange_polynomial;
    auto coefficients = lagrange_polynomial.FindLagrangePolynomial(args, vals);

    FunctionArgs nodes = args;
    args.clear();
    vals.clear();

    int no_steps = 1000;
    double delta_step = (right - left) / no_steps;
    double current_step = left;
    for (int i = 0; i < no_steps; ++i, current_step += delta_step) {
      args.push_back(current_step);
      vals.push_back(4 * current_step * current_step + exp(2 * current_step));
    }
    auto found_vals = *lagrange_polynomial.FindInterpolatedValues(args);
    
    Errors final_count_error = lagrange_polynomial.GetValuesCountError();
    double norm_error = 0.0;

    for (int i = 0; i < no_steps; ++i) {
      double cur_real_error = fabs(vals[i] - found_vals[i]) / fabs(vals[i]);
      double cur_analytic_error = fabs(std::pow(b, 4) * sin(b * x) +
          c * std::pow(d, 4) * exp(d * x));

      for (int j = 0; j < n; ++j) {
        cur_analytic_error *= args[i] - nodes[j];
      }
      cur_analytic_error /= factorial(n);
      int sign = rand() % 2;
      if (sign == 0) sign = -1;
        norm_error += sign * cur_real_error / (final_count_error[i] + cur_analytic_error);
    }
    norm_error /= no_steps;

    output_data_file << norm_error << std::endl;
  }
  output_data_file.close();
  std::cout << "Press any key to continue..." << std::endl;
  char c = getchar();
}