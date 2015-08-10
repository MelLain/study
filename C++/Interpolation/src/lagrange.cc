// This file provides functions for finding Lagrange 
// interpolation polynomial of one variable and
// counting values of it in given set of arguments.
//
// Author: Murat Apishev, great-mel@yandex.ru

#include <iostream>

#include "lagrange.h"

std::shared_ptr<Coefficients> Lagrange::FindLagrangePolynomial(const FunctionArgs& arguments,
                                                              const FunctionVals& values) {
  // check the correctness of input data
  int arg_size = arguments.size();
  int polynomial_degree = values.size();
  if (arg_size == 0 ||
      polynomial_degree == 0 ||
      arg_size != polynomial_degree) {
    std::cout << "Input Error: size(arguments) != size(values)" << std::endl;
    return nullptr;
  }

  if (polynomial_degree < 2) {
    std::cout << "Input Error: size(values) is too small (<2)" << std::endl;
    return nullptr;
  }

  if (polynomial_degree > MAX_POLYNOMIAL_DEGREE) {
    polynomial_degree = MAX_POLYNOMIAL_DEGREE;
    std::cout << "Input Warning: size(values) > MAX, first MAX elements will be used" << std::endl;
  }

  // initialize the future result
  coefficients.reset(new Coefficients);
  coefficients_count_error.clear();
  for (int i = 0; i < polynomial_degree; ++i) {
    coefficients->push_back(0.0);
    coefficients_count_error.push_back(0.0);
  }
  
  for (int cur_node = 0; cur_node < polynomial_degree; ++cur_node) {
    Coefficients current_coefficients;
    current_coefficients.push_back(1.0);
    double denominator = 1.0;
    double denominator_count_error = 0.0;
    Errors cur_coef_count_err;
    for (int i = 0; i < polynomial_degree; ++i) {
      cur_coef_count_err.push_back(0.0);
    }

    for (int node = 0; node < polynomial_degree; ++node) {
      if (cur_node != node) {
        Coefficients temp_pol;
        temp_pol.push_back(-arguments[node]);
        temp_pol.push_back(1);
        // prod one member after another
        current_coefficients = *ProdPolynomials(current_coefficients,
                                                temp_pol,
                                                &cur_coef_count_err); 
        // prod one part of denominator after another
        denominator *= (arguments[cur_node] - arguments[node]);

        denominator_count_error += DBL_EPSILON;

        // d * (a - b)
        // err --- abs error
        // e --- relative error
        // eps --- machine zero
        // e = e + eps
      }
    }
    for (int deg = 0; deg < current_coefficients.size(); ++deg) {
      current_coefficients[deg] /= denominator;

      double temp_1 = fabs(current_coefficients[deg]) * (cur_coef_count_err[deg] +
          denominator_count_error);
      double temp_2 = fabs(current_coefficients[deg]) / 2 * DBL_EPSILON;
      cur_coef_count_err[deg] = (temp_1 + temp_2) / fabs_d(current_coefficients[deg]);

      // c / d
      // err --- abs error
      // e --- relative error
      // eps --- machine zero
      // err = |c / d| * (e(c) + e(d) + eps / 2)
      // e = [ {|c / d| * (e(c) + e(d))} + {|c / d| / 2 * eps} ] / |c / d|

      // add next summand to result polynomial
      double old_value = (*coefficients)[deg];
      (*coefficients)[deg] += current_coefficients[deg] * values[cur_node];

      temp_1 = fabs(current_coefficients[deg] * values[cur_node]) * cur_coef_count_err[deg];
      temp_2 = fabs(current_coefficients[deg] * values[cur_node]) / 2 * DBL_EPSILON;
      double temp_3 = coefficients_count_error[deg] * fabs(old_value);
      double temp_4 = DBL_EPSILON * fabs((*coefficients)[deg]) / 2;
      coefficients_count_error[deg] = (temp_1 + temp_2 + temp_3 + temp_4) /
          fabs_d((*coefficients)[deg]);

      // a * x + b
      // err --- abs error
      // e --- relative error
      // eps --- machine zero
      // err = |a * x| * (e(a) + e(x) + eps / 2) + err(b) + eps / 2 * |a * x + b|
      // e = [ {|a * x| * e(a)} + {|a * x| * eps / 2} + {e(b) * |b|} + {eps * |a * x + b| / 2} ] / |a * x + b|
    }
  }

  return coefficients;
}

std::shared_ptr<FunctionVals> Lagrange::FindInterpolatedValues(const FunctionArgs& arguments) {
  std::shared_ptr<FunctionVals> values(new FunctionVals);
  values_count_error.clear();
  int no_arguments = arguments.size();
  int polynomial_degree = coefficients->size();

  // check the correctness of input data
  if (no_arguments < 1) {
    std::cout << "Input Error: size(arguments) < 1" << std::endl;
    return nullptr;
  }
  if (polynomial_degree < 1) {
    std::cout << "Input Error: polynomial degree < 0" << std::endl;
    return nullptr;
  }

  // found the values of polynomial from given arguments, 
  // one after another, using Gorner's scheme
  for (auto& arg : arguments) {
    double value = (*coefficients)[polynomial_degree - 1];
    double val_error = 0.0;
    for (int deg = coefficients->size() - 2; deg >= 0; --deg) {
      double old_value = value;
      value = value * arg + (*coefficients)[deg];

      double temp_1 = fabs(value * arg) * val_error;
      double temp_2 = fabs(old_value * arg) / 2 * DBL_EPSILON;
      double temp_3 = coefficients_count_error[deg] * fabs(old_value);
      double temp_4 = DBL_EPSILON * fabs(value) / 2;
      val_error = (temp_1 + temp_2 + temp_3 + temp_4) / fabs_d(value);

      // a * x + b
      // err --- abs error
      // e --- relative error
      // eps --- machine zero
      // err = |a * x| * (e(a) + e(x) + eps / 2) + err(b) + eps / 2 * |a * x + b|
      // e = [ {|a * x| * e(a)} + {|a * x| * eps / 2} + {e(b) * |b|} + {|a * x + b| * eps / 2} ] / |a * x + b|
    }
    values->push_back(value);
    values_count_error.push_back(val_error);
  }

  return values;
}

std::shared_ptr<Coefficients> Lagrange::ProdPolynomials(const Coefficients& pol_fst,
                                                       const Coefficients& pol_snd,
                                                       Errors* coef_err) {
  // here the full errorless of counting of pol_snd is considered
  int max_degree = pol_fst.size() + pol_snd.size() - 2;
  std::shared_ptr<Coefficients> retval(new Coefficients);
  for (int i = 0; i <= max_degree; ++i)
    retval->push_back(0.0);

  for (int i = pol_fst.size() - 1; i >= 0; --i) {
    for (int j = pol_snd.size() - 1; j >= 0; --j) {
      double old_value = (*retval)[i + j];
      (*retval)[i + j] += pol_fst[i] * pol_snd[j];

      double temp_1 = fabs(pol_fst[i] * pol_snd[j]) * (*coef_err)[i];
      double temp_2 = fabs(pol_fst[i] * pol_snd[j]) / 2 * DBL_EPSILON;
      double temp_3 = (*coef_err)[i + j] * fabs(old_value);
      double temp_4 = DBL_EPSILON * fabs((*retval)[i + j]) / 2;
      (*coef_err)[i + j] = (temp_1 + temp_2 + temp_3 + temp_4) / fabs_d((*retval)[i + j]);

      // a * x + b
      // err --- abs error
      // e --- relative error
      // eps --- machine zero
      // err = |a * x| * (e(a) + e(x) + eps / 2) + err(b) + eps / 2 * |a * x + b|
      // e = [ {|a * x| * e(a)} + {|a * x| * eps / 2} + {e(b) * |b|} + {|a * x + b| * eps / 2} ] / |a * x + b|
    }
  }

  return retval;
}

double Lagrange::fabs_d(double arg) {
  if (fabs(arg) < DBL_EPSILON) return DBL_MAX;
  return fabs(arg);
}