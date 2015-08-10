#ifndef LAGRANGE_H_
#define LAGRANGE_H_

#include <vector>
#include <memory>

// Coefficients order: ^0 -> ^n
typedef std::vector<double> Coefficients;
typedef std::vector<double> FunctionArgs;
typedef std::vector<double> FunctionVals;
typedef std::vector<double> Errors;

const int MAX_POLYNOMIAL_DEGREE = 100;

class Lagrange {
 public:
  Lagrange() {
    coefficients.reset(new Coefficients);
  }

  std::shared_ptr<Coefficients> FindLagrangePolynomial(const FunctionArgs& arguments,
                                                       const FunctionVals& values);
  std::shared_ptr<FunctionVals> FindInterpolatedValues(const FunctionArgs& arguments);

  Errors& GetValuesCountError() { return values_count_error; }
  Errors& GetCoefficientsCountError() { return coefficients_count_error; }
 private:
  std::shared_ptr<Coefficients> ProdPolynomials(const Coefficients& pol_fst,
                                                const Coefficients& pol_snd,
                                                Errors* coef_err);
  inline double fabs_d(double arg);

  Errors coefficients_count_error;
  Errors values_count_error;
  std::shared_ptr<Coefficients> coefficients;
};
#endif  // LAGRANGE_H_
