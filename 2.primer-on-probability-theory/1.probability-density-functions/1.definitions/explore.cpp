/*
exploration: write program that takes a 1d function on a finite interval [a, b], numerically integrates it (trapezoidal or simpson's rule, hand-rolled, no library), and reports whether its a valid pdf. 

test on three cases: 
  (1) uniform density: 1/(b-a)
  (2) unnormalized gaussian: exp(-x^2/2) without the 1/sqrt(2pi) factor, it should fail the integrates-to-1 test, compute what constant would fix it
  (3) function that dips negative somewhere, should fail the non-negativity test
*/

#include <iostream>
#include <cmath>
#include <string>

constexpr int N_DIVISIONS = 2000.0;

template <typename Func>
bool is_nonnegative(Func f, double a, double b) {
  for (double i = 0; i <= N_DIVISIONS; i++) {
    auto x = a + ((b - a) * (i / N_DIVISIONS));
    // std::cout << "x: " << x << std::endl;
    if (f(x) < 0) { return false; };
  }
  return true;
}

template <typename Func>
double simpson_on(Func f, double a, double b) {
  auto result = 0.0;

  for (double i = 0; i < N_DIVISIONS - 1; i += 2) {
    auto x0 = a + ((b - a) * (i / N_DIVISIONS));
    auto x1 = a + ((b - a) * ((i + 1) / N_DIVISIONS));
    auto x2 = a + ((b - a) * ((i + 2) / N_DIVISIONS));
    // std::cout << "x0: " << x0 << std::endl;
    // std::cout << "y0: " << f(x0) << std::endl;
    // std::cout << "x1: " << x1 << std::endl;
    // std::cout << "y1: " << f(x1) << std::endl;
    // std::cout << "x2: " << x2 << std::endl;
    // std::cout << "y2: " << f(x2) << std::endl;

    auto h = (b - a) / N_DIVISIONS;
    auto simpson_approx = (h / 3) * (f(x0) + 4 * f(x1) + f(x2));
    // std::cout << "h: " << h << std::endl;
    // std::cout << "integral: " << simpson_approx << std::endl;
    // std::cout << std::endl;

    result += simpson_approx;
  }

  return result;
}

template <typename Func>
std::string is_pdf(Func f, double a, double b) {
  // std::cout << simpson_on(f, a, b) << std::endl;
  // std::cout << std::abs(1.0 - simpson_on(f, a, b)) << std::endl;

  if (
    std::abs(1.0 - simpson_on(f, a, b)) < 1e-3 &&
    is_nonnegative(f, a, b)
  ) {
    return "this function is a PDF";
  } 
  return "no, this function is not a PDF";
}


int main() {
  std::ios_base::sync_with_stdio(false);

  auto uniform_a = [] (double x) { return 2.0 + x - x; };
  auto uniform_b = [] (double x) { return -2.0 + x - x; };
  auto unnormalized_gaussian = [] (double x) { return std::exp(-1 * x * x); };

  std::cout << "uniform a result: " << is_pdf(uniform_a, 0.0, 0.5) << std::endl;
  std::cout << "uniform a integral: " << simpson_on(uniform_a, 0.0, 0.5) << std::endl;

  std::cout << "uniform b result: " << is_pdf(uniform_b, 0.0, 0.5) << std::endl;
  std::cout << "uniform b integral: " << simpson_on(uniform_b, 0.0, 0.5) << std::endl;

  std::cout << "unnormalized gaussian result: " << is_pdf(unnormalized_gaussian, -15.0, 15.0) << std::endl;
  std::cout << "unnormalized gaussian integral: " << simpson_on(unnormalized_gaussian, -15.0, 15.0) << std::endl;
  std::cout << "unnormalized gaussian log: the constant that would make this pdf is: " << simpson_on(unnormalized_gaussian, -15.0, 15.0) << std::endl;

  return 0;
}