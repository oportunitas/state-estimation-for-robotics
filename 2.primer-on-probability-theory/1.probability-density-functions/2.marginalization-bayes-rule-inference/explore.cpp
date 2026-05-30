/*
exploration: 
  a robot is somewhere on a line, its position is the state 'x'.
  there's a prior p(x), a belief about position before measuring. make it a gaussian centered at 0,
  with some standard deviation 's_prior' (indicating that the position is first believed to be
  about at 0 with some uncertainty of 's_prior')
*/

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>

constexpr double N_DIVISIONS = 100000.0;

template <typename Func>
double simpson_integrate(Func f, double a, double b) {
  auto result = 0.0;

  for (double i = 0; i < N_DIVISIONS - 1; i += 2) {
    auto x0 = a + ((b - a) * (i / N_DIVISIONS));
    auto x1 = a + ((b - a) * ((i + 1) / N_DIVISIONS));
    auto x2 = a + ((b - a) * ((i + 2) / N_DIVISIONS));
    auto h = (b - a) / N_DIVISIONS;
    auto simpson_approx = (h / 3) * (f(x0) + 4 * f(x1) + f(x2));
    result += simpson_approx;
  }

  return result;
}

int main() {
  std::ios_base::sync_with_stdio(false);

  double gaussian_mean = 1.0;
  double gaussian_stdev = 10.0;

  auto gaussian = [gaussian_mean, gaussian_stdev] (double x) { 
    double c = 1 / (gaussian_stdev * std::sqrt(2 * M_PI));
    // std::cout << "c: " << c << std::endl;

    double x_part = (x - gaussian_mean) / gaussian_stdev;
    // std::cout << "x_part: " << x_part << std::endl;
    
    double exp_part = -0.5 * (x_part * x_part);
    // std::cout << "exp_part: " << exp_part << std::endl;
    
    return c * std::exp(exp_part); 
  };

  std::cout << simpson_integrate(gaussian, -10000.0, 10000.0) << std::endl;

  return 0;
}