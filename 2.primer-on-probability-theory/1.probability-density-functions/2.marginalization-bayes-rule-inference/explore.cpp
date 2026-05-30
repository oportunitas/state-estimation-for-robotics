/*
exploration: 
  a robot is somewhere on a line, its position is the state 'x'.

  there's a prior p(x), a belief about position before measuring. make it a gaussian centered at 0, with some standard deviation 's_location' (indicating that the position is first believed to be about at 0 with some uncertainty of 's_location')

  theres also a sensor model (likelihood) p(y | x). which gives the probability of a sensor reading 'y' given a true position 'x'. the sensor should return a position corrupted by gaussian noise. p(y | x) is a gaussian in x, centered at the measurement y, with the sensor's noise standard deviation of 's_sensor'. so if the sensor reads y = 2, the likelihood is a gaussian bump centered at x = 2

  task: say a measurement 'y' comes in (lets say 2). compute the posterior p(x | y) using Bayes' rule.
*/

#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>

constexpr double N_DIVISIONS = 1000000.0;

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

  auto gaussian = [] (double mean, double stdev, double x) { 
    double c = 1.0 / (stdev * std::sqrt(2 * M_PI));
    // std::cout << "c: " << c << std::endl;

    double x_part = (x - mean) / stdev;
    // std::cout << "x_part: " << x_part << std::endl;
    
    double exp_part = -0.5 * (x_part * x_part);
    // std::cout << "exp_part: " << exp_part << std::endl;
    
    return c * std::exp(exp_part); 
  };

  auto normal_gaussian = [gaussian] (double x) {
    return gaussian(0.0, 1.0, x);
  };

  auto location = [gaussian] (double x) {
    double s_location = 1.0;
    double mean_location = 0.0;
    return gaussian(mean_location, s_location, x);
  };

  auto sensor_likelihood = [gaussian] (double y, double x) {
    double s_sensor = 1.0;
    return gaussian(x, s_sensor, y);
  };

  auto true_loc_likelihood = [sensor_likelihood, location] (double x, double y) {
    auto sensor_derivative = [sensor_likelihood, location, y] (double x) {
      return sensor_likelihood(y, x) * location(x);
    };

    return (sensor_likelihood(y, x) * location(x)) / simpson_integrate(sensor_derivative, -10000.0, 10000.0);
  };

  // std::cout << sensor_likelihood(1.0, 2.0) << std::endl;

  // std::cout << simpson_integrate(normal_gaussian, -10000.0, 10000.0) << std::endl;

  std::cout << true_loc_likelihood(1.0, 2.0) << std::endl;

  return 0;
}