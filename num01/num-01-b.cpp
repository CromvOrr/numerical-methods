#include <cmath>
#include <fstream>
#include <iomanip>
#include <random>

using namespace std;

float float_derivative(float x) { return 2 * x * cos(x * x); }

double double_derivative(double x) { return 2 * x * cos(x * x); }

float float_approx_derivative(float x, float h) {
  return (sin((x + h) * (x + h)) - sin((x - h) * (x - h))) / (2 * h);
}

double double_approx_derivative(double x, double h) {
  return (sin((x + h) * (x + h)) - sin((x - h) * (x - h))) / (2 * h);
}

int main() {
  int float_samples = 120, double_samples = 60;
  float float_x = 0.2, float_generated, float_min = 1e-8;
  double double_x = 0.2, double_generated, double_min = 1e-17;

  ofstream data_float("data_float_b.dat");
  ofstream data_double("data_double_b.dat");
  random_device rd;
  mt19937 gen(rd());

  for (int i = 0; i < 8; i++) {
    uniform_real_distribution<float> float_dist(
        static_cast<float>(float_min), static_cast<float>(float_min * 10));

    for (int j = 0; j < float_samples; j++) {
      float_generated = float_dist(gen);
      data_float << float_generated << "\t"
                 << abs(float_approx_derivative(float_x, float_generated) -
                        float_derivative(float_x))
                 << endl;
    }

    float_min *= 10;
  }

  for (int i = 0; i < 17; i++) {
    uniform_real_distribution<double> double_dist(
        static_cast<double>(double_min), static_cast<double>(double_min * 10));

    for (int j = 0; j < double_samples; j++) {
      double_generated = double_dist(gen);
      data_double << double_generated << "\t"
                  << abs(double_approx_derivative(double_x, double_generated) -
                         double_derivative(double_x))
                  << endl;
    }

    double_min *= 10;
  }

  data_float.close();
  data_double.close();

  return 0;
}