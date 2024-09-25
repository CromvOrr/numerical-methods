#define _USE_MATH_DEFINES

#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

#define n 30
#define m 1500

double y(double x) // y(x)
{
  return 1.0 / (1.0 + (50.0 * (x * x)));
}

double h(double x) // h(x)
{
  return sin(cos(tan(3 * x)));
}

double axi(double i) // (a) xi
{
  return -1.0 + (2.0 * (i / (n + 1.0)));
}

double bxi(double i) // (b) xi
{
  double temp = (((2.0 * i) + 1.0) / (2.0 * (n + 1.0))) * M_PI;
  return cos(temp);
}

void generatePoints(vector<double> &points) {
  double step = 2.0 / (m - 1);
  for (double i = 0; i < m; i++) {
    points.push_back(-1.0 + (i * step));
  }
}

void Wnx() // Wn(x)
{
  ofstream results("results.txt");

  double xes[n];
  vector<double> points;
  generatePoints(points);
  for (int i = 0; i < n; i++) {
    xes[i] = bxi(i);
  }

  for (int i = 0; i < m; i++) {
    double sum = 0.0;
    for (int j = 0; j < n; j++) {
      double phi = 1.0;
      for (int k = 0; k < n; k++) {
        if (k != j) {
          phi *= (points[i] - xes[k]) / (xes[j] - xes[k]);
        }
      }
      sum += y(xes[j]) * phi;
    }
    results << points[i] << "\t" << sum << "\t" << y(points[i]) << endl;
  }

  results.close();
}

int main() {
  Wnx();

  return 0;
}