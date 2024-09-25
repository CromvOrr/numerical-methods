#define _USE_MATH_DEFINES

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

#define epsilon 1e-15
#define exactRoot asin(0.4)
#define minX 0.0
#define maxX M_PI_2

using namespace std;

double f(double x) { return sin(x) - 0.4; }

double df(double x) { return cos(x); }

double g(double x) { return (sin(x) - 0.4) * (sin(x) - 0.4); }

double dg(double x) { return 2 * (sin(x) - 0.4) * cos(x); }

double h(double x) { return ((sin(x) - 0.4)) / (2 * cos(x)); }

double dh(double x) { return 0.5 * ((((sin(x) - 0.4) * tan(x)) / cos(x)) + 1); }

void bisekcji() {
  ofstream results("bisekcji.txt");

  int i = 0;
  double a = minX, b = maxX; //* dla funkcji h zmienna b = 1.4
  double (*function)(double) = f;

  while (true) {
    double c = (a + b) / 2;
    results << i++ << "\t" << abs(c - exactRoot) << endl;

    if (function(a) * function(c) < 0)
      b = c;
    else if (function(c) * function(b) < 0)
      a = c;

    if (abs(c - exactRoot) < epsilon) {
      cout << "bisekcji: " << c << endl;
      break;
    }
  }

  results.close();
}

void falsi() {
  ofstream results("falsi.txt");

  int i = 0;
  double a = minX, b = maxX; //* dla funkcji h zmienna b = 1.4
  double (*function)(double) = f;

  while (true) {
    double c =
        ((-b * function(a)) + (a * function(b))) / (function(b) - function(a));
    results << i++ << "\t" << abs(c - exactRoot) << endl;

    if (function(a) * function(c) < 0)
      b = c;
    else if (function(c) * function(b) < 0)
      a = c;
    else
      break;

    if (abs(c - exactRoot) < epsilon) {
      cout << "falsi: " << c << endl;
      break;
    }
  }

  results.close();
}

void siecznych() {
  ofstream results("siecznych.txt");

  int i = 0;
  double xi = minX, xiMinus1 = maxX;
  double (*function)(double) = f;

  while (true) {
    double xiPlus1 = ((xi * function(xiMinus1)) - (xiMinus1 * function(xi))) /
                     (function(xiMinus1) - function(xi));
    results << i++ << "\t" << abs(xiPlus1 - exactRoot) << endl;

    if (abs(xiPlus1 - exactRoot) < epsilon) {
      cout << "siecznych: " << xiPlus1 << endl;
      break;
    }

    xiMinus1 = xi;
    xi = xiPlus1;
  }

  results.close();
}

void Newtona() {
  ofstream results("Newtona.txt");

  int i = 0;
  double xi = minX;
  double (*function)(double) = f;
  double (*derivative)(double) = df;

  while (true) {
    double xiPlus1 = xi - (function(xi) / derivative(xi));
    results << i++ << "\t" << abs(xiPlus1 - exactRoot) << endl;

    if (abs(xiPlus1 - exactRoot) < epsilon) {
      cout << "Newtona: " << xiPlus1 << endl;
      break;
    }

    xi = xiPlus1;
  }

  results.close();
}

int main() {
  cout << setprecision(150);
  cout << "znalezione pierwiastki:\n";

  bisekcji();
  falsi();
  siecznych();
  Newtona();

  return 0;
}