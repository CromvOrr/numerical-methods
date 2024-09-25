#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace Eigen;

double norm(int N, double oldX[], double newX[]) {
  double temp[N];
  for (int i = 0; i < N; i++) {
    temp[i] = oldX[i] - newX[i];
  }

  double sum = 0;
  for (int i = 0; i < N; i++) {
    sum += temp[i] * temp[i];
  }
  return sqrt(sum);
}

double normm(int N, VectorXd &oldX, double newX[]) {
  double sum = 0;
  for (int i = 0; i < N; i++) {
    sum += (oldX(i) - newX[i]) * (oldX(i) - newX[i]);
  }
  return sqrt(sum);
}

VectorXd exactSolution(int N) {
  MatrixXd A = MatrixXd::Zero(N, N);
  VectorXd b(N);
  for (int i = 0; i < N; i++) {
    if (i + 1 < N) {
      A(i, i + 1) = 1;
    }
    if (i + 2 < N) {
      A(i, i + 2) = 0.15;
    }
    A(i, i) = 3;
    if (i - 2 >= 0) {
      A(i, i - 2) = 0.15;
    }
    if (i - 1 >= 0) {
      A(i, i - 1) = 1;
    }

    b(i) = i + 1;
  }

  VectorXd x = A.fullPivLu().solve(b);
  return x;
}

int main() {
  int k = 0, N = 124;
  bool foo = false;
  VectorXd x = exactSolution(N);
  double oldX[N], newX[N];
  for (int i = 0; i < N; i++) {
    oldX[i] = static_cast<double>(10000);
  }

  ofstream data("gaussSeidel.txt");
  if (foo) {
    while (true) {
      for (int i = 0; i < N; i++) {
        newX[i] = static_cast<double>(i + 1);
      }

      for (int i = 0; i < N; i++) {
        if (i + 1 < N) {
          newX[i] -= oldX[i + 1];
        }
        if (i + 2 < N) {
          newX[i] -= 0.15 * oldX[i + 2];
        }
        if (i - 2 >= 0) {
          newX[i] -= 0.15 * oldX[i - 2];
        }
        if (i - 1 >= 0) {
          newX[i] -= oldX[i - 1];
        }

        newX[i] /= 3;
      }

      if (norm(N, oldX, newX) < 1e-13) {
        data << k << "\t" << normm(N, x, newX) << endl;
        break;
      }

      for (int i = 0; i < N; i++) {
        oldX[i] = newX[i];
      }

      data << k << "\t" << normm(N, x, newX) << endl;
      k++;
    }
  } else {
    while (true) {
      for (int i = 0; i < N; i++) {
        newX[i] = static_cast<double>(i + 1);
      }

      for (int i = 0; i < N; i++) {
        if (i + 1 < N) {
          newX[i] -= oldX[i + 1];
        }
        if (i + 2 < N) {
          newX[i] -= 0.15 * oldX[i + 2];
        }
        if (i - 2 >= 0) {
          newX[i] -= 0.15 * newX[i - 2];
        }
        if (i - 1 >= 0) {
          newX[i] -= newX[i - 1];
        }

        newX[i] /= 3;
      }

      if (norm(N, oldX, newX) < 1e-13) {
        data << k << "\t" << normm(N, x, newX) << endl;
        break;
      }

      for (int i = 0; i < N; i++) {
        oldX[i] = newX[i];
      }

      data << k << "\t" << normm(N, x, newX) << endl;
      k++;
    }
  }

  data.close();
  return 0;
}