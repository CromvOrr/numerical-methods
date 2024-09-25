#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>

#define n 4
#define epsilon 1e-16

using namespace std;
using namespace Eigen;

void setValues(double M[n][n]) {
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      M[i][j] = 0;
    }
  }

  M[0][0] = 8;
  M[0][1] = 1;
  M[1][0] = 1;
  M[1][1] = 7;
  M[1][2] = 2;
  M[2][1] = 2;
  M[2][2] = 6;
  M[2][3] = 3;
  M[3][2] = 3;
  M[3][3] = 5;
}

void multiply(double M[n][n], double V[n], double result[n]) {
  for (int i = 0; i < n; i++) {
    result[i] = 0;
    for (int j = 0; j < n; j++) {
      result[i] += M[i][j] * V[j];
    }
  }
}

double norm(double V[n]) {
  double sum = 0;
  for (int i = 0; i < n; i++) {
    sum += V[i] * V[i];
  }
  return sqrt(sum);
}

void powerMethod(double M[n][n], double p = 0.0) {
  ofstream data("powerMethodP0.txt");

  for (int i = 0; i < n; i++) {
    M[i][i] -= p;
  }

  double vector[n];

  for (int i = 0; i < n; i++) {
    vector[i] = 1;
  }
  double lambda = 1.0;

  int k = 0;
  while (true) {
    double multipled[n];
    multiply(M, vector, multipled);
    double new_lambda = norm(multipled);

    data << k << "\t" << abs(new_lambda - lambda) << endl;

    if (abs(new_lambda - lambda) < epsilon) {
      break;
    }
    lambda = new_lambda;

    for (int i = 0; i < n; i++) {
      vector[i] = multipled[i] / lambda;
    }

    k++;
  }

  data.close();

  cout << "Najwieksza wartosc wlasna: " << lambda << endl;
  cout << "Odpowiadajacy wektor wlasny:\n";
  for (int i = 0; i < n; i++) {
    cout << vector[i] << "\n";
  }
  cout << "Iteracje: " << k << endl;
  cout << "\n";
}

void setValuesEigen(MatrixXd &M) {
  M = MatrixXd::Zero(n, n);
  M(0, 0) = 8;
  M(0, 1) = 1;
  M(1, 0) = 1;
  M(1, 1) = 7;
  M(1, 2) = 2;
  M(2, 1) = 2;
  M(2, 2) = 6;
  M(2, 3) = 3;
  M(3, 2) = 3;
  M(3, 3) = 5;
}

void calcuteByEigen(MatrixXd &M) {
  SelfAdjointEigenSolver<MatrixXd> solver(M);
  auto eigenValues = solver.eigenvalues();
  auto eigenVectors = solver.eigenvectors();
  int maxIndex;
  double maxEigenValue = eigenValues.cwiseAbs().maxCoeff(&maxIndex);

  cout << "[EIGEN] Najwieksza wartosc wlasna: " << maxEigenValue << endl;
  cout << "[EIGEN] Odpowiadajacy wektor wlasny:\n"
       << eigenVectors.col(maxIndex) << endl
       << endl;
}

double underDiagonalSum(MatrixXd &M) {
  double sum = 0.0;
  for (int i = 1; i < n; i++) {
    sum += abs(M(i, i - 1));
  }
  return sum;
}

void qrMethod(MatrixXd &M) {
  ofstream dataa("qrMethod.txt");
  ofstream sum("underDiagonalSum.txt");

  const double smallEpsilon = 1e-14;
  SelfAdjointEigenSolver<MatrixXd> solver(M);
  VectorXd eigenValues = solver.eigenvalues();

  double underDiagonalSummed = 0.0;
  int i = 0, k0 = 0, k1 = 0, k2 = 0, k3 = 0;
  do {
    HouseholderQR<MatrixXd> qr(M);
    MatrixXd Q = qr.householderQ();
    MatrixXd R = qr.matrixQR().triangularView<Upper>();

    MatrixXd newM = R * Q;
    M = newM;
    VectorXd newEigenValues = M.diagonal();

    dataa << i << "\t";
    if (abs(eigenValues(3) - newEigenValues(0)) > smallEpsilon) {
      dataa << "\t" << abs(eigenValues(3) - newEigenValues(0));
    } else {
      dataa << "\t";
      k0 = 1;
    }
    if (abs(eigenValues(2) - newEigenValues(1)) > smallEpsilon) {
      dataa << "\t" << abs(eigenValues(2) - newEigenValues(1));
    } else {
      dataa << "\t";
      k1 = 1;
    }
    if (abs(eigenValues(1) - newEigenValues(2)) > smallEpsilon) {
      dataa << "\t" << abs(eigenValues(1) - newEigenValues(2));
    } else {
      dataa << "\t";
      k2 = 1;
    }
    if (abs(eigenValues(0) - newEigenValues(3)) > smallEpsilon) {
      dataa << "\t" << abs(eigenValues(0) - newEigenValues(3));
    } else {
      dataa << "\t";
      k3 = 1;
    }
    dataa << endl;

    double sumK = k0 + k1 + k2 + k3;
    underDiagonalSummed = underDiagonalSum(M);
    sum << i << "\t" << underDiagonalSummed << endl;
    if (underDiagonalSummed < 1e-15 || sumK == 4) {
      break;
    }

    i++;
  } while (true);

  cout << "Wartosci wlasne odczytane z diagonali:\n"
       << M.diagonal() << "\nIteracje: " << i << "\n\n";

  dataa.close();
  sum.close();
}

void eigenValuesEigen(MatrixXd &M) {
  SelfAdjointEigenSolver<MatrixXd> solver(M);
  auto eigenValues = solver.eigenvalues();
  cout << "[EIGEN] Wartosci wlasne odczytane z diagonali:\n"
       << eigenValues << endl;
}

int main() {
  double M[n][n];
  setValues(M);
  powerMethod(M);

  MatrixXd M2(n, n);
  setValuesEigen(M2);
  calcuteByEigen(M2);

  MatrixXd M3(n, n);
  setValuesEigen(M3);
  qrMethod(M3);

  MatrixXd M4(n, n);
  setValuesEigen(M4);
  eigenValuesEigen(M4);

  return 0;
}