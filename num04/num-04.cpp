#include <Eigen/Dense>
#include <chrono>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
using namespace chrono;
using namespace Eigen;

void setMarixAandVectorB(int sizeA, vector<vector<int>> &matrixA,
                         vector<int> &b) {
  matrixA.resize(sizeA, vector<int>(sizeA, 1));
  for (int i = 0; i < sizeA; i++) {
    matrixA[i][i] = 12;

    if (i < sizeA - 1) {
      matrixA[i][i + 1] = 8;
    }
  }

  b.resize(sizeA, 5);
}

void breakMatrixA(int sizeA, vector<vector<int>> &smallMatrixA, vector<int> &u,
                  vector<vector<int>> &vT) {
  smallMatrixA.resize(2, vector<int>(sizeA, 0));
  for (int i = 0; i < sizeA; i++) {
    if (i < sizeA - 1) {
      smallMatrixA[0][i + 1] = 7;
    }
    smallMatrixA[1][i] = 11;
  }

  u.resize(sizeA, 1);
  vT.resize(1, vector<int>(sizeA, 1));
}

void backSubstitution(int sizeA, vector<vector<int>> &smallMatrixA,
                      vector<double> &toBeCalculated,
                      vector<int> &knownVector) {
  toBeCalculated.resize(sizeA, 0.0);

  toBeCalculated[sizeA - 1] = static_cast<double>(knownVector[sizeA - 1]) /
                              static_cast<double>(smallMatrixA[1][sizeA - 1]);
  for (int i = sizeA - 2; i >= 0; i--) {
    toBeCalculated[i] = (static_cast<double>(knownVector[i]) -
                         (static_cast<double>(smallMatrixA[0][i + 1]) *
                          toBeCalculated[i + 1])) /
                        static_cast<double>(smallMatrixA[1][i]);
  }
}

void substituteIntoShermanMorrisonFormula(int sizeA, vector<double> &y,
                                          vector<double> &q,
                                          vector<vector<int>> &vT,
                                          vector<double> &z) {
  double vTq = 0.0, vTz = 0.0;
  y.resize(sizeA, 0.0);
  vector<double> vTqz(sizeA);

  for (int i = 0; i < sizeA; i++) {
    vTq += static_cast<double>(vT[0][i]) * q[i];
  }

  for (int i = 0; i < sizeA; i++) {
    vTqz[i] = vTq * z[i];
  }

  for (int i = 0; i < sizeA; i++) {
    vTz += static_cast<double>(vT[0][i]) * z[i];
  }
  vTz++;

  for (int i = 0; i < sizeA; i++) {
    vTqz[i] /= vTz;
  }

  for (int i = 0; i < sizeA; i++) {
    y[i] = q[i] - vTqz[i];
  }
}

void check(int size) {
  MatrixXd A = MatrixXd::Ones(size, size);
  A.diagonal() = VectorXd::Constant(size, 12);
  for (int i = 0; i < size - 1; i++) {
    A(i, i + 1) = 8;
  }

  VectorXd b = VectorXd::Constant(size, 5);
  VectorXd y = A.fullPivLu().solve(b);

  cout << y << endl;
}

int main() {
  int size = 80;
  vector<vector<int>> A, smallA, vT;
  vector<int> b, u;
  vector<double> q, z, y;

  setMarixAandVectorB(size, A, b);
  breakMatrixA(size, smallA, u, vT);
  backSubstitution(size, smallA, q, b);
  backSubstitution(size, smallA, z, u);
  substituteIntoShermanMorrisonFormula(size, y, q, vT, z);

  cout << "Wynik obliczony za pomoca biblioteki Eigen:\n";
  check(size);

  cout << "\nWynik obliczony za pomoca mojego programu:\n";
  for (double i : y) {
    cout << i << endl;
  }

  //* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

  ofstream dataa("dataa.dat");
  for (int i = 1; i <= 6001; i += 10) {
    //* Inicjalizacja zmiennych i ustawianie danych
    int sizee = i;
    vector<vector<int>> AA, smallAA, vTT;
    vector<int> bb, uu;
    vector<double> qq, zz, yy;
    setMarixAandVectorB(sizee, AA, bb);
    breakMatrixA(sizee, smallAA, uu, vTT);

    //* Właściwe obliczenia i pomiar czasu
    auto start = high_resolution_clock::now();
    backSubstitution(sizee, smallAA, qq, bb);
    backSubstitution(sizee, smallAA, zz, uu);
    substituteIntoShermanMorrisonFormula(sizee, yy, qq, vTT, zz);
    auto end = high_resolution_clock::now();
    auto timePassed = duration_cast<microseconds>(end - start).count() / 1000.0;

    dataa << sizee << "\t" << timePassed << endl;
  }

  dataa.close();
  return 0;
}