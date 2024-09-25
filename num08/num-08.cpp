#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using namespace std;
using namespace Eigen;

#define linesInData 100
#define linesInGData 100

void insertValues(double xValues[], double yValues[]) {
  ifstream data("data/data2023.dat");

  int i = 0;
  string line;
  while (getline(data, line)) {
    int position = line.find(",");
    string xString = line.substr(0, position);
    string yString = line.substr(position + 1, line.length());
    double x = stod(xString);
    double y = stod(yString);
    xValues[i] = x;
    yValues[i] = y;
    i++;
  }

  data.close();
}

double F(double x, double abcd[]) {
  return (abcd[0] * pow(x, 2)) + (abcd[1] * sin(x)) + (abcd[2] * cos(5 * x)) +
         (abcd[3] * exp(-x));
}

double G(double x) {
  return (2.5 * sin(x)) + (1.25 * cos(x)) + (0.625 * x) + (0.3125 * pow(x, 2));
}

double G(double x, double abcd[]) {
  return (abcd[0] * sin(x)) + (abcd[1] * cos(x)) + (abcd[2] * x) +
         (abcd[3] * pow(x, 2));
}

MatrixXd insertData(double xValues[]) {
  MatrixXd M(linesInData, 4);
  for (int i = 0; i < linesInData; i++) {
    M(i, 0) = pow(xValues[i], 2);
    M(i, 1) = sin(xValues[i]);
    M(i, 2) = cos(5 * xValues[i]);
    M(i, 3) = exp(-xValues[i]);
  }
  return M;
}

MatrixXd insertGdata(double xValues[]) {
  MatrixXd M(linesInGData, 4);
  for (int i = 0; i < linesInGData; i++) {
    M(i, 0) = sin(xValues[i]);
    M(i, 1) = cos(xValues[i]);
    M(i, 2) = xValues[i];
    M(i, 3) = pow(xValues[i], 2);
  }
  return M;
}

VectorXd SVD(MatrixXd &matrix, VectorXd &yValues) {
  JacobiSVD<MatrixXd> svd(matrix, ComputeThinU | ComputeThinV);
  VectorXd EInverse = svd.singularValues().array().inverse().matrix();
  VectorXd abcd =
      svd.matrixV() * EInverse.asDiagonal() * svd.matrixU().adjoint() * yValues;
  return abcd;
}

int main() {
  ofstream results("results.txt");

  double xValues[linesInData], yValues[linesInData];
  insertValues(xValues, yValues);

  VectorXd eigenYvalues(linesInData);
  for (int i = 0; i < linesInData; i++) {
    eigenYvalues(i) = yValues[i];
  }

  MatrixXd M = insertData(xValues);
  VectorXd eigenABCD = SVD(M, eigenYvalues);
  double abcd[4] = {eigenABCD(0), eigenABCD(1), eigenABCD(2), eigenABCD(3)};

  double xMin = 0.0;
  double xMax = 9.9;
  double step = (xMax - xMin) / (linesInData - 1);

  for (double i = xMin; i <= xMax; i += step) {
    results << i << "\t" << F(i, abcd) << endl;
  }
  cout << eigenABCD;

  results.close();

  cout << endl;
  //* -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  cout << endl;

  ofstream outputa("outputa.txt");
  ofstream output("output.txt");
  ofstream points("points.txt");

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<double> dis(-2.0, 2.0);

  double xxValues[linesInGData], yyValues[linesInGData], gValues[linesInGData];
  insertValues(xxValues, yyValues);

  double xxMin = 0.0;
  double xxMax = 9.9;
  double sstep = (xxMax - xxMin) / (linesInGData - 1);

  double d[linesInGData];
  for (int i = 0; i < linesInGData; i++) {
    uniform_real_distribution<double> dis(-2.0, 2.0);
    d[i] = dis(gen);
  }

  int ii = 0;
  for (double i = xxMin; i <= xxMax; i += sstep) {
    gValues[ii] = G(i);
    gValues[ii] += d[ii];
    points << i << "\t" << gValues[ii] << endl;
    ii++;
  }

  VectorXd eigenGvalues(linesInGData);
  for (int i = 0; i < linesInGData; i++) {
    eigenGvalues(i) = gValues[i];
  }

  MatrixXd MM = insertGdata(xxValues);
  VectorXd eeigenABCD = SVD(MM, eigenGvalues);
  double aabcd[4] = {eeigenABCD(0), eeigenABCD(1), eeigenABCD(2),
                     eeigenABCD(3)};

  for (double i = xxMin; i <= xxMax; i += sstep) {
    outputa << i << "\t" << G(i, aabcd) << endl;
    output << i << "\t" << G(i) << endl;
  }
  cout << eeigenABCD;

  outputa.close();
  output.close();
  points.close();

  return 0;
}