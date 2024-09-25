#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

bool compareFirstColumn(const string &a, const string &b) {
  return stod(a) < stod(b);
}

int main(int nargs, char *args[]) {
  ifstream input(args[1]);
  ofstream output(args[2]);

  vector<string> lines;
  string line;

  while (getline(input, line)) {
    lines.push_back(line);
  }

  sort(lines.begin(), lines.end(), compareFirstColumn);

  for (const string &sortedLine : lines) {
    output << sortedLine << "\n";
  }

  input.close();
  output.close();

  return 0;
}