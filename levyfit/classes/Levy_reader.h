#ifndef Levy_reader_H
#define Levy_reader_H

#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <cmath>
using namespace std;

class Levy_reader
{
 public:
  Levy_reader(const char* filename);
  virtual ~Levy_reader();
  double read_table(const double alpha, const double x) const;
  double getValue(const double alpha, const double x) const; // this is the proper one, with linear interpolation

 private:
  double** levy_array;
  double alpha_min;
  double alpha_max;
  double x_min;
  double x_max;
  int Nalpha;
  int Nx;
  double d_alpha;
  double d_x;

  void InitTable(const char* filename);
};

#endif // Levy_reader_H

