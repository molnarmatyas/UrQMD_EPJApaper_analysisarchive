#ifndef _levy_reader_h_
#define _levy_reader_h_

#include <my_includes.h>

class Levy_reader
{
 public:
  Levy_reader(const char* tablefile_name);
  ~Levy_reader();
  double read_table_3d(const double alpha, const double x) const;
  double read_table_1d(const double alpha, const double x) const;
  double getValue_3d(const double alpha, const double x) const; // this is the proper one, with linear interpolation
  double getValue_1d(const double alpha, const double x) const; // this is the proper one, with linear interpolation
 private:
  float** levy_3d_array;
  float** levy_1d_array;
//  float** levy_cyl_array;
  double alpha_min;
  double alpha_max;
  double x_min;
  double x_max;
  int N_alpha;
  int N_x;
  double d_alpha;
  double d_x;

  void InitTable(const char* tablefile_name);
};

#endif // _levy_reader_h_
