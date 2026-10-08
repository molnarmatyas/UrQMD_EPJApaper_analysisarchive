#include "Levy_reader.h"
using namespace std;

Levy_reader::Levy_reader(const char* filename)
{
  levy_array = NULL;
  InitTable(filename);
}

Levy_reader::~Levy_reader()
{
  if(levy_array) delete [] levy_array;
}

double Levy_reader::read_table(const double alpha, const double x) const
{
  if(!levy_array)
    return 0;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  int ix = floor((x - x_min) / d_x);
  if(ialpha >= 0 && ialpha < (Nalpha + 1))
    if(ix >= 0 && ix < (Nx + 1))
      return levy_array[ialpha][ix];
  return 0;
}

double Levy_reader::getValue(const double alpha, const double x) const
{
  if(!levy_array)
    return 0;
  int ix = floor((x - x_min) / d_x);
  double xlo = x_min + ix * d_x;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  if(ix < 0 || ix > Nx || ialpha < 0 || ialpha > Nalpha)
    return 0;
  if(ix < Nx && ialpha < Nalpha)
  {
    double xlo_lint = levy_array[ialpha][ix] * (xlo - x + d_x) / d_x + levy_array[ialpha][ix + 1] * (x - xlo) / d_x;
    double xhi_lint = levy_array[ialpha + 1][ix] * (xlo - x + d_x) / d_x + levy_array[ialpha + 1][ix + 1] * (x - xlo) /d_x;
    double alpha_linit_total = xlo_lint * (alphalo - alpha + d_alpha) / d_alpha + xhi_lint * (alpha - alphalo) / d_alpha;
    return alpha_linit_total;
  }
  if(ix == Nx && ialpha < Nalpha) // Iterpolation only in alpha
    return levy_array[ialpha][ix] * (alphalo - alpha + d_alpha) / d_alpha + levy_array[ialpha + 1][ix] * (alpha - alphalo) / d_alpha;
  if(ialpha == Nalpha && ix < Nx) // Iterpolation only in x
    return levy_array[ialpha][ix] * (xlo - x + d_x) / d_x + levy_array[ialpha][ix + 1] * (x - xlo) / d_x;
  if(ix == Nx && ialpha == Nalpha)
    return levy_array[ialpha][ix];
  return 0;
}

void Levy_reader::InitTable(const char* filename)
{
  ifstream infile;
  infile.open(filename, ios::in | ios::binary);
  if(!infile.is_open()) return;
  cout << "Initializing Levy_reader table in memory using file: " << filename << endl;
  infile.read((char*)&alpha_min, sizeof(alpha_min));  cout << "  alpha_min: " << alpha_min << endl;
  infile.read((char*)&alpha_max, sizeof(alpha_max));  cout << "  alpha_max: " << alpha_max << endl;
  infile.read((char*)&Nalpha,    sizeof(Nalpha));     cout << "  Nalpha: "    << Nalpha << endl;
  infile.read((char*)&x_min,     sizeof(x_min));      cout << "  x_min: "     << x_min << endl;
  infile.read((char*)&x_max,     sizeof(x_max));      cout << "  x_max: "     << x_max << endl;
  infile.read((char*)&Nx,        sizeof(Nx));         cout << "  Nx: "        << Nx << endl;
  d_alpha = (alpha_max - alpha_min) * 1.0 / Nalpha;
  d_x = (x_max - x_min) * 1.0 / Nx;
  cout << "    d_alpha: " << d_alpha << " d_x: " << d_x << endl;

  double value = 0;
  levy_array = new double*[Nalpha + 1];
  for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)
    levy_array[ialpha] = new double[Nx + 1];

  int ix = -1, Nx_filled = 0;
  while(Nx_filled < Nx)
  {
    infile.read((char*)&ix, sizeof(ix));
    for(int ialpha = 0; ialpha < (Nalpha + 1); ialpha++)
    {
      infile.read((char*)&value, sizeof(value));
      levy_array[ialpha][ix] = value;
    }
    Nx_filled++;
  }
  infile.close();
}

