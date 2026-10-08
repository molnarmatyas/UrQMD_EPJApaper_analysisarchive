#include <Levy_proj_reader.h>
using namespace std;

Levy_reader::Levy_reader(const char* tablefile_name)
{
  x_min = 0.;
  x_max = 0.;
  d_x = 0.;
  alpha_min = 0.;
  alpha_max = 0.; 
  d_alpha = 0.; 
  N_x = 0;
  N_alpha = 0;
  levy_3d_array = NULL;
  levy_1d_array = NULL;
//  levy_cyl_array = NULL;
  InitTable(tablefile_name);
} 

Levy_reader::~Levy_reader()
{
  if(levy_3d_array) delete [] levy_3d_array;
//  if(levy_cyl_array) delete [] levy_cyl_array;
}

double Levy_reader::read_table_3d(const double x, const double alpha) const
{
  int ialpha = floor((alpha - alpha_min) / d_alpha); 
  int ix = floor((x - x_min) / d_x); 
  if(ialpha >= 0 && ialpha <= N_alpha)
    if(ix >= 0 && ix <= N_x)
      return (double)levy_3d_array[ialpha][ix];
  return 0.;

}

double Levy_reader::read_table_1d(const double x, const double alpha) const
{
  int ialpha = floor((alpha - alpha_min) / d_alpha); 
  int ix = floor((x - x_min) / d_x); 
  if(ialpha >= 0 && ialpha <= N_alpha)
    if(ix >= 0 && ix <= N_x)
      return (double)levy_1d_array[ialpha][ix];
  return 0.;

}

double Levy_reader::getValue_3d(const double alpha, const double x) const
{
  if(x < x_min || x > x_max) return 0.;
  if(alpha < alpha_min || alpha > alpha_max) return 0.;

  int ix = floor((x - x_min) / d_x);
  double xlo = x_min + ix * d_x;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  if(ix < 0 || ix > N_x || ialpha < 0 || ialpha > N_alpha)
    return 0.;
  if(ix < N_x && ialpha < N_alpha)
  {   
    double xlo_lint = levy_3d_array[ialpha][ix] * (xlo - x + d_x) / d_x + levy_3d_array[ialpha][ix + 1] * (x - xlo) / d_x;
    double xhi_lint = levy_3d_array[ialpha + 1][ix] * (xlo - x + d_x) / d_x + levy_3d_array[ialpha + 1][ix + 1] * (x - xlo) /d_x;
    double alpha_linit_total = xlo_lint * (alphalo - alpha + d_alpha) / d_alpha + xhi_lint * (alpha - alphalo) / d_alpha;
    return alpha_linit_total; 
  }
  if(ix == N_x && ialpha < N_alpha) // Iterpolation only in alpha
    return levy_3d_array[ialpha][ix] * (alphalo - alpha + d_alpha) / d_alpha + levy_3d_array[ialpha + 1][ix] * (alpha - alphalo) / d_alpha;
  if(ialpha == N_alpha && ix < N_x) // Iterpolation only in x
    return levy_3d_array[ialpha][ix] * (xlo - x + d_x) / d_x + levy_3d_array[ialpha][ix + 1] * (x - xlo) / d_x;
  if(ix == N_x && ialpha == N_alpha)
    return levy_3d_array[ialpha][ix];
  return 0.;
}

double Levy_reader::getValue_1d(const double alpha, const double x) const
{
  if(x < x_min || x > x_max) return 0.;
  if(alpha < alpha_min || alpha > alpha_max) return 0.;

  int ix = floor((x - x_min) / d_x);
  double xlo = x_min + ix * d_x;
  int ialpha = floor((alpha - alpha_min) / d_alpha);
  double alphalo = alpha_min + ialpha * d_alpha;
  if(ix < 0 || ix > N_x || ialpha < 0 || ialpha > N_alpha)
    return 0.;
  if(ix < N_x && ialpha < N_alpha)
  {   
    double xlo_lint = levy_1d_array[ialpha][ix] * (xlo - x + d_x) / d_x + levy_1d_array[ialpha][ix + 1] * (x - xlo) / d_x;
    double xhi_lint = levy_1d_array[ialpha + 1][ix] * (xlo - x + d_x) / d_x + levy_1d_array[ialpha + 1][ix + 1] * (x - xlo) /d_x;
    double alpha_linit_total = xlo_lint * (alphalo - alpha + d_alpha) / d_alpha + xhi_lint * (alpha - alphalo) / d_alpha;
    return alpha_linit_total; 
  }
  if(ix == N_x && ialpha < N_alpha) // Iterpolation only in alpha
    return levy_1d_array[ialpha][ix] * (alphalo - alpha + d_alpha) / d_alpha + levy_1d_array[ialpha + 1][ix] * (alpha - alphalo) / d_alpha;
  if(ialpha == N_alpha && ix < N_x) // Iterpolation only in x
    return levy_1d_array[ialpha][ix] * (xlo - x + d_x) / d_x + levy_1d_array[ialpha][ix + 1] * (x - xlo) / d_x;
  if(ix == N_x && ialpha == N_alpha)
    return levy_1d_array[ialpha][ix];
  return 0.;
}

void Levy_reader::InitTable(const char* tablefile_name)
{
  size_t floatsize = sizeof(float);
  size_t intsize = sizeof(int);
  ifstream tablefile;
  tablefile.open(tablefile_name, ios::in | ios::binary);
  cout << "Initializing Levy_reader table in memory using files: " << tablefile_name << endl;
  if(!tablefile.is_open())
  {
    cout << "No table found, initialization failed." << endl;
    return;
  }

  float alpha_min_float = 0.;
  float alpha_max_float = 0.;
  float x_min_float = 0.;
  float x_max_float = 0.;

  tablefile.read((char*)&alpha_min_float, floatsize);
  tablefile.read((char*)&alpha_max_float, floatsize);
  tablefile.read((char*)&N_alpha,         intsize);
  tablefile.read((char*)&x_min_float,     floatsize);
  tablefile.read((char*)&x_max_float,     floatsize);
  tablefile.read((char*)&N_x,             intsize);   

  alpha_min = alpha_min_float;
  alpha_max = alpha_max_float;
  x_min = x_min_float;
  x_max = x_max_float;

  cout << "  alpha_min: " << alpha_min << endl;
  cout << "  alpha_max: " << alpha_max << endl;
  cout << "  N_alpha: "    << N_alpha << endl;
  cout << "  x_min: "     << x_min << endl;
  cout << "  x_max: "     << x_max << endl;
  cout << "  N_x: "        << N_x << endl;

  d_alpha = 1./N_alpha * (alpha_max - alpha_min);
  d_x = 1./N_x * (x_max - x_min);
  cout << endl;

  levy_3d_array = new float*[N_alpha + 1];
  levy_1d_array = new float*[N_alpha + 1];
  for(int ialpha = 0; ialpha < (N_alpha + 1); ialpha++)
  {
    levy_3d_array[ialpha] = new float[N_x + 1];
    levy_1d_array[ialpha] = new float[N_x + 1];
  }

  for(int ix = 0; ix < (N_x + 1); ix++)
    for(int ialpha = 0; ialpha < (N_alpha + 1); ialpha++)
    {
      float l3d_val;
      float l1d_val;
      tablefile.read((char*)&l3d_val, floatsize);
      tablefile.read((char*)&l1d_val, floatsize);
      levy_3d_array[ialpha][ix] = l3d_val;
      levy_1d_array[ialpha][ix] = l1d_val;
    }
  tablefile.close();
}

