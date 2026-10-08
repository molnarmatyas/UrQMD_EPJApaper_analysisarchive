#include <my_includes.h>
#include <complex>
#include <iostream>
using namespace std;
using std::complex;

complex <double> I(0.0, 1.0);

float levy_3d(const double x, const double alpha); // reduced r; i.e. x = r/R 
float levy_1d(const double x, const double alpha); // reduced r; i.e. x = r/R 
//float levy_cyl(const double x, const double alpha); // reduced r; i.e. x = r/R 

//const double limit = 1e-12;
//const int nsamples0 = 200;
//const int nsamples1 = 1000;
const double limit = 1e-18;
const int nsamples0 = 400;
const int nsamples1 = 1500;

const float x_min = 0.;
const float x_max = 100.;
//const int N_x = 100;
const int N_x = 10000;
const float alpha_min = 0.5;
const float alpha_max = 2.0;
//const int N_alpha = 15;
const int N_alpha = 200;

const size_t floatsize = sizeof(float); 
const size_t intsize = sizeof(int); 

int main(int argc, char** argv)
{
  if(argc != 1)
    return 0;
  ofstream tablefile;
  tablefile.open("levy_values_moreprecise.dat", ios::out | ios::binary); 

  tablefile.write((char*)&alpha_min , floatsize);
  tablefile.write((char*)&alpha_max , floatsize);
  tablefile.write((char*)&N_alpha   , intsize);
  tablefile.write((char*)&x_min , floatsize);
  tablefile.write((char*)&x_max , floatsize);
  tablefile.write((char*)&N_x   , intsize);

  for(int ix = 0; ix < (N_x + 1); ix++)
  {
//    if(ix % (N_x / 100) == 0)
      cerr << ix << " ";
    double x =  x_min + (x_max - x_min) * ix * 1. / N_x;
    for(int ialpha = 0; ialpha < (N_alpha + 1); ialpha++)
    {
      double alpha =  alpha_min + (alpha_max - alpha_min) * ialpha * 1. / N_alpha;
      float levy_3d_value = levy_3d(x, alpha);
      float levy_1d_value = levy_1d(x, alpha);
//      float levy_cyl_value = levy_sph(x, alpha);
      tablefile.write((char*)&levy_3d_value, floatsize);
      tablefile.write((char*)&levy_1d_value, floatsize);
//      tablefile.write((char*)&levy_cyl_value, floatsize);
    }
  }

  tablefile.close();
  return 0;
}

float levy_3d(const double x, const double alpha)
{
   if(x < 1e-16)
     return (float)(1./2./SQR(M_PI)/alpha*pow(2., 3./alpha) * tgamma(3./alpha));
   double rmax = (alpha < 1.) ? 0.002 + 4.* SQR(alpha - 0.5) : 1.5 + 13.* (alpha - 1.);
   if((1e-16 < x) && (x < rmax))
   {
     double qmin = 0.;
     double qmax = pow(-2. * log(limit), 1./alpha);
     double dq1 = 2. * M_PI / x * 1. / nsamples0;
     double dq2 = pow(2., 1./alpha) * 1. / nsamples0;
     double dq = (dq1 > dq2) ? dq2 : dq1;
     double F = 0.;
     for(double qi = qmin; qi < qmax; qi += dq)
       F += qi * sin(qi * x) * exp(-0.5 * pow(qi, alpha));
     double c = 1. / (2. * SQR(M_PI) * x) * dq;
     return (float)(c * F);
   }
   else
   {
     double tmax = 3.*50./2./sqrt(2.)/x; // 50 stands for e^{-50}, thiswas investigated on 05/18/2023.
//    double dt = 0.2 * M_PI / r * 1. / nsamples1;
     double dt = 2. * M_PI / x * 1. / nsamples1;
     double ti = dt;
     double X = SQR(alpha) * pow(ti, alpha - 1.) / (2. * x);
//    double X = SQR(alpha) * Q;
//    double up = r *ti *(sqrt(SQR(X) + 2.) * (SQR(X) + 2.) / 3. + Q - X- 1./3. * pow(X, 3));
     complex <double> F(0.,0.);
     while(ti < tmax)
     {
       double fi = sqrt(X*X + 2.) - X;
       double dfi = (alpha-1.) * (pow(alpha,4) /4. /SQR(x) * pow(ti,(2.* alpha - 3.)) * 1. / sqrt(X*X + 2.) - X / ti);
       complex<double> eiphi = exp(I * fi);
       complex<double> eiaphi = exp(I * alpha * fi);
       F += SQR(eiphi) * (1. + ti * I * dfi) * ti * exp( I * ti * x *eiphi - 0.5 * pow(ti, alpha) * eiaphi);
       ti += dt;
       X = SQR(alpha) * pow(ti, alpha - 1.) / (2. * x);
//      up = r *ti *(sqrt(SQR(X) + 2) *(SQR(X) + 2) / 3. + Q - X - 1./3.* pow(X, 3));
     }
     double c = 1./(2. * SQR(M_PI) * x) * dt;
     return (float)(c * imag(F));
   }

} 

float levy_1d(const double x, const double alpha)
{
   if(x < 1e-16)
     return (float)(1./M_PI/alpha * pow(2., 1./alpha) * tgamma(1./alpha));
   double rmax = (alpha < 1.) ? 0.002 + 4.* SQR(alpha - 0.5) : 1.5 + 13.* (alpha - 1.);
   if((1e-16 < x) && (x < rmax))
   {
     double qmin = 0.;
     double qmax = pow(-2. * log(limit), 1./alpha);
     double dq1 = 2. * M_PI / x * 1. / nsamples0;
     double dq2 = pow(2., 1./alpha) * 1. / nsamples0;
     double dq = (dq1 > dq2) ? dq2 : dq1;
     double F = 0.;
     for(double qi = qmin; qi < qmax; qi += dq)
       F += cos(qi * x) * exp(-0.5 * pow(qi, alpha));
     double c = 1. / M_PI * dq;
     return (float)(c * F);
   }
   else
   {
     double tmax = 3.*50./2./sqrt(2.)/x; // 50 stands for e^{-50}, thiswas investigated on 05/18/2023.
//    double dt = 0.2 * M_PI / r * 1. / nsamples1;
     double dt = 2. * M_PI / x * 1. / nsamples1;
     double ti = dt;
     double X = SQR(alpha) * pow(ti, alpha - 1.) / (2. * x);
//    double X = SQR(alpha) * Q;
//    double up = r *ti *(sqrt(SQR(X) + 2.) * (SQR(X) + 2.) / 3. + Q - X- 1./3. * pow(X, 3));
     complex <double> F(0.,0.);
     while(ti < tmax)
     {
       double fi = sqrt(X*X + 2.) - X;
       double dfi = (alpha-1.) * (pow(alpha,4) /4. /SQR(x) * pow(ti,(2.* alpha - 3.)) * 1. / sqrt(X*X + 2.) - X / ti);
       complex<double> eiphi = exp(I * fi);
       complex<double> eiaphi = exp(I * alpha * fi);
       F += eiphi * (1. + ti * I * dfi) * exp( I * ti * x *eiphi - 0.5 * pow(ti, alpha) * eiaphi);
       ti += dt;
       X = SQR(alpha) * pow(ti, alpha - 1.) / (2. * x);
//      up = r *ti *(sqrt(SQR(X) + 2) *(SQR(X) + 2) / 3. + Q - X - 1./3.* pow(X, 3));
     }
     double c = 1./M_PI * dt;
     return (float)(c * real(F));
   }
} 

//float levy_cyl(const double x, const double alpha)
//{
//  if(x < 1e-10)
//    return 1./2./SQR(M_PI)/alpha*pow(2., 3./alpha) * Gamma(3./alpha); // EZT MÉG ELLENŐRIZZÜK!!!!!
//  double qmin = 0.;
//  double qmax = pow(-2. * log(limit), 1./alpha);
//  double dq1 = 2. * M_PI / x * 1. / nsamples;
//  double dq2 = pow(2., 1./alpha) * 1. / nsamples;
//  double dq = (dq1 > dq2) ? dq2 : dq1;
//  double F = 0.;
//  for(double qi = qmin; qi < qmax; qi += dq)
//    F += qi * cyl_bessel_j(0., qi * x) * exp(-0.5 * pow(qi, alpha));
//  double c = 1. / (2. * SQR(M_PI) * x) * dq; 
//  return (float)(c * F);
//} 
