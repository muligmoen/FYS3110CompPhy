#include <cmath>

const int MAXIT = 10;
const double EPS = 3e-14;

double gammln(double);

void gauss_laguerre(double *x_return, double *w_return,
                    const int n, const double alf)
{
  double *x = new double[n+1];
  double *w = new double[n+1];
  double p2, p3, pp, z;

  for (int i=1;i<=n;i++) {
    if (i == 1) {
      z=(1.0+alf)*(3.0+0.92*alf)/(1.0+2.4*n+1.8*alf);
    } else if (i == 2) {
      z += (15.0+6.25*alf)/(1.0+0.9*alf+2.5*n);
    } else {
      double ai=i-2;
      z += ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/
              (1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
    }
    for (int its=0;its<MAXIT;its++) {
      double p1=1.0;
      //double p2=0.0;
      p2 = 0.0;
      for (int j=1;j<=n;j++) {
        p3=p2;
        p2=p1;
        p1=((2*j-1+alf-z)*p2-(j-1+alf)*p3)/j;
      }
      pp=(n*p1-(n+alf)*p2)/z;
      double z1 = z;
      z = z1 - p1/pp;
      if (std::abs(z-z1) <= EPS) break;
    }
    //int N = i-1;
    //x[N] = z;
    //w[N] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
    x[i] = z;
    w[i] = -exp(gammln(alf+n)-gammln((double)n))/(pp*n*p2);
  }
  for (int ii=0; ii<n; ii++)
  {
    x_return[ii] = x[ii+1];
    w_return[ii] = w[ii+1];
  }
  delete[] x;
  delete[] w;
}


double gammln( double xx)
{
  const static double cof[6]={76.18009172947146,-86.50532032941677,
                        24.01409824083091,-1.231739572450155,
                        0.1208650973866179e-2,-0.5395239384953e-5};

  double y = xx;
  double x = xx;
  double tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  double ser=1.000000000190015;
  for (int j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}