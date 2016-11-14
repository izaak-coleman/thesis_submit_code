#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mathutil.h"

double z(unsigned int n1, unsigned int n2, double u)
{
  double mu, sigma;
  mu = (n1*n2)/2.0;
  sigma = sqrt((n1*n2)*(n1 + n2 + 1.0)/12.0);
  return fabs(u - mu)/sigma;
}

#define A1	0.31938153
#define A2	-0.356563782
#define A3	1.781477937
#define A4	-1.821255978
#define A5	1.330274429
#define B	0.2316419
#define PI	3.1415927
#define ONE_DIVIDED_BY_SQRT_2PI 0.39894228
#define THRESHOLD 0.05

double exactP(double z)
{
  double t = 1.0 / (1.0 + B * z);
  double t2 = t*t;
  double t3 = t2*t;
  double t4 = t3*t;
  double t5 = t4*t;
  double sum = A1*t + A2 *t2 + A3*t3 + A4*t4 + A5*t5;

  return ONE_DIVIDED_BY_SQRT_2PI * exp((-1.0)*z*z/2.0) * sum;
}

double z2p(double z)
{
  double p;

  p = 0.5 * (1.0 - sqrt(1 - exp(-2.0*z*z/PI)));

  if (p  >= THRESHOLD) return p;
  else return exactP(z);
}
