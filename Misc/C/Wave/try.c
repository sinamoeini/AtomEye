#include <math.h>
#include "nr.h"


/* Legendre polynomial */
float P(int l,float x)
{
  float i,p0,p1,p2;
  
  if (l==0) return(1);
  if (l==1) return(x);
  
  p0 = 1;
  p1 = x;
  for (i=1;i<l;i++)
    {
      p2 = ((2*i+1)*x*p1-i*p0)/(i+1);
      p0 = p1;
      p1 = p2;
    }
  return(p1);
}

void main()
{
  float a,b,c,d;
  int l;
  
  for (l=0;l<=100;l++) printf ("%f\n",P(l,1.0));
}
