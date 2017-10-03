/***********************************************************/
/*  TWO BODY QUANTUM SCATTERING USING PARTIAL WAVE METHOD  */
/*  FOR CENTRAL POTENTIAL ONLY        LI JU. MARCH.28      */
/***********************************************************/

#include <math.h>
#include <stdio.h>
#include "nr.h"
/* The Numerical Recipes Package header */

#define PI    3.141592654
#define INTERGRATION 2E5
#define FNAME "cs.out"   /* ouput filename                             */
#define ANGSTEP (PI/500) /* Finery of cross section output             */

double Width,Height,K;
/* The Width of the potential barrier         */
/* The reduced height of the potential barrier     */
/* The incoming wave vector                        */

double W(float x)
{
  return(10000*exp(-x*x));
}

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

int main()
{ 
  FILE * OUTPUT;
  double x,y0,y1,y2,G;
  float j1,j2,n1,n2,a1,a2;
  double k0,k1,k2; 
  double X1,X2,Y1,Y2,f,fi;
  int l,OK,FIRST;
  double WIDTH1, WIDTH2,LMAX,XSTEP;
  double delta[1000];

/* x is the reduced radial length = kr */
/* l is the angular momentum number    */
/* G,j1,j2,n1,n2 are the parameters in */
/* determing phase shift using interpolation */ 
/* delta[LMAX+1] are the matrix which store  */
/* phase shifts for each angular momentum state */
/* W(x) is the reduced central potential        */
/* X1,X2,Y1,Y2 are the points ready for interpoation */
/* y2,y1,y0 are the steps in Numerov intergral method */

  Width  = 4;
  K      = 1;

  WIDTH1 = K*Width;
  XSTEP  = WIDTH1/INTERGRATION;
  WIDTH2 = WIDTH1+XSTEP*1000;
  LMAX   = WIDTH1*2+5;
  OUTPUT = fopen(FNAME,"w");
  for (l=0;l<=LMAX;l++)
    {
      y0 = 0;
      y1 = 1E-300;
      printf ("%E\n",y1);
      k0 = 0;
      x  = XSTEP;
      k1 = (1-l*(l+1)/x/x-W(x)/K/K)*XSTEP*XSTEP/12.;
      
      OK=FIRST=1;
      while (OK)
	{
	  x += XSTEP;
	  k2 = (1-l*(l+1)/x/x-W(x)/K/K)*XSTEP*XSTEP/12.;
	  y2 = (2*(1-5*k1)*y1-(1+k0)*y0)/(1+k2);
	  
	  if ((x>WIDTH1)&&FIRST)
	    {
	      X1 = x;
	      Y1 = y2;
	      sphbes (l,X1,&j1,&n1,&a1,&a2);
	      FIRST = 0;
	    }
	  
	  if (x>WIDTH2)
	    {
	      X2 = x;
	      Y2 = y2;
	      sphbes (l,X2,&j2,&n2,&a1,&a2);
	      G = X1*Y2/X2/Y1;
	      printf ("%lf %lf\n",y1,y2);
              delta[l] = atan((G*j1-j2)/(G*n1-n2));
	      OK = 0;
	    }

	  k0 = k1;
	  y0 = y1;
	  k1 = k2;
	  y1 = y2;
	}            
      printf ("Partial wave No.%d intergrated over. delta = %f\n",l,delta[l]);
    }
  
  /* Differential cross section */
  for (a1=0;a1<PI;a1+=ANGSTEP)
    {
      f = 0.;
      fi = 0.;
      for (l=0;l<=LMAX;l++)
	{
	  f += sin(delta[l])*cos(delta[l])*(2*l+1)*P(l,cos(a1));
	  fi += sin(delta[l])*sin(delta[l])*(2*l+1)*P(l,cos(a1));
	}
      fprintf (OUTPUT,"%f %f \n",a1,(f*f+fi*fi)/K/K);
    }
  
  /* total cross-section */
  a2 = 0;
  for (l=0;l<=LMAX;l++)
    a2 += sin(delta[l])*sin(delta[l])*(2*l+1);
  printf ("The total cross section is %f \n",4*PI*a2/K/K);
  
  fclose(OUTPUT);
  return(1);
}
