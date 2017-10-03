/*********************************************************/
/* One-dimensional time dependent Schoedinger Equation.  */
/* The scattering of a wavepack on a square barrier.     */
/* Use Crank-Nicholson/Caylay algorithm...               */
/* Stable, Norm Conserving.     Li Ju. March.28,1995     */
/*********************************************************/

#include <math.h>
#include <stdio.h>

#define MESH 1000    /* # of points from 0-1 */
#define TSTEP 2E-6   /* Time inteval         */

float buffer[6*MESH]; 
float *a2,*b2,*a1,*b1,*a0,*b0;

/* buffer is where wavefunctions are kept */
/* (a0,b0)/(a1,b1) are the real and imaginery */
/* part of wave function */



/*
** Potential function
*/
float V(int i)
{
  if (abs(i-500)<32) return (100);
  else return(0);
}

/*
** Initialize the wave pack.
*/
void InitWave()
{
  int i;
  float module;
  for (i=1;i<MESH-1;i++)
    {
      module = exp(-((i-250)/40)*((i-250)/40)/2.);
/* The wavepack is Guassian, character width = 40 */
/* It's also moving,  wave vector = 0.05          */

      a0[i] = module * sin(0.05*i/MESH);
      b0[i] = module * cos(0.05*i/MESH);
    }
  a0[0]=a0[MESH-1]=b0[0]=b0[MESH-1]=0;
  a1[0]=a1[MESH-1]=b1[0]=b1[MESH-1]=0;
  return;
}


/*
** The explicit operator part
** of the norm-conserving Cayley propagator.
*/
void Explicit()
{
 int i,w;
 for (i=1;i<MESH-1;i++)
  {
   w=(a0[i-1]-2*a0[i]+a0[i+1])*MESH*MESH-V(i)*a0[i];
   v=(b0[i-1]-2*b0[i]+b0[i+1])*MESH*MESH-V(i)*b0[i];
   a1[i] = a0[i]+v*TSTEP;
   b1[i] = b0[i]-w*TSTEP;
 }
 return();
}


/* 
** The implicit operator part 
** of the norm-conserving Cayley propagator.
*/
void Implicit()
{
 float alpha[MESH+2],ialpha[MESH+2];
 float beta[MESH+2],ibeta[MESH+2];
 float gamma[MESH+2],igamma[MESH+2];

 int i,N=MESH-2;

 alpha[N-1]=ialpha[N-1]=0;
 beta[N-1] = a1[N];
 ibeta[N-1] = b1[N];

 for (i=N-2;i>=0;i--)
  {
    


 }
 return();
}


int main()
{
 long Steps;
 float Width,Height;
/* Total time step,Barrier width,height */
 
 scanf ("%ld %f %f",&Steps,&Width,&Height);
 a0 = buffer;
 b0 = a0+MESH;
 a1 = b0+MESH;
 b1 = a1+MESH;
 a2 = b1+MESH;
 b2 = a2+MESH;

 InitWave();  

 



