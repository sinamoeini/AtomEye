/**********************************************************************/
/*  Monte Carlo Simulation   Oct.6,1994  Li Ju                        */
/**********************************************************************/
#include <stdio.h>
#include <math.h>

/*****************************/
/* set constants             */
/*****************************/
#define AVO 6.0225E+23
#define BOLZ 1.38054
#define PI 3.14159265

/***** Krypton ********/
#define SIGMA 3.405
#define EPSI 120.0
#define WTMOL 39.948
#define RCUT 3.0
/**********************/

long NP; /*number of particles*/
long STEPS; /*number of steps*/
float DR; /*reduced density*/
float TR; /*reduced temperature*/
float ALPHA;/*random walk distance*/

float * X;/* position of particles*/
float * BUFFER;
float LR; /* reduce length of a unit cell*/
long NC; /* number of LRs per side   */
float LENGTH; /*length of one side of the bulk*/
double U_all;/* total configuration energy of the system */
long seed; /* random seed */

float frand()
{
 long k=5701;
 long j=75017;
 long m=266927;
 long buf=m*(m/seed)+k*seed*seed+k;
 seed=buf-buf/j*j;
/*printf("%f\n",((float)seed)/j);*/
 return(((float)seed)/j);
 }

/* LR is the reduced length of a unit cell */
/* NC is the number of LRs per side        */
void FCC(float * X, long NP, float LR,long NC)
{
 long i,j,k,l,m=0;
 float *x=X;
 float *y=x+NP;
 float *z=y+NP;
 x[0]=y[0]=z[0]=z[1]=y[2]=x[3]=0;
 x[1]=y[1]=x[2]=z[2]=y[3]=z[3]=LR/2;
 for(i=0;i<NC;i++)
  for(j=0;j<NC;j++)
    for(k=0;k<NC;k++)
    {
     for (l=0;l<=3;l++)
      {
       x[m+l]=x[l]+LR*i;
       y[m+l]=y[l]+LR*j; 
       z[m+l]=z[l]+LR*k;
       }
      m+=4;
     }
  return;
  }  


/*************************************************/
/*  computing potential energy of the current    */
/*  configuration using Lennard-Jones cut off    */
/*  at r=RCUT(defined as constant).              */ 
/*************************************************/
double U(float * X,long NP)
{
 long i,j;
 double u=0;
 double r2,r6,r12,R2=RCUT*RCUT;
 float *x=X;
 float *y=x+NP;
 float *z=y+NP;
 float dx,dy,dz,side=LENGTH;
 float side2=side/2;
 for (i=NP-1;i>=1;i--)
  for (j=i-1;j>=0;j--)
   {
    dx=fabs(x[i]-x[j]);
    dy=fabs(y[i]-y[j]);
    dz=fabs(z[i]-z[j]);
    dx=dx>side2?side-dx:dx;
    dy=dy>side2?side-dy:dy;
    dz=dz>side2?side-dz:dz;
    r2=dx*dx+dy*dy+dz*dz;      
    if (r2<R2)
     {
      r6=r2*r2*r2;
      r12=r6*r6;
      u+=4.*(1/r12-1/r6);
      }
    }
  return (u);
}

/**********************/
/* system initialize  */
/**********************/
void init()
{
 LR=cbrt(4/DR);
 NC=rint(cbrt(NP/4)+0.1);
 LENGTH=LR*NC;
 FCC(X,NP,LR/20*19,NC); /* Not let the particle stick together */
 U_all=U(X,NP);
 seed=37771;
 }

/*************************/
/* if rejected,return 1  */
/*************************/  
int Monte_Carlo()
{
 float *x=X;
 float *y=x+NP;
 float *z=y+NP;
 float *x1=BUFFER;
 float *y1=x1+NP;
 float *z1=y1+NP;
 double UU,dU; 
 long i; 
 for(i=0;i<NP;i++)
  {
   x1[i]=x[i]+ALPHA*frand(); 
   y1[i]=y[i]+ALPHA*frand(); 
   z1[i]=z[i]+ALPHA*frand();
   if ((x1[i]<0)||(x1[i]>LENGTH)||(y1[i]<0)||(y1[i]>LENGTH)||(z1[i]<0)||(z1[i]>LENGTH))
    {
     while (x1[i]<0) x1[i]+=LENGTH;
     while (x1[i]>LENGTH) x1[i]-=LENGTH;
     while (y1[i]<0) y1[i]+=LENGTH;
     while (y1[i]>LENGTH) y1[i]-=LENGTH;
     while (z1[i]<0) z1[i]+=LENGTH;
     while (z1[i]>LENGTH) z1[i]-=LENGTH;
    }
   }
 UU=U(x1,NP);
 dU=UU-U_all;
 if ((dU>0)&&((frand()+1)/2>exp(-1*TR*dU))) return(1) ; /* rejected */
 y1=X;
 X=BUFFER;
 BUFFER=y1;
 U_all=UU;
 return(0);
}

void main()
{ 
  char item[100];
  long i;
  long rejection=0;
  FILE * all;
  printf("Input particle number:\n");
  scanf("%ld",&NP);
  printf("Input reduced density:\n");
  scanf("%f",&DR);
  printf("Input reduced temperature:\n");
  scanf("%f",&TR);
  printf("Input sampling steps:\n");
  scanf("%ld",&STEPS);
  printf("Input ALPHA:\n");
  scanf("%f",&ALPHA);
  if (!(all=fopen("all.out","w+"))) 
   {
    printf ("Unable to creat file all.out...quit.\n");
    exit(0);
   }
  if (!(X=(float *)malloc(sizeof(float)*NP*6)))
   {
    printf("Unable to open such a big memory\n");
    exit(1);
    }
   BUFFER=X+3*NP;
   init();
   for (i=1;i<=STEPS;i++)
    {
     rejection+=Monte_Carlo();
     printf("%6d %13.7f\n",(int)i,(float)U_all);
/*     if (!fwrite(item,strlen(item),1,all))
      {
       printf("Disk full...exit.\n");
       exit(2);
       }*/
     }
   fclose(all);  
   printf("\n\nRejection=%ld\n",rejection);
   return;
  } 
 




