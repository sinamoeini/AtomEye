/**********************************************************************/
/*  Monte Carlo Simulation   Oct.30,1994  Li Ju                       */
/*  This version move one particle at a time.                          */
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
#define RCUT 4.0
/**********************/

/***** Control *******/
#define TURN 20      /* print on file every TURN */
#define GR_DR 0.025  /* the increment of GR sampling */
/********************/
#define RAN_K 5701   /* RAN_? are the constants used in producing randon number */
#define RAN_J 3612
#define RAN_M 566927
/*********************/

long NP; /*number of particles*/
long STEPS; /*number of steps*/
long EQUILIBRIUM; /* After which system should be in equilibrium */
float DR; /*reduced density*/
float TR; /*reduced temperature*/
float ALPHA;/*random walk distance*/

float * X;/* position of particles*/
float * DX;/* diplacement from original place */
long * GR;  /* the G(r) */   
long  GR_STEP; /* the sampling units=LENGTH2/GR_DR */ 
float LR; /* reduced length of a unit cell*/
long NC; /* number of LRs per side   */
double LENGTH; /*length of one side of the bulk*/
double LENGTH2;/*half of LENGTH */
double MLENGTH2;/* minus LENGTH2 */
double U_all;/* total configuration energy of the system */
float seed; /* random seed */


/****************************************************/
/*     function declaration area                    */
/****************************************************/
float frand();
void init();
void FCC(float * X, long NP, float LR,long NC);
float Monte_Carlo();
void PRINT_OUT(long step, FILE * all);
/****************************************************/


/*****************************************/
/* random float number generator 0.- 1.  */
/*****************************************/
float frand()
{
 long irand=(RAN_J*((int)(seed*RAN_M))+RAN_K)%RAN_M;
 return(seed=((float)irand+0.5)/RAN_M);
 }

/* LR is the reduced length of a unit cell */
/* NC is the number of LRs per side        */
void FCC(float * X, long NP, float LR,long NC)
{
 long i,j,k,l,m=0;
 float *x=X;
 float *y=x+NP;
 float *z=y+NP;
 x[0]=y[0]=z[0]=z[1]=y[2]=x[3]=LR/4;  /* keep the atoms from the walls */
 x[1]=y[1]=x[2]=z[2]=y[3]=z[3]=LR*3/4;
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

/**********************/
/* system initialize  */
/**********************/
void init()
{
 long i;
 LR=cbrt(4/DR);
 NC=rint(cbrt(NP/4)+0.1);
 LENGTH=LR*NC;
 LENGTH2=LENGTH/2;
 MLENGTH2=-1*LENGTH2;
 GR_STEP=LENGTH2/GR_DR+1;
 GR=(long *)malloc(sizeof(long)*GR_STEP);
 for (i=0;i<GR_STEP;i++) GR[i]=0; 
 for (i=0;i<3*NP;i++) DX[i]=0;
 FCC(X,NP,LR,NC);
 seed=0.37771;
 }


/* return the accept rate of this round  */
/* Don't laugh at the clumsyness of the  */
/* program, decision tree is enforced to */
/* accelerate.                           */

float Monte_Carlo()
{
 float *x=X;
 float *y=x+NP;
 float *z=y+NP;
 float xx,yy,zz;
 float Dx,Dy,Dz,dx,dy,dz;
 float try_x,try_y,try_z;
 double r2,r6,r12,R2=RCUT*RCUT;
 double U,UU,UUU,UUUU;
 double dU;
 long i,j,accept=0;

 U_all=0; 
 for(i=0;i<NP;i++)
  {
   Dx=ALPHA*(frand()-0.5); 
   Dy=ALPHA*(frand()-0.5); 
   Dz=ALPHA*(frand()-0.5);
   xx=x[i];
   yy=y[i];
   zz=z[i];
   try_x=xx+Dx;
   try_y=yy+Dy;
   try_z=zz+Dz;
   if (try_x<0) try_x+=LENGTH;
   if (try_x>LENGTH) try_x-=LENGTH;
   if (try_y<0) try_y+=LENGTH;
   if (try_y>LENGTH) try_y-=LENGTH;
   if (try_z<0) try_z+=LENGTH;
   if (try_z>LENGTH) try_z-=LENGTH;
   U=UU=0.;

  for (j=i-1;j>=0;j--)
   {
    dx=x[j]-xx;
    dy=y[j]-yy;
    dz=z[j]-zz;
    if (dx<MLENGTH2) dx+=LENGTH;
    if (dy<MLENGTH2) dy+=LENGTH;
    if (dz<MLENGTH2) dz+=LENGTH;
    if (dx>LENGTH2) dx=LENGTH-dx;
    if (dy>LENGTH2) dy=LENGTH-dy;
    if (dz>LENGTH2) dz=LENGTH-dz;
    r2=dx*dx+dy*dy+dz*dz;
    if(r2<R2)
    {      
     r6=1/(r2*r2*r2);
     U+=r6*(r6-1.);
     }
    }
  UUU=U;
 
  for (j=NP-1;j>i;j--)
   {
    dx=x[j]-xx;
    dy=y[j]-yy;
    dz=z[j]-zz;
    if (dx<MLENGTH2) dx+=LENGTH;
    if (dy<MLENGTH2) dy+=LENGTH;
    if (dz<MLENGTH2) dz+=LENGTH;
    if (dx>LENGTH2) dx=LENGTH-dx;
    if (dy>LENGTH2) dy=LENGTH-dy;
    if (dz>LENGTH2) dz=LENGTH-dz;
    r2=dx*dx+dy*dy+dz*dz;
    if(r2<R2)
    {      
     r6=1/(r2*r2*r2);
     U+=r6*(r6-1.);
     }
    }

  for (j=i-1;j>=0;j--)
   {
    dx=x[j]-try_x;
    dy=y[j]-try_y;
    dz=z[j]-try_z;
    if (dx<MLENGTH2) dx+=LENGTH;
    if (dy<MLENGTH2) dy+=LENGTH;
    if (dz<MLENGTH2) dz+=LENGTH;
    if (dx>LENGTH2) dx=LENGTH-dx;
    if (dy>LENGTH2) dy=LENGTH-dy;
    if (dz>LENGTH2) dz=LENGTH-dz;
    r2=dx*dx+dy*dy+dz*dz;
    if(r2<R2)
    {      
     r6=1/(r2*r2*r2);
     UU+=r6*(r6-1.);
     }
    }
    UUUU=UU;
 
  for (j=NP-1;j>i;j--)
   {
    dx=x[j]-try_x;
    dy=y[j]-try_y;
    dz=z[j]-try_z;
    if (dx<MLENGTH2) dx+=LENGTH;
    if (dy<MLENGTH2) dy+=LENGTH;
    if (dz<MLENGTH2) dz+=LENGTH;
    if (dx>LENGTH2) dx=LENGTH-dx;
    if (dy>LENGTH2) dy=LENGTH-dy;
    if (dz>LENGTH2) dz=LENGTH-dz;
    r2=dx*dx+dy*dy+dz*dz;
    if(r2<R2)
    {      
     r6=1/(r2*r2*r2);
     UU+=r6*(r6-1.);
     }
    }

 dU=4.*(UU-U);
 if (dU<0)
  {
   DX[i]+=Dx;
   DX[i+NP]+=Dy;
   DX[i+NP+NP]+=Dz;
   x[i]=try_x;
   y[i]=try_y;
   z[i]=try_z;
   U_all+=4.*UUUU;
   accept+=1;
   }
  else
   {
   if (frand()<exp(-1*dU/TR))
    {
     DX[i]+=Dx;
     DX[i+NP]+=Dy;
     DX[i+NP+NP]+=Dz;
     x[i]=try_x;
     y[i]=try_y;
     z[i]=try_z;
     U_all+=4.*UUUU; 
     accept+=1;
     }
     else U_all+=4.*UUU;
    }
  }
  return ((float)accept/NP);
}


/***************************************************/
/* This is the function to compute the properties  */
/* of the system, and then print out.              */
/***************************************************/
void Print_out(long step, FILE * all)
 {
  char item[150];
  double MSD=0;
  double Pressure,Virial=0;
  double totalx=0;
  double totaly=0;
  double totalz=0;
  long i,j;
  float *x=X;
  float *y=x+NP;
  float *z=y+NP;
  double r2,r6;
  float R2=RCUT*RCUT;
  double R3=LENGTH2*LENGTH2;
  float dx,dy,dz;
  for (i=0;i<NP;i++)
   {
    totalx+=DX[i];
    totaly+=DX[i+NP];
    totalz+=DX[i+NP+NP];
    }
  totalx/=NP;
  totaly/=NP;
  totalz/=NP;
  for (i=0;i<NP;i++)
   {
    DX[i]-=totalx;
    DX[i+NP]-=totaly;
    DX[i+NP+NP]-=totalz;
    MSD+=DX[i]*DX[i]+DX[i+NP]*DX[i+NP]+DX[i+NP+NP]*DX[i+NP+NP];
    } 
  MSD/=NP;
  for (i=NP-1;i>=1;i--)
   for (j=i-1;j>=0;j--)
    {
     dx=x[i]-x[j];
     dy=y[i]-y[j];
     dz=z[i]-z[j];
     if (dx<MLENGTH2) dx+=LENGTH;
     if (dy<MLENGTH2) dy+=LENGTH;
     if (dz<MLENGTH2) dz+=LENGTH;
     if (dx>LENGTH2) dx=LENGTH-dx;
     if (dy>LENGTH2) dy=LENGTH-dy;
     if (dz>LENGTH2) dz=LENGTH-dz;
     r2=dx*dx+dy*dy+dz*dz;
     if ((r2<R3)&&(step>EQUILIBRIUM)) GR[(int)(sqrt(r2)/GR_DR)]++; 
     if(r2<R2)
      {      
       r6=1/(r2*r2*r2);
       Virial+=r6*(r6+r6-1);
       }
      }
  Pressure=TR*DR+DR*Virial*8/NP;
  sprintf(item,"%7ld  %f  %f  %f\n",step,(float)U_all,(float)MSD, (float)Pressure);
  printf(item);
  if (!fwrite(item,strlen(item),1,all))
   {
    printf("Disk full...exit.\n");
    exit(4);
    }
  return;
 }


void main()
{ 
  char item[100];
  long i,j;
  double acceptance=0;
  FILE * all; /* the step,configuration energy,MSD,pressure */
  FILE * gr;  /* the G(r) */
  printf("Input particle number:\n");
  scanf("%ld",&NP);
  printf("Input reduced density:\n");
  scanf("%f",&DR);
  printf("Input reduced temperature:\n");
  scanf("%f",&TR);
  printf("Input sampling steps:\n");
  scanf("%ld",&STEPS);
  printf("After how many steps to equilibrate?\n");
  scanf("%ld",&EQUILIBRIUM);
  printf("Input ALPHA:\n");
  scanf("%f",&ALPHA);
  if (!(all=fopen("all.out","w+"))) 
   {
    printf ("Unable to creat file all.out...quit.\n");
    exit(5);
   }
  if (!(gr=fopen("gr.out","w+"))) 
   {
    printf ("Unable to creat file gr.out...quit.\n");
    exit(1);
   }
  if (!(X=(float *)malloc(sizeof(float)*NP*6)))
   {
    printf("Unable to open such a big memory\n");
    exit(2);
    }
   DX=X+3*NP;
   init();
   for (i=1;i<=STEPS;i++)
    {
     if (i%TURN==1) Print_out(i,all);
     acceptance+=Monte_Carlo();
     }
   fclose(all);
   for (i=1;i<GR_STEP;i++)
     {
      sprintf(item,"%7ld  %f  %f\n",i,i*GR_DR,(float)GR[i]/i/i);
      printf(item);
      if (!fwrite(item,strlen(item),1,gr))
       {
        printf("Disk full...exit\n");
        exit(3);
        }
      }
   fclose(gr);
   printf("NC=%d\n",NC);  
   printf("\n\nAcceptance rate = %f\n",acceptance/STEPS);
   return;                                     
  } 



 






