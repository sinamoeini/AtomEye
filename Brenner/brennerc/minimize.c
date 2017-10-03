/*
 This code was released into the public domain by Peter McCluskey on 7/11/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/
/*     This used to have an "EMIN" argument, but that argument was not
 *     used, and then we got compiler errors when other routines called
 *    minimize with no argument. 
 */

#include "brenner.h"
#include <math.h>
#include <stdlib.h>

static int nmin,idmin[MAX_ATOMS+1];

static void
conmin(BrennerMainInfo *, int n, double *x, double *f, double *g,
       int *ifun, int *iter, double eeps, int *nflag, int mxfun, double *wr,
       int iout, int mdim, FILE *idev, double acc, int nmeth);

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef min
#define min(a,b) ((a) <= (b) ? (a) : (b))
#endif
#ifndef max
#define max(a,b) ((a) >= (b) ? (a) : (b))
#endif

void minimize(BrennerMainInfo *info)
{
  double *rr = (double *)malloc((3*info->num_atms + 1)*sizeof(double));
  double *g = (double *)malloc((3*info->num_atms + 1)*sizeof(double));
  double *wr = (double *)malloc((15*info->num_atms + 3)*sizeof(double));
  double eeps = 1.0e-7;
  double acc = 1.0e-20;
  double f;
  int mxfun = 1500;
  int i, n;
  int mdim, nmeth;
  int ifun, iter, nflag;
  int iout = 1;

  for(i = 1; i <= 3*info->num_atms; ++i)
  {
    rr[i] = 0;
    g[i] = 0;
  }
  for(i = 1; i <= 15*info->num_atms+2; ++i)
    wr[i] = 0;
  /*
   * Set up number of atoms and their id's to minimize
   */
  nmin = 0;
  for(i = 1; i <= info->num_atms; ++i)
  {
    if(info->atm_num[i].thermostated)
      idmin[++nmin] = i;
  }

  info->lchk = 1;
  find_atom_energies(info);
  info->lchk = choose_lj(info);
  f = info->system_energy;
  printf(" INITIAL ENERGY IS %f\n",f);

  for(i = 1; i <= nmin; ++i)
  {
    int id = idmin[i];
    rr[i] = info->atm_num[id].coord.x;
    g[i] = -info->atm_num[id].force.x;
    rr[nmin+i] = info->atm_num[id].coord.y;
    g[nmin+i] = -info->atm_num[id].force.y;
    rr[2*nmin+i] = info->atm_num[id].coord.z;
    g[2*nmin+i] = -info->atm_num[id].force.z;
  }
  n = 3*nmin;

  mdim = 5*n+2;
  nmeth = 0;
  /* IF (NMETH.EQ.1) MDIM = N*(N+7)/2 */

  conmin(info, n, rr, &f, g, &ifun, &iter, eeps, &nflag, mxfun,
	 wr, iout, mdim, stdout, acc ,nmeth);
  printf("nflag= %d\n",nflag);

  for(i = 1; i <= nmin; ++i)
  {
    int id = idmin[i];
#ifndef INFINITE_CUBE
    rr[i] -= info->cube[0]*floor(rr[i]/info->cube[0] + 0.5);
    rr[nmin+i] -= info->cube[1]*floor(rr[nmin+i]/info->cube[1] + 0.5);
    rr[2*nmin+i] -= info->cube[2]*floor(rr[2*nmin+i]/info->cube[2] + 0.5);
#endif
    info->atm_num[id].coord.x = rr[i];
    info->atm_num[id].coord.y = rr[nmin+i];
    info->atm_num[id].coord.z = rr[2*nmin+i];
    printf("minimize %2d %2d %f %f %f\n",i,id,rr[i],rr[nmin+i],rr[2*nmin+i]);
  }
  free(wr);
  free(rr);
  free(g);
}

static void
calcfg(BrennerMainInfo *info, double *rr, double *f, double *g)
{
  int i;
  BrenAtom *atm_num = info->atm_num;
  for(i = 1; i <= nmin; ++i)
  {
    int id = idmin[i];
    atm_num[id].coord.x = rr[i];
    atm_num[id].coord.y = rr[nmin+i];
    atm_num[id].coord.z = rr[2*nmin+i];
  }

  info->lchk = 1;
  info->system_energy = 0.0;
  find_atom_energies(info);
  info->lchk = choose_lj(info);
  *f = info->system_energy;

  for(i = 1; i <= nmin; ++i)
  {
    int id = idmin[i];
    g[i] = -atm_num[id].force.x;
    if(0)printf("g[%2d] %2d %11.9f\n",i,id,g[i]);
    g[nmin+i] = -atm_num[id].force.y;
    g[2*nmin+i] = -atm_num[id].force.z;
  }
}

static void
conmin(BrennerMainInfo *info, int n, double *x, double *f, double *g,
       int *ifun, int *iter, double eeps, int *nflag, int mxfun, double *wr,
       int iout, int mdim, FILE *idev, double acc, int nmeth)
{
/*
* PURPOSE:    SUBROUTINE CONMIN MINIMIZES AN UNCONSTRAINED NONLINEAR
*             SCALAR VALUED FUNCTION OF A VECTOR VARIABLE X
*             EITHER BY THE BFGS VARIABLE METRIC ALGORITHM OR ON A
*             BEALE RESTARTED CONJUGATE GRADIENT ALGORITHM.
*
* USAGE:      CALL CONMIN(N,X,F,G,IFUN,ITER,EEPS,NFLAG,MXFUN,W,
*             IOUT,MDIM,IDEV,ACC,NMETH)
*
* PARAMETERS: N      THE NUMBER OF VARIABLES IN THE FUNCTION TO
*                    BE MINIMIZED.
*             X      THE VECTOR CONTAINING THE CURRENT ESTIMATE TO
*                    THE MINIMIZER. ON ENTRY TO CONMIN,X MUST CONTAIN
*                    AN INITIAL ESTIMATE SUPPLIED BY THE USER.
*                    ON EXITING,X WILL HOLD THE BEST ESTIMATE TO THE
*                    MINIMIZER OBTAINED BY CONMIN. X MUST BE DOUBLE
*                    PRECISIONED AND DIMENSIONED N.
*             F      ON EXITING FROM CONMIN,F WILL CONTAIN THE LOWEST
*                    VALUE OF THE OBJECT FUNCTION OBTAINED.
*                    F IS DOUBLE PRECISIONED.
*             G      ON EXITING FROM CONMIN,G WILL CONTAIN THE
*                    ELEMENTS OF THE GRADIENT OF F EVALUATED AT THE
*                    POINT CONTAINED IN X. G MUST BE DOUBLE
*                    PRECISIONED AND DIMENSIONED N.
*             IFUN   UPON EXITING FROM CONMIN,IFUN CONTAINS THE
*                    NUMBER OF TIMES THE FUNCTION AND GRADIENT
*                    HAVE BEEN EVALUATED.
*             ITER   UPON EXITING FROM CONMIN,ITER CONTAINS THE
*                    TOTAL NUMBER OF SEARCH DIRECTIONS CALCULATED
*                    TO OBTAIN THE CURRENT ESTIMATE TO THE MINIMIZER.
*             EEPS    EEPS IS THE USER SUPPLIED CONVERGENCE PARAMETER.
*                    CONVERGENCE OCCURS WHEN THE NORM OF THE GRADIENT
*                    IS LESS THAN OR EQUAL TO EEPS TIMES THE MAXIMUM
*                    OF ONE AND THE NORM OF THE VECTOR X. EEPS
*                    MUST BE DOUBLE PRECISIONED.
*             NFLAG  UPON EXITING FROM CONMIN,NFLAG STATES WHICH
*                    CONDITION CAUSED THE EXIT.
*                    IF NFLAG=0, THE ALGORITHM HAS CONVERGED.
*                    IF NFLAG=1, THE MAXIMUM NUMBER OF FUNCTION
*                       EVALUATIONS HAVE BEEN USED.
*                    IF NFLAG=2, THE LINEAR SEARCH HAS FAILED TO
*                       IMPROVE THE FUNCTION VALUE. THIS IS THE
*                       USUAL EXIT IF EITHER THE FUNCTION OR THE
*                       GRADIENT IS INCORRECTLY CODED.
*                    IF NFLAG=3, THE SEARCH VECTOR WAS NOT
*                       A DESCENT DIRECTION. THIS CAN ONLY BE CAUSED
*                       BY ROUNDOFF,AND MAY SUGGEST THAT THE
*                       CONVERGENCE CRITERION IS TOO STRICT.
*             MXFUN  MXFUN IS THE USER SUPPLIED MAXIMUM NUMBER OF
*                    FUNCTION AND GRADIENT CALLS THAT CONMIN WILL
*                    BE ALLOWED TO MAKE.
*             W      W IS A VECTOR OF WORKING STORAGE.IF NMETH=0,
*                    W MUST BE DIMENSIONED 5*N+2. IF NMETH=1,
*                    W MUST BE DIMENSIONED N*(N+7)/2. IN BOTH CASES,
*                    W MUST BE DOUBLE PRECISIONED.
*             IOUT   IOUT IS A USER  SUPPLIED OUTPUT PARAMETER.
*                    IF IOUT = 0, THERE IS NO PRINTED OUTPUT FROM
*                    CONMIN. IF IOUT > 0,THE VALUE OF F AND THE
*                    NORM OF THE GRADIENT SQUARED,AS WELL AS ITER
*                    AND IFUN,ARE WRITTEN EVERY IOUT ITERATIONS.
*             MDIM   MDIM IS THE USER SUPPLIED DIMENSION OF THE
*                    VECTOR W. IF NMETH=0,MDIM=5*N+2. IF NMETH=1,
*                    MDIM=N*(N+7)/2.
*             IDEV   IDEV IS THE USER SUPPLIED NUMBER OF THE OUTPUT
*                    DEVICE ON WHICH OUTPUT IS TO BE WRITTEN WHEN
*                    IOUT>0.
*             ACC    ACC IS A USER SUPPLIED ESTIMATE OF MACHINE
*                    ACCURACY. A LINEAR SEARCH IS UNSUCCESSFULLY
*                    TERMINATED WHEN THE NORM OF THE STEP SIZE
*                    BECOMES SMALLER THAN ACC. IN PRACTICE,
*                    ACC=10.D-20 HAS PROVED SATISFACTORY. ACC IS
*                    DOUBLE PRECISIONED.
*             NMETH  NMETH IS THE USER SUPPLIED VARIABLE WHICH
*                    CHOOSES THE METHOD OF OPTIMIZATION. IF
*                    NMETH=0,A CONJUGATE GRADIENT METHOD IS
*                    USED. IF NMETH=1, THE BFGS METHOD IS USED.
*
* REMARKS:    IN ADDITION TO THE SPECIFIED VALUES IN THE ABOVE
*             ARGUMENT LIST, THE USER MUST SUPPLY A SUBROUTINE
*             CALCFG WHICH CALCULATES THE FUNCTION AND GRADIENT AT
*             X AND PLACES THEM IN F AND G(1),...,G(N) RESPECTIVELY.
*             THE SUBROUTINE MUST HAVE THE FORM:
*                    SUBROUTINE CALCFG(N,X,F,G)
*                    DOUBLE PRECISION X(N),G(N),F
*
*             AN EXAMPLE SUBROUTINE FOR THE ROSENBROCK FUNCTION IS:
*
*                    SUBROUTINE CALCFG(N,X,F,G)
*                    DOUBLE PRECISION X(N),G(N),F,T1,T2
*                    T1=X(2)-X(1)*X(1)
*                    T2=1.0-X(1)
*                    F=100.0*T1*T1+T2*T2
*                    G(1)=-400.0*T1*X(1)-2.0*T2
*                    G(2)=200.0*T1
*                    RETURN
*                    END
*      DOUBLE PRECISION F,FP,FMIN,ALPHA,AT,AP,GSQ,DG,DG1
*      DOUBLE PRECISION DP,STEP,ACC,DAL,U1,U2,U3,U4,EEPS
*      DOUBLE PRECISION XSQ,RTST,DSQRT,DMIN1,DMAX1,DABS
*      LOGICAL RSW
*/
  int i, j, ii, ij;
  int nx = 0, ng = 0;
  int nry = 0, nrd = 0, ncons = 0, ncons1 = 0, ncons2 = 0, nrst = 0;
  int ncalls = 0, nxpi = 0, ngpi = 0, nrdpi = 0, nrypi = 0, ngpj = 0;
  int rsw;
  int ioutk = 1;
  double dg1 = 0, xsq = 0, gsq = 0;
  double fmin = 0, alpha = 0, dal = 0;
  double dg = 0, ap = 0, fp = 0, dp = 0, at = 0;
  double step = 0, rtst = 0;
  double u1 = 0, u2 = 0, u3 = 0, u4 = 0;
  int cntxx = 0;
  printf("N IN CONMIN= %d\n",n);
  /*
   * INITIALIZE ITER,IFUN,NFLAG,AND IOUTK,WHICH COUNTS OUTPUT ITERATIONS.
   */
  *iter = 0;
  *ifun = 0;
  *nflag = 0;
/*
* SET PARAMETERS TO EXTRACT VECTORS FROM W.
* WR(I) HOLDS THE SEARCH VECTOR,WR(NX+I) HOLDS THE BEST CURRENT
* ESTIMATE TO THE MINIMIZER,AND WR(NG+I) HOLDS THE GRADIENT
* AT THE BEST CURRENT ESTIMATE.
*/
  nx = n;
  ng = nx + n;
/*
* TEST WHICH METHOD IS BEING USED.
* IF NMETH=0, WR(NRY+I) HOLDS THE RESTART Y VECTOR AND
* WR(NRD+I) HOLDS THE RESTART SEARCH VECTOR.
*/
  if(nmeth != 1)
  {
    nry = ng+n;
    nrd = nry+n;
    ncons = 5*n;
    ncons1 = ncons+1;
    ncons2 = ncons+2;
  }
  else
  {
/*
* IF NMETH=1,WR(NCONS+I) HOLDS THE APPROXIMATE INVERSE HESSIAN.
*/
    ncons = 3*n;
  }
/*
* CALCULATE THE FUNCTION AND GRADIENT AT THE INITIAL
* POINT AND INITIALIZE NRST,WHICH IS USED TO DETERMINE
* WHETHER A BEALE RESTART IS BEING DONE. NRST=N MEANS THAT THIS
* ITERATION IS A RESTART ITERATION. INITIALIZE RSW,WHICH INDICATES
* THAT THE CURRENT SEARCH DIRECTION IS A GRADIENT DIRECTION.
*/
 label20:
  calcfg(info,x,f,g);
  ++*ifun;
  nrst = n;
  rsw = TRUE;
/*
* CALCULATE THE INITIAL SEARCH DIRECTION , THE NORM OF X SQUARED,
* AND THE NORM OF G SQUARED. DG1 IS THE CURRENT DIRECTIONAL
* DERIVATIVE,WHILE XSQ AND GSQ ARE THE SQUARED NORMS.
*/
  dg1 = 0;
  xsq = 0;
  for(i = 1; i <= n; ++i)
  {
    wr[i] = -g[i];
    xsq += x[i]*x[i];
    dg1 -= g[i]*g[i];
  }
  gsq = -dg1;
  printf("gsq.a %7.3f\n", gsq);
/*
* TEST IF THE INITIAL POINT IS THE MINIMIZER.
*/
  if(gsq <= eeps*eeps*max(1.0,xsq)) return;
/*
* BEGIN THE MAJOR ITERATION LOOP. NCALLS IS USED TO GUARANTEE THAT
* AT LEAST TWO POINTS HAVE BEEN TRIED WHEN NMETH=0. FMIN IS THE
* CURRENT FUNCTION VALUE.
*/
 label40:
  fmin = *f;
  ncalls = *ifun;
/*
* IF OUTPUT IS DESIRED,TEST IF THIS IS THE CORRECT ITERATION
* AND IF SO, WRITE OUTPUT.
*/
      if(iout)
      {
	if(!ioutk)
	  fprintf(idev,
		  "  ITER %5d FNCTN CALLS %6d F =%15.8f  G-SQUARED = %15.8f\n",
		  *iter, *ifun, fmin, gsq);
	if(++ioutk == iout) ioutk = 0;
      }
/*
* BEGIN LINEAR SEARCH. ALPHA IS THE STEPLENGTH.
* SET ALPHA TO THE NONRESTART CONJUGATE GRADIENT ALPHA.
*/
    alpha *= dg/dg1;
/*
* IF NMETH=1 OR A RESTART HAS BEEN PERFORMED, SET ALPHA=1.0.
*/
    if(nrst == 1 || nmeth == 1) alpha = 1.0;
/*
* IF IT IS THE FIRST ITERATION, SET ALPHA=1.0/DSQRT(GSQ),
* WHICH SCALES THE INITIAL SEARCH VECTOR TO UNITY.
*/
    if(rsw) alpha = 1.0/sqrt(gsq);
/*	
* THE LINEAR SEARCH FITS A CUBIC TO F AND DAL, THE FUNCTION AND ITS
* DERIVATIVE AT ALPHA, AND TO FP AND DP,THE FUNCTION
* AND DERIVATIVE AT THE PREVIOUS TRIAL POINT AP.
* INITIALIZE AP ,FP,AND DP.
*/
    ap = 0;
    fp = fmin;
    dp = dg1;
/*
* SAVE THE CURRENT DERIVATIVE TO SCALE THE NEXT SEARCH VECTOR.
*/
    dg = dg1;
/*
* UPDATE THE ITERATION.
*/
    ++*iter;
/*
* CALCULATE THE CURRENT STEPLENGTH  AND STORE THE CURRENT X AND G.
*/
    step = 0;
    for(i = 1; i <= n; ++i)
    {
        step += wr[i]*wr[i];
        nxpi = nx+i;
        ngpi = ng+i;
        wr[nxpi] = x[i];
	wr[ngpi] = g[i];
    }
    step = sqrt(step);
/*
* BEGIN THE LINEAR SEARCH ITERATION.
* TEST FOR FAILURE OF THE LINEAR SEARCH.
*/
 label80:
    if(alpha*step <= acc)
    {
      /*
       * TEST IF DIRECTION IS A GRADIENT DIRECTION.
       */
      if(!rsw) goto label20;
      *nflag = 2;
      return;
    }
/*
* CALCULATE THE TRIAL POINT.
*/
    for(i = 1; i <= n; ++i)
    {
      nxpi = nx+i;
      x[i] = wr[nxpi]+alpha*wr[i];
    }
    if(0)printf("trial x[1] %12.10f %12.10f %12.10f %12.10f\n",
	   x[1],wr[nx+1],alpha,wr[1]);
/*
* EVALUATE THE FUNCTION AT THE TRIAL POINT.
*/
    calcfg(info,x,f,g);
/*
* TEST IF THE MAXIMUM NUMBER OF FUNCTION CALLS HAVE BEEN USED.
*/
    ++*ifun;
    if(*ifun > mxfun)
    {
      *nflag = 1;
      return;
    }
/*
* COMPUTE THE DERIVATIVE OF F AT ALPHA.
*/
    dal = 0;
    for(i = 1; i <= n; ++i)
      dal += g[i]*wr[i];
/*
* TEST WHETHER THE NEW POINT HAS A NEGATIVE SLOPE BUT A HIGHER
* FUNCTION VALUE THAN ALPHA = 0. IF THIS IS THE CASE,THE SEARCH
* HAS PASSED THROUGH A LOCAL MAX AND IS HEADING FOR A DISTANT LOCAL
* MINIMUM.
*/
    if(0) printf("f %f fmin %f dal %f alpha %f dg %f\n",*f,fmin,dal, alpha, dg);
    if(*f > fmin && dal < 0) goto label160;
/*
* IF NOT, TEST WHETHER THE STEPLENGTH CRITERIA HAVE BEEN MET.
*/
    if(!(*f > (fmin+.0001*alpha*dg) || fabs(dal/dg) > .9))
    {
      /*
       * IF THEY HAVE BEEN MET, TEST IF TWO POINTS HAVE BEEN TRIED
       * IF NMETH = 0 AND IF THE TRUE LINE MINIMUM HAS NOT BEEN FOUND.
       */
      if(!((*ifun-ncalls) <= 1 && fabs(dal/dg) > eeps && nmeth == 0))
      {
	goto label170;
      }
    }
/*
* A NEW POINT MUST BE TRIED. USE CUBIC INTERPOLATION TO FIND
* THE TRIAL POINT AT.
*/
    u1 = dp+dal-3.0*(fp-*f)/(ap-alpha);
    u2 = u1*u1-dp*dal;
    if(u2 < 0) u2 = 0;
    else u2 = sqrt(u2);
    at = alpha-(alpha-ap)*(dal+u2-u1)/(dal-dp+2.*u2);
/*
* TEST WHETHER THE LINE MINIMUM HAS BEEN BRACKETED.
*/
    if(0) printf("dal %f dp %f at %f ap %f alpha %f\n",dal,dp, at,ap,alpha);
    if((dal/dp) <= 0)
    {
      /*
       * THE MINIMUM HAS BEEN BRACKETED. TEST WHETHER THE TRIAL POINT LIES
       * SUFFICIENTLY WITHIN THE BRACKETED INTERVAL.
       * IF IT DOES NOT, CHOOSE AT AS THE MIDPOINT OF THE INTERVAL.
       */
      if(at < (1.01*min(alpha,ap)) || at > (.99*max(alpha,ap)))
	at = (alpha+ap)/2.0;
      goto label150;
    }
/*
* THE MINIMUM HAS NOT BEEN BRACKETED. TEST IF BOTH POINTS ARE
* GREATER THAN THE MINIMUM AND THE TRIAL POINT IS SUFFICIENTLY
* SMALLER THAN EITHER.
*/
   if(!(dal > 0.0 && 0.0 < at && at < (.99*min(ap,alpha))))
   {
     /*
      * TEST IF BOTH POINTS ARE LESS THAN THE MINIMUM AND THE TRIAL POINT
      * IS SUFFICIENTLY LARGE.
      */
     if(!(dal <= 0.0 && at > (1.01*max(ap,alpha))))
     {
       /*
	* IF THE TRIAL POINT IS TOO SMALL,DOUBLE THE LARGEST PRIOR POINT.
	*/
       if(dal <= 0) at = 2.0*max(ap,alpha);
       /*
	* IF THE TRIAL POINT IS TOO LARGE, HALVE THE SMALLEST PRIOR POINT.
	*/
       if(dal > 0) at = min(ap,alpha)/2.0;
     }
    }
/*
* SET AP = ALPHA, ALPHA = AT,AND CONTINUE SEARCH.
*/
 label150:
   ap = alpha;
   fp = *f;
   dp = dal;
   alpha = at;
   goto label80;
/*
* A RELATIVE MAX HAS BEEN PASSED.REDUCE ALPHA AND RESTART THE SEARCH.
*/
 label160:
   alpha /= 3;
   ap = 0;
   fp = fmin;
   dp = dg;
   goto label80;
/*
* THE LINE SEARCH HAS CONVERGED. TEST FOR CONVERGENCE OF THE ALGORITHM.
*/
label170:
   gsq = 0.0;
   xsq = 0.0;
   for(i = 1; i <= n; ++i)
   {
     gsq += g[i]*g[i];
     xsq += x[i]*x[i];
   }
   ++cntxx;
   if(0) printf("gsq.c %10.6f %d\n",gsq, cntxx);
   if(gsq <= eeps*eeps*max(1.0,xsq)) return;
/*
* SEARCH CONTINUES. SET WR[i] = ALPHA*WR[i],THE FULL STEP VECTOR.
*/
   for(i = 1; i <= n; ++i)
     wr[i] *= alpha;
/*
* COMPUTE THE NEW SEARCH VECTOR. FIRST TEST WHETHER A
* CONJUGATE GRADIENT OR A VARIABLE METRIC VECTOR IS USED.
*/
   if(nmeth == 1) goto label330;
/*
* CONJUGATE GRADIENT UPDATE SECTION.
* TEST IF A POWELL RESTART IS INDICATED.
*/
   rtst = 0;
   for(i = 1; i <= n; ++i)
   {
     ngpi = ng+i;
     rtst += g[i]*wr[ngpi];
   }
   if(fabs(rtst/gsq) > 0.2) nrst = n;
/*
* IF A RESTART IS INDICATED, SAVE THE CURRENT D AND Y
* AS THE BEALE RESTART VECTORS AND SAVE D'Y AND Y'Y
* IN WR(NCONS+1) AND WR(NCONS+2).
*/
   if(nrst == n)
   {
     wr[ncons+1] = 0;
     wr[ncons+2] = 0;
     for(i = 1; i <= n; ++i)
     {
       nrdpi = nrd+i;
       nrypi = nry+i;
       ngpi = ng+i;
       wr[nrypi] = g[i]-wr[ngpi];
       wr[nrdpi] = wr[i];
       wr[ncons1] += wr[nrypi]*wr[nrypi];
       wr[ncons2] += wr[i]*wr[nrypi];
     }
   }
/*
* CALCULATE  THE RESTART HESSIAN TIMES THE CURRENT GRADIENT.
*/
   u1 = 0.0;
   u2 = 0.0;
   for(i = 1; i <= n; ++i)
   {
     nrdpi = nrd+i;
     nrypi = nry+i;
     u1 -= wr[nrdpi]*g[i]/wr[ncons1];
     u2 += wr[nrdpi]*g[i]*2./wr[ncons2]-wr[nrypi]*g[i]/wr[ncons1];
   }
   u3 = wr[ncons2]/wr[ncons1];
   for(i = 1; i <= n; ++i)
   {
     nxpi = nx+i;
     nrdpi = nrd+i;
     nrypi = nry+i;
     wr[nxpi] = -u3*g[i]-u1*wr[nrypi]-u2*wr[nrdpi];
   }
/*
* IF THIS IS A RESTART ITERATION,WR[NX+I] CONTAINS THE NEW SEARCH
* VECTOR.
*/
   if(nrst == n) goto label300;
/*
* NOT A RESTART ITERATION. CALCULATE THE RESTART HESSIAN
* TIMES THE CURRENT Y.
*/
label250:
   u1 = 0.;
   u2 = 0.;
   u3 = 0.;
   u4 = 0.;
   for(i = 1; i <= n; ++i)
   {
     ngpi = ng+i;
     nrdpi = nrd+i;
     nrypi = nry+i;
     u1 -= (g[i]-wr[ngpi])*wr[nrdpi]/wr[ncons1];
     u2 += -(g[i]-wr[ngpi])*wr[nrypi]/wr[ncons1]
       + 2.0*wr[nrdpi]*(g[i]-wr[ngpi])/wr[ncons2];
     u3 += wr[i]*(g[i]-wr[ngpi]);
   }
   step = 0;
   for(i = 1; i <= n; ++i)
   {
     ngpi = ng+i;
     nrdpi = nrd+i;
     nrypi = nry+i;
     step = (wr[ncons2]/wr[ncons1])*(g[i]-wr[ngpi])
       +u1*wr[nrypi]+u2*wr[nrdpi];
     u4 += step*(g[i]-wr[ngpi]);
     wr[ngpi] = step;
   }
/*
* CALCULATE THE DOUBLY UPDATED HESSIAN TIMES THE CURRENT
* GRADIENT TO OBTAIN THE SEARCH VECTOR.
*/
   u1 = 0.0;
   u2 = 0.0;
   for(i = 1; i <= n; ++i)
   {
     u1 -= wr[i]*g[i]/u3;
     ngpi = ng+i;
     u2 += (1.0+u4/u3)*wr[i]*g[i]/u3-wr[ngpi]*g[i]/u3;
   }
   for(i = 1; i <= n; ++i)
   {
     ngpi = ng+i;
     nxpi = nx+i;
     wr[nxpi] -= u1*wr[ngpi] + u2*wr[i];
   }
/*
* CALCULATE THE DERIVATIVE ALONG THE NEW SEARCH VECTOR.
*/
label300:
   dg1 = 0;
   for(i = 1; i <= n; ++i)
   {
     nxpi = nx+i;
     wr[i] = wr[nxpi];
     dg1 += wr[i]*g[i];
   }
/*
* IF THE NEW DIRECTION IS NOT A DESCENT DIRECTION,STOP.
*/
   if (dg1 <= 0)
   {
     /*
      * UPDATE NRST TO ASSURE AT LEAST ONE RESTART EVERY N ITERATIONS.
      */
     if(nrst == n) nrst = 0;
     ++nrst;
     rsw = FALSE;
     goto label40;
   }
/*
* ROUNDOFF HAS PRODUCED A BAD DIRECTION.
*/
label320:
   *nflag = 3;
   return;
/*
* A VARIABLE METRIC ALGORITM IS BEING USED. CALCULATE Y AND D'Y.
*/
label330:
   u1 = 0.0;
   for(i = 1; i <= n; ++i)
   {
     ngpi = ng+i;
     wr[ngpi] = g[i]-wr[ngpi];
     u1 += wr[i]*wr[ngpi];
   }
/*
* IF RSW = .TRUE.,SET UP THE INITIAL SCALED APPROXIMATE HESSIAN.
*/
   if(!rsw) goto label380;
/*
* CALCULATE Y'Y.
*/
   u2 = 0;
   for(i = 1; i <= n; ++i)
   {
     ngpi = ng+i;
     u2 += wr[ngpi]*wr[ngpi];
   }
/*
* CALCULATE THE INITIAL HESSIAN AS H = (P'Y/Y'Y)*I
* AND THE INITIAL U2 = Y'HY AND WR(NX+I) = HY.
*/
   ij = 1;
   u3 = u1/u2;
   for(i = 1; i <= n; ++i)
   {
     for(j = i; j <= n; ++j)
     {
       ncons1 = ncons+ij;
       /* WRITE(6,*) 'NCONS1 =  ',NCONS1 */
       wr[ncons1] = 0.0;
       if(i == j) wr[ncons1] = u3;
       ++ij;
     }
     nxpi = nx+i;
     ngpi = ng+i;
     wr[nxpi] *= u3;
   }
   u2 *= u3;
   goto label430;
/*
* CALCULATE WR(NX+I) = HY AND U2 = Y'HY.
*/
label380:
   u2 = 0.0;
   for(i = 1; i <= n; ++i)
   {
     u3 = 0.0;
     ij = i;
     if(i != 1)
     {
       int ii = i-1;
       for(j = 1; j <= ii; ++j)
       {
	 ngpj = ng+j;
	 ncons1 = ncons+ij;
	 u3 += wr[ncons1]*wr[ngpj];
	 ij += n-j;
       }
     }
     for(j = 1; j <= n; ++j)
     {
       ncons1 = ncons+ij;
       ngpj = ng+j;
       u3 += wr[ncons1]*wr[ngpj];
       ++ij;
     }
     ngpi = ng+i;
     u2 += u3*wr[ngpi];
     nxpi = nx+i;
     wr[nxpi] = u3;
   }
/*
* CALCULATE THE UPDATED APPROXIMATE HESSIAN.
*/
label430:
   u4 = 1.0+u2/u1;
   for(i = 1; i <= n; ++i)
   {
     nxpi = nx+i;
     ngpi = ng+i;
     wr[ngpi] = u4*wr[i]-wr[nxpi];
   }
   ij = 1;
   for(i = 1; i <= n; ++i)
   {
     nxpi = nx+i;
     u3 = wr[i]/u1;
     u4 = wr[nxpi]/u1;
     for(j = i; j <= n; ++j)
     {
       ncons1 = ncons+ij;
       ngpj = ng+j;
       wr[ncons1] = wr[ncons1]+u3*wr[ngpj]-u4*wr[j];
       ++ij;
     }
   }
/*
* CALCULATE THE NEW SEARCH DIRECTION WR[i] = -HG AND ITS DERIVATIVE.
*/
   dg1 = 0.0;
   for(i = 1; i <= n; ++i)
   {
     u3 = 0.0;
     ij = i;
     if(i != 1)
     {
       ii = i-1;
       for(j = 1; j <= ii; ++j)
       {
	 ncons1 = ncons+ij;
	 u3 -= wr[ncons1]*g[j];
	 ij += n-j;
       }
     }
     for(j = i; j <= n; ++j)
     {
       ncons1 = ncons+ij;
       u3 -= wr[ncons1]*g[j];
       ++ij;
     }
     dg1 += u3*g[i];
     wr[i] = u3;
   }
/*
* TEST FOR A DOWNHILL DIRECTION.
*/
   if(dg1 > 0) goto label320;
   rsw = FALSE;
   goto label40;
}
