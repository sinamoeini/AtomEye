/********************************************************************/
/* Ewald V2.0:                                                      */
/*                                                                  */
/* Calculate the total energy, force, stress tensor and force       */
/* constant matrix (in k-space) of general 3D triclinic PBC rigid   */
/* ion system with zero total charge. The result is the limit of    */
/* (infinitely) many replicas stacked in a nearly spherical cluster */
/* and surrounded by conductive medium outside, thus the cluster    */
/* has no surface charge distribution and there is no uniform       */
/* macroscopic electric field; it is the only self-consistent       */
/* set-up for both static and MD simulations in PBC. For details    */
/* see reference "Doc/moldy.ps", and websites                       */
/*                                                                  */
/*   http://www.ee.duke.edu/~ayt/ewaldpaper/                        */
/*   http://www.keele.ac.uk/depts/ph/tc/cph_res/dynamo.html         */
/*   http://www.earth.ox.ac.uk/%7Ekeith/moldy.html                  */
/*   http://www.fos.su.se/physical/sasha/md_prog.html               */
/*                                                                  */
/*                                                Aug. 26, 1998     */
/*                                    Developed by Ju Li at MIT     */
/********************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI acos(-1.0)

/* number of particles */
static int NP;
/* real space positions and charges */
static double *x,*y,*z,*q;
/* Ewald summation parameters */
static double my_accuracy,alpha,rcut,kcut; 
/* H matrix */
static double h11,h12,h13,h21,h22,h23,h31,h32,h33; 
/* G matrix: H*G^T = 2\pi I */
static double g11,g12,g13,g21,g22,g23,g31,g32,g33,volume;
/* Distance from origin to three unit    */
/* cell surfaces and reciprocal surfaces */
static double rd1,rd2,rd3,kd1,kd2,kd3;
/* real space lattices to be be summed */
static int max_n1,max_n2,max_n3,max_xvec,num_xvec,*nvec;
static double *xvec;
/* reciprocal lattices to be be summed */
static int max_g1,max_g2,max_g3,max_gvec,num_gvec,*gvec; 
#ifdef _TABULATE_ERFC
static double derr,*err,*cli;
#endif
#ifdef _TABULATE_SINE
static double *sinn,*cosn;
#endif

/* Our convention is that (h11,h12,h13) (in f77 index)   */
/* is the (x,y,z) component of the first cell edge, etc. */ 
void init_cell (double H[3][3])
{ 
    /* Because H is assumed to be coming from Fortran */
    /* we have to do a little conversion */
    h11 = H[0][0];  /* 1st element */
    h12 = H[1][0];  /* 4th element */
    h13 = H[2][0];  /* 7th element */
    h21 = H[0][1];  /* 2nd element */
    h22 = H[1][1];  /* 5th element */
    h23 = H[2][1];  /* 8th element */
    h31 = H[0][2];  /* 3rd element */
    h32 = H[1][2];  /* 6th element */
    h33 = H[2][2];  /* 9th element */
    /* Get the reciprocal lattice vectors */
    g11 = h22*h33 - h23*h32;
    g22 = h33*h11 - h31*h13;
    g33 = h11*h22 - h12*h21;
    g12 = h23*h31 - h21*h33;
    g23 = h31*h12 - h32*h11;
    g31 = h12*h23 - h13*h22;
    g13 = h21*h32 - h31*h22;
    g21 = h32*h13 - h12*h33;
    g32 = h13*h21 - h23*h11;
    volume = h11*g11 + h12*g12 + h13*g13;
    /* the shortest distance to respective surfaces */
    rd1 = fabs(volume)/sqrt(g11*g11+g12*g12+g13*g13);
    rd2 = fabs(volume)/sqrt(g21*g21+g22*g22+g23*g23);
    rd3 = fabs(volume)/sqrt(g31*g31+g32*g32+g33*g33);
    /* reciprocal lattice vectors */
    g11 *= 2*PI/volume;
    g12 *= 2*PI/volume;
    g13 *= 2*PI/volume;
    g21 *= 2*PI/volume;
    g22 *= 2*PI/volume;
    g23 *= 2*PI/volume;
    g31 *= 2*PI/volume;
    g32 *= 2*PI/volume;
    g33 *= 2*PI/volume;
    /* the shortest distance to respective reciprocal surfaces */
    kd1 = 2*PI/sqrt(h11*h11+h12*h12+h13*h13);
    kd2 = 2*PI/sqrt(h21*h21+h22*h22+h23*h23);
    kd3 = 2*PI/sqrt(h31*h31+h32*h32+h33*h33);
    volume = fabs(volume);
    return;
}  /* end init_cell() */


/* the time ratio of evaluating a single term */
/* in the real and reciprocal space summation */
#define TRTF 5.5
/* empirical constants of Ju Li to the Fincham */
/* formulas to achieve the acclaimed accuracy  */
#define RCUT_COEFF 1.2
#define KCUT_COEFF RCUT_COEFF

void init_ewald (int *number_particles,
		 double H[3][3], double *accuracy)
{
    int i,j,k,max_err;
    double xx,yy,zz,r2,maxr;

    NP = *number_particles;
    my_accuracy = *accuracy;
    init_cell(H);
    
    /* Set the parameters alpha and cutoffs based on
       the formula by Fincham CCP5 38, p17 (1993) */ 
    alpha = pow(NP*PI*PI*PI*TRTF/volume/volume,1./6.);
    rcut = sqrt(-log(my_accuracy))/alpha*RCUT_COEFF;
    kcut = 2*alpha*sqrt(-log(my_accuracy))*KCUT_COEFF;

    /* calculate the maximum separation between */
    /* two particles inside the unit cell: must */
    /* be one of the eight corners. */ 
    maxr = 0.;
    for (i=-1; i<=1; i+=2)
	for (j=-1; j<=1; j+=2)
	    for (k=-1; k<=1; k+=2)
	    {
		xx = i*h11 + j*h21 + k*h31;
		yy = i*h12 + j*h22 + k*h32;
		zz = i*h13 + j*h23 + k*h33;
		r2 = xx*xx + yy*yy + zz*zz;
		if (r2>maxr*maxr) maxr=sqrt(r2);
	    }
    
    /* construct a list of important real-space cells */
    maxr += rcut;
    max_xvec = ceil(4.*PI*maxr*maxr*maxr/3./volume*1.2/16)*16;
    nvec = (int *) malloc(3*max_xvec*sizeof(int));
    xvec = (double *) malloc(3*max_xvec*sizeof(double));
    max_n1 = ceil(rcut/rd1);
    max_n2 = ceil(rcut/rd2);
    max_n3 = ceil(rcut/rd3);
    /* first record is the bare unit cell */
    num_xvec = 3;
    xvec[0] = 0.;
    xvec[1] = 0.;
    xvec[2] = 0.;
    for (i=-max_n1; i<=max_n1; i++)
	for (j=-max_n2; j<=max_n2; j++)
	    for (k=-max_n3; k<=max_n3; k++)
		if (!((i==0)&&(j==0)&&(k==0)))
		{
		    xx = i*h11 + j*h21 + k*h31;
		    yy = i*h12 + j*h22 + k*h32;
		    zz = i*h13 + j*h23 + k*h33;
		    r2 = xx*xx + yy*yy + zz*zz;
		    /* only these cells are possible to */
		    /* have an interaction within rcut  */
		    if (r2<maxr*maxr)
		    {
			num_xvec += 3;
			if (num_xvec >= 3*max_xvec)
			{
			    printf ("init_ewald(): max_xvec reached.\n");
			    exit(1);
			}
			nvec[num_xvec-3] = i;
			nvec[num_xvec-2] = j;
			nvec[num_xvec-1] = k;
		    }
		}

    /* construct a list of necessary reciprocal k-points */
    max_gvec = ceil(2.*PI*kcut*kcut*kcut/3./
		    (8*PI*PI*PI/volume)*1.2/16)*16;
    gvec = (int *) malloc(3*max_gvec*sizeof(int));
    max_g1 = ceil(kcut/kd1);
    max_g2 = ceil(kcut/kd2);
    max_g3 = ceil(kcut/kd3);
    /* first record is G=0, which has no inversion image */
    num_gvec = 3;
    gvec[0] = 0;
    gvec[1] = 0;
    gvec[2] = 0;
    /* There are inversion symmetry in energy,  */
    /* force, stress, but not in force constant */
    /* matrix calculations in a general system  */ 
    for (k=0; k<=max_g3; k++)
	for (j=-max_g2; j<=max_g2; j++)
	    for (i=-max_g1; i<=max_g1; i++)
		if ((k>0)||(j>0)||((j==0)&&(i>0)))
		{
		    xx = i*g11 + j*g21 + k*g31;
		    yy = i*g12 + j*g22 + k*g32;
		    zz = i*g13 + j*g23 + k*g33;
		    r2 = xx*xx + yy*yy + zz*zz;
		    if (r2 < kcut*kcut)
		    {
			num_gvec += 3;
			if (num_gvec >= 3*max_gvec)
			{
			    printf ("init_ewald(): max_gvec reached.\n");
			    exit(1);
			}
			gvec[num_gvec-3] = i;
			gvec[num_gvec-2] = j;
			gvec[num_gvec-1] = k;
		    }
		}
    
    /* allocate real space positions and charges */
    x = (double *) malloc(NP*sizeof(double));
    y = (double *) malloc(NP*sizeof(double));
    z = (double *) malloc(NP*sizeof(double));
    q = (double *) malloc(NP*sizeof(double));

#ifdef _TABULATE_ERFC
    /* tabulate the error function and its derivative:     */
    /* because we will do expansion to the 2nd order, the  */
    /* interval should be proportional to (accuracy)^(1/3) */
    derr = pow(*accuracy,1./3.)/2./_TABULATE_ERFC;
    /* default _TABULATE_ERFC = 1. */
    max_err = ceil(alpha*rcut/derr);
    err = (double *)malloc(max_err*sizeof(double));
    cli = (double *)malloc(max_err*sizeof(double));
    for (i=0; i<max_err; i++)
    {
	err[i] = erfc(i*derr);
	cli[i] = 2./sqrt(PI)*exp(-i*i*derr*derr);
    }
#endif

#ifdef _TABULATE_SINE
    /* massive tabulation of sine and cosine values */
    sinn = (double *) malloc(NP*(2*max_g1+1)*(2*max_g2+1)
			     *(2*max_g3+1)*sizeof(double));
    cosn = (double *) malloc(NP*(2*max_g1+1)*(2*max_g2+1)
			     *(2*max_g3+1)*sizeof(double));
#endif

#ifdef _PRINT_EWALD
  printf ("\nAlpha = %f  Rcut = %f  Kcut = %f\n",
	  alpha, rcut, kcut);
  printf ("Max_n1 = %d  Max_n2 = %d  Max_n3 = %d => %d L-points\n",
	  max_n1, max_n2, max_n3, num_xvec/3);
  printf ("Max_g1 = %d  Max_g2 = %d  Max_g3 = %d => %d G-points\n\n",
	  max_g1, max_g2, max_g3, num_gvec/3); 
#endif
  
  return;
} /* end init_ewald() */


void exit_ewald()
{
    free(gvec); free(xvec); free(nvec);
    free(x); free(y); free(z); free(q); 
#ifdef _TABULATE_ERFC 
    free(err); free(cli);
#endif
#ifdef _TABULATE_SINE 
    free(sinn); free(cosn); 
#endif
    return;
} /* end exit_ewald() */


void ewald (charge, s1, s2, s3, H, pote, fx, fy, fz, stress)
/* charge of particles from 1 to NP */
double *charge;
/* reduced coordinates of particles in [-0.5, 0.5] */
double *s1, *s2, *s3;
/* H matrix of the cell, assumed to be passed from fortran */
/* code (column order): in there H(1,1),H(1,2),H(1,3) are  */
/* the xyz components of the first edge, etc. */
double H[3][3]; 
/* Coulomb energy per cell of the system */
double *pote;
/* Coulomb force on particles */
double *fx, *fy, *fz;
/* stress tensor due to Coulomb interactions */
double stress[3][3]; 
{
    double dx,dy,dz,qsum=0.;
    double rx,ry,rz,r2,r,product,xx,margin,factor,gactor;
    double qij,ff;
    double kx,ky,kz,k2,cc,cossum,sinsum,ak;
    double sinx,cosx,siny,cosy,sinz,cosz;
    int i,j,k,l,n;
    int i000,ij00,ijk0,ijkl,jkl;

    /* make it work for Parrinello-Rahman MD */
    init_cell(H);
    alpha = pow(NP*PI*PI*PI*TRTF/volume/volume,1./6.);
    rcut = sqrt(-log(my_accuracy))/alpha*RCUT_COEFF;
    kcut = 2*alpha*sqrt(-log(my_accuracy))*KCUT_COEFF;
    /* the formulas have nice scaling with volume, so */
    /* generated L and G lattices need not be revised */
    /* unless large shape changes happen during MD:   */
    for (n=0; n<num_xvec; n+=3)
    {
	xvec[n]   = nvec[n]*h11 + nvec[n+1]*h21 + nvec[n+2]*h31;
	xvec[n+1] = nvec[n]*h12 + nvec[n+1]*h22 + nvec[n+2]*h32;
	xvec[n+2] = nvec[n]*h13 + nvec[n+1]*h23 + nvec[n+2]*h33;
    }
    /* even the erfc() table has the same range as before */
    
    for (i=0; i<NP; i++)
    {
	fx[i] = 0.;
	fy[i] = 0.;
	fz[i] = 0.;
	while ((s1[i]>0.5)||(s1[i]<-0.5)||
	       (s2[i]>0.5)||(s2[i]<-0.5)||
	       (s3[i]>0.5)||(s3[i]<-0.5))
	{
	    printf("ewald(): reduced coordinates > 0.5.\n");
	    printf("may lead to inaccurate results.\n");
	    /* exit(1); */
	    if (s1[i]<-0.5) s1[i]++;
	    if (s1[i]>0.5)  s1[i]--;
	    if (s2[i]<-0.5) s2[i]++;
	    if (s2[i]>0.5)  s2[i]--;
	    if (s3[i]<-0.5) s3[i]++;
	    if (s3[i]>0.5)  s3[i]--;
	    printf("reduced coordinates modified.\n");
	}
	x[i] = s1[i]*h11 + s2[i]*h21 + s3[i]*h31;
	y[i] = s1[i]*h12 + s2[i]*h22 + s3[i]*h32;
	z[i] = s1[i]*h13 + s2[i]*h23 + s3[i]*h33;
	qsum += charge[i];
    }
    
    if (fabs(qsum/NP) > 0.0001)
	printf("\nwarning from ewald():"\
	       "significant net charge in the system.\n");
    
    for (i=0; i<3; i++) 
	for (j=0; j<3; j++) 
	    stress[i][j] = 0.;
  
    *pote = 0.;
    for (i=0; i<NP; i++)
    {
	q[i] = charge[i] - qsum/NP;
	/* the self-energy */
	*pote -= alpha*q[i]*q[i]/sqrt(PI);
    } 

    /* Do the real space summation */
    for (i=0; i<NP; i++)
	for (j=i; j<NP; j++)
	{
	    dx = x[i] - x[j];
	    dy = y[i] - y[j];
	    dz = z[i] - z[j];
	    qij = q[i] * q[j];
	    for (n=0; n<num_xvec; n+=3)
	    {
		rx = dx + xvec[n];
		ry = dy + xvec[n+1];
		rz = dz + xvec[n+2];
		r2 = rx*rx + ry*ry + rz*rz;
		if ((r2>0.)&&(r2<rcut*rcut))
		{
		    r = sqrt(r2);
		    product = alpha * r;
#ifdef _TABULATE_ERFC
		    k = floor(product/derr);
		    xx = k * derr;
		    margin = product - xx;
		    gactor = err[k] - cli[k] * margin
			* (1-xx*margin);
#else
		    gactor = erfc(product);
#endif
		    /* energy */
		    *pote += (i==j)? gactor*qij/r*0.5:
			gactor*qij/r;
		    
		    /* force */
#ifdef _TABULATE_ERFC
		    margin = product*product - xx*xx;
		    factor = gactor + product * cli[k] *
			(1. - margin * (1 - margin * 0.5));
#else
		    factor = gactor + 2 * product *
			exp(-product*product) / sqrt(PI);
#endif
		    ff = factor * qij / r2 / r;
		    fx[i] += ff * rx;
		    fy[i] += ff * ry;
		    fz[i] += ff * rz; 
		    fx[j] -= ff * rx;
		    fy[j] -= ff * ry;
		    fz[j] -= ff * rz;
		    
		    /* stress */
		    stress[0][0] += (i==j)?ff*rx*rx*0.5:ff*rx*rx;
		    stress[1][1] += (i==j)?ff*ry*ry*0.5:ff*ry*ry;
		    stress[2][2] += (i==j)?ff*rz*rz*0.5:ff*rz*rz;
		    stress[0][1] += (i==j)?ff*rx*ry*0.5:ff*rx*ry;
		    stress[0][2] += (i==j)?ff*rx*rz*0.5:ff*rx*rz;
		    stress[1][2] += (i==j)?ff*ry*rz*0.5:ff*ry*rz;
		}
	    }
	}

    /* Do the reciprocal space summation */

#ifdef _TABULATE_SINE
    /* bootstrap the sine and cosine arrays */
    for (i=0; i<NP; i++)
    {
	sinx = sin(s1[i]*2*PI);
	cosx = cos(s1[i]*2*PI);
	siny = sin(s2[i]*2*PI);
	cosy = cos(s2[i]*2*PI);
	sinz = sin(s3[i]*2*PI);
	cosz = cos(s3[i]*2*PI);

	/* corner */
	i000 = i*(2*max_g1+1)*(2*max_g2+1)*(2*max_g3+1);
	sinn[i000] = sin(2*PI*(-max_g1*s1[i]
			       -max_g2*s2[i]-max_g3*s3[i]));
	cosn[i000] = cos(2*PI*(-max_g1*s1[i]
			       -max_g2*s2[i]-max_g3*s3[i]));

	/* corner -> line */
	for (j=1; j<=2*max_g1; j++)
	{
	    ij00 = i000 + j*(2*max_g2+1)*(2*max_g3+1);
	    sinn[ij00] =
		sinn[ij00-(2*max_g2+1)*(2*max_g3+1)]*cosx + 
		cosn[ij00-(2*max_g2+1)*(2*max_g3+1)]*sinx;
	    cosn[ij00] =
		cosn[ij00-(2*max_g2+1)*(2*max_g3+1)]*cosx -
		sinn[ij00-(2*max_g2+1)*(2*max_g3+1)]*sinx;
	}
	
	/* line -> surface */
	for (j=0; j<=2*max_g1; j++)
	{
	    ij00 = i000 + j*(2*max_g2+1)*(2*max_g3+1);
	    for (k=1; k<=2*max_g2; k++)
	    {
		ijk0 = ij00 + k*(2*max_g3+1);
		sinn[ijk0] =
		    sinn[ijk0-(2*max_g3+1)]*cosy +
		    cosn[ijk0-(2*max_g3+1)]*siny;
		cosn[ijk0] =
		    cosn[ijk0-(2*max_g3+1)]*cosy -
		    sinn[ijk0-(2*max_g3+1)]*siny;
	    }
	}

	/* surface -> cube */
	for (j=0; j<=2*max_g1; j++)
	{
	    ij00 = i000 + j*(2*max_g2+1)*(2*max_g3+1);
	    for (k=0; k<=2*max_g2; k++)
	    {
		ijk0 = ij00 + k*(2*max_g3+1);
		for (l=1; l<=2*max_g3; l++)
		{
		    ijkl = ijk0 + l;
		    sinn[ijkl] =
			sinn[ijkl-1]*cosz +
			cosn[ijkl-1]*sinz;
		    cosn[ijkl] =
			cosn[ijkl-1]*cosz -
			sinn[ijkl-1]*sinz;
		}
	    }
	}
    }  
#endif
    
    /* in the summation, omit G = 0 */
    for (n=3; n<num_gvec; n+=3)
    {
	cossum = 0.;
	sinsum = 0.;
	
#ifdef _TABULATE_SINE
	jkl = ((gvec[n]+max_g1)*(2*max_g2+1)+ 
	        gvec[n+1]+max_g2)*(2*max_g3+1) +
	        gvec[n+2]+max_g3;
#endif

	for (i=0; i<NP; i++)
	{
#ifdef _TABULATE_SINE
	    ijkl = i*(2*max_g1+1)*(2*max_g2+1)
		*(2*max_g3+1) + jkl;
	    cossum += q[i]*cosn[ijkl];
	    sinsum += q[i]*sinn[ijkl];
#else
	    cossum += q[i]*cos(2*PI*(gvec[n]*s1[i] + 
				     gvec[n+1]*s2[i] +
				     gvec[n+2]*s3[i]));
	    sinsum += q[i]*sin(2*PI*(gvec[n]*s1[i] +
				     gvec[n+1]*s2[i] +
				     gvec[n+2]*s3[i]));
#endif
	}

	kx = gvec[n]*g11+gvec[n+1]*g21+gvec[n+2]*g31;
	ky = gvec[n]*g12+gvec[n+1]*g22+gvec[n+2]*g32;
	kz = gvec[n]*g13+gvec[n+1]*g23+gvec[n+2]*g33;
	k2 = kx*kx + ky*ky + kz*kz;
	/* with inversion symmetry, each G represents two */
	ak = 4.*PI*exp(-k2/alpha/alpha/4.)/k2/volume;

	/* energy */
	cc = ak*(cossum*cossum+sinsum*sinsum);
	*pote += cc;

	/* stress */
	stress[0][0] += cc*(1.-(2./k2+1./alpha/alpha*0.5)*kx*kx);
	stress[1][1] += cc*(1.-(2./k2+1./alpha/alpha*0.5)*ky*ky);
	stress[2][2] += cc*(1.-(2./k2+1./alpha/alpha*0.5)*kz*kz);
	stress[0][1] -= cc*(2./k2+1./alpha/alpha*0.5)*kx*ky;
	stress[0][2] -= cc*(2./k2+1./alpha/alpha*0.5)*kx*kz;
	stress[1][2] -= cc*(2./k2+1./alpha/alpha*0.5)*ky*kz;

	/* force */
	for (i=0; i<NP; i++)
	{
	    
#ifdef _TABULATE_SINE             
	    ijkl = i*(2*max_g1+1)*(2*max_g2+1)
		*(2*max_g3+1) + jkl;
	    cc = 2.*ak*q[i]*(cossum*sinn[ijkl]-
			     sinsum*cosn[ijkl]);
#else
	    cc = 2*ak*q[i]*(cossum*sin(2*PI*(gvec[n]*s1[i] +
					     gvec[n+1]*s2[i] +
					     gvec[n+2]*s3[i])) -
			    sinsum*cos(2*PI*(gvec[n]*s1[i] +
					     gvec[n+1]*s2[i] +
					     gvec[n+2]*s3[i])));
#endif
	    fx[i] += kx*cc;
	    fy[i] += ky*cc;
	    fz[i] += kz*cc;
	}
    }

    stress[1][0] = stress[0][1];
    stress[2][0] = stress[0][2];
    stress[2][1] = stress[1][2];
    for (i=0; i<3; i++) 
	for (j=0; j<3; j++) 
	    stress[i][j] /= volume;
    
    return;    
} /* end ewald() */


/****************************************************/
/* see <Theory of Lattice Dynamics> ed. Maradudin,  */
/* Solid State Phys. Academic Press, 1978, Page 209 */
/****************************************************/

/* last term of (6.2.47) with removal of divergence
   in Page 210: writes into a 3x3 complex matrix H */   
void calc_ewald_dynamical_H  
(double fk[3], int m, int n, double H[3][6])
{ 
    int i,j,l;
    double rx[3],r2,r,product;
    double xx,margin,factor,gactor;
    double phase,a1,a2,cc;
    for (i=0; i<3; i++)
	for (j=0; j<6; j++) 
	    H[i][j] = 0.;
    /* m stands for \kappa, n stands for \kappa^\prime  */
    /* l is actually l^\prime in (6.2.47), the real l=0 */
    for (l=0; l<num_xvec; l+=3)
    {
	rx[0] = x[m] - x[n] + xvec[l];
	rx[1] = y[m] - y[n] + xvec[l+1];
	rx[2] = z[m] - z[n] + xvec[l+2];
	phase = -fk[0]*rx[0]-fk[1]*rx[1]-fk[2]*rx[2];
	
	if ((l==0)&&(m==n))
	{ 
	    /* C convention is row order */ 
	    H[0][0] += 4./3./sqrt(PI);
	    H[1][2] += 4./3./sqrt(PI);
	    H[2][4] += 4./3./sqrt(PI);
	}
	else
	{
	    r2 = rx[0]*rx[0]+rx[1]*rx[1]+rx[2]*rx[2];
	    if (r2 < rcut*rcut)
	    {
		r = sqrt(r2);
		/* P corresponds to our alpha^2 */
		product = alpha * r; 
		a1 = 3./(product*product*product)
		    *erfc(product) + 2/sqrt(PI)*
		    (3./(product*product)+2.)
		    *exp(-product*product);
		a2 = 1./(product*product*product)
		    *erfc(product) + 2/sqrt(PI)*
		    exp(-product*product)/
		    (product*product);
		for (i=0; i<3; i++)
		    for (j=0; j<3; j++)
		    {
			cc = rx[i]*rx[j]/r2*a1;
			if (i==j) cc -= a2;
			H[i][2*j] += cc*cos(phase);
			H[i][2*j+1] += cc*sin(phase);
		    }
	    }
	}
    }
    for (i=0; i<3; i++)
	for (j=0; j<6; j++) 
	    H[i][j] *= alpha*alpha*alpha;
    return;
} /* end calc_ewald_dynamical_H() */


/* short-ranged response */
void calc_ewald_dynamical_Q
(double fk[3], int m, int n, double Q[3][6])
{
    int i,j,l,p;
    double fk2 = fk[0]*fk[0]+fk[1]*fk[1]+fk[2]*fk[2];
    double tk[3],tk2,phase,cc,sincc,coscc;
    double H[3][6];

    /* first term of (6.2.47) */
    for (i=0; i<3; i++)
	for (j=0; j<3; j++) 
	{
	    Q[i][j*2] = -4.*PI/volume*fk[i]*fk[j]/
		fk2*(exp(-fk2/(alpha*alpha)/4.)-1.);
	    Q[i][j*2+1] = 0.;
	}
    /* is regular at k=0 */

    /* the third term of (6.2.47) */
    calc_ewald_dynamical_H (fk, m, n, H);
    for (i=0; i<3; i++)
	for (j=0; j<6; j++) 
	    Q[i][j] += H[i][j];

    /* the second term of (6.2.47), with G=0 */
    /* excluded, and symmetry cannot be used */
    for (l=3; l<num_gvec; l+=3)
	for (p=-1;p<=1;p+=2)
	{
	    tk[0] = fk[0] + p * 
		(gvec[l]*g11+gvec[l+1]*g21+gvec[l+2]*g31);
	    tk[1] = fk[1] + p * 
		(gvec[l]*g12+gvec[l+1]*g22+gvec[l+2]*g32);
	    tk[2] = fk[2] + p * 
		(gvec[l]*g13+gvec[l+1]*g23+gvec[l+2]*g33);
	    tk2 = tk[0]*tk[0] + tk[1]*tk[1] + tk[2]*tk[2];
	    phase = p * (gvec[l]  *(g11*(x[m]-x[n])+
				    g12*(y[m]-y[n])+
				    g13*(z[m]-z[n])) +
			 gvec[l+1]*(g21*(x[m]-x[n])+
				    g22*(y[m]-y[n])+
				    g23*(z[m]-z[n])) +
			 gvec[l+2]*(g31*(x[m]-x[n])+
				    g32*(y[m]-y[n])+
				    g33*(z[m]-z[n])));
	    cc = -4.*PI/volume/tk2*exp(-tk2/alpha/alpha/4.);
	    sincc = sin(phase) * cc;
	    coscc = cos(phase) * cc;
	    for (i=0; i<3; i++)
		for (j=0; j<3; j++)
		{
		    Q[i][2*j]   += coscc * tk[i] * tk[j];
		    Q[i][2*j+1] += sincc * tk[i] * tk[j];
		}
	}
    
    /* multiply by respective charges: */
    for (i=0; i<3; i++)
	for (j=0; j<6; j++) 
	    Q[i][j] *= q[m] * q[n];
    
    return;
} /* end calc_ewald_dynamical_Q() */


#define TINY 1e-12
/*******************************************************/
/* Force constant matrix in k-space; must call ewald() */
/* first. fk[3] is (kx,ky,kz) in Cartesian frame. Phi  */
/* should be declared "double complex phi(L,*)" in     */
/* Fortran with L>=3*NP.                               */
/*******************************************************/
void ewald_dynamical(double fk[3], double *phi, int *L)
{
    int m,n,p,i,j,l; 
    /* particle indices [0,NP-1], standing for kappa,   */
    /* kappa prime, and kappa double prime respectively */
    double zero_k[3] = {TINY,0.,0.};
    double Q[3][6], total[3][6];
    double fk2 = fk[0]*fk[0] + fk[1]*fk[1] + fk[2]*fk[2];
  
    if (fk2 < TINY*TINY)
    { 
	printf ("\nwarning from ewald_dynamical():\n");
	printf ("|k| is too small\n");
	printf ("(kx=%e, ky=0, kz=0) is assumed.\n\n",TINY);
	fk[0] = TINY;
	fk[1] = fk[2] = 0.;
	fk2 = fk[0]*fk[0] + fk[1]*fk[1] + fk[2]*fk[2];
    }
    
    /* Initialize the incoming matrix; in f77, */
    /* complex number is stored in the simple  */
    /* way of (real,imag) in memory            */
    for (i=0; i<3*NP; i++)
	for (j=0; j<6*NP; j++)
	    phi[i*2*(*L)+j] = 0.;

    /* evaluating the upper half is enough: */
    for (m=0; m<NP; m++)
	for (n=m; n<NP; n++)
	{
	    for (i=0; i<3; i++)
		for (j=0; j<6; j++)
		    total[i][j] = 0.;
	    
	    /* the second term of (6.2.46) */ 
	    if (m==n)
		for (p=0; p<NP; p++)
		{
		    calc_ewald_dynamical_Q
			(zero_k, m, p, Q);
		    for (i=0; i<3; i++)
			for (j=0; j<6; j++)
			    total[i][j] += Q[i][j];
		}
	    
	    /* the third term of (6.2.46) */ 
	    calc_ewald_dynamical_Q (fk, m, n, Q);
	    for (i=0; i<3; i++)
		for (j=0; j<6; j++)
		    total[i][j] -= Q[i][j];
	    
	    /* the fourth term of (6.2.46) */ 
	    /* macroscopic electric field  */
	    /* contribution: see (6.2.27)  */ 
	    for (i=0; i<3; i++)
		for (j=0; j<3; j++)
		    total[i][2*j] += q[m] * q[n] *
			4.*PI/volume*fk[i]*fk[j]/fk2;
	    
	    /* write in Phi, make it Hermitian: */
	    for (i=0; i<3; i++)
		for (j=0; j<3; j++)
		{
		    phi[(n*3+j)*2*(*L)+(m*3+i)*2]
			= total[i][2*j]; 
		    phi[(n*3+j)*2*(*L)+(m*3+i)*2+1]
			= total[i][2*j+1];
		    phi[(m*3+i)*2*(*L)+(n*3+j)*2]
			= total[i][2*j]; 
		    phi[(m*3+i)*2*(*L)+(n*3+j)*2+1]
			= -total[i][2*j+1];
		}
	}
    return;
} /* end ewald_dynamical() */


/* Make the subroutines accessible to Fortran */
void _init_cell (double H[3][3])
{init_cell(H);}

void init_cell_ (double H[3][3])
{init_cell(H);}
    
void _init_ewald
(int *number_particles,double H[3][3],double *accuracy)
{init_ewald(number_particles, H, accuracy);}

void init_ewald_
(int *number_particles,double H[3][3],double *accuracy)
{init_ewald(number_particles, H, accuracy);}

void _exit_ewald()
{exit_ewald();} 

void exit_ewald_()
{exit_ewald();} 

void _ewald
(double *charge, double *s1, double *s2,
 double *s3, double H[3][3], double *pote,
 double *fx, double *fy, double *fz, double stress[3][3])
{ewald(charge,s1,s2,s3,H,pote,fx,fy,fz,stress);}

void ewald_
(double *charge, double *s1, double *s2,
 double *s3, double H[3][3], double *pote,
 double *fx, double *fy, double *fz, double stress[3][3])
{ewald(charge,s1,s2,s3,H,pote,fx,fy,fz,stress);}

void _ewald_dynamical (double fk[3], double *phi, int *L)
{ewald_dynamical(fk, phi, L);}
    
void ewald_dynamical_ (double fk[3], double *phi, int *L)
{ewald_dynamical(fk, phi, L);}


/*
  Notes:
  
  1. Phi((m-1)*3+i,(n-1)*3+j) (m,n=1..NP, i,j=1..3, in f77 notation)
  is the force constant matrix of ideal rigid-ion system: to get the
  dynamical matrix we have to divide each term by sqrt(mass_m*mass_n),
  before we do the diagonalization. Because we are yet going to add 
  the short ranged interaction part, it is not done here.
  
  2. It is important to note that this matrix is only complete for
  rigid-ion model. In general if we write E_c = E_c(.. ,q_m,q_n,
  .. ,x_m,x_n, ..), that is, the long ranged Coulomb interaction is
  modelled "as if" by putting point charges at sites, this force
  matrix is complete only if {q_m} doesn't vary with x_n's. There are
  two extra terms entering the second derivative and unfortunately
  they are both long-ranged. This "dynamic charge transfer" is the
  reason for the discrepency between Z^* and e_T^*, described in
  Harrison <Elec. Struc. & Prop. Solids>, Page 220. In fact, it is
  easy to see that the long range effect (change of forces on other
  atoms) of infinitesimal charge transfer occuring in a finite region,
  due to the infinitesimal displacement of one atom, is equivalent to
  that of the displacement itself, both by creating a small dipole.
  For a tetrahedral system or systems with higher symmetry, that
  dipole created by charge transfer must also be in the direction of
  the displacement, thus defining an effective e_T^* is completely
  plausible.

  3. For lattice waves of displacing one type of atoms sinusoidally
  along certain k, we can add all the small dipoles together using
  linear superposition, which will create a polarization field with
  the same wave-vector k. The problem simplifies as k->0: if the
  displacements are perpendicular to k (transverse direction), there
  will be no spatial charge accumulation (\rho = -\nabla \cdot P) and
  so there will be no macroscopic electric field; if the displacements
  are in logitudinal direction then there will be spatial charge
  accumulation and corresponding excited electric field, but that
  field will induce further polarization of the electrons (charge
  re-distributed again) and so it will be screened. Thus the actual
  polarization field induced by logitudinal displacements is divided
  by the electronic dielectric constant \epsilon(\omega=0,k=0), and so
  are the forces on individual atoms. So the end the LO-TO splitting
  at k=0 is proportional to (e_T^*)^2 / \epsilon.
  
*/
