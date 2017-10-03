/*******************************************************************/
/* Simple Ewald V1.0:                                              */
/*                                                                 */
/* Calculate the electrostatic potential coefficients of general   */
/* 3D triclinic PBC rigid-ion systems with zero total charge. The  */
/* result is for the limit of (infinitely) many replicas stacked   */
/* in a nearly spherical cluster and surrounded by metallic sheath */
/* outside, so the cluster has no surface charge distribution and  */
/* there is no uniform macroscopic electric field; it is the only  */
/* self-consistent set-up for both static and MD calculations in   */
/* PBC. For details see reference "Doc/moldy.ps" and websites      */
/*                                                                 */
/*   http://www.ee.duke.edu/~ayt/ewaldpaper/                       */
/*   http://www.keele.ac.uk/depts/ph/tc/cph_res/dynamo.html        */
/*   http://www.earth.ox.ac.uk/%7Ekeith/moldy.html                 */
/*   http://www.fos.su.se/physical/sasha/md_prog.html              */
/*                                                                 */
/*                                                 Aug.26, 1998    */
/*                                    Developed by Ju Li at MIT    */
/*******************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define PI acos(-1.0)

/* calculate the electrostatic potential coefficients */
/* (dimensionless) which depends only on structure.   */
/* pote_matrix[] should be larger than NP*NP.         */
void simple_ewald
(double *h11,double *h12,double *h13,
 double *h21,double *h22,double *h23,
 double *h31,double *h32,double *h33,
 int *NP,double *s1,double *s2,double *s3,
 double *accuracy, double *pote_matrix)
{
    /* Our convention is that (h11,h12,h13) is the    */
    /* (x,y,z) component of the first cell edge, etc, */
    /* and so is (g11,g12,g13): H*G^T = 2\pi I.       */
    double g11,g12,g13,g21,g22,g23,g31,g32,g33,volume;
    /* Distances from origin to the three */
    /* real and reciprocal unit surfaces  */ 
    double rd1, rd2, rd3, kd1, kd2, kd3;
    double alpha, rcut, kcut;
    int i, j, ij_index;
    int n1, n2, n3, max_n1, max_n2, max_n3;
    int k1, k2, k3, max_k1, max_k2, max_k3;
    double ds1, ds2, ds3, dx, dy, dz, r2, r;
    
    /* Get the reciprocal lattice vectors */
    g11 = *h22**h33 - *h23**h32;
    g22 = *h33**h11 - *h31**h13;
    g33 = *h11**h22 - *h12**h21;
    g12 = *h23**h31 - *h21**h33;
    g23 = *h31**h12 - *h32**h11;
    g31 = *h12**h23 - *h13**h22;
    g13 = *h21**h32 - *h31**h22;
    g21 = *h32**h13 - *h12**h33;
    g32 = *h13**h21 - *h23**h11;
    volume = *h11*g11 + *h12*g12 + *h13*g13;
    /* the shortest distances to respective surfaces */
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
    kd1 = 2*PI/sqrt(*h11**h11+*h12**h12+*h13**h13);
    kd2 = 2*PI/sqrt(*h21**h21+*h22**h22+*h23**h23);
    kd3 = 2*PI/sqrt(*h31**h31+*h32**h32+*h33**h33);
    volume = fabs(volume);
    
    /* Set the parameters alpha and cutoffs based on
       the formula by Fincham CCP5 38, p17 (1993) */ 
    alpha = pow(*NP*PI*PI*PI*5.5/volume/volume,1./6.);
    rcut = sqrt(-log(*accuracy))/alpha*1.2;
    kcut = 2*alpha*sqrt(-log(*accuracy))*1.2;
    
    max_n1 = ceil(rcut/rd1);
    max_n2 = ceil(rcut/rd2);
    max_n3 = ceil(rcut/rd3);
    max_k1 = ceil(kcut/kd1);
    max_k2 = ceil(kcut/kd2);
    max_k3 = ceil(kcut/kd3);
    
    /* make the particles inside [-0.5,0.5] so the  */
    /* maximum unit cell extension scheme will work */
    for (i=0; i<*NP; i++)
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
    
    for (i=0; i<*NP; i++)
	for (j=i; j<*NP; j++)
	{
	    ds1 = s1[j] - s1[i];
	    ds2 = s2[j] - s2[i];
	    ds3 = s3[j] - s3[i];
	    /* pointer in the electrostatic */
	    /* potential coefficients array */
	    ij_index = i*(*NP) + j;
	    pote_matrix [ij_index] = 0.;
	    
	    /* real space lattice sum */
	    for (n1=-max_n1; n1<=max_n1; n1++)
		for (n2=-max_n2; n2<=max_n2; n2++)
		    for (n3=-max_n3; n3<=max_n3; n3++)
		    {
			dx = (n1+ds1)**h11+(n2+ds2)**h21+(n3+ds3)**h31;
			dy = (n1+ds1)**h12+(n2+ds2)**h22+(n3+ds3)**h32;
			dz = (n1+ds1)**h13+(n2+ds2)**h23+(n3+ds3)**h33;
			r2 = dx*dx + dy*dy + dz*dz;
			if (r2 < rcut*rcut)
			{
			    r = sqrt(r2);
			    pote_matrix [ij_index] += 
				(r!=0.)?erfc(alpha*r)/r:-2.*alpha/sqrt(PI);
			}
		    }

	    /* reciprocal space lattice sum */
	    for (k1=-max_k1; k1<=max_k1; k1++)
		for (k2=-max_k2; k2<=max_k2; k2++)
		    for (k3=0; k3<=max_k3; k3++)
			if (!((k1==0)&&(k2==0)&&(k3==0)))
			{
			    dx = k1*g11 + k2*g21 + k3*g31;
			    dy = k1*g12 + k2*g22 + k3*g32;
			    dz = k1*g13 + k2*g23 + k3*g33;
			    r2 = dx*dx + dy*dy + dz*dz;
			    /* use inversion symmetry */
			    if (r2 <= kcut*kcut)
				pote_matrix [ij_index] += (k3>0)?
				    8.*PI*exp(-r2/alpha/alpha/4.)/r2/volume
				    *cos(2*PI*(k1*ds1+k2*ds2+k3*ds3)):
				4.*PI*exp(-r2/alpha/alpha/4.)/r2/volume
				    *cos(2*PI*(k1*ds1+k2*ds2+k3*ds3));
			}
	    
	    /* symmetrize the coefficient matrix */
	    pote_matrix [j*(*NP)+i] = pote_matrix [ij_index];
	}
    return;
} /* end simple_ewald() */

/******************************************************/
/* Calculate the electrostatic potential at particle  */
/* i after the coefficient matrix has been evaluated: */
/* *scale is the dimensional prefactor. Notice that   */
/* if Q is electron occupation, then -eU is what is   */
/* returned, which is OK if used consistently.        */
/******************************************************/
double ewald_potential
(int *j, double *scale, int *NP, double *Q, double *pote_matrix)
{
    int i;
    double totalcharge=0.,potential=0.;
    /* renormalize the total charge to 0 */
    for (i=0; i<*NP; i++) totalcharge += Q[i];
    /* notice that we expect *j to start from 1 (f77 convention) */
    for (i=0; i<*NP; i++) potential += pote_matrix [(*j-1)**NP+i] * 
			      (*scale) * (Q[i]-totalcharge/(*NP));
    return (potential);
} /* end ewald_potential() */


/* Calculate the total energy after the   */
/* coefficient matrix has been evaluated. */
double total_ewald_energy
(double *scale, int *NP, double *Q, double *pote_matrix)
{
    int i,j;
    double totalcharge=0.,tote=0.;
    for (i=0; i<*NP; i++) totalcharge += Q[i]; 
    for (j=0; j<*NP; j++)
	for (i=0; i<*NP; i++)
	    tote += (Q[i]-totalcharge/(*NP)) *
		(Q[j]-totalcharge/(*NP)) *
		(*scale) * pote_matrix[i**NP+j] / 2.;
    return (tote);
} /* end total_ewald_energy() */

/* Make C subroutines accessible to all kinds of Fortran compilers */
void _simple_ewald
(double *h11,double *h12,double *h13,
 double *h21,double *h22,double *h23,
 double *h31,double *h32,double *h33,
 int *NP,double *s1,double *s2,double *s3,
 double *accuracy, double *pote_matrix)
{simple_ewald(h11,h12,h13,h21,h22,h23,h31,h32,h33,
	      NP,s1,s2,s3,accuracy,pote_matrix);}
void simple_ewald_
(double *h11,double *h12,double *h13,
 double *h21,double *h22,double *h23,
 double *h31,double *h32,double *h33,
 int *NP,double *s1,double *s2,double *s3,
 double *accuracy, double *pote_matrix)
{simple_ewald(h11,h12,h13,h21,h22,h23,h31,h32,h33,
	      NP,s1,s2,s3,accuracy,pote_matrix);}

double _ewald_potential
(int *j,double *scale,int *NP,double *Q,double *pote_matrix)
{return(ewald_potential(j,scale,NP,Q,pote_matrix));}
double ewald_potential_
(int *j,double *scale,int *NP,double *Q,double *pote_matrix)
{return(ewald_potential(j,scale,NP,Q,pote_matrix));}

double _total_ewald_energy
(double *scale, int *NP, double *Q, double *pote_matrix)
{return(total_ewald_energy(scale,NP,Q,pote_matrix));}
double total_ewald_energy_
(double *scale, int *NP, double *Q, double *pote_matrix)
{return(total_ewald_energy(scale,NP,Q,pote_matrix));}
