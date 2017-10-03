/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
*/

#include <stdlib.h>
#include <string.h>
#ifdef __GNUC__
/* g2c.h comes with g77 and provides typedefs for integer and doublereal
 * needed by c_interface.h; if not available, guess at the typedefs
 */
#include "g2c.h"
#else
typedef long int integer;
typedef double doublereal;
#endif
#include "c_interface.h"

#include "brenfort.h"

int
set_cell_threshold(int n)
{
  integer nf = n;
  return setcel_(&nf);
}

BrennerMainInfo *
alloc_bren(int *kflag)
{
	BrennerMainInfo *info = (BrennerMainInfo*)malloc(sizeof(BrennerMainInfo));
	memset(info, 0, sizeof(*info));
	info->temperature = 300;
	info->num_atms = 0;
	info->steps = 100000;
	info->timestep = 0.5;
	info->starttime = 0;
	info->lchk = 1;
	info->num_thermostated_atoms = 0;
	info->head[0] = 0;
	info->xmass[1] = 12.0;
	info->xmass[2] = 1.0;
	info->xmass[3] = 28.0;
	info->xmass[4] = 72.0;
	info->nxmol = 500;
	info->RLL = 0.5;
	info->PSEED = 0.6955;
	info->kt[HYDROGEN] = 2;
	info->kt[CARBON] = 1;
	info->kt[SILICON] = 3;
	info->kt[GERMANIUM] = 4;
	info->kt2[2] = HYDROGEN;
	info->kt2[1] = CARBON;
	info->kt2[3] = SILICON;
	info->kt2[4] = GERMANIUM;
	info->param_file_dir = "";
	return info;
}

static void
read_d_files(const BrennerMainInfo *info)
{
	int i, j, k, l, m, n;
	integer i_f, j_f, k_f, l_f, m_f, n_f;
	int itd, i2d, i3d;
	FILE *fd;
	char filename[512];
	sprintf(filename, "%s%s", info->param_file_dir, "inter2d_iv.d");
	if ((fd = fopen(filename, "r")) == NULL)
	{
		fprintf(stderr, "opening %s\n", filename);
		exit(-1);
	}
	fscanf(fd, "%d\n",&i2d);

	/*Read bicubic spline coefficients*/
	while (1) {
		double xhh;
		doublereal xhh_dr;
		if(fscanf(fd, "%d %d %d %lf\n",&i,&j,&k,&xhh) != 4)
		{
		  fprintf(stderr, "error reading inter2d_iv.d\n");
		  exit(-1);
		}
		if (i <= 0) 
			break;
		xhh_dr = xhh;
		i_f = i;
		j_f = j;
		k_f = k;
		setxh_(&i_f, &j_f, &k_f, &xhh_dr);
	}
	for(k = 1; k < 73; ++k)
	{
	    fscanf(fd, "%d %d %d\n",&i, &l, &m);
	    i_f = i;
	    l_f = l;
	    m_f = m;
	    for(j = 1; j < 16; j += 4)
	    {
		double clm_tmp[4];
		doublereal clm_f[4];
		if(fscanf(fd, "%lf %lf %lf %lf\n", &clm_tmp[0], &clm_tmp[1],
			  &clm_tmp[2], &clm_tmp[3]) != 4)
		{
		  fprintf(stderr, "error reading inter2d_iv.d clm %d\n", j);
		  exit(-1);
		}
		j_f = j;
		clm_f[0] = clm_tmp[0];
		clm_f[1] = clm_tmp[1];
		clm_f[2] = clm_tmp[2];
		clm_f[3] = clm_tmp[3];
		setclm_(&i_f, &l_f, &m_f, &j_f, &clm_f[0], &clm_f[1],
			&clm_f[2], &clm_f[3]);
	    }
	}
	fclose(fd);

/* Read tricubic spline coefficients*/

	for(i = 0; i < 3; ++i)
	{
	  static const char *inter_file[] = {"inter3d_iv_new.d",
					     "inter3d_ch.d",
					     "inter3d_h.d" };
	  FILE *fd;
	  sprintf(filename, "%s%s", info->param_file_dir, inter_file[i]);
	  if ((fd = fopen(filename, "r")) == NULL)
	  {
	    fprintf(stderr,"opening %s failed", filename);
	    exit(-1);
	  }
	  fscanf(fd, "%d\n",&i3d);
	  while(1) {
		if(fscanf(fd, "%d %d %d\n",&l,&m,&n) != 3)
		  break;

		l_f = l;
		m_f = m;
		n_f = n;
		for(j = 1; j <= 64; ++j) {
			double clmn1;
			doublereal clmnf;
			integer i1 = i+1;
			if(!fscanf(fd,"%lf\n", &clmn1))
			{
			  fprintf(stderr, "parse error in file %s, %d %d %d %d\n",
				  inter_file[i], l, m, n, i);
			  exit(-1);
			}
			j_f = j;
			clmnf = clmn1;
			setcl2_(&i1, &l_f, &m_f, &n_f, &j_f, &clmnf);
		}
	  }
	  fclose(fd);
	}
/* Read tricubic spline coefficients for torsional potential*/
	
	sprintf(filename, "%s%s", info->param_file_dir, "inter3dtors.d");
	if ((fd = fopen(filename, "r")) == NULL)
	{
	  fprintf(stderr, "opening %s failed\n", filename);
	  exit(-1);
	}
	fscanf(fd, "%d\n",&itd);
	for(j = 1; j < 109; ++j) {
		double t;
		doublereal t_f;
		if(fscanf(fd, "%d %d %d\n",&l,&m,&n) != 3)
		{
		  fprintf(stderr, "read error in inter3dtors.d %d\n", j);
		  exit(-1);
		}
		l_f = l;
		m_f = m;
		n_f = n;
		for(i = 1; i < 64; i += 3)
		{
			double t1,t2,t3;
			if(fscanf(fd, "%lf %lf %lf\n", &t1, &t2, &t3) != 3)
			{
			  fprintf(stderr,
				  "error reading inter3dtors.d @ %d, %d\n",
				  j, i);
			  exit(-1);
			}
			i_f = i;
			t_f = t1;
			sett_(&l_f, &m_f, &n_f, &i_f, &t_f);
			i_f = i + 1;
			t_f = t2;
			sett_(&l_f, &m_f, &n_f, &i_f, &t_f);
			i_f = i + 2;
			t_f = t3;
			sett_(&l_f, &m_f, &n_f, &i_f, &t_f);
		}
		fscanf(fd, "%lf\n",&t);
		i_f = i;
		t_f = t;
		sett_(&l_f, &m_f, &n_f, &i_f, &t_f);
	}
	fclose(fd);

	if((itd != i2d) || (itd != i3d))
	{
            fprintf(stderr, "Incompatable potential types %d %d %d\n",
		    itd, i2d, i3d);
	    exit(-1);
	}
}

void
init_bren(BrennerMainInfo *info, int kflag, int tight_binding)
{
/*    Tim Freeman's comments from the Fortran code:
*     KFLAG says what sort of work to do.
*     It is read from input.d, third number on the first line.
*     See cases in subroutine thermos in thermostats.f. 
*     -1 is friction plus random force.
*     1 is a Berndstein thermostat.
*     2 means to set all the velocities to zero after each timestep.
*     3 is an Evans-Hoover thermostat.
*     4 has no significance whatsoever; I'm hoping that will give me
*       simple physics.
*     5 also calls BERE to do a Berndstein thermostat.
*     6 means to do energy minimization.
*     8 also has significance, but I don't know what it means yet.
*       vscale sets kflag to 8 temporarily then restores it to its old
*       value.
*/

      integer kflag_ir = kflag;
      integer ipot = (tight_binding ? 2 : 1);
      doublereal pseed_dr = info->PSEED;
      doublereal rll_dr = info->RLL;
      doublereal temperature_dr = info->temperature;
      int i;
      doublereal t0 = info->starttime;
      doublereal delta0 = info->timestep;
      doublereal cube0[3];
      integer num_atms = info->num_atms;
      doublereal *pos0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *vel0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *acc0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *r30 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *r40 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      integer *itr0 = (integer *)malloc(num_atms*sizeof(integer));
      integer *atm_no = (integer *)malloc(num_atms*sizeof(integer));

      /* init input data */
      inputd_(&kflag_ir, &pseed_dr, &rll_dr, &temperature_dr, &ipot);

      for(i = 0; i < info->num_atms; ++i)
      {
	pos0[3*i] = info->atm_num[i].coord.x;
	pos0[3*i+1] = info->atm_num[i].coord.y;
	pos0[3*i+2] = info->atm_num[i].coord.z;
	vel0[3*i] = info->atm_num[i].velocity.x;
	vel0[3*i+1] = info->atm_num[i].velocity.y;
	vel0[3*i+2] = info->atm_num[i].velocity.z;
	acc0[3*i] = info->atm_num[i].accel.x;
	acc0[3*i+1] = info->atm_num[i].accel.y;
	acc0[3*i+2] = info->atm_num[i].accel.z;
	r30[3*i] = info->atm_num[i].dx3dt3.x;
	r30[3*i+1] = info->atm_num[i].dx3dt3.y;
	r30[3*i+2] = info->atm_num[i].dx3dt3.z;
	r40[3*i] = 0;
	r40[3*i+1] = 0;
	r40[3*i+2] = 0;
#if 0
	if(info->atm_num[i].thermostated)
	  ++info->num_thermostated_atoms;
#endif
	itr0[i] = (!info->atm_num[i].movable ? 2 : info->atm_num[i].thermostated);
	atm_no[i] = info->atm_num[i].type;
      }
#ifdef INFINITE_CUBE
      cube0[0] = 200;		/* arbitrary; this program probably should */
      cube0[1] = 200;		/* not be used with INFINITE_CUBE */
      cube0[2] = 200;
      fprintf(stderr, "warning - Fortran version has not been tested with INFINITE_CUBE\n");
#else
      cube0[0] = info->cube[0];
      cube0[1] = info->cube[1];
      cube0[2] = info->cube[2];
#endif
      coordd_(pos0, vel0, acc0, r30, r40, itr0, atm_no,
	      &num_atms, &t0, &delta0, cube0);
      free(pos0);
      free(vel0);
      free(acc0);
      free(r30);
      free(r40);
      free(itr0);
      free(atm_no);

      /* setup potential parameters  */
      read_d_files(info);
      csetpp_();
      /* setup Langevin parameters  */
      setgle_();
      /* initialize random number generator */
      setran_();

      if(kflag == 6)
      {
/**********************
* Energy minimization *
**********************/
	   csetmi_();
           minimize_();
      }
      else 
/**********************
* Dynamic Simulation  *
**********************/
      {
	csetmd_();
      }
}

void
update_info(BrennerMainInfo *info)
{
      int num_atms = info->num_atms;
      int i;
      doublereal *pos0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *vel0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *acc0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *r30 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *r40 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *rnp0 = (doublereal *)malloc(3*num_atms*sizeof(doublereal));
      doublereal *energy = (doublereal *)malloc(num_atms*sizeof(doublereal));
      doublereal syse;

      getd_(pos0, vel0, acc0, r30, r40, rnp0, energy, &syse);

      for(i = 0; i < info->num_atms; ++i)
      {
	/*	info->atm_num[i].energy = energy[i]; */
	info->atm_num[i].coord.x = pos0[3*i];
	info->atm_num[i].coord.y = pos0[3*i+1];
	info->atm_num[i].coord.z = pos0[3*i+2];
	info->atm_num[i].velocity.x = vel0[3*i];
	info->atm_num[i].velocity.y = vel0[3*i+1];
	info->atm_num[i].velocity.z = vel0[3*i+2];
	info->atm_num[i].accel.x = acc0[3*i];
	info->atm_num[i].accel.y = acc0[3*i+1];
	info->atm_num[i].accel.z = acc0[3*i+2];
	info->atm_num[i].dx3dt3.x = r30[3*i];
	info->atm_num[i].dx3dt3.y = r30[3*i+1];
	info->atm_num[i].dx3dt3.z = r30[3*i+2];
	info->atm_num[i].force.x = rnp0[3*i];
	info->atm_num[i].force.y = rnp0[3*i+1];
	info->atm_num[i].force.z = rnp0[3*i+2];

	/* # bonds? */
      }
      info->system_energy = syse;
      free(pos0);
      free(vel0);
      free(acc0);
      free(r30);
      free(r40);
      free(rnp0);
      free(energy);
}

void
bren_1_step(BrennerMainInfo *info, int kflag)
{
	/* predictor */
	cpred_();
	/* calculate energy and forces */
	model_();
	/* apply thermostats */
	thermos_();
	/* corrector */
	ccorr_();
	if(kflag == 5) bere_();
	if(info)
	  update_info(info);
}
