/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
*/

/*
 * This is a sample program to demonstrate the C Brenner potential library.
 */

#include "brenner.h"
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>

#ifndef __SAFE_FGETS_H__
#include "safe_fgets.h"
#endif

#ifdef MCHECK
#ifndef __mcheck_h__
#include <mcheck.h>
#define __mcheck_h__
#endif
#endif

void
init_coord(BrennerMainInfo *info)
{
  int itr; /* thermostat flag */
  BrenAtom *atm_num = info->atm_num;
  // FIXME We're probably vulnerable to buffer overrun errors on buf.  Tim
  // Freeman 13 Sep 2000.
  char buf[128];
  int np;
  int i, k, natom;
  int dummy;
  FILE *fp = fopen("coord.d", "r");
  if(!fp)
  {
    fprintf(stderr, "can't read coord.d\n");
    my_exit(-1);
  }
  
  safe_fgets(info->head, sizeof(info->head), fp,
	     "header line from coord.d");
  i = strlen(info->head);
  /* Throw out the tailing newline, if there is anything there. */
  if(i)
    info->head[i-1] = 0;
  safe_fgets (buf, sizeof(buf), fp, "number of atoms from coord.d");
  sscanf(buf, "%d", &np);
  safe_fgets(buf, sizeof(buf), fp, "start time and timestep from coord.d");
  {
    double starttime, timestep;
    if (sscanf(buf, "%lf %lf", &starttime, &timestep) != 2) {
      fprintf (stderr, "Failed to parse start time and timestep from coord.d, "
	       "this line:\n%s", buf);
      my_exit (-1);
    }
    info->starttime = starttime;
    info->timestep = timestep;
}
  printf ("Start time is %f, timestep is %f.\n",
	  info->starttime, info->timestep);
  safe_fgets(buf, sizeof(buf), fp, "periodic cube from coord.d");
  {
    double c0, c1, c2;
    if (sscanf(buf, "%lf %lf %lf", &c0, &c1, &c2) != 3) {
      fprintf (stderr, "Failed to parse periodic cube from coord.d, "
	       "this line:\n%s", buf);
      my_exit (-1);
    }
#ifndef INFINITE_CUBE
    info->cube[0] = c0;
    info->cube[1] = c1;
    info->cube[2] = c2;
#endif
  }  
  if(np > MAX_ATOMS)
  {
    fprintf(stderr, "np= %d greater than npmax= %d\n", np, MAX_ATOMS);
    my_exit(-1);
  }
  for(i = 0; i < np; ++i)
  {
    if(!fgets(buf, sizeof(buf), fp))
    {
      fprintf(stderr, "warning - expected %d atoms in coord.d, found %d\n",
	      np, i);
      break;
    }
    sscanf(buf, "%d", &k);
    --k;
    assert(k >= 0 && k < MAX_ATOMS);
    {
      double x, y, z;
      if (sscanf(buf, "%d %d %lf %lf %lf %d", &dummy, &natom, &x, &y, &z,
		 &itr) != 6) {
	fprintf (stderr,
		 "Failed to parse this coordinate line from coord.d:\n%s",
		 buf);
	my_exit (-1);
      }
      atm_num[k].coord.x = x;
      atm_num[k].coord.y = y;
      atm_num[k].coord.z = z;
    }
#ifndef INFINITE_CUBE
    if(atm_num[k].coord.x > info->cube[0]/2)
      atm_num[k].coord.x -= info->cube[0];
    if(atm_num[k].coord.y > info->cube[1]/2)
      atm_num[k].coord.y -= info->cube[1];
    if(atm_num[k].coord.z > info->cube[2]/2)
      atm_num[k].coord.z -= info->cube[2];
#endif
    atm_num[k].prev_coord.x = atm_num[k].coord.x;
    atm_num[k].prev_coord.y = atm_num[k].coord.y;
    atm_num[k].prev_coord.z = atm_num[k].coord.z;
    atm_num[k].ktype = info->kt[natom];
    ++info->noa[info->kt[natom]];

    atm_num[k].number = k;
    atm_num[k].type = natom;
    atm_num[k].mass = info->xmass[info->kt[natom]];

    atm_num[k].thermostated = 0;
    if(itr == 2)
      atm_num[k].movable = 0;
    else
    {
      atm_num[k].movable = 1;
      if(itr == 1)
	atm_num[k].thermostated = 1;
    }
    ++info->num_atms;
  }
  for(i = 0; i < np; ++i)
  {
    if(!fgets(buf, sizeof(buf), fp))
    {
      fprintf(stderr, "warning - expected %d atoms in coord.d velocity, found %d\n",
	      np, i);
      break;
    }
    sscanf(buf, "%d", &k);
    --k;
    assert(k >= 0 && k < MAX_ATOMS);
    {
      double x, y, z;
      if(sscanf(buf, "%d %lf %lf %lf", &dummy, &x, &y, &z) != 4)
	{
	  fprintf(stderr, "Error parsing velocity %d, %s\n", k, buf);
	  my_exit(-1);
	}
      atm_num[k].velocity.x = x;
      atm_num[k].velocity.y = y;
      atm_num[k].velocity.z = z;
    }
    atm_num[k].velocity.x *= info->timestep;
    atm_num[k].velocity.y *= info->timestep;
    atm_num[k].velocity.z *= info->timestep;
  }
  for(i = 0; i < np; ++i)
  {
    if(!fgets(buf, sizeof(buf), fp))
    {
      fprintf(stderr, "warning - expected %d atoms in coord.d accel, found %d\n",
	      np, i);
      break;
    }
    sscanf(buf, "%d", &k);
    --k;
    assert(k >= 0 && k < MAX_ATOMS);
    {
      double x, y, z;
      if(sscanf(buf, "%d %lf %lf %lf", &dummy, &x, &y, &z) != 4)
	{
	  fprintf(stderr, "Error parsing accel %d, %s\n", k, buf);
	  my_exit(-1);
	}
      atm_num[k].accel.x = x;
      atm_num[k].accel.y = y;
      atm_num[k].accel.z = z;
    }
  }
  for(i = 0; i < np; ++i)
  {
    if(!fgets(buf, sizeof(buf), fp))
    {
      fprintf(stderr, "warning - expected %d atoms in coord.d dx3dt3, found %d\n",
	      np, i);
      break;
    }
    sscanf(buf, "%d", &k);
    --k;
    assert(k >= 0 && k < MAX_ATOMS);
    {
      double x, y, z;
      if(sscanf(buf, "%d %lf %lf %lf", &dummy, &x, &y, &z) != 4)
	{
	  fprintf(stderr, "Error parsing dx3dt3 %d, %s\n", k, buf);
	  my_exit(-1);
	}
      atm_num[k].dx3dt3.x = x;
      atm_num[k].dx3dt3.y = y;
      atm_num[k].dx3dt3.z = z;
    }
  }
  printf("read %d atoms\n", info->num_atms);
  fclose(fp);
}

static void
append_file(FILE *fp, const BrennerMainInfo *info, Float TTIME)
{
      int i;
      int atm;
      int NRA = 0, NLA = 0;
      fprintf(fp, "%s\n", info->head);
      fprintf(fp, "%6d%6d%6d%6d\n", info->num_atms, 3, NRA,NLA);
      fprintf(fp, "%20.11e%20.11e\n", TTIME, info->timestep);
#ifdef INFINITE_CUBE
      fprintf(fp, "%20.11e%20.11e%20.11e\n", 1000.0, 1000.0, 1000.0);
#else
      fprintf(fp, "%20.11e%20.11e%20.11e\n", info->cube[0],
	      info->cube[1], info->cube[2]);
#endif
      for (atm = 0; atm < info->num_atms; atm++)
      {
           fprintf(fp, "%5d%5d%20.11e%20.11e%20.11e\n",
		   atm+1, info->kt2[info->atm_num[atm].ktype],
		   info->atm_num[atm].coord.x,
		   info->atm_num[atm].coord.y, info->atm_num[atm].coord.z);
      }
}

int
main(int argc, char **argv)
{
	int i, j;
	const char *coord_file = NULL;
	FILE *fp;
	int kflag;
	BrennerMainInfo *info = alloc_bren(&kflag);

#ifdef MCHECK
    mcheck (0);
    // printf ("Initialized mcheck.\n");
#endif

	init_coord(info);

	while ((i = (int) getopt(argc, argv, "?0c:s:t:k:m:v:zR:L:")) != -1)
	{
	  extern char *optarg;
	  switch(i)
	  {
	  case '?':
	    printf("usage: %s [-s <#steps>] [-t <temperature(K)>] [-c <name of output file>]\n", argv[0]);
	    exit(-1);
	  case 'c':
	    coord_file = strdup(optarg);
	    break;
	  case 's':
	    info->steps = atoi(optarg);
	    if (info->steps == 0) {
		printf("Invalid entry: '%s', setting steps to 20\n", optarg);
		info->steps = 20;
	    }
	    break;
	  case 'L':
	    info->min_lj_range_tree = atoi(optarg);
	    break;
	  case 'R':
	    info->min_range_tree = atoi(optarg);
	    break;
	  case 'm':
	    info->nxmol = atoi(optarg);
	    break;
	  case 't':
	    info->temperature = atof(optarg);
	    break;
	  case 'k':
	    kflag = atoi(optarg);
	    break;
	  case 'z':
	    info->squeeze = 1;
	    break;
	  case 'v':
	    info->volume_scale_dir = atoi(optarg);
	    if(info->volume_scale_dir < -1 || info->volume_scale_dir > 2)
	    {
	      fprintf(stderr,
		      "warning, invalid volume_scale_dir (%s), changed to -1\n",
		      optarg);
	      info->volume_scale_dir = -1;
	    }
	  case '0':
	    info->zero_com_velocity = 1;
	    break;
	  }
	}
	if(kflag == 6 && !coord_file)
	  printf("Warning: -k 6 isn't of much use without -c\n");
	fp = fopen("xmol.d", "w");
	if(!fp)
	{
	  fprintf(stderr, "cant write to xmol.d\n");
	  my_exit(-1);
	}
	init_bren(info, kflag, 0);
	if(kflag != 6)
      for (j=0; j < info->steps; ) {
        bren_1_step(info, kflag);
        j++;
        if(!(j % info->nxmol)) {
          Float now = info->starttime + j*info->timestep;
          append_file(fp, info, now);
          printf("Time = %.3f fs ", now);
          printf("Total energy = %.5f eV \n", info->system_energy);
        }
      }
	fclose(fp);
	if(coord_file) {
      FILE *fp1 = fopen(coord_file, "w");
      info->starttime += info->steps * info->timestep;
	  write_coord_file(fp1, info);
	  fclose(fp1);
	}
	free(info);
	return 0;
}

void my_exit(int code)
{
  exit(code);
}
