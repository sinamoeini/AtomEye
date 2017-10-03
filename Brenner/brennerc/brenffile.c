/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
*/

/*
 * This is a sample program to demonstrate the C interface to the Fortran
 * Brenner potential library.
 */

#include "brenfort.h"
#include <stdlib.h>
#include <assert.h>

static void
init_coord(BrennerMainInfo *info)
{
  int itr; /* thermostat flag */
  BrenAtom *atm_num = info->atm_num;
  char buf[128];
  int np;
  int i, k, natom;
  int dummy;
  FILE *fp = fopen("coord.d", "r");
  if(!fp)
  {
    fprintf(stderr, "can't read coord.d\n");
    exit(-1);
  }
  fgets(info->head, sizeof(info->head), fp);
  i = strlen(info->head);
  if(i)
    info->head[i-1] = 0;
  fgets(buf, sizeof(buf), fp);
  sscanf(buf, "%d", &np);
  fgets(buf, sizeof(buf), fp);
  sscanf(buf, "%f %f", &info->starttime, &info->timestep);
  fgets(buf, sizeof(buf), fp);
#ifndef INFINITE_CUBE
     		/* this program probably should */
     		/* not be used with INFINITE_CUBE */
  if(info->cube[0] < 1.e-6)
    sscanf(buf, "%f %f %f", &info->cube[0], &info->cube[1], &info->cube[2]);
#endif

  if(np > MAX_ATOMS)
  {
    fprintf(stderr, "np= %d greater than npmax= %d\n", np, MAX_ATOMS);
    exit(-1);
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
    sscanf(buf, "%d %d %lf %lf %lf %d", &dummy, &natom,
	   &atm_num[k].coord.x,
	   &atm_num[k].coord.y, &atm_num[k].coord.z, &itr);
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
    if(sscanf(buf, "%d %f %f %f", &dummy, &atm_num[k].velocity.x,
	   &atm_num[k].velocity.y, &atm_num[k].velocity.z) != 4)
    {
      fprintf(stderr, "Error parsing velocity %d, %s\n", k, buf);
      exit(-1);
    }
#if 0
    atm_num[k].velocity.x *= info->timestep;
    atm_num[k].velocity.y *= info->timestep;
    atm_num[k].velocity.z *= info->timestep;
#endif
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
    if(sscanf(buf, "%d %f %f %f", &dummy, &atm_num[k].accel.x,
	   &atm_num[k].accel.y, &atm_num[k].accel.z) != 4)
    {
      fprintf(stderr, "Error parsing accel %d, %s\n", k, buf);
      exit(-1);
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
    if(sscanf(buf, "%d %f %f %f", &dummy, &atm_num[k].dx3dt3.x,
	   &atm_num[k].dx3dt3.y, &atm_num[k].dx3dt3.z) != 4)
    {
      fprintf(stderr, "Error parsing dx3dt3 %d, %s\n", k, buf);
      exit(-1);
    }
  }
  printf("read %d atoms\n", info->num_atms);
  fclose(fp);
}

static void
append_file(FILE *fp, const BrennerMainInfo *info, double TTIME)
{
      int i;
      int atm;
      int NRA = 0, NLA = 0;
      fprintf(fp, "%s\n", info->head);
      fprintf(fp, "%6d%6d%6d%6d\n", info->num_atms, 3, NRA,NLA);
      fprintf(fp, "%20.11e%20.11e\n", TTIME, info->timestep);
#ifndef INFINITE_CUBE
      fprintf(fp, "%20.11e%20.11e%20.11e\n", info->cube[0],
	      info->cube[1], info->cube[2]);
#else
      fprintf(fp, "%20.11e%20.11e%20.11e\n", 200.0, 200.0, 200.0);
#endif
      for (atm = 0; atm < info->num_atms; atm++)
      {
           fprintf(fp, "%5d%5d%20.11e%20.11e%20.11e\n",
		   atm+1, info->kt2[info->atm_num[atm].ktype],
		   info->atm_num[atm].coord.x,
		   info->atm_num[atm].coord.y, info->atm_num[atm].coord.z);
      }
}

int main(int argc, char **argv)
{
      int lstep;
      int kflag = (argc < 2 ? 1 : atoi(argv[1]));
      FILE *fp;
      BrennerMainInfo *info = alloc_bren(&kflag);
      info->steps = (argc < 3 ? 100000 : atoi(argv[2]));
#ifndef INFINITE_CUBE
      if(argc > 3)
      {
	int cube_dim = set_cell_threshold(atoi(argv[3]));
	info->cube[0] = info->cube[1] = info->cube[2] = cube_dim;
      }
#endif
      init_coord(info);
      init_bren(info, kflag, 0);
      if(kflag == 6)
      {
	update_info(info);
	printf("final energy %f\n", info->system_energy);
	return 0;
      }

      fp = fopen("xmol.d", "w");
      if(!fp)
      {
	fprintf(stderr, "cant write to xmol.d\n");
	exit(-1);
      }

      if(info->nxmol > info->steps && info->steps >= 10)
	info->nxmol = info->steps/10;
      for(lstep = 0; lstep < info->steps; ++lstep)
      {
	  bren_1_step(info, kflag);
	  /* write out position file to be post converted to xmol format */
	  if(!((lstep+1) % info->nxmol))
	    append_file(fp, info, info->starttime + lstep+info->timestep);
	  printf("Time = %.3f fs ", lstep*info->timestep);
	  printf("Total energy = %.5f eV \n", info->system_energy);
      }
      
      return 0;
}
