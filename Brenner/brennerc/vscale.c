/*
 This code was released into the public domain by Peter McCluskey on 8/18/2000.
 It is adapted from code released into the public domain by Donald Brenner on
 2/29/2000.
*/

#include "brenner.h"

/* scale volume according to potential energy. */

/* I don't know how to adapt this routine to the infinite cube case, so instead
   I'll just chop it out.  Tim Freeman 26 Aug 2000. */

#ifndef INFINITE_CUBE
void vscale(BrennerMainInfo *info)
{
  vector rtemp[2000], all, ey, ckeep;
  int i;
  /* isc is scaling direction */
  int isc = info->volume_scale_dir;
  Float cube2[3];
  Float escale, scale;
  ckeep.x = info->cube[0];
  ckeep.y = info->cube[1];
  ckeep.z = info->cube[2];
  for(i = 0; i < info->num_atms; ++i)
    {
      rtemp[i].x = info->atm_num[i].coord.x;
      rtemp[i].y = info->atm_num[i].coord.y;
      rtemp[i].z = info->atm_num[i].coord.z;
    }
  /*
    kkeep = kflag;
    mkeep = maxkb;
    maxkb = 1;
    kflag = 8;
  */
  find_atom_energies(info);
  
  escale = info->system_energy;
  scale = -0.0010;
  /* change here to control scaling direction */
  while(1)
    {
      cube2[isc] = info->cube[isc]/2.0;
      if(isc == 0)
        {
          info->cube[isc] = ckeep.x*(scale + 1.0);
          for(i = 0; i < info->num_atms; ++i)
            {
              info->atm_num[i].coord.x = rtemp[i].x * (scale + 1.0);
            }
        }
      else if(isc == 1)
        {
          info->cube[isc] = ckeep.y*(scale + 1.0);
          for(i = 0; i < info->num_atms; ++i)
            {
              info->atm_num[i].coord.y = rtemp[i].y * (scale + 1.0);
            }
        }
      else if(isc == 2)
        {
          info->cube[isc] = ckeep.z*(scale + 1.0);
          for(i = 0; i < info->num_atms; ++i)
            {
              info->atm_num[i].coord.z = rtemp[i].z * (scale + 1.0);
            }
        }
      
      find_atom_energies(info);
      
      printf("scale,escale,tote: %f %f %f\n",scale,escale,info->system_energy);
      if(info->system_energy > escale)
        {
          if(scale == -0.0010) /* the Fortran code had scale.eq.0.-0010d0,
                                  which was equivalent to scale == -10;
                                  which was never true */
            scale = -scale;
          else
            {
              scale = -scale/2.0;
              if(fabs(scale) < 0.0002) break;
            }
        }
      else
        {
          escale = info->system_energy;
          ckeep.x = info->cube[0];
          ckeep.y = info->cube[1];
          ckeep.z = info->cube[2];
          for(i = 0; i < info->num_atms; ++i)
            {
              rtemp[i].x = info->atm_num[i].coord.x;
              rtemp[i].y = info->atm_num[i].coord.y;
              rtemp[i].z = info->atm_num[i].coord.z;
            }
        }
    }
  info->cube[0] = ckeep.x;
  info->cube[1] = ckeep.y;
  info->cube[2] = ckeep.z;
  cube2[0] = info->cube[0]/2.0;
  cube2[1] = info->cube[1]/2.0;
  cube2[2] = info->cube[2]/2.0;
  for(i = 0; i < info->num_atms; ++i)
    {
      info->atm_num[i].coord.x = rtemp[i].x;
      info->atm_num[i].coord.y = rtemp[i].y;
      info->atm_num[i].coord.z = rtemp[i].z;
    }
  /* kflag = kkeep; */
  /*    maxkb = mkeep */
  /*      write(42,*) cube(1)*cube(2)*cube(3) , escale */
}
#endif
