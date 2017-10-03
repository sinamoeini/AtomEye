/* hc_bren.c - Copyright (c) 1998 Zyvex LLC.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    or its derived works must display the following acknowledgement:
 * 	This product includes software developed by Zyvex LLC.
 * 
 * This software is provided "as is" and any express or implied warranties,
 * including, but not limited to, the implied warranties of merchantability
 * or fitness for any particular purpose are disclaimed. In no event shall
 * Zyvex LLC be liable for any direct, indirect, incidental, special,
 * exemplary, or consequential damages (including, but not limited to,
 * procurement of substitute goods or services; loss of use, data, or
 * profits; or business interruption) however caused and on any theory of
 * liability, whether in contract, strict liability, or tort (including
 * negligence or otherwise) arising in any way out of the use of this
 * software, even if advised of the possibility of such damage.
 */

/*
 Code written by John Michelsen to plug into Hyperchem
 */

#include "brenner.h"
#ifdef USE_HYPERCHEM
#define MAX_MOLECULES 20

typedef struct { double x, y, z; } HC_COORD;

#include "Hc.h"
#include "hsv.h"

void init_hyperchem(void)
{
	if (!LoadHAPI("hapi.dll")) {
		printf("Couldn't find hapi.dll\n");
		exit(1);
	}
	if (!hcConnect("")) {
		printf("Couldn't find HyperChem!\n");
		exit(1);
	}
}

void quit_hyperchem(void)
{
	hcDisconnect();
}

void quit_quickly(char *reason)
{
	printf("Quitting in the middle because '%s'!\n", reason);
	quit_hyperchem();
	exit(1);
}

void offer_cancel(int offer)
{
	hcSetInt(cancel_menu, offer);
}

int should_I_cancel(void)
{
	/* If the menu's still there, the user didn't hit cancel yet. */
	return(!hcGetInt(cancel_menu));
}

HC_COORD hc_coords[MAX_ATOMS];
int hc_types[MAX_ATOMS];
struct {int atom, molecule;} atms_selected[MAX_ATOMS];
int atms_in_molecule[MAX_MOLECULES];

void read_movable(BrenAtom *atm_num, int num_atms)
{
	int num_selected, i, atm;

	for (atm = 0; atm < num_atms; atm++) {
		atm_num[atm].movable = 1;
	}
	num_selected = hcGetInt(selected_atom_count);
	/*printf("There's %d atom[s] selected.\n", num_selected);*/
	hcGetIntVec(selected_atom, (int *)atms_selected, num_selected*2);
	for (i = 0; i < num_selected; i++) {
		atm = atms_selected[i].atom-1;
		/* If it's the 2nd atom in the 3rd molecule, we need to add
			the number of atoms in the 2nd and 1st molecule to get the
			real 'atm' number. */
		while (atms_selected[i].molecule > 1) {
			atms_selected[i].molecule--;
			atm += atms_in_molecule[atms_selected[i].molecule-1];
		}
		atm_num[atm].movable = 0;
	}
}

int read_atom_info(BrenAtom *atm_num)
{
	int atm, atm_total, num_molecules, m;
	int itr = 0;
	int num_atms;
	
	num_molecules = hcGetInt(molecule_count);
	if (num_molecules > MAX_MOLECULES) {
		quit_quickly("Too many molecules");
	}
	hcGetIntVec(atom_count, atms_in_molecule, num_molecules);
	atm_total = 0;
	for (m = 0; m < num_molecules; m++) {
		atm_total += atms_in_molecule[m];
	}
	if (atm_total > MAX_ATOMS) {
		quit_quickly("Too many atoms");
	}
	num_atms=atm_total;      

	hcGetRealArr(coordinates, (double *)hc_coords, atm_total*3);
	hcGetIntArr(atomic_number, hc_types, atm_total);
	for (atm = 0; atm < atm_total; atm++) {
		atm_num[atm].number = atm;		
		if (hc_types[atm]==1)
			atm_num[atm].type = 1;
		else
			atm_num[atm].type = hc_types[atm]*2;    /* type = atomic mass */
		atm_num[atm].coord.x = hc_coords[atm].x;
		atm_num[atm].coord.y = hc_coords[atm].y;
		atm_num[atm].coord.z = hc_coords[atm].z;

		/* natom and itr need to be implemented!! */
		atm_num[atm].ktype = kt[natom];
		++noa[kt[natom]];

		atm_num[atm].mass = xmass[kt[natom]];

		atm_num[atm].thermostated = 0;
		if(itr == 2)
		  atm_num[atm].movable = 0;
		else
		{
		  atm_num[atm].movable = 1;
		  if(itr == 1)
		    atm_num[atm].thermostated = 1;
		}
	}
	read_movable(atm_num, num_atms);
	return num_atms;
}

/* Checks to see whether the set of movable atoms has changed.  Then,
     for each atom that's not movable by this program, it re-loads
	 the coordinates from HyperChem in case the user has moved it.
	 */
void check_for_user_actions(BrenAtom *atm_num, int num_atms)
{
	int atm;

	read_movable(atm_num, num_atms);
	hcGetRealArr(coordinates, (double *)hc_coords, num_atms*3);
	for (atm = 0; atm < num_atms; atm++) {
		if (!atm_num[atm].movable) {
			atm_num[atm].coord.x = hc_coords[atm].x;
		    atm_num[atm].coord.y = hc_coords[atm].y;
		    atm_num[atm].coord.z = hc_coords[atm].z;
		}
	}
}

void write_atom_info(const BrenAtom *atm_num, int num_atms)
{
	int atm;

	for (atm = 0; atm < num_atms; atm++) {
		  hc_coords[atm].x = atm_num[atm].coord.x;
	      hc_coords[atm].y = atm_num[atm].coord.y;
	      hc_coords[atm].z = atm_num[atm].coord.z;
	}
	hcSetRealArr(coordinates, (double *)hc_coords, num_atms*3);
}


#include <unistd.h>
#include <getopt.h>

int
main(int argc, char **argv)
{
	int i, j;
	BrennerMainInfo *info = (BrennerMainInfo*)malloc(sizeof(BrennerMainInfo));

	initbrenner1(info);
	init_hyperchem();

	num_atoms = read_atom_info(atm_num);
	initialize_atoms(info->cube);

	initbrenner2(info);
	offer_cancel(1);
	
	for (j=0; j < info->steps; j++) {
		perform1brennerstep(info);
		
		check_for_user_actions(atm_num, num_atms);
		write_atom_info(atm_num, num_atms);
		printf("Time = %.3f ps ", j*TIMESTEP);
		printf("Total energy = %.5f eV \n", system_energy);
		if (should_I_cancel()) {
			quit_quickly("You told me to");
		}
	}
	offer_cancel(0);

	quit_hyperchem();
	printf("%d steps %d atoms most in a box %d\n",
	       j, num_atms, most_atoms_in_a_box);
	free(info);
	return 0;
}

void my_exit(int code)
{
  exit(code);
}

#endif /* USE_HYPERCHEM */
