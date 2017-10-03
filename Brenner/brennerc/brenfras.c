/*
 This code was released into the public domain by Peter McCluskey on 6/2/2000.
*/

/*
 * This is a sample program to show how to combine the Fortran Brenner
 * potential library with a RasMol display.
 */

#include "brenfort.h"
#include <unistd.h>
#include <getopt.h>
#include <rasapi.h>
#include <tokens.h>

#define max(x,y) ((x) > (y) ? (x) : (y))

static Molecule __far *Database;

static void read_atom_info(char *filename, BrennerMainInfo *info, double);
static void RefreshRasmolDisplay(BrennerMainInfo *info);

#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(aptr=group->alist;aptr;aptr=aptr->anext)


static void
setDisplayMode(ViewStruct *bitmap)
{
  int tokens[32];
  double dtokens[2];
  bitmap->display_mode = CPKTok;
  bitmap->display_arg1 = NumberTok;
  bitmap->display_arg2 = 0.475;

  SetRadiusValue(120,Database,bitmap);
  EnableWireframe(CylinderFlag,40,Database, bitmap);
  SetRibbonStatus(False,0,0,Database);
  DisableBackbone(Database);

  tokens[0] = bitmap->display_mode;
  tokens[1] = bitmap->display_arg1;
  tokens[2] = 0;
  dtokens[0] = bitmap->display_arg2;
  dtokens[1] = -1;
  ExecuteParsedCommand(Database, bitmap, tokens, dtokens);
  ReDrawFlag |= RFColour;
}

static void
read_atom_info(char *filename, BrennerMainInfo *info, double cube_init)
{
	register Chain __far *chain;
	register Group __far *group;
	register Atom __far *aptr;
	extern Long OrigCX, OrigCY, OrigCZ;
	BrenAtom *atm_num = info->atm_num;
	int auto_cube = (cube_init == 0);
	const double tiny_a = 1.e-13;
	const double tiny_v = 1.e-4;

	Database = CreateMolGroup();
	RasMolInitFile(Database, NULL, filename);

	ForEachAtom {
		int atm = info->num_atms;
		int natom;
		double t;
		if(info->num_atms >= MAX_ATOMS)
		{
		  fprintf(stderr, "too many atoms %d\n", info->num_atms);
		  my_exit(-1);
		}
		atm_num[atm].number = atm;
		natom = aptr->elemno;
		atm_num[atm].coord.x = (aptr->xorg + OrigCX) * 4 / 1000.0;
		atm_num[atm].coord.y = (aptr->yorg + OrigCY) * 4 / 1000.0;
		atm_num[atm].coord.z = (aptr->zorg + OrigCZ) * -4 / 1000.0;

		atm_num[atm].velocity.x = (rand()/(double)RAND_MAX)*2*tiny_v-tiny_v;
		atm_num[atm].velocity.y = (rand()/(double)RAND_MAX)*2*tiny_v-tiny_v;
		atm_num[atm].velocity.z = (rand()/(double)RAND_MAX)*2*tiny_v-tiny_v;
		atm_num[atm].accel.x = (rand()/(double)RAND_MAX)*2*tiny_a-tiny_a;
		atm_num[atm].accel.y = (rand()/(double)RAND_MAX)*2*tiny_a-tiny_a;
		atm_num[atm].accel.z = (rand()/(double)RAND_MAX)*2*tiny_a-tiny_a;
		
		atm_num[atm].dx3dt3.x = (rand()/(double)RAND_MAX)*2*tiny_a-tiny_a;
		atm_num[atm].dx3dt3.y = (rand()/(double)RAND_MAX)*2*tiny_a-tiny_a;
		atm_num[atm].dx3dt3.z = (rand()/(double)RAND_MAX)*2*tiny_a-tiny_a;

		t = max(atm_num[atm].coord.x, atm_num[atm].coord.y);
		t = max(t, atm_num[atm].coord.z);
		cube_init = max(cube_init, 2*t);

		atm_num[atm].ktype = info->kt[natom];
		if(!atm_num[atm].ktype)
		  fprintf(stderr, "warning - unknown type %d, atm number %d\n",
			  natom, atm);
		++info->noa[info->kt[natom]];

		atm_num[atm].type = natom;
		atm_num[atm].mass = info->xmass[info->kt[natom]];

		atm_num[atm].movable = 1;
		atm_num[atm].thermostated = 1;

		++info->num_atms;
	}
	info->cube[0] = info->cube[1] = info->cube[2] = cube_init;
#if 1		/* disable this for wireframe display */
	setDisplayMode(GetScreenView());
#endif
}

static void
RefreshRasmolDisplay(BrennerMainInfo *info)
{
    ViewStruct *bitmap = GetScreenView();
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    extern Long OrigCX, OrigCY, OrigCZ;
    int i = 0;
    int detect_bonds = False;
    BrenAtom *atm_num = info->atm_num;

    /*    ReDrawFlag |= RFInitial; */
    ForEachAtom {
	aptr->xorg = (int)(atm_num[i].coord.x * 1000 / 4) - OrigCX;
	aptr->yorg = (int)(atm_num[i].coord.y * 1000 / 4) - OrigCY;
	aptr->zorg = (int)(atm_num[i].coord.z * 1000 / -4) - OrigCZ;
	++i;
    }
#if 0
    if( InfoBondCount < (MainAtomCount+HetaAtomCount)-InfoChainCount
	&& detect_bonds)
    {   if( MainAtomCount+HetaAtomCount > 255 )
        {   CreateMoleculeBonds(False,False,Database);
        } else CreateMoleculeBonds(False,True,Database);
    }

    /* Explicit Hydrogen Bonds! */
    if( InfoHBondCount > 0 )
        SetHBondStatus(True,True,0,Database);

    InitialTransform(Database);
    EnableWireframe(0,0,Database, bitmap); /* restore prior values */
#endif

    VoxelsClean = False;
    ApplyTransform(Database, bitmap);

    ReDrawFlag |= RFRefresh | RFApply;
    if(0 /*!update_only */)
    {
        SetRadiusValue(120,Database,bitmap);
	EnableWireframe(CylinderFlag,40,Database, bitmap);
	SetRibbonStatus(False,0,0,Database);
	DisableBackbone(Database);
    }
    /*    else repeatDisplayMode(bitmap); */
    if(0) CPKColourAttrib(Database);
    RefreshScreen(Database);
}

int
main(int argc, char **argv)
{
	int i, j;
	int kflag = 1;
	const char *coord_file = NULL;
	double delta = 0.5;
	double cube_init = 0;
	BrennerMainInfo *info = alloc_bren(&kflag);

	while ((i = (int) getopt(argc, argv, "?k:d:s:t:")) != -1)
	{
	  extern char *optarg;
	  switch(i)
	  {
	  case '?':
	    printf("useage: %s [-s <#steps>] [-t <temperature(K)>] [-k kflag] [-d <size of cube>] input_file", argv[0]);
	    exit(-1);
	  case 'd':
	    cube_init = atof(optarg);
	    break;
	  case 'k':
	    kflag = atoi(optarg);
	    break;
	  case 's':
	    info->steps = atoi(optarg);
	    if (info->steps == 0) {
		printf("Warning - you asked for no steps\n");
	    }
	    break;
	  case 't':
	    info->temperature = atof(optarg);
	    break;
	  }
	}

	if(argc < 2)
	{
	  fprintf(stderr, "error - no filename specified\n");
	  exit(-1);
	}
	read_atom_info(argv[argc-1], info, cube_init);

	init_bren(info, kflag, 0);
	for (j = 0; j < info->steps; ++j)
	{
		bren_1_step(info, kflag);
		
		if(!(j % 10))
		  RefreshRasmolDisplay(info);
		printf("Time = %.3f fs ", j*info->timestep);
		printf("Total energy = %.5f eV \n", info->system_energy);
	}
	free(info);
	RasMolExit();
	return 0;
}

void my_exit(int code)
{
  RasMolExit();
  exit(code);
}
