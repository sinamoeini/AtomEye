/************************************************/
/* C2P Mar 19 1999 <liju99@mit.edu>             */
/* Convert configuration file of Ju Li' format  */
/* to .pdb file which can be viewed by Rasmol.  */
/* cc -o c2p c2p.c; c2p -f config.amorphous.864 */
/************************************************/

#include <stdio.h>
#include <math.h>
#define MAX_STRLEN 128
#define TP_STRLEN 2
#define MAX_TP 65536
#define STATS_FILE_HEADER "stats.c2p"
#define MACHINE_EPS 2.220446e-16
#define TINY (MACHINE_EPS*0)
#define equal(a,b) (fabs((a)-(b))<=TINY)

struct ATOMTYPE
{
    int occurrence; /* number of atoms of such type */
    char name[TP_STRLEN+1]; /* type name */
    double mass; /* atomic mass */
};

/* see if a file exists, and if so ask */
/* the user whether to overwrite it.   */
int freetowrite (char filename[])
{
    FILE *file;
    char c;
    file = fopen (filename, "r");
    if (file != NULL) 
    {
	fclose (file);
	printf ("file \"%s\" exists, overwrite (y/n)? ",
		filename);
        c = getc(stdin);
	if (c!='\n') while(getc(stdin)!='\n');
        if (!((c=='y')||(c=='Y'))) return(0);
    }
    return(1);
} /* end freetowrite() */

double matinv(double A[3][3], double B[3][3])
{
    double D11,D22,D33,D12,D21,D13,D31,D32,D23,determinant;
    D11=A[1][1]*A[2][2]-A[1][2]*A[2][1];
    D22=A[2][2]*A[0][0]-A[2][0]*A[0][2];
    D33=A[0][0]*A[1][1]-A[0][1]*A[1][0];
    D12=A[1][2]*A[2][0]-A[1][0]*A[2][2];
    D23=A[2][0]*A[0][1]-A[2][1]*A[0][0];
    D31=A[0][1]*A[1][2]-A[0][2]*A[1][1];
    D13=A[1][0]*A[2][1]-A[2][0]*A[1][1];
    D21=A[2][1]*A[0][2]-A[0][1]*A[2][2];
    D32=A[0][2]*A[1][0]-A[1][2]*A[0][0];
    determinant=A[0][0]*D11+A[0][1]*D12+A[0][2]*D13;
    B[0][0]=D11/determinant;
    B[1][1]=D22/determinant;
    B[2][2]=D33/determinant;
    B[0][1]=D21/determinant;
    B[1][2]=D32/determinant;
    B[2][0]=D13/determinant;
    B[1][0]=D12/determinant;
    B[2][1]=D23/determinant;
    B[0][2]=D31/determinant;
    return (determinant);
} /* end matinv() */


int read_nc_line (FILE *in)
{
    char c;
    int i, nc = -1;
    c = (char) getc(in);
    if (c == 'n')
	fscanf (in, "c(%d)  = %d\n\n", &i, &nc);
    else if (c != '\n')
	ungetc ((int)c, in);
    return (nc);
} /* end read_nc_line() */


#define max(x,y) ((x)>(y)?(x):(y))

int main (int argc, char *argv[])
{
    int i, j, k, np, nc[3], ntp, to_stdout = 0,
	force_write = 0, save_stats = 0;
    char buf[MAX_STRLEN], read[MAX_STRLEN],
	write[2*MAX_STRLEN], tp[TP_STRLEN+1],
	tps[MAX_TP][TP_STRLEN+1], *r, *p, *q;
    struct ATOMTYPE at[MAX_TP];
    double H[3][3], HI[3][3], mass, s[3], x[3], shift[3]={0.,0.,0.};
    FILE *in, *out, *stats;

    for (k=1; k<max(argc,2); k++)
	if (argc == 1)
	{ /* no input argument */
	    printf ("Convert config file of JL's format: ");
	    for (r=read; r<read+MAX_STRLEN-1; r++)
		if ((*r=(char)getc(stdin))=='\n') break;
	    *r = (char)0;
	}
	else if (argv[k][0] == '-')
	    switch (*(argv[k]+1))
	    {
	    case 'f': /* forced write like % mv -f */
		force_write = 1;
		break;
	    case 'o':
		to_stdout = 1;
		break;	
	    case 's':
		save_stats = 1;
		break;
	    case 'h':
		printf ("-------------------------------------------\n"
			"Convert configuration file of Ju Li' format\n"
			"to .pdb file which can be viewed by Rasmol.\n"
			"-------------------------------------------\n"
			"%% c2p -f config.amorphous config.xtal\n"
			"-f: forced overwrite existent .pdb file\n"
			"-o: send converted .pdb to stdout\n"
			"-s: save stats to \"%s.*\" instead of stdout\n"
			"-x 0.01: shift observer sx by 0.01\n"
			"-y -0.03: shift observer sy by -0.03\n"
			"-z 0.02: shift observer sz by 0.02\n",
			STATS_FILE_HEADER);
		break;
	    case 'x':
		sscanf (argv[++k], "%lf", shift);
		break;
	    case 'y':
		sscanf (argv[++k], "%lf", shift+1);
		break;
	    case 'z':
		sscanf (argv[++k], "%lf", shift+2);
		break;
	    default:
		printf ("c2p: unknown option \"-%c\".\n", *(argv[k]+1));
		exit(1);
	    }
	else sprintf (read, "%s", argv[k]);
    /* my handling of string input is more */
    /* robust than scanf(): */
    for (r=read; (*r==' ')||(*r=='\t'); r++);
    for (q=r; *q!=(char)0; q++)
	if ((*q==' ')||(*q=='\t')) *q=(char)0;
    /* open that configuration file */
    if ( (in=fopen(r,"r")) == NULL )
    {
	printf("c2p: cannot open config file \"%s\".\n", r);
	exit(1);
    }
    fscanf (in, "Number of particles = %d\n\n", &np);
    fscanf (in, "H(1,1) = %lf A\n", &H[0][0]);
    fscanf (in, "H(1,2) = %lf A\n", &H[0][1]);
    fscanf (in, "H(1,3) = %lf A\n", &H[0][2]);
    nc[0] = read_nc_line (in);
    fscanf (in, "H(2,1) = %lf A\n", &H[1][0]);
    fscanf (in, "H(2,2) = %lf A\n", &H[1][1]);
    fscanf (in, "H(2,3) = %lf A\n", &H[1][2]);
    nc[1] = read_nc_line (in);
    fscanf (in, "H(3,1) = %lf A\n", &H[2][0]);
    fscanf (in, "H(3,2) = %lf A\n", &H[2][1]);
    fscanf (in, "H(3,3) = %lf A\n", &H[2][2]);
    nc[2] = read_nc_line (in);
    /* descriptor line */
    fscanf (in, "TP Mass(amu)     sx        sy        sz");
    /* absorb rest of the stuff on this line */
    while (((buf[0]=(char)getc(in))!='\n')&&(buf[0]!=(char)EOF));
    if ( !to_stdout )
    {
	sprintf (write, "%s.pdb", r);
	if (argc == 1)
	{ /* interactive mode */
	    printf("Write to (default=\"%s\"): ", write);
	    q = write + strlen(write);
	    for (p=q; p<write+2*MAX_STRLEN-5; p++)
		if ((*p=(char)getc(stdin))=='\n') break;
	    if (p==q) q = write;
	    *p = (char)0;
	    /* add file suffix if without */
	    if (strcasecmp(p-4,".pdb")) sprintf(p,".pdb");
	}
	else q = write;
	for (; (*q==' ')||(*q=='\t'); q++);
	if ( force_write || freetowrite(q) ) out = fopen (q, "w");
	else exit(1);
    }
    else out = stdout;
    /* .pdb header: */
    fprintf (out, "HEADER    Converted from \"%s\" by c2p.\n", r);
    /* H matrix in footnote -1 and 0: */
    fprintf (out,
	     "FTNOTE   1  Supercell Parallelepiped (H) in Angstroms:\n"
	     "FTNOTE   2  a = (%.4f %.4f %.4f)\n"
	     "FTNOTE   3  b = (%.4f %.4f %.4f)\n"
	     "FTNOTE   4  c = (%.4f %.4f %.4f)\n",
	     H[0][0], H[0][1], H[0][2],
	     H[1][0], H[1][1], H[1][2],
	     H[2][0], H[2][1], H[2][2]);
    /* crystal dimensions: */
    fprintf (out, "CRYST1 %8.3f %8.3f %8.3f "
	     "90.000 90.000 90.000 P 1           1\n",
	     H[0][0], H[1][1], H[2][2]); 
    for (ntp=i=0; i<np; i++)
    {  /* atom types and positions: */
	fscanf (in, "%s %lf %lf %lf %lf",
		tp, &mass, &s[0], &s[1], &s[2]);
	s[0] -= shift[0];
	s[1] -= shift[1];
	s[2] -= shift[2];
	while (s[0]<0) s[0]++;
	while (s[0]>1) s[0]--;
	while (s[1]<0) s[1]++;
	while (s[1]>1) s[1]--;
	while (s[2]<0) s[2]++;
	while (s[2]>1) s[2]--;
	for (j=0; j<ntp; j++)
	    if ( (!strcmp(tp, at[j].name)) &&
		 equal(mass,at[j].mass) )
	    {
		at[j].occurrence++;
		break;
	    }
	if (j==ntp) /* found a new atom type */
	    if (ntp<MAX_TP)
	    {
		sprintf(at[ntp].name, "%s", tp);
		at[ntp].mass = mass;
		at[ntp++].occurrence = 1;
	    }
	    else
	    {
		printf ("v: MAX_TP = %d exceeded.\n", MAX_TP);
		exit(1);
	    }
	while (((buf[0]=(char)getc(in))!='\n')&&(buf[0]!=(char)EOF));
	if (buf[0] == (char)EOF)
	{
	    printf ("c2p: premature end of config file \"%s\":", r);
	    printf ("claimed number of particles = %d, actual = %d\n",
		    np, i+1);
	    goto finish;
	}
	x[0] = s[0]*H[0][0]+s[1]*H[1][0]+s[2]*H[2][0];
	x[1] = s[0]*H[0][1]+s[1]*H[1][1]+s[2]*H[2][1];
	x[2] = s[0]*H[0][2]+s[1]*H[1][2]+s[2]*H[2][2];
	fprintf (out, "ATOM  %5d %2s          %2s    "
		 "%8.3f%8.3f%8.3f\n",
		 i+1, tp, " 1", x[0], x[1], x[2]);
    }
finish: /* print statistics */
    close(in);
    if ( save_stats )
    { /* look for true file name without "/" */
	for (p=r+strlen(r)-1; p>=r; p--)
	    if (*p=='/') break;
	sprintf (buf, "%s.%s", STATS_FILE_HEADER, p+1);
	stats = fopen (buf, "w");
    }
    else stats = stdout;
    if ( (!to_stdout) || save_stats )
    {
	fprintf (stats, "\nIn \"%s\", we found\n", r);
	fprintf (stats, "--------------------------------------\n");
	fprintf (stats, "    Occurrences  Type   Mass\n");
	for (j=0; j<ntp; j++)
	    fprintf (stats, "      %7d     %2s   %.3f\n",
		     at[j].occurrence, at[j].name, at[j].mass);
	    fprintf (stats, "--------------------------------------\n");
	    if (!to_stdout) fprintf (stats, "saved to \"%s\".\n\n", q);
	    close(stats);
	    if (out!=NULL) close(out);
    }
    return(0);
} /* end main() */
