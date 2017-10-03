/****************************************/
/*  Grab a portion of a current file    */
/*  and save plots in a Matlab script:  */
/*  % grab <filename> <length> <start>  */
/****************************************/

#include <stdio.h>
#include <string.h>

#define GM_FILE         "/home/Moon/g.m"
#define MAX_STRLEN      128
#define DEFAULT_LENGTH  (1024*24)
#define DEFAULT_START   0
#define DELTA_FIRST 

int main (int argc, char *argv[])
{
    int i;
    unsigned long size, length, start;
    double delta, J;
    char fn[MAX_STRLEN]={0}, fo[MAX_STRLEN]={0};
    FILE *current, *out;
    
    strncpy (fn, (argc<2)?"jx.out":argv[1], MAX_STRLEN-1);
    length = (argc<3)?DEFAULT_LENGTH:atoi(argv[2]);
    start  = (argc<4)?DEFAULT_START:atoi(argv[3]);
    printf ("J file = \"%s\",\n", fn);
    printf ("start = %ld, length = %ld.\n", start, length);
    while ((current=fopen(fn,"r")) == NULL)
    {
	printf ("cannot open \"%s\", input another: ", fn);
	scanf ("%s\n", fn);
    }
    fseek (current, 0L, SEEK_END);
    size = ftell(current)/sizeof(double);
#ifdef DELTA_FIRST
    /* the first double precision number is delta in ps */
    fseek (current, 0L, SEEK_SET);
    fread (&delta, sizeof(double), 1, current);
    /* original start is based the first J data */
    start++;
#else
    /* the last double precision number is delta in ps */
    fseek (current, -sizeof(double), SEEK_END);
    fread (&delta, sizeof(double), 1, current);
#endif    
    printf ("delta = %.4f ps, J data size = %ld.\n",
	    delta, size-1);
    if (start+length>size-1) length = size-1-start;
    fseek (current, sizeof(double)*start, SEEK_SET);
    
    strncpy (fo, GM_FILE, MAX_STRLEN-1);
    while ((out=fopen(fo,"w+")) == NULL)
    {
	printf ("invalid output file \"%s\", input another: ", fo);
	scanf ("%s\n", fo);
    }
    fprintf (out, "clf; delta = %f; J = [\n", delta);
    for (i=0; i<length; i++)
    {
	fread (&J, sizeof(double), 1, current);
	fprintf (out, "%f\n", J);
    }
    fprintf(out, "]; u=[%d %d]*delta; t=(%d:%d)*delta;\n",
	    start, start+length-1, start, start+length-1);
    fprintf(out, "plot(t,J); hold on; plot(u, [0 0],'r--');\n");
    fprintf(out, "xlabel('t [ps]'); ylabel('J [arbitrary unit]');\n");
    fprintf(out, "v = axis; axis([u v(3) v(4)]);\n");
    fprintf(out, "tt=0; ss = input('Average over interval [ps]: ', 's');\n");
    fprintf(out, "while ~isempty(ss)\n");
    fprintf(out, " filt_ps = str2num(ss);\n");
    fprintf(out, " n = ceil(filt_ps/delta);\n");
    fprintf(out, " a = filtfilt(ones(1,n)/n, [1], J);\n");
    fprintf(out, " if (tt==0) tt=1; \n");
    fprintf(out, " else clf;\n");
    fprintf(out, "  plot(t,J); hold on; plot(u, [0 0],'r--');\n");
    fprintf(out, "  xlabel('t [ps]'); ylabel('J [arbitrary unit]');\n");
    fprintf(out, "  axis([u v(3) v(4)]);\n");
    fprintf(out, " end;\n");
    fprintf(out, " plot (t, a, 'w');\n");
    fprintf(out, " title (['Smoother line is J''s running average "
	    "over \\Delta = ' num2str(filt_ps) ' ps (\\delta = ' "
	    "num2str(delta) ' ps).']);\n");
    fprintf(out, " ss = input('Average over interval [ps]: ', 's');\n");
    fprintf(out, "end;\n");
    fclose(out);
    printf ("%.2f ps of current saved on \"%s\".\n",
	    length*delta, fo);
    return(0);
} /* end main() */

/* % cc -o grab grab.c;  grab */
