/************************************************************/
/* Fetch a numbered line (plus following lines) from a file */
/************************************************************/
/* cc -o lth lth.c */
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char*argv[])
{
    FILE *file;
    int c, i, l, m;
    if (argc<3)
    {
	printf ("Usage: lth <filename> <line-number> [+follow]\n");
	return(1);
    }
    file = fopen(argv[1], "r");
    if (file == NULL)
    {
	printf ("file \"%s\" does not exist.\n", argv[1]);
	return(1);
    }
    l = atoi(argv[2]);
    if (argc==4) m = atoi(argv[3]); else m=1;
    for (i=1; i<l; i++)
    {
	while (((c=getc(file))!=EOF)&&(c!='\n'));
	if (c==EOF)
	{
	    printf ("it seems that file \"%s\" does not contain %d lines, try\n",
		    argv[1], l);
	    printf ("%% wc -l %s\n", argv[1]);
	    return(1);
	}
    }
    for (i=1; i<=m; i++)
    {
	while (((c=getc(file))!=EOF)&&(c!='\n'))
	    putc (c, stdout);
	if (c==EOF) return(1);
	printf ("\n");
    }
    fclose(file);
    return (0);
}
