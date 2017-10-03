/**********************************************************/
/* Print out file between line number X and line number Y */
/*                                                        */
/* cc -static headtail.c -o $BIN_PATH/headtail            */
/* $BIN_PATH/headtail /etc/passwd 3 4                     */
/**********************************************************/

#include <stdio.h>

#define MAXCHAR (1024*128)

int main(int argc, char *argv[])
{
    FILE *fp;
    char buffer[MAXCHAR], *p;
    int line, start, finis;
    
    if (argc == 3)
    {
        start = atoi(argv[2]);
        finis = start;
    }
    else if (argc == 4)
    {
        start = atoi(argv[2]);
        finis = atoi(argv[3]);
    }
    else
    {
        printf ("\nPurpose: print out file between "
                "line number X and line number Y.\n\n");
        printf ("Usage: %s filename X\n"
                "       %s filename X Y\n\n", argv[0], argv[0]);
        return (1);
    }

    if ( (fp = fopen(argv[1],"r")) == NULL )
        err(1, "Cannot read file");

    for (line=1; fgets(buffer, MAXCHAR, fp); line++)
    {
        if (line > finis) break;
        if (line >= start)
        {
            for (p=buffer; *p!=0; p++);
            if (*(p-1)!='\n')
            {
                fprintf (stderr, "MAXCHAR = %d exceeded\n", MAXCHAR);
                return(1);
            }
            else *(p-1) = 0;
            puts (buffer);
        }
    }
    fclose(fp);
    return 0;
}
