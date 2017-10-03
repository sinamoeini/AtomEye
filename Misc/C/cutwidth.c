/**********************************************************/
/* Print out file between line number X and line number Y */
/*                                                        */
/* cc -static cutwidth.c -o $BIN_PATH/cutwidth            */
/* $BIN_PATH/cutwidth /etc/passwd 10                      */
/**********************************************************/

#include <stdio.h>

#define MAXCHAR (1024*128)

int main(int argc, char *argv[])
{
    FILE *fp;
    char buffer[MAXCHAR], *p;
    int width;

    if (argc != 3)
    {
        printf ("\nPurpose: print out file in 1 to width.\n\n");
        printf ("Usage: %s filename width\n\n", argv[0]);
        return (1);
    }

    width = atoi(argv[2]);
    if ( (fp = fopen(argv[1],"r")) == NULL )
        err(1, "Cannot read file");

    while (fgets(buffer, MAXCHAR, fp))
    {
        for (p=buffer; *p!=0; p++);
        if (*(p-1)!='\n')
        {
            fprintf (stderr, "MAXCHAR = %d exceeded\n", MAXCHAR);
            return(1);
        }
        else *(p-1) = 0;
        buffer[width]=0;
        puts (buffer);
    }
    fclose(fp);
    fflush(stdout);
    return 0;
}
