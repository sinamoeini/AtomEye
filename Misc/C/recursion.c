/***************/
/* recursion.c */
/***************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <unistd.h>
#include <time.h>
#include <sys/time.h>
#include <sys/types.h>

#define MAXDATA  100
#define MAXLEVEL 100

int Level;
struct List
{
    int size;
    int data[MAXDATA];
} D[] = {
    {2, {1,2}},
    {2, {4,5}},
    {3, {7,8,9}} };

int dumpdata (int *father, int level)
{
    int j,k;
    if (level == Level-1)
        for (j=0; j<D[level].size; j++)
        {
            father[level] = D[level].data[j];
            for (k=0; k<Level; k++)
                printf ("%d ", father[k]);
            printf ("\n");
        }
    else
        for (j=0; j<D[level].size; j++)
        {
            father[level] = D[level].data[j];
            dumpdata (father, level+1);
        }
}

int main (int argc, char *argv[])
{
    int i,j;
    int father[MAXLEVEL];
    Level = sizeof(D) / sizeof(struct List);
    for (i=0; i<Level; i++)
    {
        for (j=0; j<D[i].size; j++)
            printf ("%d ", D[i].data[j]);
        printf ("\n");
    }
    printf ("********\n");
    dumpdata (father, 0);
    return (0);
} /* end main() */


/* cc -c recursion.c; cc recursion.o -lm -o recursion; recursion */
