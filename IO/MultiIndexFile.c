/****************************************/
/* libIO:                               */
/*                                      */
/* Strings, parsers & files             */
/*                                      */
/* Dec.12, 1999  Ju Li <liju99@mit.edu> */
/****************************************/

#include "IO.h"

/*******************************************************************/
/* Multi-Index file support: a multi-index file starts with a line */
/* of integers like "3 1 2", followed by one-integer lines like    */
/* 789                                                             */
/* 38                                                              */
/* 9800                                                            */
/* 560                                                             */
/* 2345                                                            */
/* 456                                                             */
/* of which the first three integers belong to the first category, */
/* 560 belongs to the second category, and 2345, 456 belong to the */
/* last. Usually these integers are the indices of a certain array */
/* that will undergo some (here, three kinds of) operations.       */
/*******************************************************************/


void MultiIndexFree (MultiIndex *m)
{
    if ( ISNULL(m) ) return;
    Free (m->size);
    Free (m->cumulant);
    Free (m->index);
    m->ncat = 0;
    return;
} /* end MultiIndexFree() */


/* return total number of indices read */
int MultiIndexRead (char *filename, MultiIndex *m)
{
    FILE *fp;
    int i, j;
    TermString separator;
    fp = rOpen(filename);
    MultiIndexFree (m);
    REALLOC_ZERO (MultiIndexRead, m->cumulant, 1, int);
    while ( fscanf(fp, "%d%[^0-9]", &j, separator) != EOF )
    {
        m->ncat++;
        REALLOC (MultiIndexRead, m->size, m->ncat, int);
        m->size[m->ncat-1] = j;
        REALLOC (MultiIndexRead, m->cumulant, m->ncat+1, int);
        m->cumulant[m->ncat] = m->cumulant[m->ncat-1] + m->size[m->ncat-1];
        if (strchr(separator, '\n')) break;
    }
    REALLOC (MultiIndexRead, m->index, m->cumulant[m->ncat], int);
    for (i=0; i<m->ncat; i++)
        for (j=m->cumulant[i]; j<m->cumulant[i+1]; j++)
            if ( fscanf(fp, "%d", m->index+j) == EOF )
                pe ("%s ends prematurely, it does not contain %d indices.\n",
                    filename, m->cumulant[m->ncat]);
    return(m->cumulant[m->ncat]);
} /* end MultiIndexRead() */


void MultiIndexSave (MultiIndex *m, char *filename)
{
    FILE *fp;
    int i, j;
    fp = wOpen(filename);
    if (m->ncat == 0) return;
    fprintf (fp, "%d", m->size[0]);
    for (i=1; i<m->ncat; i++)
        fprintf (fp, " %d", m->size[i]);
    fprintf (fp, "\n");
    for (i=0; i<m->ncat; i++)
        for (j=m->cumulant[i]; j<m->cumulant[i+1]; j++)
            fprintf (fp, "%d\n", m->index[j]);
    return;
} /* end MultiIndexSave() */


#ifdef _MultiIndexRead_TEST
/*
cat<<EOF>/tmp/a
 3    1 2  
789
38
9800
560
2345
456
EOF
*/
int main (int argc, char *argv[])
{
    MultiIndex m[1]={0};
    MultiIndexRead ("/tmp/a", m);
    MultiIndexSave (m, "/tmp/b");
    return (0);
}
#endif /* _MultiIndexRead_TEST */
