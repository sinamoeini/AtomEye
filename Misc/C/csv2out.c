/*************/
/* csv2out.c */
/*************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define HEADERLINES 1
/* number of lines at the beginning of the csv file to skip */
#define CUSIP_FNAME "CUSIP.out"
/* row vector of company CUSIP strings */
#define DATES_FNAME "DATES.out"
/* column vector of date strings */
#define MAXLINESIZE (1024*16)
/* maximum number of chars in a string, including the terminating 0 */
#define MAXNUMCUSIP (1024*16)
/* maximum number of equities in the file */
#define MAXNUMDATES (1024*16)
/* maximum number of dates in the file */
#define MAXNUMPROPERTY 256
/* maximum number of property entries in one line */
#define PROG_REPORT_FREQ  100000
/* report progress of loading entries from csv file */
#define DOUBLE_PRECISION_INFINITY (1e300)
/* impossible number to represent NaN */
#define CUSIP_RADIX  43
/* 'Z' - '0' + 1 */
#define CUSIP_BASE   3418801
/* 43^4 = 3418801 */
#define EPS_REMAINING  (5e-2)
/* Accuracy left in the double precision after */
/* representing last digit of 3418801.3418801  */
#define TOLERANCE    (EPS_REMAINING / CUSIP_BASE)


static int quadchar_to_int (char a3, char a2, char a1, char a0)
{
    return( (int)
            (a0-'0') +
            CUSIP_RADIX *
            ( (a1-'0') +
              CUSIP_RADIX *
              ( (a2-'0') +
                CUSIP_RADIX *
                (a3-'0') ) ) );
} /* end quadchar_to_int() */


static void int_to_quadchar (int id, char *a3, char *a2, char *a1, char *a0)
{
    int residue;
    residue = id % CUSIP_RADIX;
    *a0 = residue + '0';
    id = (id - residue) / CUSIP_RADIX;
    residue = id % CUSIP_RADIX;
    *a1 = residue + '0';
    id = (id - residue) / CUSIP_RADIX;
    residue = id % CUSIP_RADIX;
    *a2 = residue + '0';
    id = (id - residue) / CUSIP_RADIX;
    residue = id % CUSIP_RADIX;
    *a3 = residue + '0';
    return;
} /* end int_to_quadchar() */


static double twoint_to_double (int i1, int i2)
{
    return( i1 + (double) i2 / CUSIP_BASE );
} /* end twoint_to_double() */


static void double_to_twoint (double id, int *i1, int *i2)
{
    *i1 = (int) id;
    *i2 = (int) (CUSIP_BASE * (id - *i1) + EPS_REMAINING);
    return;
} /* end double_to_twoint() */


/* "01/31/1995" -> unique double */
static double DATESstring_to_double (char *DATESstring)
{
    return( twoint_to_double (
        quadchar_to_int ( DATESstring[6],
                          DATESstring[7], +
                          DATESstring[8],
                          DATESstring[9] ),
        quadchar_to_int ( DATESstring[0],
                          DATESstring[1], +
                          DATESstring[3],
                          DATESstring[4] ) ) );
} /* end DATESstring_to_double() */


/* unique double -> "01/31/1995" */
static char *double_to_DATESstring (double id)
{
    static char DATESstring[11]={'0','0','/','0','0','/','0','0','0','0',0};
    static int i1, i2;
    double_to_twoint (id, &i1, &i2);
    int_to_quadchar( i1,
                     &DATESstring[6],
                     &DATESstring[7], 
                     &DATESstring[8],
                     &DATESstring[9] );
    int_to_quadchar( i2,
                     &DATESstring[0],
                     &DATESstring[1], 
                     &DATESstring[3],
                     &DATESstring[4] );
    return( DATESstring );
} /* end double_to_DATESstring() */


#ifdef _DATESstring_to_double
#define TRY_STRING "01/31/1995"
int main (int argc, char *argv[])
{
    double id = DATESstring_to_double(TRY_STRING);
    printf ("\"%s\" -> %.15g\n", TRY_STRING, id);
    printf ("%.15g -> \"%s\"\n", id, double_to_DATESstring(id));
    return (0);
}
/* cc -D_DATESstring_to_double -D_NO_main csv2out.c -o csv2out; csv2out */
#endif /* _DATESstring_to_double */


/* "40049W30" -> unique double */
static double CUSIPstring_to_double (char *CUSIPstring)
{
    return( twoint_to_double (
        quadchar_to_int ( CUSIPstring[0],
                          CUSIPstring[1], 
                          CUSIPstring[2],
                          CUSIPstring[3] ),
        quadchar_to_int ( CUSIPstring[4],
                          CUSIPstring[5], 
                          CUSIPstring[6],
                          CUSIPstring[7] ) ) );
} /* end CUSIPstring_to_double() */


/* unique double -> "40049W30" */
static char *double_to_CUSIPstring (double id)
{
    static char CUSIPstring[9]={0};
    static int i1, i2;
    double_to_twoint (id, &i1, &i2);
    int_to_quadchar( i1,
                     &CUSIPstring[0],
                     &CUSIPstring[1], 
                     &CUSIPstring[2],
                     &CUSIPstring[3] );
    int_to_quadchar( i2,
                     &CUSIPstring[4],
                     &CUSIPstring[5], 
                     &CUSIPstring[6],
                     &CUSIPstring[7] );
    return( CUSIPstring );
} /* end double_to_CUSIPstring() */

#ifdef _CUSIPstring_to_double
#define TRY_STRING "WZADGHD9"
int main (int argc, char *argv[])
{
    double id = CUSIPstring_to_double(TRY_STRING);
    printf ("\"%s\" -> %.15g\n", TRY_STRING, id);
    printf ("%.15g -> \"%s\"\n", id, double_to_CUSIPstring(id));
    return (0);
}
/* cc -D_CUSIPstring_to_double -D_NO_main csv2out.c -o csv2out; csv2out */
#endif /* _CUSIPstring_to_double */


#define SWAP(x,y,tmp) ((tmp)=(x), (x)=(y), (y)=(tmp))
/* qsort_num_recipes() */
#define QUICKSORT_M  7
static void qsort_num_recipes (int N, double x[], int idx[])
{
    register int i, j, tmp;
    int idxt, ir=N, k, l=1, jstack=0;
    static int istack[128];
    register double a;
    for (i=0; i<N; i++) idx[i] = i;
    idx--;
    for (;;)
    {
	if (ir-l < QUICKSORT_M)
	{
	    for (j=l+1; j<=ir; j++)
	    {
		idxt = idx[j];
		a = x[idxt];
		for (i=j-1; i>=l; i--)
		{
		    if (x[idx[i]] <= a) break;
		    idx[i+1]=idx[i];
		}
		idx[i+1]=idxt;
	    }
	    if (jstack == 0) break;
	    ir = istack[jstack--];
	    l = istack[jstack--];
	}
	else
	{
	    k = (l+ir) >> 1;
	    SWAP(idx[k],idx[l+1],tmp);
	    if (x[idx[l]] > x[idx[ir]]) SWAP(idx[l],idx[ir],tmp);
	    if (x[idx[l+1]] > x[idx[ir]]) SWAP(idx[l+1],idx[ir],tmp);
	    if (x[idx[l]] > x[idx[l+1]]) SWAP(idx[l],idx[l+1],tmp);
	    i = l+1;
	    j = ir;
	    idxt = idx[l+1];
	    a = x[idxt];
	    for (;;)
	    {
		do i++; while (x[idx[i]] < a);
		do j--; while (x[idx[j]] > a);
		if (j < i) break;
		SWAP (idx[i],idx[j],tmp);
	    }
	    idx[l+1]=idx[j];
	    idx[j]=idxt;
	    jstack += 2;
	    if (ir-i+1 >= j-l)
	    {
		istack[jstack] = ir;
		istack[jstack-1] = i;
		ir = j-1;
	    }
	    else
	    {
		istack[jstack] = j-1;
		istack[jstack-1] = l;
		l = i;
	    }
	}
    }
}
#undef QUICKSORT_M
/* end qsort_num_recipes() */

#ifdef _qsort_num_recipes
int main (int argc, char *argv[])
{
    double x[] = { 10, 10, 4, 7, 10, 54354, 34 };
    int i, N = sizeof(x) / sizeof(double);
    int *idx = (int *) malloc(N * sizeof(int));
    qsort_num_recipes (N, x, idx);
    for (i=0; i<N; i++) printf ("%.15g\n", x[idx[i]]);
    free (idx);
    return (0);
}
/* cc -D_qsort_num_recipes -D_NO_main csv2out.c -o csv2out; csv2out */
#endif /* _qsort_num_recipes */


/* Sort stack from small to big and delete identical */
/* entries; return new stack length. */
int reorganize (int num, double *stack)
{
    int i, newnum, *idx;
    double *buf;
    idx = (int *)  malloc(num * sizeof(int));
    buf = (double *) malloc(num * sizeof(double));
    qsort_num_recipes (num, stack, idx);
    for (i=0; i<num; i++) buf[i] = stack[idx[i]];
    free (idx);
    stack[0] = buf[0];
    for (newnum=1,i=1; i<num; i++)
        if ( buf[i] != stack[newnum-1] )
            stack[newnum++] = buf[i];
    free (buf);
    return (newnum);
} /* end reorganize() */

#ifdef _reorganize
int main (int argc, char *argv[])
{
    double x[] = { 10, 10, 4, 2, 3, 7, 4, 10, 54354, 34, 0 };
    int i, N = sizeof(x) / sizeof(double);
    N = reorganize (N, x);
    for (i=0; i<N; i++) printf ("%.15g\n", x[i]);
    return (0);
}
/* cc -D_reorganize -D_NO_main csv2out.c -o csv2out; csv2out */
#endif /* _reorganize */


/* locate index of value in stack[], which is sorted and non-repeating */
int find_in_sorted_stack (double value, int num, double *stack)
{
    int i=0, j=num-1, k=(i+j)/2;
    while ( (value > stack[k] + TOLERANCE) ||
            (value < stack[k] - TOLERANCE) )
    {
        if (value < stack[k]) j = k-1;
        else i = k+1;
        k = (i+j) / 2;
    }
    return (k);
} /* end find_in_sorted_stack() */

#ifdef _find_in_sorted_stack
#define GUESS 8
int main (int argc, char *argv[])
{
    double x[] = { 4, 7, 8, 9, 10, 34, 543, 54324, 2454255 };
    printf ("%d %d\n", GUESS,
            find_in_sorted_stack (x[GUESS], sizeof(x) / sizeof(double), x));
    return (0);
}
/* cc -D_find_in_sorted_stack -D_NO_main csv2out.c -o csv2out; csv2out */
#endif /* _find_in_sorted_stack */


static int num_property;
static char *property_names[MAXNUMPROPERTY];

/* fill in property_names[] */
void get_property_names (char *csv_fname)
{
    int i;
    FILE *fp_in;
    char buffer[MAXLINESIZE], *p, *q, c;
    if ( (fp_in = fopen(csv_fname,"r")) == NULL )
    {
        printf ("get_property_names: csv file name \"%s\" is invalid.\n",
                csv_fname);
        exit(1);
    }
    for (i=0; i<HEADERLINES; i++) fgets(buffer, MAXLINESIZE, fp_in);
    printf ("\n%d header lines read in.\n", HEADERLINES);
    for (p=buffer-1,i=0; i<7; i++) for (p++; *p!=','; p++);
    num_property = 0;
    while (1)
    {
        for (q=p+1; (*q!=',') && (*q!='\n') && (*q!=0); q++);
        c = *q;
        *q = 0;
        if (*(q-1)=='\015') *(q-1) = 0;
        property_names[num_property] = (char *) malloc(q-p);
        strcpy (property_names[num_property], p+1);
        p = q;
        num_property++;
        if (num_property == MAXNUMPROPERTY)
        {
            printf ("get_property_names: overflow -\n"
                    "please increase MAXNUMPROPERTY and recompile.\n");
            exit(1);
        }
        if (c != ',') break;
    }
    fclose (fp_in);
    return;
} /* end get_property_names() */
/* cc -O3 csv2out.c -o csv2out; csv2out *.csv */


#ifndef _NO_main
int main (int argc, char *argv[])
{
    int i, j, property_idx;
    char *csv_fname, property_fname[MAXLINESIZE], buffer[MAXLINESIZE];
    char *p, *q;
    FILE *fp_in, *fp_out;
    double CUSIP_stack[MAXNUMCUSIP], DATES_stack[MAXNUMDATES], CUSIP, DATES;
    int num_CUSIP = 0, num_DATES = 0, num_entry_lines = 0, num_entries = 0;
    double *property, value;

    if (argc == 2)
    {
        get_property_names (argv[1]);
        printf ("\nproperty_idx\tproperty_name\n"
                "------------    -------------\n");
        for (i=0; i<num_property; i++)
            printf ("     %d   \t%s\n", i+1, property_names[i]);
        printf ("\n");
        return (1);
    }
    else if (argc != 3) 
    {
        printf ("\nPurpose: parse a csv file and save a certain property\n"
                "         as a matrix in which each column is time history\n"
                "         of a certain company. The company designations\n"
                "         will be stored in \"%s\", the time strings\n"
                "         will be stored in \"%s\". Inside the matrix\n"
                "         file, unavailable entries are stored as NaN.\n\n",
                CUSIP_FNAME, DATES_FNAME);
        printf ("Usage: First run\n\n"
                "> %s csv_file_name\n\n"
                "to list property indices. Then run,\n\n"
                "> %s csv_file_name property_idx\n\n"
                "to get the specific property matrix.\n\n",
                argv[0], argv[0]);
        return (1);
    }
    csv_fname = argv[1];
    get_property_names (csv_fname);
    property_idx = atoi(argv[2]);
    if ( (property_idx < 1) || (property_idx > num_property ) )
    {
        printf ("property_idx must be between 1 and %d, it is now %d.\n",
                num_property, property_idx);
        return (1);
    }
    if ( (fp_in = fopen(csv_fname,"r")) == NULL )
    {
        printf ("csv file name \"%s\" is invalid.\n", csv_fname);
        return (1);
    }

    printf ("\nFirst sweep:\n");
    for (i=0; i<HEADERLINES; i++) fgets(buffer, MAXLINESIZE, fp_in);

    while (fgets(buffer, MAXLINESIZE, fp_in) != NULL)
    {
        for (p=buffer-1,i=0; i<2; i++)
        {
            for (p++; *p!=','; p++);
            if (p[1]=='"') for (p=p+2; *p!='"'; p++);
        }
        p[9] = 0;
        CUSIP_stack[num_CUSIP++] = CUSIPstring_to_double(p+1);
        p += 9;
        for (i=0; i<3; i++)
        {
            for (p++; *p!=','; p++);
            if (p[1]=='"') for (p=p+2; *p!='"'; p++);
        }
        p[11] = 0;
        DATES_stack[num_DATES++] = DATESstring_to_double(p+1);

        if (num_CUSIP == MAXNUMCUSIP)
            num_CUSIP = reorganize (num_CUSIP, CUSIP_stack);
        if (num_CUSIP == MAXNUMCUSIP)
        {
            printf ("overflow - please increase MAXNUMCUSIP and recompile.\n");
            return(1);
        }

        if (num_DATES == MAXNUMDATES)
            num_DATES = reorganize (num_DATES, DATES_stack);
        if (num_DATES == MAXNUMDATES)
        {
            printf ("overflow - please increase MAXNUMDATES and recompile.\n");
            return(1);
        }

        num_entry_lines++;
        if (num_entry_lines % PROG_REPORT_FREQ == 0)
            printf ("%d entry lines read in.\n", num_entry_lines);
    }

    num_CUSIP = reorganize (num_CUSIP, CUSIP_stack);
    if ( (fp_out = fopen(CUSIP_FNAME,"w")) == NULL )
    {
        printf ("Output file \"%s\" cannot be opened.\n", CUSIP_FNAME);
        return (1);
    }
    for (i=0; i<num_CUSIP; i++)
        fprintf (fp_out, "%s ", double_to_CUSIPstring(CUSIP_stack[i]));
    fprintf (fp_out, "\n");
    fclose (fp_out);
    printf ("\n%d sorted CUSIPs saved on \"%s\"\n", num_CUSIP, CUSIP_FNAME);

    num_DATES = reorganize (num_DATES, DATES_stack);
    if ( (fp_out = fopen(DATES_FNAME,"w")) == NULL )
    {
        printf ("Output file \"%s\" cannot be opened.\n", DATES_FNAME);
        return (1);
    }
    for (i=0; i<num_DATES; i++)
        fprintf (fp_out, "%s\n", double_to_DATESstring(DATES_stack[i]));
    fclose(fp_out);
    printf ("%d sorted dates saved on \"%s\"\n\n", num_DATES, DATES_FNAME);

    if ( (property = (double *)malloc(num_DATES*num_CUSIP*sizeof(double)))
         == NULL )
    {
        printf ("Cannot allocate %d x %d matrix, please free some memory.\n",
                num_DATES, num_CUSIP);
        return (1);
    }
    for (i=0; i<num_DATES*num_CUSIP; i++)
        property[i] = DOUBLE_PRECISION_INFINITY;

    rewind (fp_in);
    property_idx--;
    sprintf (property_fname, "%s.out", property_names[property_idx]);
    if ( (fp_out = fopen(property_fname,"w")) == NULL )
    {
        printf ("Output file \"%s\" cannot be opened.\n", property_fname);
        return (1);
    }

    printf ("Second sweep:\n");
    for (i=0; i<HEADERLINES; i++) fgets(buffer, MAXLINESIZE, fp_in);
    num_entry_lines = num_entries = 0;
    while (fgets(buffer, MAXLINESIZE, fp_in) != NULL)
    {
        for (p=buffer-1,i=0; i<2; i++)
        {
            for (p++; *p!=','; p++);
            if (p[1]=='"') for (p=p+2; *p!='"'; p++);
        }
        p[9] = 0;
        CUSIP = CUSIPstring_to_double(p+1);
        p += 9;
        for (i=0; i<3; i++)
        {
            for (p++; *p!=','; p++);
            if (p[1]=='"') for (p=p+2; *p!='"'; p++);
        }
        p[11] = 0;
        DATES = DATESstring_to_double(p+1);
        p = p + 11;
        for (i=0; i<property_idx; i++) for (p++; *p!=','; p++);
        for (q=p+1; (*q!=',') && (*q!='\n') && (*q!=0); q++);
        if (q > p+1)
        {
            *q = 0;
            property[
                find_in_sorted_stack(DATES,num_DATES,DATES_stack) * num_CUSIP +
                find_in_sorted_stack(CUSIP,num_CUSIP,CUSIP_stack) ] =
                atof(p+1);
            num_entries ++;
        }
        num_entry_lines++;
        if (num_entry_lines % PROG_REPORT_FREQ == 0)
            printf ("%d entry lines read in.\n", num_entry_lines);
    }

    for (j=0; j<num_DATES; j++)
    {
        for (i=0; i<num_CUSIP; i++)
        {
            value = property[ j * num_CUSIP + i ];
            if (value == DOUBLE_PRECISION_INFINITY) fprintf (fp_out, "NaN ");
            else fprintf (fp_out, "%.16g ", value);
        }
        fprintf (fp_out, "\n");
    }
    fclose (fp_out);
    printf ("\n%d x %d matrix saved on \"%s\",\n"
            "%d entry lines (%.1f%%), %d entries (%.1f%%).\n\n",
            num_DATES, num_CUSIP, property_fname,
            num_entry_lines, 100.*num_entry_lines/num_DATES/num_CUSIP,
            num_entries,     100.*num_entries/num_DATES/num_CUSIP);

    fclose(fp_in);
    return (0);
} /* end main() */
/* cc -O3 csv2out.c -o csv2out; csv2out *.csv 1 */
#endif
