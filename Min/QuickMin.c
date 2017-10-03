
#include "Min.h"

/********************************************************/
/* QuickMin minimizers:                                 */
/* A. Christensen, Int. J. Modern Phys. C 16 (2005) 193 */
/********************************************************/

double QuickMin_delta = 0.01;


/* Collective version of QuickMin */
double QuickMinCollective (double (*potential) (double *x, double *gradient),
                           int N, double *x)
{
    register int i;
    double *gradient, *displacement;
    double P, gradient_square_sum, step_square_sum;
    MALLOC (QuickMinCollective, gradient, N, double);
    CALLOC (QuickMinCollective, displacement, N, double);
    for (CG_state.number_of_evaluations=0;
         CG_state.number_of_evaluations<CG_option.max_potential_evaluations;
         CG_state.number_of_evaluations++)
    {
        CG_state.potential_min = potential(x, gradient);
        gradient_square_sum = 0;
        step_square_sum = 0;
        P = 0;
        for (i=N; i--;)
        {
            gradient_square_sum += SQUARE(gradient[i]);
            step_square_sum += SQUARE(displacement[i]);
            P -= gradient[i] * displacement[i];
        }
        if ( (gradient_square_sum < CG_option.gradient_square_sum_tolerance) &&
             (step_square_sum < CG_option.step_square_sum_tolerance) ) break;
        if (P > 0)
            for (i=N; i--;)
                displacement[i] = -gradient[i] *
                    (QuickMin_delta + P/gradient_square_sum);
        else
            for (i=N; i--;)
                displacement[i] = -gradient[i] * QuickMin_delta;
        for (i=N; i--;) x[i] += displacement[i];
    }
    free (gradient);
    free (displacement);
    if (CG_state.number_of_evaluations == CG_option.max_potential_evaluations)
    {
        CG_state.error = 1;
        CG_state.error_mesg =
            "QuickMinCollective: ITERATION WAS TERMINATED BECAUSE POTENTIAL\n"
            "EVALUATION LIMIT CG_option.max_potential_evaluations WAS "
            "EXCEEDED.\n";
    }
    else CG_state.error = 0;
    return(CG_state.potential_min);
} /* end QuickMinCollective() */


/* Entry-by-entry version of QuickMin */
double QuickMinIndividual (double (*potential) (double *x, double *gradient),
                           int N, double *x)
{
    register int i;
    double *gradient, *displacement;
    double P, gradient_square_sum, step_square_sum;
    MALLOC (QuickMinIndividual, gradient, N, double);
    CALLOC (QuickMinIndividual, displacement, N, double);
    for (CG_state.number_of_evaluations=0;
         CG_state.number_of_evaluations<CG_option.max_potential_evaluations;
         CG_state.number_of_evaluations++)
    {
        CG_state.potential_min = potential(x, gradient);
        gradient_square_sum = 0;
        step_square_sum = 0;
        for (i=N; i--;)
        {
            gradient_square_sum += SQUARE(gradient[i]);
            step_square_sum += SQUARE(displacement[i]);
        }
        if ( (gradient_square_sum < CG_option.gradient_square_sum_tolerance) &&
             (step_square_sum < CG_option.step_square_sum_tolerance) ) break;
        for (i=N; i--;)
        {
            P = -gradient[i] * displacement[i];
            if (P > 0)
                displacement[i] += -gradient[i] * QuickMin_delta;
            else
                displacement[i] = -gradient[i] * QuickMin_delta;
        }
        for (i=N; i--;) x[i] += displacement[i];
    }
    free (gradient);
    free (displacement);
    if (CG_state.number_of_evaluations == CG_option.max_potential_evaluations)
    {
        CG_state.error = 1;
        CG_state.error_mesg =
            "QuickMinIndividual: ITERATION WAS TERMINATED BECAUSE POTENTIAL\n"
            "EVALUATION LIMIT CG_option.max_potential_evaluations WAS "
            "EXCEEDED.\n";
    }
    else CG_state.error = 0;
    return(CG_state.potential_min);
} /* end QuickMinIndividual() */


#ifdef _minimizer_contest_TEST
#define N 10000
double xanswer[N];
double test_potential (double *x, double *gradient)
{
    int i;
    double value;
    for (value=i=0; i<N; i++)
    {
        value += SQUARE(x[i] - xanswer[i]) / (i+1);
        gradient[i] = 2 * (x[i] - xanswer[i]) / (i+1);
    }
    return (value);
} /* end test_potential() */

int main (int argc, char *argv[])
{
    register int i,j;
    double xstart[N], xmin[N], value, gradient[N], g2, xdiff2;
    struct List
    {
        char *name;
        double (*minimizer) (double (*potential) (double *x, double *gradient),
                             int n, double *x);
    } funs[] =
          {{"IMSL zxcgr", &zxcgr},
           {"QuickMinCollective", &QuickMinCollective},
           {"QuickMinIndividual", &QuickMinIndividual},
           {"FIRE",               &FIREminimize},};
    Vfrandom (N, xanswer);
    Vfrandom (N, xstart);
    for (i=0; i<sizeof(funs)/sizeof(struct List); i++)
    {
        VEQV(N, xstart, xmin);
        printf ("Testing %s...", funs[i].name);
        value = funs[i].minimizer (test_potential, N, xmin);
        if (CG_state.error) pe(CG_state.error_mesg);
        printf (" %d evaluations\n", CG_state.number_of_evaluations);
        if (value!=test_potential(xmin,gradient))
            printf("warning: value incongruence of %e vs %e\n",
                   value, test_potential(xmin,gradient));
        VLENGTH2 (N, gradient, j, g2);
        VDIFF2 (N, xmin, xanswer, j, xdiff2);
        printf ("|g|^2=%e, |Dx|^2=%e, e=%e\n", g2, xdiff2, value);
    }
    return (0);
}
#endif /* _minimizer_contest_TEST */
