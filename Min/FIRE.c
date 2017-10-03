
#include "Min.h"

/*****************************************************/
/* FIRE minimizer: Bitzek, Koskinen, Gähler, Moseler */
/* and Gumbsch, Phys. Rev. Lett. 97 (2006) 170201.   */
/*****************************************************/

FIRE_Option FIRE_option = {
    5,    /* Nmin */
    1.1,  /* finc */
    0.5,  /* fdec */
    0.1,  /* alphastart */
    0.99, /* falpha */
    0.5   /* timestepmax */
};


double FIREminimize (double (*potential) (double *x, double *gradient),
                     int N, double *x)
{
    register int i;
    int iter_since_negative = 0;
    double *gradient, *velocity;
    double P, gradient_square_sum, velocity_square_sum;
    double timestep, alpha;
  
    timestep = FIRE_option.timestepmax / 10; 
    alpha = FIRE_option.alphastart; 

    MALLOC (FIREminimize, gradient, N, double);
    CALLOC (FIREminimize, velocity, N, double);

    for (CG_state.number_of_evaluations=0;
         CG_state.number_of_evaluations<CG_option.max_potential_evaluations;
         CG_state.number_of_evaluations++)
    {
        CG_state.potential_min = potential(x, gradient);
        gradient_square_sum = 0;
        velocity_square_sum = 0;
        P = 0;
        for (i=N; i--;)
        {
            velocity[i] -= gradient[i] * timestep;
            P -= gradient[i] * velocity[i];
            gradient_square_sum += SQUARE(gradient[i]);
            velocity_square_sum += SQUARE(velocity[i]);
        }
        if ( (gradient_square_sum < CG_option.gradient_square_sum_tolerance) &&
             (velocity_square_sum * SQUARE(timestep) <
              CG_option.step_square_sum_tolerance) ) break;
        if (P <= 0)
        {
            timestep *= FIRE_option.fdec;
            for (i=N; i--;) velocity[i] = 0;
            alpha = FIRE_option.alphastart;
            iter_since_negative = 0;
        }
        else
        {
            for (i=N; i--;)
                velocity[i] = (1-alpha) * velocity[i] - alpha * gradient[i] *
                    sqrt(velocity_square_sum/gradient_square_sum);
            iter_since_negative++;
            if (iter_since_negative > FIRE_option.Nmin)
            {
                timestep = MIN( timestep*FIRE_option.finc,
                                FIRE_option.timestepmax );
                alpha *= FIRE_option.falpha;
            }
        }
        for (i=N; i--;) x[i] += velocity[i] * timestep;
    }
    free (gradient);
    free (velocity);

    if (CG_state.number_of_evaluations == CG_option.max_potential_evaluations)
    {
        CG_state.error = 1;
        CG_state.error_mesg =
            "FIREminimize: ITERATION WAS TERMINATED BECAUSE POTENTIAL\n"
            "EVALUATION LIMIT CG_option.max_potential_evaluations WAS "
            "EXCEEDED.\n";
    }
    else CG_state.error = 0;
    return(CG_state.potential_min);
} /* end FIREminimize() */
