#ifndef BRENFORT_H
#define BRENFORT_H

#include "brenner.h"

/* Fortran subroutines called from C other than those in c_interface.f
 (those get automatically prototyped in c_interface.h): */

/* setup Langevin parameters  */
void setgle_(void);
/* initialize random number generator */
void setran_(void);
/* predictor */
void cpred_(void);
/* calculate energy and forces */
void model_(void);
/* apply thermostats */
void thermos_(void);
/* corrector */
void ccorr_(void);

void bere_(void);
/* Energy minimization */
void minimize_(void);

#endif /* BRENFORT_H */
