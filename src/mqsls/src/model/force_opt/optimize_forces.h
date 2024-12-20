/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: optimize_forces.h
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024-12-19 20:42:18
 */

#ifndef OPTIMIZE_FORCES_H
#define OPTIMIZE_FORCES_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void optimize_forces(const double center[3],
                            const double initial_guess[9],
                            const double T_min[3], const double T_max[3],
                            const double psi[3], double result[9],
                            double *radius, double *exitflag);

extern void optimize_forces_initialize(void);

extern void optimize_forces_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for optimize_forces.h
 *
 * [EOF]
 */
