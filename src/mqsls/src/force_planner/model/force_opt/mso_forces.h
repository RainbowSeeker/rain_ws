/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mso_forces.h
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2025-01-02 16:21:33
 */

#ifndef MSO_FORCES_H
#define MSO_FORCES_H

/* Include Files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void mso_forces(const double center[3], const double T_min[3],
                       const double T_max[3], double best_result[9],
                       double *best_radius, double *best_exitflag);

extern void mso_forces_initialize(void);

extern void mso_forces_terminate(void);

#ifdef __cplusplus
}
#endif

#endif
/*
 * File trailer for mso_forces.h
 *
 * [EOF]
 */
