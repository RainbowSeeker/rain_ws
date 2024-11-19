/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: norm_FqL2z8vo.c
 *
 * Code generated for Simulink model 'Formation_FMS'.
 *
 * Model version                  : 1.259
 * Simulink Coder version         : 9.8 (R2022b) 13-May-2022
 * C/C++ source code generated on : Fri Aug 16 11:40:06 2024
 */

#include "rtwtypes.h"
#include "norm_FqL2z8vo.h"
#include <math.h>

/* Function for Chart: '<Root>/FMS State Machine' */
real32_T norm_FqL2z8vo(const real32_T x[2])
{
  real32_T absxk;
  real32_T scale;
  real32_T t;
  real32_T y;
  scale = 1.29246971E-26F;
  absxk = (real32_T)fabs(x[0]);
  if (absxk > 1.29246971E-26F) {
    y = 1.0F;
    scale = absxk;
  } else {
    t = absxk / 1.29246971E-26F;
    y = t * t;
  }

  absxk = (real32_T)fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0F;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * (real32_T)sqrt(y);
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
