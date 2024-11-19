/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: rt_atan2f_snf.c
 *
 * Code generated for Simulink model 'Formation_FMS'.
 *
 * Model version                  : 1.259
 * Simulink Coder version         : 9.8 (R2022b) 13-May-2022
 * C/C++ source code generated on : Fri Aug 16 11:40:06 2024
 */

#include "rtwtypes.h"
#include "rt_atan2f_snf.h"
#include "rt_nonfinite.h"
#include <math.h>
#include "rt_defines.h"
#include "rtGetNaN.h"

real32_T rt_atan2f_snf(real32_T u0, real32_T u1)
{
  int32_T tmp;
  int32_T tmp_0;
  real32_T y;
  if (rtIsNaNF(u0) || rtIsNaNF(u1)) {
    y = (rtNaNF);
  } else if (rtIsInfF(u0) && rtIsInfF(u1)) {
    if (u0 > 0.0F) {
      tmp = 1;
    } else {
      tmp = -1;
    }

    if (u1 > 0.0F) {
      tmp_0 = 1;
    } else {
      tmp_0 = -1;
    }

    y = (real32_T)atan2((real32_T)tmp, (real32_T)tmp_0);
  } else if (u1 == 0.0F) {
    if (u0 > 0.0F) {
      y = RT_PIF / 2.0F;
    } else if (u0 < 0.0F) {
      y = -(RT_PIF / 2.0F);
    } else {
      y = 0.0F;
    }
  } else {
    y = (real32_T)atan2(u0, u1);
  }

  return y;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
