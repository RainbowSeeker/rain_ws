//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: control_3dof.cpp
//
// Code generated for Simulink model 'control_3dof'.
//
// Model version                  : 1.701
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Fri Dec 20 15:10:40 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "control_3dof.h"
#include "rtwtypes.h"
#include <cmath>
#include <emmintrin.h>
#include <cstring>
#include "cmath"
#include "limits"

// Exported block parameters
struct_VO89MlSMmM5pi10TnqAIZ CONTROL_PARAM{
  1.0,
  1.5,
  1.0,
  0.0,
  20.0,

  { 3.0, 5.0, 10.0 },

  { 0.0, 20.0, 20.0 },
  0.2,
  1.0,
  0.5,
  1.0,
  0.5
} ;                                    // Variable: CONTROL_PARAM
                                          //  Referenced by:
                                          //    '<S4>/Force Saturation && Disturbution'
                                          //    '<S4>/Position 2nd ESO'
                                          //    '<S6>/Parallel Control'
                                          //    '<S7>/Vertical Control'
                                          //    '<S7>/Gain'
                                          //    '<S10>/Parallel Control'
                                          //    '<S11>/Vertical Control'
                                          //    '<S11>/Gain'
                                          //    '<S14>/Parallel Control'
                                          //    '<S15>/Vertical Control'
                                          //    '<S15>/Gain'
                                          //    '<S19>/KP'
                                          //    '<S19>/KV'


extern "C"
{
  real_T rtNaN { -std::numeric_limits<real_T>::quiet_NaN() };

  real_T rtInf { std::numeric_limits<real_T>::infinity() };

  real_T rtMinusInf { -std::numeric_limits<real_T>::infinity() };

  real32_T rtNaNF { -std::numeric_limits<real32_T>::quiet_NaN() };

  real32_T rtInfF { std::numeric_limits<real32_T>::infinity() };

  real32_T rtMinusInfF { -std::numeric_limits<real32_T>::infinity() };
}

extern "C"
{
  // Return rtNaN needed by the generated code.
  static real_T rtGetNaN(void)
  {
    return rtNaN;
  }

  // Return rtNaNF needed by the generated code.
  static real32_T rtGetNaNF(void)
  {
    return rtNaNF;
  }
}

extern "C"
{
  // Test if value is infinite
  static boolean_T rtIsInf(real_T value)
  {
    return std::isinf(value);
  }

  // Test if single-precision value is infinite
  static boolean_T rtIsInfF(real32_T value)
  {
    return std::isinf(value);
  }

  // Test if value is not a number
  static boolean_T rtIsNaN(real_T value)
  {
    return std::isnan(value);
  }

  // Test if single-precision value is not a number
  static boolean_T rtIsNaNF(real32_T value)
  {
    return std::isnan(value);
  }
}

//
// Output and update for atomic system:
//    '<S6>/Parallel Control'
//    '<S10>/Parallel Control'
//    '<S14>/Parallel Control'
//
void control_3dof::control_3dof_ParallelControl(const real_T rtu_q[3], const
  real_T rtu_w[3], const real_T rtu_vforce[3], const real_T rtu_ai_sp[3], real_T
  rty_f_parallel[3], real_T rtp_li, real_T rtp_mi,
  B_ParallelControl_control_3do_T *localB)
{
  localB->scale = 3.3121686421112381E-170;
  localB->absxk = std::abs(rtu_w[0]);
  if (localB->absxk > 3.3121686421112381E-170) {
    localB->a = 1.0;
    localB->scale = localB->absxk;
  } else {
    localB->t = localB->absxk / 3.3121686421112381E-170;
    localB->a = localB->t * localB->t;
  }

  localB->absxk = std::abs(rtu_w[1]);
  if (localB->absxk > localB->scale) {
    localB->t = localB->scale / localB->absxk;
    localB->a = localB->a * localB->t * localB->t + 1.0;
    localB->scale = localB->absxk;
  } else {
    localB->t = localB->absxk / localB->scale;
    localB->a += localB->t * localB->t;
  }

  localB->absxk = std::abs(rtu_w[2]);
  if (localB->absxk > localB->scale) {
    localB->t = localB->scale / localB->absxk;
    localB->a = localB->a * localB->t * localB->t + 1.0;
    localB->scale = localB->absxk;
  } else {
    localB->t = localB->absxk / localB->scale;
    localB->a += localB->t * localB->t;
  }

  localB->scale *= std::sqrt(localB->a);
  localB->scale = rtp_mi * rtp_li * (localB->scale * localB->scale);
  for (localB->i = 0; localB->i < 3; localB->i++) {
    localB->rtu_q[3 * localB->i] = rtu_q[0] * rtu_q[localB->i];
    localB->rtp_mi[3 * localB->i] = rtp_mi * rtu_q[0] * rtu_q[localB->i];
    localB->rtu_q_tmp = 3 * localB->i + 1;
    localB->rtu_q[localB->rtu_q_tmp] = rtu_q[1] * rtu_q[localB->i];
    localB->rtp_mi[localB->rtu_q_tmp] = rtp_mi * rtu_q[1] * rtu_q[localB->i];
    localB->rtu_q_tmp = 3 * localB->i + 2;
    localB->rtu_q[localB->rtu_q_tmp] = rtu_q[2] * rtu_q[localB->i];
    localB->rtp_mi[localB->rtu_q_tmp] = rtp_mi * rtu_q[2] * rtu_q[localB->i];
  }

  localB->absxk = rtu_vforce[1];
  localB->t = rtu_vforce[0];
  localB->a = rtu_vforce[2];
  localB->rtu_ai_sp = rtu_ai_sp[1];
  localB->rtu_ai_sp_m = rtu_ai_sp[0];
  localB->rtu_ai_sp_c = rtu_ai_sp[2];
  for (localB->i = 0; localB->i <= 0; localB->i += 2) {
    __m128d tmp;
    __m128d tmp_0;
    __m128d tmp_1;
    __m128d tmp_2;
    __m128d tmp_3;
    __m128d tmp_4;
    __m128d tmp_5;
    tmp = _mm_loadu_pd(&localB->rtu_q[localB->i + 3]);
    tmp_0 = _mm_loadu_pd(&localB->rtu_q[localB->i]);
    tmp_1 = _mm_loadu_pd(&localB->rtu_q[localB->i + 6]);
    tmp_2 = _mm_loadu_pd(&rtu_q[localB->i]);
    tmp_3 = _mm_loadu_pd(&localB->rtp_mi[localB->i + 3]);
    tmp_4 = _mm_loadu_pd(&localB->rtp_mi[localB->i]);
    tmp_5 = _mm_loadu_pd(&localB->rtp_mi[localB->i + 6]);
    _mm_storeu_pd(&rty_f_parallel[localB->i], _mm_add_pd(_mm_add_pd(_mm_add_pd
      (_mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(localB->absxk)), _mm_mul_pd(tmp_0,
      _mm_set1_pd(localB->t))), _mm_mul_pd(tmp_1, _mm_set1_pd(localB->a))),
      _mm_mul_pd(_mm_set1_pd(localB->scale), tmp_2)), _mm_add_pd(_mm_add_pd
      (_mm_mul_pd(tmp_3, _mm_set1_pd(localB->rtu_ai_sp)), _mm_mul_pd(tmp_4,
      _mm_set1_pd(localB->rtu_ai_sp_m))), _mm_mul_pd(tmp_5, _mm_set1_pd
      (localB->rtu_ai_sp_c)))));
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    rty_f_parallel[localB->i] = (((localB->rtu_q[localB->i + 3] * localB->absxk
      + localB->rtu_q[localB->i] * localB->t) + localB->rtu_q[localB->i + 6] *
      localB->a) + localB->scale * rtu_q[localB->i]) + ((localB->rtp_mi
      [localB->i + 3] * localB->rtu_ai_sp + localB->rtp_mi[localB->i] *
      localB->rtu_ai_sp_m) + localB->rtp_mi[localB->i + 6] * localB->rtu_ai_sp_c);
  }
}

//
// Output and update for atomic system:
//    '<S7>/Vertical Control'
//    '<S11>/Vertical Control'
//    '<S15>/Vertical Control'
//
void control_3dof::control_3dof_VerticalControl(const real_T rtu_q_sp[3], const
  real_T rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
  rty_f_vertical[3], real_T rty_dot_err[3], real_T rty_state[12], real_T rtp_Kq,
  real_T rtp_Kw, real_T rtp_li, real_T rtp_mi, B_VerticalControl_control_3do_T
  *localB)
{
  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;
  localB->Sqi[0] = 0.0;
  localB->Sqi[3] = -rtu_q[2];
  localB->Sqi[6] = rtu_q[1];
  localB->Sqi[1] = rtu_q[2];
  localB->Sqi[4] = 0.0;
  localB->Sqi[7] = -rtu_q[0];
  localB->Sqi[2] = -rtu_q[1];
  localB->Sqi[5] = rtu_q[0];
  localB->Sqi[8] = 0.0;
  _mm_storeu_pd(&localB->e_q[0], _mm_sub_pd(_mm_mul_pd(_mm_set_pd(rtu_q[0],
    rtu_q_sp[1]), _mm_set_pd(rtu_q_sp[2], rtu_q[2])), _mm_mul_pd(_mm_set_pd
    (rtu_q_sp[0], rtu_q[1]), _mm_set_pd(rtu_q[2], rtu_q_sp[2]))));
  localB->e_q[2] = rtu_q_sp[0] * rtu_q[1] - rtu_q[0] * rtu_q_sp[1];
  for (localB->i = 0; localB->i <= 0; localB->i += 2) {
    tmp_5 = _mm_loadu_pd(&localB->e_q[localB->i]);
    _mm_storeu_pd(&localB->w_sp_tmp[localB->i], _mm_mul_pd(_mm_set1_pd(-rtp_Kq),
      tmp_5));
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    localB->w_sp_tmp[localB->i] = -rtp_Kq * localB->e_q[localB->i];
  }

  for (localB->i_m = 0; localB->i_m < 3; localB->i_m++) {
    for (localB->i = 0; localB->i <= 0; localB->i += 2) {
      tmp_5 = _mm_loadu_pd(&localB->Sqi[localB->i + 3]);
      tmp_3 = _mm_loadu_pd(&localB->Sqi[localB->i]);
      tmp_4 = _mm_loadu_pd(&localB->Sqi[localB->i + 6]);
      _mm_storeu_pd(&localB->Sqi_m[localB->i + 3 * localB->i_m], _mm_add_pd
                    (_mm_add_pd(_mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_m + 1]), tmp_5), _mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_m]), tmp_3)), _mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_m + 2]), tmp_4)));
    }

    for (localB->i = 2; localB->i < 3; localB->i++) {
      localB->Sqi_m[localB->i + 3 * localB->i_m] = (localB->Sqi[3 * localB->i_m
        + 1] * localB->Sqi[localB->i + 3] + localB->Sqi[3 * localB->i_m] *
        localB->Sqi[localB->i]) + localB->Sqi[3 * localB->i_m + 2] * localB->
        Sqi[localB->i + 6];
    }
  }

  localB->a = rtp_mi * rtp_li;
  localB->rtu_q = 0.0;
  for (localB->i = 0; localB->i < 3; localB->i++) {
    localB->e_w_tmp_tmp = (localB->Sqi_m[localB->i + 3] * localB->w_sp_tmp[1] +
      localB->Sqi_m[localB->i] * localB->w_sp_tmp[0]) + localB->Sqi_m[localB->i
      + 6] * localB->w_sp_tmp[2];
    localB->e_w_tmp[localB->i] = localB->e_w_tmp_tmp;
    localB->e_w[localB->i] = rtu_w[localB->i] + localB->e_w_tmp_tmp;
    localB->rtu_q += rtu_q[localB->i] * localB->w_sp_tmp[localB->i];
  }

  _mm_storeu_pd(&localB->b_a[0], _mm_mul_pd(_mm_sub_pd(_mm_mul_pd(_mm_set_pd
    (rtu_q_sp[0], localB->w_sp_tmp[1]), _mm_set_pd(localB->w_sp_tmp[2],
    rtu_q_sp[2])), _mm_mul_pd(_mm_set_pd(localB->w_sp_tmp[0], rtu_q_sp[1]),
    _mm_set_pd(rtu_q_sp[2], localB->w_sp_tmp[2]))), _mm_set1_pd(localB->rtu_q)));
  localB->b_a[2] = (localB->w_sp_tmp[0] * rtu_q_sp[1] - rtu_q_sp[0] *
                    localB->w_sp_tmp[1]) * localB->rtu_q;
  for (localB->i_m = 0; localB->i_m < 3; localB->i_m++) {
    localB->w_sp_tmp_c[localB->i_m] = ((localB->w_sp_tmp[localB->i_m] - rtp_Kw *
      localB->e_w[localB->i_m]) - localB->e_w_tmp[localB->i_m]) - localB->
      b_a[localB->i_m];
    for (localB->i = 0; localB->i < 3; localB->i++) {
      localB->Sqi_m[localB->i_m + 3 * localB->i] = (localB->Sqi[localB->i_m + 3]
        * rtp_mi * localB->Sqi[3 * localB->i + 1] + rtp_mi * localB->Sqi
        [localB->i_m] * localB->Sqi[3 * localB->i]) + localB->Sqi[localB->i_m +
        6] * rtp_mi * localB->Sqi[3 * localB->i + 2];
    }
  }

  localB->rtu_q = -1.0 / rtp_mi / rtp_li;
  tmp_5 = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(2.0), _mm_loadu_pd(&localB->e_w[0])),
                     _mm_loadu_pd(&localB->e_q[0]));
  _mm_storeu_pd(&localB->dv[0], tmp_5);
  localB->e_w_tmp_tmp = localB->dv[0];
  localB->e_q_idx_1 = localB->dv[1];
  localB->e_q_idx_2 = 2.0 * localB->e_w[2] + localB->e_q[2];
  localB->w_sp_tmp_p = localB->w_sp_tmp_c[1];
  localB->w_sp_tmp_cv = localB->w_sp_tmp_c[0];
  localB->w_sp_tmp_f = localB->w_sp_tmp_c[2];
  localB->rtu_ai_sp = rtu_ai_sp[1];
  localB->rtu_ai_sp_g = rtu_ai_sp[0];
  localB->rtu_ai_sp_g1 = rtu_ai_sp[2];
  for (localB->i = 0; localB->i <= 0; localB->i += 2) {
    __m128d tmp;
    __m128d tmp_0;
    __m128d tmp_1;
    __m128d tmp_2;
    tmp_5 = _mm_loadu_pd(&localB->Sqi[localB->i + 3]);
    tmp_3 = _mm_set1_pd(localB->a);
    tmp_4 = _mm_loadu_pd(&localB->Sqi[localB->i]);
    tmp = _mm_loadu_pd(&localB->Sqi[localB->i + 6]);
    tmp_0 = _mm_loadu_pd(&localB->Sqi_m[localB->i + 3]);
    tmp_1 = _mm_loadu_pd(&localB->Sqi_m[localB->i]);
    tmp_2 = _mm_loadu_pd(&localB->Sqi_m[localB->i + 6]);
    _mm_storeu_pd(&rty_f_vertical[localB->i], _mm_sub_pd(_mm_add_pd(_mm_add_pd
      (_mm_mul_pd(_mm_mul_pd(tmp_5, tmp_3), _mm_set1_pd(localB->w_sp_tmp_p)),
       _mm_mul_pd(_mm_mul_pd(tmp_3, tmp_4), _mm_set1_pd(localB->w_sp_tmp_cv))),
      _mm_mul_pd(_mm_mul_pd(tmp, tmp_3), _mm_set1_pd(localB->w_sp_tmp_f))),
      _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_0, _mm_set1_pd(localB->rtu_ai_sp)),
      _mm_mul_pd(tmp_1, _mm_set1_pd(localB->rtu_ai_sp_g))), _mm_mul_pd(tmp_2,
      _mm_set1_pd(localB->rtu_ai_sp_g1)))));
    tmp_3 = _mm_set1_pd(localB->rtu_q);
    _mm_storeu_pd(&rty_dot_err[localB->i], _mm_add_pd(_mm_add_pd(_mm_mul_pd
      (_mm_mul_pd(tmp_5, tmp_3), _mm_set1_pd(localB->e_q_idx_1)), _mm_mul_pd
      (_mm_mul_pd(tmp_3, tmp_4), _mm_set1_pd(localB->e_w_tmp_tmp))), _mm_mul_pd
      (_mm_mul_pd(tmp, tmp_3), _mm_set1_pd(localB->e_q_idx_2))));
    tmp_5 = _mm_loadu_pd(&rtu_q_sp[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i], tmp_5);
    tmp_5 = _mm_loadu_pd(&localB->e_q[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i + 3], tmp_5);
    tmp_5 = _mm_loadu_pd(&localB->w_sp_tmp[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i + 6], tmp_5);
    tmp_5 = _mm_loadu_pd(&localB->e_w[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i + 9], tmp_5);
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    localB->Sqi_k = localB->Sqi[localB->i + 3];
    localB->Sqi_c = localB->Sqi[localB->i];
    localB->Sqi_b = localB->Sqi[localB->i + 6];
    rty_f_vertical[localB->i] = ((localB->Sqi_k * localB->a * localB->w_sp_tmp_p
      + localB->a * localB->Sqi_c * localB->w_sp_tmp_cv) + localB->Sqi_b *
      localB->a * localB->w_sp_tmp_f) - ((localB->Sqi_m[localB->i + 3] *
      localB->rtu_ai_sp + localB->Sqi_m[localB->i] * localB->rtu_ai_sp_g) +
      localB->Sqi_m[localB->i + 6] * localB->rtu_ai_sp_g1);
    rty_dot_err[localB->i] = (localB->Sqi_k * localB->rtu_q * localB->e_q_idx_1
      + localB->rtu_q * localB->Sqi_c * localB->e_w_tmp_tmp) + localB->Sqi_b *
      localB->rtu_q * localB->e_q_idx_2;
    rty_state[localB->i] = rtu_q_sp[localB->i];
    rty_state[localB->i + 3] = localB->e_q[localB->i];
    rty_state[localB->i + 6] = localB->w_sp_tmp[localB->i];
    rty_state[localB->i + 9] = localB->e_w[localB->i];
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
real_T control_3dof::control_3dof_norm(const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = std::abs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = std::abs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = std::abs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * std::sqrt(y);
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
real_T control_3dof::control_3dof_xzlangeM(const real_T x[9])
{
  real_T y;
  int32_T k;
  boolean_T exitg1;
  y = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 9)) {
    real_T absxk;
    absxk = std::abs(x[k]);
    if (std::isnan(absxk)) {
      y = (rtNaN);
      exitg1 = true;
    } else {
      if (absxk > y) {
        y = absxk;
      }

      k++;
    }
  }

  return y;
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xzlascl(real_T cfrom, real_T cto, int32_T m,
  int32_T n, real_T A[9], int32_T iA0, int32_T lda)
{
  control_3dof_B.cfromc_n = cfrom;
  control_3dof_B.ctoc_p = cto;
  control_3dof_B.notdone_o = true;
  while (control_3dof_B.notdone_o) {
    control_3dof_B.cfrom1_l = control_3dof_B.cfromc_n * 2.0041683600089728E-292;
    control_3dof_B.cto1_j = control_3dof_B.ctoc_p / 4.9896007738368E+291;
    if ((std::abs(control_3dof_B.cfrom1_l) > std::abs(control_3dof_B.ctoc_p)) &&
        (control_3dof_B.ctoc_p != 0.0)) {
      control_3dof_B.mul_d = 2.0041683600089728E-292;
      control_3dof_B.cfromc_n = control_3dof_B.cfrom1_l;
    } else if (std::abs(control_3dof_B.cto1_j) > std::abs
               (control_3dof_B.cfromc_n)) {
      control_3dof_B.mul_d = 4.9896007738368E+291;
      control_3dof_B.ctoc_p = control_3dof_B.cto1_j;
    } else {
      control_3dof_B.mul_d = control_3dof_B.ctoc_p / control_3dof_B.cfromc_n;
      control_3dof_B.notdone_o = false;
    }

    for (control_3dof_B.j_bn = 0; control_3dof_B.j_bn < n; control_3dof_B.j_bn++)
    {
      control_3dof_B.offset_h = (control_3dof_B.j_bn * lda + iA0) - 2;
      control_3dof_B.scalarLB_e = (m / 2) << 1;
      control_3dof_B.vectorUB_bj = control_3dof_B.scalarLB_e - 2;
      for (control_3dof_B.b_i_d = 0; control_3dof_B.b_i_d <=
           control_3dof_B.vectorUB_bj; control_3dof_B.b_i_d += 2) {
        __m128d tmp;
        control_3dof_B.i1 = (control_3dof_B.b_i_d + control_3dof_B.offset_h) + 1;
        tmp = _mm_loadu_pd(&A[control_3dof_B.i1]);
        _mm_storeu_pd(&A[control_3dof_B.i1], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.mul_d)));
      }

      for (control_3dof_B.b_i_d = control_3dof_B.scalarLB_e;
           control_3dof_B.b_i_d < m; control_3dof_B.b_i_d++) {
        control_3dof_B.i1 = (control_3dof_B.b_i_d + control_3dof_B.offset_h) + 1;
        A[control_3dof_B.i1] *= control_3dof_B.mul_d;
      }
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
real_T control_3dof::control_3dof_xnrm2(int32_T n, const real_T x[9], int32_T
  ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      control_3dof_B.scale_l = 3.3121686421112381E-170;
      control_3dof_B.kend_n = ix0 + n;
      for (control_3dof_B.k_i = ix0; control_3dof_B.k_i < control_3dof_B.kend_n;
           control_3dof_B.k_i++) {
        control_3dof_B.absxk_d = std::abs(x[control_3dof_B.k_i - 1]);
        if (control_3dof_B.absxk_d > control_3dof_B.scale_l) {
          control_3dof_B.t_d = control_3dof_B.scale_l / control_3dof_B.absxk_d;
          y = y * control_3dof_B.t_d * control_3dof_B.t_d + 1.0;
          control_3dof_B.scale_l = control_3dof_B.absxk_d;
        } else {
          control_3dof_B.t_d = control_3dof_B.absxk_d / control_3dof_B.scale_l;
          y += control_3dof_B.t_d * control_3dof_B.t_d;
        }
      }

      y = control_3dof_B.scale_l * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
real_T control_3dof::control_3dof_xdotc(int32_T n, const real_T x[9], int32_T
  ix0, const real_T y[9], int32_T iy0)
{
  real_T d;
  d = 0.0;
  if (n >= 1) {
    for (int32_T k{0}; k < n; k++) {
      d += x[(ix0 + k) - 1] * y[(iy0 + k) - 1];
    }
  }

  return d;
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xaxpy(int32_T n, real_T a, int32_T ix0, real_T
  y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (int32_T k{0}; k < n; k++) {
      int32_T tmp;
      tmp = (iy0 + k) - 1;
      y[tmp] += y[(ix0 + k) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
real_T control_3dof::control_3dof_xnrm2_n(int32_T n, const real_T x[3], int32_T
  ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      control_3dof_B.scale_g = 3.3121686421112381E-170;
      control_3dof_B.kend = ix0 + n;
      for (control_3dof_B.k_o = ix0; control_3dof_B.k_o < control_3dof_B.kend;
           control_3dof_B.k_o++) {
        control_3dof_B.absxk = std::abs(x[control_3dof_B.k_o - 1]);
        if (control_3dof_B.absxk > control_3dof_B.scale_g) {
          control_3dof_B.t = control_3dof_B.scale_g / control_3dof_B.absxk;
          y = y * control_3dof_B.t * control_3dof_B.t + 1.0;
          control_3dof_B.scale_g = control_3dof_B.absxk;
        } else {
          control_3dof_B.t = control_3dof_B.absxk / control_3dof_B.scale_g;
          y += control_3dof_B.t * control_3dof_B.t;
        }
      }

      y = control_3dof_B.scale_g * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xaxpy_f(int32_T n, real_T a, const real_T x[9],
  int32_T ix0, real_T y[3], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_ju = (n / 2) << 1;
    control_3dof_B.vectorUB_j = control_3dof_B.scalarLB_ju - 2;
    for (control_3dof_B.k_a = 0; control_3dof_B.k_a <= control_3dof_B.vectorUB_j;
         control_3dof_B.k_a += 2) {
      __m128d tmp;
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_a) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i3]);
      _mm_storeu_pd(&y[control_3dof_B.i3], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k_a) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k_a = control_3dof_B.scalarLB_ju; control_3dof_B.k_a < n;
         control_3dof_B.k_a++) {
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_a) - 1;
      y[control_3dof_B.i3] += x[(ix0 + control_3dof_B.k_a) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3],
  int32_T ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_j = (n / 2) << 1;
    control_3dof_B.vectorUB_f = control_3dof_B.scalarLB_j - 2;
    for (control_3dof_B.k = 0; control_3dof_B.k <= control_3dof_B.vectorUB_f;
         control_3dof_B.k += 2) {
      __m128d tmp;
      control_3dof_B.i2 = (iy0 + control_3dof_B.k) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i2]);
      _mm_storeu_pd(&y[control_3dof_B.i2], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k = control_3dof_B.scalarLB_j; control_3dof_B.k < n;
         control_3dof_B.k++) {
      control_3dof_B.i2 = (iy0 + control_3dof_B.k) - 1;
      y[control_3dof_B.i2] += x[(ix0 + control_3dof_B.k) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xzlascl_l(real_T cfrom, real_T cto, int32_T m,
  int32_T n, real_T A[3], int32_T iA0, int32_T lda)
{
  control_3dof_B.cfromc = cfrom;
  control_3dof_B.ctoc = cto;
  control_3dof_B.notdone = true;
  while (control_3dof_B.notdone) {
    control_3dof_B.cfrom1 = control_3dof_B.cfromc * 2.0041683600089728E-292;
    control_3dof_B.cto1 = control_3dof_B.ctoc / 4.9896007738368E+291;
    if ((std::abs(control_3dof_B.cfrom1) > std::abs(control_3dof_B.ctoc)) &&
        (control_3dof_B.ctoc != 0.0)) {
      control_3dof_B.mul = 2.0041683600089728E-292;
      control_3dof_B.cfromc = control_3dof_B.cfrom1;
    } else if (std::abs(control_3dof_B.cto1) > std::abs(control_3dof_B.cfromc))
    {
      control_3dof_B.mul = 4.9896007738368E+291;
      control_3dof_B.ctoc = control_3dof_B.cto1;
    } else {
      control_3dof_B.mul = control_3dof_B.ctoc / control_3dof_B.cfromc;
      control_3dof_B.notdone = false;
    }

    for (control_3dof_B.j_b = 0; control_3dof_B.j_b < n; control_3dof_B.j_b++) {
      control_3dof_B.offset = (control_3dof_B.j_b * lda + iA0) - 2;
      control_3dof_B.scalarLB = (m / 2) << 1;
      control_3dof_B.vectorUB_b = control_3dof_B.scalarLB - 2;
      for (control_3dof_B.b_i_n = 0; control_3dof_B.b_i_n <=
           control_3dof_B.vectorUB_b; control_3dof_B.b_i_n += 2) {
        __m128d tmp;
        control_3dof_B.i_ln = (control_3dof_B.b_i_n + control_3dof_B.offset) + 1;
        tmp = _mm_loadu_pd(&A[control_3dof_B.i_ln]);
        _mm_storeu_pd(&A[control_3dof_B.i_ln], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.mul)));
      }

      for (control_3dof_B.b_i_n = control_3dof_B.scalarLB; control_3dof_B.b_i_n <
           m; control_3dof_B.b_i_n++) {
        control_3dof_B.i_ln = (control_3dof_B.b_i_n + control_3dof_B.offset) + 1;
        A[control_3dof_B.i_ln] *= control_3dof_B.mul;
      }
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xswap(real_T x[9], int32_T ix0, int32_T iy0)
{
  real_T temp;
  temp = x[ix0 - 1];
  x[ix0 - 1] = x[iy0 - 1];
  x[iy0 - 1] = temp;
  temp = x[ix0];
  x[ix0] = x[iy0];
  x[iy0] = temp;
  temp = x[ix0 + 1];
  x[ix0 + 1] = x[iy0 + 1];
  x[iy0 + 1] = temp;
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xrotg(real_T *a, real_T *b, real_T *c, real_T *s)
{
  control_3dof_B.roe = *b;
  control_3dof_B.absa = std::abs(*a);
  control_3dof_B.absb = std::abs(*b);
  if (control_3dof_B.absa > control_3dof_B.absb) {
    control_3dof_B.roe = *a;
  }

  control_3dof_B.scale = control_3dof_B.absa + control_3dof_B.absb;
  if (control_3dof_B.scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    control_3dof_B.ads = control_3dof_B.absa / control_3dof_B.scale;
    control_3dof_B.bds = control_3dof_B.absb / control_3dof_B.scale;
    control_3dof_B.scale *= std::sqrt(control_3dof_B.ads * control_3dof_B.ads +
      control_3dof_B.bds * control_3dof_B.bds);
    if (control_3dof_B.roe < 0.0) {
      control_3dof_B.scale = -control_3dof_B.scale;
    }

    *c = *a / control_3dof_B.scale;
    *s = *b / control_3dof_B.scale;
    if (control_3dof_B.absa > control_3dof_B.absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = control_3dof_B.scale;
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xrot(real_T x[9], int32_T ix0, int32_T iy0,
  real_T c, real_T s)
{
  control_3dof_B.temp = x[iy0 - 1];
  control_3dof_B.temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = control_3dof_B.temp * c - control_3dof_B.temp_tmp * s;
  x[ix0 - 1] = control_3dof_B.temp_tmp * c + control_3dof_B.temp * s;
  control_3dof_B.temp = x[ix0] * c + x[iy0] * s;
  x[iy0] = x[iy0] * c - x[ix0] * s;
  x[ix0] = control_3dof_B.temp;
  control_3dof_B.temp = x[iy0 + 1];
  control_3dof_B.temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = control_3dof_B.temp * c - control_3dof_B.temp_tmp * s;
  x[ix0 + 1] = control_3dof_B.temp_tmp * c + control_3dof_B.temp * s;
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_svd(const real_T A[9], real_T U[9], real_T s[3],
  real_T V[9])
{
  __m128d tmp;
  boolean_T exitg1;
  control_3dof_B.b_s[0] = 0.0;
  control_3dof_B.e[0] = 0.0;
  control_3dof_B.work[0] = 0.0;
  control_3dof_B.b_s[1] = 0.0;
  control_3dof_B.e[1] = 0.0;
  control_3dof_B.work[1] = 0.0;
  control_3dof_B.b_s[2] = 0.0;
  control_3dof_B.e[2] = 0.0;
  control_3dof_B.work[2] = 0.0;
  for (control_3dof_B.qjj = 0; control_3dof_B.qjj < 9; control_3dof_B.qjj++) {
    control_3dof_B.b_A[control_3dof_B.qjj] = A[control_3dof_B.qjj];
    U[control_3dof_B.qjj] = 0.0;
    control_3dof_B.Vf[control_3dof_B.qjj] = 0.0;
  }

  control_3dof_B.doscale = false;
  control_3dof_B.anrm = control_3dof_xzlangeM(A);
  control_3dof_B.cscale = control_3dof_B.anrm;
  if ((control_3dof_B.anrm > 0.0) && (control_3dof_B.anrm <
       6.7178761075670888E-139)) {
    control_3dof_B.doscale = true;
    control_3dof_B.cscale = 6.7178761075670888E-139;
    control_3dof_xzlascl(control_3dof_B.anrm, control_3dof_B.cscale, 3, 3,
                         control_3dof_B.b_A, 1, 3);
  } else if (control_3dof_B.anrm > 1.4885657073574029E+138) {
    control_3dof_B.doscale = true;
    control_3dof_B.cscale = 1.4885657073574029E+138;
    control_3dof_xzlascl(control_3dof_B.anrm, control_3dof_B.cscale, 3, 3,
                         control_3dof_B.b_A, 1, 3);
  }

  for (control_3dof_B.m = 0; control_3dof_B.m < 2; control_3dof_B.m++) {
    control_3dof_B.qp1 = control_3dof_B.m + 2;
    control_3dof_B.qq_tmp = 3 * control_3dof_B.m + control_3dof_B.m;
    control_3dof_B.qq = control_3dof_B.qq_tmp + 1;
    control_3dof_B.apply_transform = false;
    control_3dof_B.nrm = control_3dof_xnrm2(3 - control_3dof_B.m,
      control_3dof_B.b_A, control_3dof_B.qq_tmp + 1);
    if (control_3dof_B.nrm > 0.0) {
      control_3dof_B.apply_transform = true;
      if (control_3dof_B.b_A[control_3dof_B.qq_tmp] < 0.0) {
        control_3dof_B.nrm = -control_3dof_B.nrm;
      }

      control_3dof_B.b_s[control_3dof_B.m] = control_3dof_B.nrm;
      if (std::abs(control_3dof_B.nrm) >= 1.0020841800044864E-292) {
        control_3dof_B.nrm = 1.0 / control_3dof_B.nrm;
        control_3dof_B.d = (control_3dof_B.qq_tmp - control_3dof_B.m) + 3;
        control_3dof_B.qjj = ((((control_3dof_B.d - control_3dof_B.qq_tmp) / 2) <<
          1) + control_3dof_B.qq_tmp) + 1;
        control_3dof_B.kase = control_3dof_B.qjj - 2;
        for (control_3dof_B.e_k = control_3dof_B.qq; control_3dof_B.e_k <=
             control_3dof_B.kase; control_3dof_B.e_k += 2) {
          tmp = _mm_loadu_pd(&control_3dof_B.b_A[control_3dof_B.e_k - 1]);
          _mm_storeu_pd(&control_3dof_B.b_A[control_3dof_B.e_k - 1], _mm_mul_pd
                        (tmp, _mm_set1_pd(control_3dof_B.nrm)));
        }

        for (control_3dof_B.e_k = control_3dof_B.qjj; control_3dof_B.e_k <=
             control_3dof_B.d; control_3dof_B.e_k++) {
          control_3dof_B.b_A[control_3dof_B.e_k - 1] *= control_3dof_B.nrm;
        }
      } else {
        control_3dof_B.d = (control_3dof_B.qq_tmp - control_3dof_B.m) + 3;
        control_3dof_B.qjj = ((((control_3dof_B.d - control_3dof_B.qq_tmp) / 2) <<
          1) + control_3dof_B.qq_tmp) + 1;
        control_3dof_B.kase = control_3dof_B.qjj - 2;
        for (control_3dof_B.e_k = control_3dof_B.qq; control_3dof_B.e_k <=
             control_3dof_B.kase; control_3dof_B.e_k += 2) {
          tmp = _mm_loadu_pd(&control_3dof_B.b_A[control_3dof_B.e_k - 1]);
          _mm_storeu_pd(&control_3dof_B.b_A[control_3dof_B.e_k - 1], _mm_div_pd
                        (tmp, _mm_set1_pd(control_3dof_B.b_s[control_3dof_B.m])));
        }

        for (control_3dof_B.e_k = control_3dof_B.qjj; control_3dof_B.e_k <=
             control_3dof_B.d; control_3dof_B.e_k++) {
          control_3dof_B.b_A[control_3dof_B.e_k - 1] /=
            control_3dof_B.b_s[control_3dof_B.m];
        }
      }

      control_3dof_B.b_A[control_3dof_B.qq_tmp]++;
      control_3dof_B.b_s[control_3dof_B.m] =
        -control_3dof_B.b_s[control_3dof_B.m];
    } else {
      control_3dof_B.b_s[control_3dof_B.m] = 0.0;
    }

    for (control_3dof_B.kase = control_3dof_B.qp1; control_3dof_B.kase < 4;
         control_3dof_B.kase++) {
      control_3dof_B.qjj = (control_3dof_B.kase - 1) * 3 + control_3dof_B.m;
      if (control_3dof_B.apply_transform) {
        control_3dof_xaxpy(3 - control_3dof_B.m, -(control_3dof_xdotc(3 -
          control_3dof_B.m, control_3dof_B.b_A, control_3dof_B.qq_tmp + 1,
          control_3dof_B.b_A, control_3dof_B.qjj + 1) /
          control_3dof_B.b_A[control_3dof_B.qq_tmp]), control_3dof_B.qq_tmp + 1,
                           control_3dof_B.b_A, control_3dof_B.qjj + 1);
      }

      control_3dof_B.e[control_3dof_B.kase - 1] =
        control_3dof_B.b_A[control_3dof_B.qjj];
    }

    for (control_3dof_B.kase = control_3dof_B.m + 1; control_3dof_B.kase < 4;
         control_3dof_B.kase++) {
      control_3dof_B.qjj = (3 * control_3dof_B.m + control_3dof_B.kase) - 1;
      U[control_3dof_B.qjj] = control_3dof_B.b_A[control_3dof_B.qjj];
    }

    if (control_3dof_B.m + 1 <= 1) {
      control_3dof_B.nrm = control_3dof_xnrm2_n(2, control_3dof_B.e, 2);
      if (control_3dof_B.nrm == 0.0) {
        control_3dof_B.e[0] = 0.0;
      } else {
        if (control_3dof_B.e[1] < 0.0) {
          control_3dof_B.e[0] = -control_3dof_B.nrm;
        } else {
          control_3dof_B.e[0] = control_3dof_B.nrm;
        }

        control_3dof_B.nrm = control_3dof_B.e[0];
        if (std::abs(control_3dof_B.e[0]) >= 1.0020841800044864E-292) {
          control_3dof_B.nrm = 1.0 / control_3dof_B.e[0];
          control_3dof_B.qjj = ((((2 - control_3dof_B.m) / 2) << 1) +
                                control_3dof_B.m) + 2;
          control_3dof_B.kase = control_3dof_B.qjj - 2;
          for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq <=
               control_3dof_B.kase; control_3dof_B.qq += 2) {
            tmp = _mm_loadu_pd(&control_3dof_B.e[control_3dof_B.qq - 1]);
            _mm_storeu_pd(&control_3dof_B.e[control_3dof_B.qq - 1], _mm_mul_pd
                          (tmp, _mm_set1_pd(control_3dof_B.nrm)));
          }

          for (control_3dof_B.qq = control_3dof_B.qjj; control_3dof_B.qq < 4;
               control_3dof_B.qq++) {
            control_3dof_B.e[control_3dof_B.qq - 1] *= control_3dof_B.nrm;
          }
        } else {
          control_3dof_B.qjj = ((((2 - control_3dof_B.m) / 2) << 1) +
                                control_3dof_B.m) + 2;
          control_3dof_B.kase = control_3dof_B.qjj - 2;
          for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq <=
               control_3dof_B.kase; control_3dof_B.qq += 2) {
            tmp = _mm_loadu_pd(&control_3dof_B.e[control_3dof_B.qq - 1]);
            _mm_storeu_pd(&control_3dof_B.e[control_3dof_B.qq - 1], _mm_div_pd
                          (tmp, _mm_set1_pd(control_3dof_B.nrm)));
          }

          for (control_3dof_B.qq = control_3dof_B.qjj; control_3dof_B.qq < 4;
               control_3dof_B.qq++) {
            control_3dof_B.e[control_3dof_B.qq - 1] /= control_3dof_B.nrm;
          }
        }

        control_3dof_B.e[1]++;
        control_3dof_B.e[0] = -control_3dof_B.e[0];
        for (control_3dof_B.qjj = control_3dof_B.qp1; control_3dof_B.qjj < 4;
             control_3dof_B.qjj++) {
          control_3dof_B.work[control_3dof_B.qjj - 1] = 0.0;
        }

        for (control_3dof_B.qjj = control_3dof_B.qp1; control_3dof_B.qjj < 4;
             control_3dof_B.qjj++) {
          control_3dof_xaxpy_f(2, control_3dof_B.e[control_3dof_B.qjj - 1],
                               control_3dof_B.b_A, 3 * (control_3dof_B.qjj - 1)
                               + 2, control_3dof_B.work, 2);
        }

        for (control_3dof_B.qjj = control_3dof_B.qp1; control_3dof_B.qjj < 4;
             control_3dof_B.qjj++) {
          control_3dof_xaxpy_fh(2, -control_3dof_B.e[control_3dof_B.qjj - 1] /
                                control_3dof_B.e[1], control_3dof_B.work, 2,
                                control_3dof_B.b_A, 3 * (control_3dof_B.qjj - 1)
                                + 2);
        }
      }

      for (control_3dof_B.qjj = control_3dof_B.qp1; control_3dof_B.qjj < 4;
           control_3dof_B.qjj++) {
        control_3dof_B.Vf[control_3dof_B.qjj - 1] =
          control_3dof_B.e[control_3dof_B.qjj - 1];
      }
    }
  }

  control_3dof_B.m = 1;
  control_3dof_B.b_s[2] = control_3dof_B.b_A[8];
  control_3dof_B.e[1] = control_3dof_B.b_A[7];
  control_3dof_B.e[2] = 0.0;
  U[6] = 0.0;
  U[7] = 0.0;
  U[8] = 1.0;
  for (control_3dof_B.qp1 = 1; control_3dof_B.qp1 >= 0; control_3dof_B.qp1--) {
    control_3dof_B.qq = 3 * control_3dof_B.qp1 + control_3dof_B.qp1;
    if (control_3dof_B.b_s[control_3dof_B.qp1] != 0.0) {
      for (control_3dof_B.kase = control_3dof_B.qp1 + 2; control_3dof_B.kase < 4;
           control_3dof_B.kase++) {
        control_3dof_B.qjj = ((control_3dof_B.kase - 1) * 3 + control_3dof_B.qp1)
          + 1;
        control_3dof_xaxpy(3 - control_3dof_B.qp1, -(control_3dof_xdotc(3 -
          control_3dof_B.qp1, U, control_3dof_B.qq + 1, U, control_3dof_B.qjj) /
          U[control_3dof_B.qq]), control_3dof_B.qq + 1, U, control_3dof_B.qjj);
      }

      for (control_3dof_B.kase = control_3dof_B.qp1 + 1; control_3dof_B.kase < 4;
           control_3dof_B.kase++) {
        control_3dof_B.qjj = (3 * control_3dof_B.qp1 + control_3dof_B.kase) - 1;
        U[control_3dof_B.qjj] = -U[control_3dof_B.qjj];
      }

      U[control_3dof_B.qq]++;
      if (control_3dof_B.qp1 - 1 >= 0) {
        U[3 * control_3dof_B.qp1] = 0.0;
      }
    } else {
      U[3 * control_3dof_B.qp1] = 0.0;
      U[3 * control_3dof_B.qp1 + 1] = 0.0;
      U[3 * control_3dof_B.qp1 + 2] = 0.0;
      U[control_3dof_B.qq] = 1.0;
    }
  }

  for (control_3dof_B.qp1 = 2; control_3dof_B.qp1 >= 0; control_3dof_B.qp1--) {
    if ((control_3dof_B.qp1 + 1 <= 1) && (control_3dof_B.e[0] != 0.0)) {
      control_3dof_xaxpy(2, -(control_3dof_xdotc(2, control_3dof_B.Vf, 2,
        control_3dof_B.Vf, 5) / control_3dof_B.Vf[1]), 2, control_3dof_B.Vf, 5);
      control_3dof_xaxpy(2, -(control_3dof_xdotc(2, control_3dof_B.Vf, 2,
        control_3dof_B.Vf, 8) / control_3dof_B.Vf[1]), 2, control_3dof_B.Vf, 8);
    }

    control_3dof_B.Vf[3 * control_3dof_B.qp1] = 0.0;
    control_3dof_B.Vf[3 * control_3dof_B.qp1 + 1] = 0.0;
    control_3dof_B.Vf[3 * control_3dof_B.qp1 + 2] = 0.0;
    control_3dof_B.Vf[control_3dof_B.qp1 + 3 * control_3dof_B.qp1] = 1.0;
  }

  control_3dof_B.qp1 = 0;
  control_3dof_B.nrm = 0.0;
  for (control_3dof_B.qq = 0; control_3dof_B.qq < 3; control_3dof_B.qq++) {
    control_3dof_B.r = control_3dof_B.b_s[control_3dof_B.qq];
    if (control_3dof_B.r != 0.0) {
      control_3dof_B.rt = std::abs(control_3dof_B.r);
      control_3dof_B.r /= control_3dof_B.rt;
      control_3dof_B.b_s[control_3dof_B.qq] = control_3dof_B.rt;
      if (control_3dof_B.qq + 1 < 3) {
        control_3dof_B.e[control_3dof_B.qq] /= control_3dof_B.r;
      }

      control_3dof_B.qq_tmp = 3 * control_3dof_B.qq + 1;
      control_3dof_B.qjj = 2 + control_3dof_B.qq_tmp;
      control_3dof_B.kase = control_3dof_B.qq_tmp;
      for (control_3dof_B.d = control_3dof_B.qq_tmp; control_3dof_B.d <=
           control_3dof_B.kase; control_3dof_B.d += 2) {
        tmp = _mm_loadu_pd(&U[control_3dof_B.d - 1]);
        _mm_storeu_pd(&U[control_3dof_B.d - 1], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.r)));
      }

      for (control_3dof_B.d = control_3dof_B.qjj; control_3dof_B.d <=
           control_3dof_B.qq_tmp + 2; control_3dof_B.d++) {
        U[control_3dof_B.d - 1] *= control_3dof_B.r;
      }
    }

    if (control_3dof_B.qq + 1 < 3) {
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq];
      if (control_3dof_B.r != 0.0) {
        control_3dof_B.rt = std::abs(control_3dof_B.r);
        control_3dof_B.r = control_3dof_B.rt / control_3dof_B.r;
        control_3dof_B.e[control_3dof_B.qq] = control_3dof_B.rt;
        control_3dof_B.b_s[control_3dof_B.qq + 1] *= control_3dof_B.r;
        control_3dof_B.qq_tmp = (control_3dof_B.qq + 1) * 3 + 1;
        control_3dof_B.qjj = 2 + control_3dof_B.qq_tmp;
        control_3dof_B.kase = control_3dof_B.qq_tmp;
        for (control_3dof_B.d = control_3dof_B.qq_tmp; control_3dof_B.d <=
             control_3dof_B.kase; control_3dof_B.d += 2) {
          tmp = _mm_loadu_pd(&control_3dof_B.Vf[control_3dof_B.d - 1]);
          _mm_storeu_pd(&control_3dof_B.Vf[control_3dof_B.d - 1], _mm_mul_pd(tmp,
            _mm_set1_pd(control_3dof_B.r)));
        }

        for (control_3dof_B.d = control_3dof_B.qjj; control_3dof_B.d <=
             control_3dof_B.qq_tmp + 2; control_3dof_B.d++) {
          control_3dof_B.Vf[control_3dof_B.d - 1] *= control_3dof_B.r;
        }
      }
    }

    control_3dof_B.nrm = std::fmax(control_3dof_B.nrm, std::fmax(std::abs
      (control_3dof_B.b_s[control_3dof_B.qq]), std::abs
      (control_3dof_B.e[control_3dof_B.qq])));
  }

  while ((control_3dof_B.m + 2 > 0) && (control_3dof_B.qp1 < 75)) {
    control_3dof_B.qq = control_3dof_B.m + 1;
    exitg1 = false;
    while (!(exitg1 || (control_3dof_B.qq == 0))) {
      control_3dof_B.rt = std::abs(control_3dof_B.e[control_3dof_B.qq - 1]);
      if (control_3dof_B.rt <= (std::abs(control_3dof_B.b_s[control_3dof_B.qq -
            1]) + std::abs(control_3dof_B.b_s[control_3dof_B.qq])) *
          2.2204460492503131E-16) {
        control_3dof_B.e[control_3dof_B.qq - 1] = 0.0;
        exitg1 = true;
      } else if ((control_3dof_B.rt <= 1.0020841800044864E-292) ||
                 ((control_3dof_B.qp1 > 20) && (control_3dof_B.rt <=
                   2.2204460492503131E-16 * control_3dof_B.nrm))) {
        control_3dof_B.e[control_3dof_B.qq - 1] = 0.0;
        exitg1 = true;
      } else {
        control_3dof_B.qq--;
      }
    }

    if (control_3dof_B.m + 1 == control_3dof_B.qq) {
      control_3dof_B.kase = 4;
    } else {
      control_3dof_B.qjj = control_3dof_B.m + 2;
      control_3dof_B.kase = control_3dof_B.m + 2;
      exitg1 = false;
      while ((!exitg1) && (control_3dof_B.kase >= control_3dof_B.qq)) {
        control_3dof_B.qjj = control_3dof_B.kase;
        if (control_3dof_B.kase == control_3dof_B.qq) {
          exitg1 = true;
        } else {
          control_3dof_B.rt = 0.0;
          if (control_3dof_B.kase < control_3dof_B.m + 2) {
            control_3dof_B.rt = std::abs(control_3dof_B.e[control_3dof_B.kase -
              1]);
          }

          if (control_3dof_B.kase > control_3dof_B.qq + 1) {
            control_3dof_B.rt += std::abs(control_3dof_B.e[control_3dof_B.kase -
              2]);
          }

          control_3dof_B.r = std::abs(control_3dof_B.b_s[control_3dof_B.kase - 1]);
          if ((control_3dof_B.r <= 2.2204460492503131E-16 * control_3dof_B.rt) ||
              (control_3dof_B.r <= 1.0020841800044864E-292)) {
            control_3dof_B.b_s[control_3dof_B.kase - 1] = 0.0;
            exitg1 = true;
          } else {
            control_3dof_B.kase--;
          }
        }
      }

      if (control_3dof_B.qjj == control_3dof_B.qq) {
        control_3dof_B.kase = 3;
      } else if (control_3dof_B.m + 2 == control_3dof_B.qjj) {
        control_3dof_B.kase = 1;
      } else {
        control_3dof_B.kase = 2;
        control_3dof_B.qq = control_3dof_B.qjj;
      }
    }

    switch (control_3dof_B.kase) {
     case 1:
      control_3dof_B.rt = control_3dof_B.e[control_3dof_B.m];
      control_3dof_B.e[control_3dof_B.m] = 0.0;
      for (control_3dof_B.qjj = control_3dof_B.m + 1; control_3dof_B.qjj >=
           control_3dof_B.qq + 1; control_3dof_B.qjj--) {
        control_3dof_xrotg(&control_3dof_B.b_s[control_3dof_B.qjj - 1],
                           &control_3dof_B.rt, &control_3dof_B.r,
                           &control_3dof_B.smm1);
        if (control_3dof_B.qjj > control_3dof_B.qq + 1) {
          control_3dof_B.rt = -control_3dof_B.smm1 * control_3dof_B.e[0];
          control_3dof_B.e[0] *= control_3dof_B.r;
        }

        control_3dof_xrot(control_3dof_B.Vf, 3 * (control_3dof_B.qjj - 1) + 1, 3
                          * (control_3dof_B.m + 1) + 1, control_3dof_B.r,
                          control_3dof_B.smm1);
      }
      break;

     case 2:
      control_3dof_B.rt = control_3dof_B.e[control_3dof_B.qq - 1];
      control_3dof_B.e[control_3dof_B.qq - 1] = 0.0;
      for (control_3dof_B.qjj = control_3dof_B.qq + 1; control_3dof_B.qjj <=
           control_3dof_B.m + 2; control_3dof_B.qjj++) {
        control_3dof_xrotg(&control_3dof_B.b_s[control_3dof_B.qjj - 1],
                           &control_3dof_B.rt, &control_3dof_B.smm1,
                           &control_3dof_B.c);
        control_3dof_B.r = control_3dof_B.e[control_3dof_B.qjj - 1];
        control_3dof_B.rt = -control_3dof_B.c * control_3dof_B.r;
        control_3dof_B.e[control_3dof_B.qjj - 1] = control_3dof_B.r *
          control_3dof_B.smm1;
        control_3dof_xrot(U, 3 * (control_3dof_B.qjj - 1) + 1, 3 *
                          (control_3dof_B.qq - 1) + 1, control_3dof_B.smm1,
                          control_3dof_B.c);
      }
      break;

     case 3:
      control_3dof_B.rt = control_3dof_B.b_s[control_3dof_B.m + 1];
      control_3dof_B.r = std::fmax(std::fmax(std::fmax(std::fmax(std::abs
        (control_3dof_B.rt), std::abs(control_3dof_B.b_s[control_3dof_B.m])),
        std::abs(control_3dof_B.e[control_3dof_B.m])), std::abs
        (control_3dof_B.b_s[control_3dof_B.qq])), std::abs
        (control_3dof_B.e[control_3dof_B.qq]));
      tmp = _mm_set1_pd(control_3dof_B.r);
      _mm_storeu_pd(&control_3dof_B.dv3[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.m], control_3dof_B.rt), tmp));
      control_3dof_B.rt = control_3dof_B.dv3[0];
      control_3dof_B.smm1 = control_3dof_B.dv3[1];
      _mm_storeu_pd(&control_3dof_B.dv3[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.qq],
         control_3dof_B.e[control_3dof_B.m]), tmp));
      control_3dof_B.smm1 = ((control_3dof_B.smm1 + control_3dof_B.rt) *
        (control_3dof_B.smm1 - control_3dof_B.rt) + control_3dof_B.dv3[0] *
        control_3dof_B.dv3[0]) / 2.0;
      control_3dof_B.c = control_3dof_B.rt * control_3dof_B.dv3[0];
      control_3dof_B.c *= control_3dof_B.c;
      if ((control_3dof_B.smm1 != 0.0) || (control_3dof_B.c != 0.0)) {
        control_3dof_B.shift = std::sqrt(control_3dof_B.smm1 *
          control_3dof_B.smm1 + control_3dof_B.c);
        if (control_3dof_B.smm1 < 0.0) {
          control_3dof_B.shift = -control_3dof_B.shift;
        }

        control_3dof_B.shift = control_3dof_B.c / (control_3dof_B.smm1 +
          control_3dof_B.shift);
      } else {
        control_3dof_B.shift = 0.0;
      }

      control_3dof_B.rt = (control_3dof_B.dv3[1] + control_3dof_B.rt) *
        (control_3dof_B.dv3[1] - control_3dof_B.rt) + control_3dof_B.shift;
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq] / control_3dof_B.r *
        control_3dof_B.dv3[1];
      for (control_3dof_B.qq_tmp = control_3dof_B.qq + 1; control_3dof_B.qq_tmp <=
           control_3dof_B.m + 1; control_3dof_B.qq_tmp++) {
        control_3dof_xrotg(&control_3dof_B.rt, &control_3dof_B.r,
                           &control_3dof_B.smm1, &control_3dof_B.c);
        if (control_3dof_B.qq_tmp > control_3dof_B.qq + 1) {
          control_3dof_B.e[0] = control_3dof_B.rt;
        }

        control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq_tmp - 1];
        control_3dof_B.shift = control_3dof_B.b_s[control_3dof_B.qq_tmp - 1];
        control_3dof_B.e[control_3dof_B.qq_tmp - 1] = control_3dof_B.r *
          control_3dof_B.smm1 - control_3dof_B.shift * control_3dof_B.c;
        control_3dof_B.rt = control_3dof_B.c *
          control_3dof_B.b_s[control_3dof_B.qq_tmp];
        control_3dof_B.b_s[control_3dof_B.qq_tmp] *= control_3dof_B.smm1;
        control_3dof_B.qjj = (control_3dof_B.qq_tmp - 1) * 3 + 1;
        control_3dof_B.kase = 3 * control_3dof_B.qq_tmp + 1;
        control_3dof_xrot(control_3dof_B.Vf, control_3dof_B.qjj,
                          control_3dof_B.kase, control_3dof_B.smm1,
                          control_3dof_B.c);
        control_3dof_B.b_s[control_3dof_B.qq_tmp - 1] = control_3dof_B.shift *
          control_3dof_B.smm1 + control_3dof_B.r * control_3dof_B.c;
        control_3dof_xrotg(&control_3dof_B.b_s[control_3dof_B.qq_tmp - 1],
                           &control_3dof_B.rt, &control_3dof_B.smm1,
                           &control_3dof_B.c);
        control_3dof_B.shift = control_3dof_B.e[control_3dof_B.qq_tmp - 1];
        control_3dof_B.rt = control_3dof_B.shift * control_3dof_B.smm1 +
          control_3dof_B.c * control_3dof_B.b_s[control_3dof_B.qq_tmp];
        control_3dof_B.b_s[control_3dof_B.qq_tmp] = control_3dof_B.shift *
          -control_3dof_B.c + control_3dof_B.smm1 *
          control_3dof_B.b_s[control_3dof_B.qq_tmp];
        control_3dof_B.r = control_3dof_B.c *
          control_3dof_B.e[control_3dof_B.qq_tmp];
        control_3dof_B.e[control_3dof_B.qq_tmp] *= control_3dof_B.smm1;
        control_3dof_xrot(U, control_3dof_B.qjj, control_3dof_B.kase,
                          control_3dof_B.smm1, control_3dof_B.c);
      }

      control_3dof_B.e[control_3dof_B.m] = control_3dof_B.rt;
      control_3dof_B.qp1++;
      break;

     default:
      if (control_3dof_B.b_s[control_3dof_B.qq] < 0.0) {
        control_3dof_B.b_s[control_3dof_B.qq] =
          -control_3dof_B.b_s[control_3dof_B.qq];
        control_3dof_B.qp1 = 3 * control_3dof_B.qq + 1;
        control_3dof_B.qjj = 2 + control_3dof_B.qp1;
        control_3dof_B.kase = control_3dof_B.qp1;
        for (control_3dof_B.qq_tmp = control_3dof_B.qp1; control_3dof_B.qq_tmp <=
             control_3dof_B.kase; control_3dof_B.qq_tmp += 2) {
          tmp = _mm_loadu_pd(&control_3dof_B.Vf[control_3dof_B.qq_tmp - 1]);
          _mm_storeu_pd(&control_3dof_B.Vf[control_3dof_B.qq_tmp - 1],
                        _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
        }

        for (control_3dof_B.qq_tmp = control_3dof_B.qjj; control_3dof_B.qq_tmp <=
             control_3dof_B.qp1 + 2; control_3dof_B.qq_tmp++) {
          control_3dof_B.Vf[control_3dof_B.qq_tmp - 1] =
            -control_3dof_B.Vf[control_3dof_B.qq_tmp - 1];
        }
      }

      control_3dof_B.qp1 = control_3dof_B.qq + 1;
      while ((control_3dof_B.qq + 1 < 3) &&
             (control_3dof_B.b_s[control_3dof_B.qq] <
              control_3dof_B.b_s[control_3dof_B.qp1])) {
        control_3dof_B.rt = control_3dof_B.b_s[control_3dof_B.qq];
        control_3dof_B.b_s[control_3dof_B.qq] =
          control_3dof_B.b_s[control_3dof_B.qp1];
        control_3dof_B.b_s[control_3dof_B.qp1] = control_3dof_B.rt;
        control_3dof_B.qjj = 3 * control_3dof_B.qq + 1;
        control_3dof_B.kase = (control_3dof_B.qq + 1) * 3 + 1;
        control_3dof_xswap(control_3dof_B.Vf, control_3dof_B.qjj,
                           control_3dof_B.kase);
        control_3dof_xswap(U, control_3dof_B.qjj, control_3dof_B.kase);
        control_3dof_B.qq = control_3dof_B.qp1;
        control_3dof_B.qp1++;
      }

      control_3dof_B.qp1 = 0;
      control_3dof_B.m--;
      break;
    }
  }

  s[0] = control_3dof_B.b_s[0];
  s[1] = control_3dof_B.b_s[1];
  s[2] = control_3dof_B.b_s[2];
  if (control_3dof_B.doscale) {
    control_3dof_xzlascl_l(control_3dof_B.cscale, control_3dof_B.anrm, 3, 1, s,
      1, 3);
  }

  for (control_3dof_B.m = 0; control_3dof_B.m < 3; control_3dof_B.m++) {
    V[3 * control_3dof_B.m] = control_3dof_B.Vf[3 * control_3dof_B.m];
    control_3dof_B.qp1 = 3 * control_3dof_B.m + 1;
    V[control_3dof_B.qp1] = control_3dof_B.Vf[control_3dof_B.qp1];
    control_3dof_B.qp1 = 3 * control_3dof_B.m + 2;
    V[control_3dof_B.qp1] = control_3dof_B.Vf[control_3dof_B.qp1];
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_pinv(const real_T A[9], real_T X[9])
{
  __m128d tmp;
  __m128d tmp_0;
  boolean_T exitg1;
  control_3dof_B.p = true;
  for (control_3dof_B.r_o = 0; control_3dof_B.r_o < 9; control_3dof_B.r_o++) {
    X[control_3dof_B.r_o] = 0.0;
    if (control_3dof_B.p) {
      control_3dof_B.absx = A[control_3dof_B.r_o];
      if (std::isinf(control_3dof_B.absx) || std::isnan(control_3dof_B.absx)) {
        control_3dof_B.p = false;
      }
    }
  }

  if (!control_3dof_B.p) {
    for (control_3dof_B.r_o = 0; control_3dof_B.r_o < 9; control_3dof_B.r_o++) {
      X[control_3dof_B.r_o] = (rtNaN);
    }
  } else {
    control_3dof_svd(A, control_3dof_B.U, control_3dof_B.s, control_3dof_B.V);
    control_3dof_B.absx = std::abs(control_3dof_B.s[0]);
    if (std::isinf(control_3dof_B.absx) || std::isnan(control_3dof_B.absx)) {
      control_3dof_B.absx = (rtNaN);
    } else if (control_3dof_B.absx < 4.4501477170144028E-308) {
      control_3dof_B.absx = 4.94065645841247E-324;
    } else {
      std::frexp(control_3dof_B.absx, &control_3dof_B.exponent);
      control_3dof_B.absx = std::ldexp(1.0, control_3dof_B.exponent - 53);
    }

    control_3dof_B.absx *= 3.0;
    control_3dof_B.r_o = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.r_o < 3)) {
      if (std::isinf(control_3dof_B.s[control_3dof_B.r_o]) || std::isnan
          (control_3dof_B.s[control_3dof_B.r_o])) {
        control_3dof_B.absx = 1.7976931348623157E+308;
        exitg1 = true;
      } else {
        control_3dof_B.r_o++;
      }
    }

    control_3dof_B.r_o = -1;
    control_3dof_B.exponent = 0;
    while ((control_3dof_B.exponent < 3) &&
           (control_3dof_B.s[control_3dof_B.exponent] > control_3dof_B.absx)) {
      control_3dof_B.r_o++;
      control_3dof_B.exponent++;
    }

    if (control_3dof_B.r_o + 1 > 0) {
      control_3dof_B.vcol = 1;
      for (control_3dof_B.j = 0; control_3dof_B.j <= control_3dof_B.r_o;
           control_3dof_B.j++) {
        control_3dof_B.absx = 1.0 / control_3dof_B.s[control_3dof_B.j];
        control_3dof_B.exponent = 2 + control_3dof_B.vcol;
        control_3dof_B.vectorUB = control_3dof_B.vcol;
        for (control_3dof_B.ar = control_3dof_B.vcol; control_3dof_B.ar <=
             control_3dof_B.vectorUB; control_3dof_B.ar += 2) {
          tmp_0 = _mm_loadu_pd(&control_3dof_B.V[control_3dof_B.ar - 1]);
          _mm_storeu_pd(&control_3dof_B.V[control_3dof_B.ar - 1], _mm_mul_pd
                        (tmp_0, _mm_set1_pd(control_3dof_B.absx)));
        }

        for (control_3dof_B.ar = control_3dof_B.exponent; control_3dof_B.ar <=
             control_3dof_B.vcol + 2; control_3dof_B.ar++) {
          control_3dof_B.V[control_3dof_B.ar - 1] *= control_3dof_B.absx;
        }

        control_3dof_B.vcol += 3;
      }

      control_3dof_B.vcol = 0;
      for (control_3dof_B.exponent = 0; control_3dof_B.exponent <= 6;
           control_3dof_B.exponent += 3) {
        for (control_3dof_B.vectorUB = control_3dof_B.exponent + 1;
             control_3dof_B.vectorUB <= control_3dof_B.exponent + 3;
             control_3dof_B.vectorUB++) {
          X[control_3dof_B.vectorUB - 1] = 0.0;
        }
      }

      for (control_3dof_B.j = 0; control_3dof_B.j <= 6; control_3dof_B.j += 3) {
        control_3dof_B.ar = -1;
        control_3dof_B.vcol++;
        control_3dof_B.b = 3 * control_3dof_B.r_o + control_3dof_B.vcol;
        for (control_3dof_B.ib = control_3dof_B.vcol; control_3dof_B.ib <=
             control_3dof_B.b; control_3dof_B.ib += 3) {
          control_3dof_B.exponent = control_3dof_B.j + 3;
          control_3dof_B.vectorUB = control_3dof_B.j + 1;
          for (control_3dof_B.b_ic = control_3dof_B.j + 1; control_3dof_B.b_ic <=
               control_3dof_B.vectorUB; control_3dof_B.b_ic += 2) {
            tmp_0 = _mm_loadu_pd(&control_3dof_B.V[(control_3dof_B.ar +
              control_3dof_B.b_ic) - control_3dof_B.j]);
            tmp = _mm_loadu_pd(&X[control_3dof_B.b_ic - 1]);
            _mm_storeu_pd(&X[control_3dof_B.b_ic - 1], _mm_add_pd(_mm_mul_pd
              (tmp_0, _mm_set1_pd(control_3dof_B.U[control_3dof_B.ib - 1])), tmp));
          }

          for (control_3dof_B.b_ic = control_3dof_B.exponent;
               control_3dof_B.b_ic <= control_3dof_B.j + 3; control_3dof_B.b_ic
               ++) {
            X[control_3dof_B.b_ic - 1] += control_3dof_B.V[(control_3dof_B.ar +
              control_3dof_B.b_ic) - control_3dof_B.j] *
              control_3dof_B.U[control_3dof_B.ib - 1];
          }

          control_3dof_B.ar += 3;
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_TCISolver(const real_T F_sp[3], const real_T
  F_fixed[3], const real_T Q[9], real_T b_min_tension, real_T b_max_tension,
  real_T result[3], real_T F_actual[3])
{
  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  int32_T exitg1;
  boolean_T exitg2;
  control_3dof_B.found[0] = false;
  result[0] = 0.0;
  control_3dof_B.F_d_idx_0 = F_sp[0];
  control_3dof_B.found[1] = false;
  result[1] = 0.0;
  control_3dof_B.F_d_idx_1 = F_sp[1];
  control_3dof_B.found[2] = false;
  result[2] = 0.0;
  control_3dof_B.F_d_idx_2 = F_sp[2];
  std::memcpy(&control_3dof_B.Q_new_m[0], &Q[0], 9U * sizeof(real_T));
  do {
    exitg1 = 0;
    control_3dof_pinv(control_3dof_B.Q_new_m, control_3dof_B.dv);
    for (control_3dof_B.i_l = 0; control_3dof_B.i_l <= 0; control_3dof_B.i_l +=
         2) {
      tmp = _mm_loadu_pd(&control_3dof_B.dv[control_3dof_B.i_l + 3]);
      tmp_0 = _mm_loadu_pd(&control_3dof_B.dv[control_3dof_B.i_l]);
      tmp_1 = _mm_loadu_pd(&control_3dof_B.dv[control_3dof_B.i_l + 6]);
      _mm_storeu_pd(&control_3dof_B.vu_c[control_3dof_B.i_l], _mm_add_pd
                    (_mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd
        (control_3dof_B.F_d_idx_1)), _mm_mul_pd(tmp_0, _mm_set1_pd
        (control_3dof_B.F_d_idx_0))), _mm_mul_pd(tmp_1, _mm_set1_pd
        (control_3dof_B.F_d_idx_2))));
    }

    for (control_3dof_B.i_l = 2; control_3dof_B.i_l < 3; control_3dof_B.i_l++) {
      control_3dof_B.vu_c[control_3dof_B.i_l] =
        (control_3dof_B.dv[control_3dof_B.i_l + 3] * control_3dof_B.F_d_idx_1 +
         control_3dof_B.dv[control_3dof_B.i_l] * control_3dof_B.F_d_idx_0) +
        control_3dof_B.dv[control_3dof_B.i_l + 6] * control_3dof_B.F_d_idx_2;
    }

    control_3dof_B.is_valid = true;
    control_3dof_B.i_l = 0;
    control_3dof_B.b_i = 0;
    exitg2 = false;
    while ((!exitg2) && (control_3dof_B.b_i < 3)) {
      control_3dof_B.i_l = control_3dof_B.b_i;
      if (control_3dof_B.found[control_3dof_B.b_i]) {
        control_3dof_B.b_i++;
      } else if (-control_3dof_B.vu_c[control_3dof_B.b_i] < b_min_tension -
                 1.0E-5) {
        control_3dof_B.vu_c[control_3dof_B.b_i] = -b_min_tension;
        control_3dof_B.is_valid = false;
        exitg2 = true;
      } else if (-control_3dof_B.vu_c[control_3dof_B.b_i] > b_max_tension +
                 1.0E-5) {
        control_3dof_B.vu_c[control_3dof_B.b_i] = -b_max_tension;
        control_3dof_B.is_valid = false;
        exitg2 = true;
      } else {
        control_3dof_B.b_i++;
      }
    }

    if (control_3dof_B.is_valid) {
      exitg1 = 1;
    } else {
      control_3dof_B.found[control_3dof_B.i_l] = true;
      result[control_3dof_B.i_l] = control_3dof_B.vu_c[control_3dof_B.i_l];
      control_3dof_B.F_d_idx_0 = F_fixed[0];
      control_3dof_B.F_d_idx_1 = F_fixed[1];
      control_3dof_B.F_d_idx_2 = F_fixed[2];
      for (control_3dof_B.b_i = 0; control_3dof_B.b_i < 3; control_3dof_B.b_i++)
      {
        if (control_3dof_B.found[control_3dof_B.b_i]) {
          control_3dof_B.result = result[control_3dof_B.b_i];
          _mm_storeu_pd(&control_3dof_B.dv2[0], _mm_sub_pd(_mm_set_pd
            (control_3dof_B.F_d_idx_1, control_3dof_B.F_d_idx_0), _mm_mul_pd
            (_mm_loadu_pd(&Q[3 * control_3dof_B.b_i]), _mm_set1_pd
             (control_3dof_B.result))));
          control_3dof_B.F_d_idx_0 = control_3dof_B.dv2[0];
          control_3dof_B.F_d_idx_1 = control_3dof_B.dv2[1];
          control_3dof_B.F_d_idx_2 -= Q[3 * control_3dof_B.b_i + 2] *
            control_3dof_B.result;
        }

        control_3dof_B.Q_new_m[control_3dof_B.b_i + 3 * control_3dof_B.i_l] =
          -(F_sp[control_3dof_B.b_i] - F_fixed[control_3dof_B.b_i]);
      }
    }
  } while (exitg1 == 0);

  if (!control_3dof_B.found[0]) {
    result[0] = control_3dof_B.vu_c[0];
  }

  if (!control_3dof_B.found[1]) {
    result[1] = control_3dof_B.vu_c[1];
  }

  if (!control_3dof_B.found[2]) {
    result[2] = control_3dof_B.vu_c[2];
  }

  control_3dof_B.result = result[1];
  control_3dof_B.F_d_idx_0 = result[0];
  control_3dof_B.F_d_idx_1 = result[2];
  for (control_3dof_B.i_l = 0; control_3dof_B.i_l <= 0; control_3dof_B.i_l += 2)
  {
    _mm_storeu_pd(&F_actual[control_3dof_B.i_l], _mm_add_pd(_mm_add_pd
      (_mm_mul_pd(_mm_loadu_pd(&Q[control_3dof_B.i_l + 3]), _mm_set1_pd
                  (control_3dof_B.result)), _mm_mul_pd(_mm_loadu_pd
      (&Q[control_3dof_B.i_l]), _mm_set1_pd(control_3dof_B.F_d_idx_0))),
      _mm_mul_pd(_mm_loadu_pd(&Q[control_3dof_B.i_l + 6]), _mm_set1_pd
                 (control_3dof_B.F_d_idx_1))));
  }

  for (control_3dof_B.i_l = 2; control_3dof_B.i_l < 3; control_3dof_B.i_l++) {
    F_actual[control_3dof_B.i_l] = (Q[control_3dof_B.i_l + 3] *
      control_3dof_B.result + Q[control_3dof_B.i_l] * control_3dof_B.F_d_idx_0)
      + Q[control_3dof_B.i_l + 6] * control_3dof_B.F_d_idx_1;
  }
}

// Model step function
void control_3dof::step()
{
  __m128d tmp;
  __m128d tmp_0;
  boolean_T exitg1;

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[0] = control_3dof_U.Dir_sp.q_sp1[0];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[3] = control_3dof_U.Dir_sp.q_sp2[0];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[6] = control_3dof_U.Dir_sp.q_sp3[0];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[1] = control_3dof_U.Dir_sp.q_sp1[1];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[4] = control_3dof_U.Dir_sp.q_sp2[1];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[7] = control_3dof_U.Dir_sp.q_sp3[1];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[2] = control_3dof_U.Dir_sp.q_sp1[2];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[5] = control_3dof_U.Dir_sp.q_sp2[2];

  // SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Dir_sp'

  control_3dof_B.Q_new[8] = control_3dof_U.Dir_sp.q_sp3[2];

  // RelationalOperator: '<S21>/FixPt Relational Operator' incorporates:
  //   Inport: '<Root>/Dir_sp'
  //   UnitDelay: '<S21>/Delay Input1'
  //
  //  Block description for '<S21>/Delay Input1':
  //
  //   Store in Global RAM

  control_3dof_B.FixPtRelationalOperator = (control_3dof_U.Dir_sp.timestamp !=
    control_3dof_DW.DelayInput1_DSTATE);
  for (control_3dof_B.i = 0; control_3dof_B.i < 9; control_3dof_B.i++) {
    // Delay: '<S20>/Delay' incorporates:
    //   Concatenate: '<S20>/Vector Concatenate'
    //   Switch: '<S20>/Switch'

    if (control_3dof_DW.icLoad) {
      control_3dof_DW.Delay_DSTATE[control_3dof_B.i] =
        control_3dof_B.Q_new[control_3dof_B.i];
    }

    // End of Delay: '<S20>/Delay'

    // Switch: '<S20>/Switch' incorporates:
    //   Concatenate: '<S20>/Vector Concatenate'

    if (control_3dof_B.FixPtRelationalOperator) {
      control_3dof_DW.Delay_DSTATE[control_3dof_B.i] =
        control_3dof_B.Q_new[control_3dof_B.i];
    }

    // End of Switch: '<S20>/Switch'
  }

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'

  control_3dof_DW.obj.K_p[0] = CONTROL_PARAM.ESO_PL[0];
  control_3dof_DW.obj.K_v[0] = CONTROL_PARAM.ESO_VL[0];
  control_3dof_DW.obj.K_p[1] = CONTROL_PARAM.ESO_PL[1];
  control_3dof_DW.obj.K_v[1] = CONTROL_PARAM.ESO_VL[1];
  control_3dof_DW.obj.K_p[2] = CONTROL_PARAM.ESO_PL[2];
  control_3dof_DW.obj.K_v[2] = CONTROL_PARAM.ESO_VL[2];
  if (!control_3dof_DW.obj.is_init) {
    control_3dof_DW.obj.is_init = true;
    control_3dof_DW.obj.Z_1[0] = control_3dof_U.Payload_Out.pL[0];
    control_3dof_DW.obj.Z_2[0] = control_3dof_U.Payload_Out.vL[0];
    control_3dof_DW.obj.Z_3[0] = 0.0;
    control_3dof_DW.obj.Z_1[1] = control_3dof_U.Payload_Out.pL[1];
    control_3dof_DW.obj.Z_2[1] = control_3dof_U.Payload_Out.vL[1];
    control_3dof_DW.obj.Z_3[1] = 0.0;
    control_3dof_DW.obj.Z_1[2] = control_3dof_U.Payload_Out.pL[2];
    control_3dof_DW.obj.Z_2[2] = control_3dof_U.Payload_Out.vL[2];
    control_3dof_DW.obj.Z_3[2] = 0.0;
  }

  // Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'

  //  Extended state observer
  //  Update states
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_sub_pd(_mm_set_pd
    (control_3dof_U.Payload_Out.vL[0], control_3dof_U.Payload_Out.pL[0]),
    _mm_set_pd(control_3dof_DW.obj.Z_2[0], control_3dof_DW.obj.Z_1[0])));

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_DW.obj.Z_1[0] += (control_3dof_DW.obj.K_p[0] *
    control_3dof_B.dv1[0] + control_3dof_DW.obj.Z_2[0]) * 0.02;
  control_3dof_DW.obj.Z_2[0] += (((control_3dof_DW.UnitDelay_DSTATE[0] +
    control_3dof_DW.obj.Z_3[0]) + control_3dof_B.dv1[0] *
    control_3dof_DW.obj.K_p[1]) + control_3dof_DW.obj.K_v[1] *
    control_3dof_B.dv1[1]) * 0.02;
  control_3dof_DW.obj.Z_3[0] += (control_3dof_B.dv1[0] *
    control_3dof_DW.obj.K_p[2] + control_3dof_B.dv1[1] *
    control_3dof_DW.obj.K_v[2]) * 0.02;

  // Gain: '<S19>/KP' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Sum: '<S19>/Sum'
  //
  control_3dof_B.u0 = (control_3dof_U.Traj_sp.pos_sp[0] -
                       control_3dof_U.Payload_Out.pL[0]) * CONTROL_PARAM.KP;

  // Saturate: '<S19>/vel_limit'
  if (control_3dof_B.u0 > 5.0) {
    control_3dof_B.u0 = 5.0;
  } else if (control_3dof_B.u0 < -5.0) {
    control_3dof_B.u0 = -5.0;
  }

  // Sum: '<S19>/Sum1' incorporates:
  //   Gain: '<S19>/KV'
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Saturate: '<S19>/vel_limit'
  //   Sum: '<S19>/Sum2'
  //   Sum: '<S19>/Sum3'
  //
  control_3dof_B.u0 = ((control_3dof_B.u0 + control_3dof_U.Traj_sp.vel_sp[0]) -
                       control_3dof_U.Payload_Out.vL[0]) * CONTROL_PARAM.KV +
    control_3dof_U.Traj_sp.acc_ff[0];

  // Saturate: '<S19>/acc_limit'
  if (control_3dof_B.u0 > 3.0) {
    control_3dof_B.u0 = 3.0;
  } else if (control_3dof_B.u0 < -3.0) {
    control_3dof_B.u0 = -3.0;
  }

  // Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_sub_pd(_mm_loadu_pd
    (&control_3dof_U.Payload_Out.pL[0]), _mm_set_pd(control_3dof_DW.obj.Z_1[1],
    control_3dof_U.Payload_Out.p_1[0])));

  // MATLAB Function: '<S4>/Force Saturation && Disturbution'
  control_3dof_B.e_1[0] = control_3dof_B.dv1[0];

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_B.e_2_b = control_3dof_U.Payload_Out.vL[1] -
    control_3dof_DW.obj.Z_2[1];
  control_3dof_DW.obj.Z_1[1] += (control_3dof_DW.obj.K_p[0] *
    control_3dof_B.dv1[1] + control_3dof_DW.obj.Z_2[1]) * 0.02;
  control_3dof_DW.obj.Z_2[1] += (((control_3dof_DW.UnitDelay_DSTATE[1] +
    control_3dof_DW.obj.Z_3[1]) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.dv1[1]) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_b) *
    0.02;
  control_3dof_DW.obj.Z_3[1] += (control_3dof_B.dv1[1] *
    control_3dof_DW.obj.K_p[2] + control_3dof_DW.obj.K_v[2] *
    control_3dof_B.e_2_b) * 0.02;

  // Gain: '<S19>/KP' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Sum: '<S19>/Sum'
  //
  control_3dof_B.u0_p = (control_3dof_U.Traj_sp.pos_sp[1] -
    control_3dof_U.Payload_Out.pL[1]) * CONTROL_PARAM.KP;

  // Saturate: '<S19>/vel_limit'
  if (control_3dof_B.u0_p > 5.0) {
    control_3dof_B.u0_p = 5.0;
  } else if (control_3dof_B.u0_p < -5.0) {
    control_3dof_B.u0_p = -5.0;
  }

  // Sum: '<S19>/Sum1' incorporates:
  //   Gain: '<S19>/KV'
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Saturate: '<S19>/vel_limit'
  //   Sum: '<S19>/Sum2'
  //   Sum: '<S19>/Sum3'
  //
  control_3dof_B.u0_p = ((control_3dof_B.u0_p + control_3dof_U.Traj_sp.vel_sp[1])
    - control_3dof_U.Payload_Out.vL[1]) * CONTROL_PARAM.KV +
    control_3dof_U.Traj_sp.acc_ff[1];

  // Saturate: '<S19>/acc_limit'
  if (control_3dof_B.u0_p > 3.0) {
    control_3dof_B.u0_p = 3.0;
  } else if (control_3dof_B.u0_p < -3.0) {
    control_3dof_B.u0_p = -3.0;
  }

  // Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_sub_pd(_mm_loadu_pd
    (&control_3dof_U.Payload_Out.pL[1]), _mm_set_pd(control_3dof_DW.obj.Z_1[2],
    control_3dof_U.Payload_Out.p_1[1])));

  // MATLAB Function: '<S4>/Force Saturation && Disturbution'
  control_3dof_B.e_1[1] = control_3dof_B.dv1[0];

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_B.e_2_b = control_3dof_U.Payload_Out.vL[2] -
    control_3dof_DW.obj.Z_2[2];
  control_3dof_DW.obj.Z_1[2] += (control_3dof_DW.obj.K_p[0] *
    control_3dof_B.dv1[1] + control_3dof_DW.obj.Z_2[2]) * 0.02;
  control_3dof_DW.obj.Z_2[2] += ((((control_3dof_DW.UnitDelay_DSTATE[2] +
    control_3dof_DW.obj.Z_3[2]) + 9.8) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.dv1[1]) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_b) *
    0.02;
  control_3dof_DW.obj.Z_3[2] += (control_3dof_B.dv1[1] *
    control_3dof_DW.obj.K_p[2] + control_3dof_DW.obj.K_v[2] *
    control_3dof_B.e_2_b) * 0.02;

  // Gain: '<S19>/KP' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Sum: '<S19>/Sum'
  //
  control_3dof_B.e_2_b = (control_3dof_U.Traj_sp.pos_sp[2] -
    control_3dof_U.Payload_Out.pL[2]) * CONTROL_PARAM.KP;

  // Saturate: '<S19>/vel_limit'
  if (control_3dof_B.e_2_b > 5.0) {
    control_3dof_B.e_2_b = 5.0;
  } else if (control_3dof_B.e_2_b < -5.0) {
    control_3dof_B.e_2_b = -5.0;
  }

  // Sum: '<S19>/Sum1' incorporates:
  //   Gain: '<S19>/KV'
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Saturate: '<S19>/vel_limit'
  //   Sum: '<S19>/Sum2'
  //   Sum: '<S19>/Sum3'
  //
  control_3dof_B.e_2_b = ((control_3dof_B.e_2_b + control_3dof_U.Traj_sp.vel_sp
    [2]) - control_3dof_U.Payload_Out.vL[2]) * CONTROL_PARAM.KV +
    control_3dof_U.Traj_sp.acc_ff[2];

  // Saturate: '<S19>/acc_limit'
  if (control_3dof_B.e_2_b > 3.0) {
    control_3dof_B.e_2_b = 3.0;
  } else if (control_3dof_B.e_2_b < -3.0) {
    control_3dof_B.e_2_b = -3.0;
  }

  // MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //
  control_3dof_B.e_1[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_1[2];
  control_3dof_B.UnitDelay_DSTATE_tmp = control_3dof_norm(control_3dof_B.e_1);
  control_3dof_B.e_1[0] /= control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.vu[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_1[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_2[0];
  control_3dof_B.e_1[1] = control_3dof_B.dv1[0] /
    control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.vu[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_1[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_2[1];

  // BusCreator generated from: '<S4>/Force Saturation && Disturbution' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_div_pd(_mm_sub_pd(_mm_set_pd
    (control_3dof_U.Payload_Out.vL[2], control_3dof_U.Payload_Out.pL[2]),
    _mm_set_pd(control_3dof_U.Payload_Out.v_1[2],
               control_3dof_U.Payload_Out.p_1[2])), _mm_set_pd
    (CONTROL_PARAM.CABLE_LEN, control_3dof_B.UnitDelay_DSTATE_tmp)));

  // MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Constant: '<S4>/-g'
  //   Inport: '<Root>/Payload_Out'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Saturate: '<S19>/acc_limit'
  //   Sum: '<S4>/Sum'
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_B.e_1[2] = control_3dof_B.dv1[0];
  control_3dof_B.q_2[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_2[2];
  _mm_storeu_pd(&control_3dof_B.e_2[0], _mm_sub_pd(_mm_mul_pd(_mm_set_pd
    (control_3dof_B.vu[0], control_3dof_B.e_1[1]), _mm_set_pd
    (control_3dof_B.dv1[0], control_3dof_B.dv1[1])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.e_1[0], control_3dof_B.vu[1]), _mm_set_pd
    (control_3dof_B.dv1[1], control_3dof_B.dv1[0]))));
  control_3dof_B.e_2[2] = control_3dof_B.e_1[0] * control_3dof_B.vu[1] -
    control_3dof_B.vu[0] * control_3dof_B.e_1[1];
  control_3dof_B.UnitDelay_DSTATE_tmp = control_3dof_norm(control_3dof_B.q_2);
  control_3dof_B.q_2[0] /= control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.vu[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_2[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_3[0];
  control_3dof_B.q_2[1] /= control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.vu[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_2[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_3[1];
  control_3dof_B.q_2[2] /= control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.vu[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_2[2]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_3[2];
  _mm_storeu_pd(&control_3dof_B.rtb_control_in2_w[0], _mm_sub_pd(_mm_mul_pd
    (_mm_set_pd(control_3dof_B.vu[0], control_3dof_B.q_2[1]), _mm_set_pd
     (control_3dof_B.q_2[2], control_3dof_B.vu[2])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_2[0], control_3dof_B.vu[1]), _mm_set_pd(control_3dof_B.vu
    [2], control_3dof_B.q_2[2]))));
  control_3dof_B.rtb_control_in2_w[2] = control_3dof_B.q_2[0] *
    control_3dof_B.vu[1] - control_3dof_B.vu[0] * control_3dof_B.q_2[1];
  control_3dof_B.UnitDelay_DSTATE_tmp = control_3dof_norm(control_3dof_B.q_3);
  control_3dof_B.q_3_k = control_3dof_B.q_3[0] /
    control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.q_3[0] = control_3dof_B.q_3_k;
  control_3dof_B.vu[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_3[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[0] = control_3dof_B.e_1[0];
  control_3dof_B.Q_new[3] = control_3dof_B.q_2[0];
  control_3dof_B.Q_new[6] = control_3dof_B.q_3_k;
  control_3dof_B.F_trim_c = (0.0 - control_3dof_DW.obj.Z_3[0]) *
    CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.F_trim[0] = control_3dof_B.F_trim_c;
  control_3dof_B.q_3_k = CONTROL_PARAM.MASS_LOAD * control_3dof_B.u0;
  control_3dof_B.F_way[0] = control_3dof_B.q_3_k;
  control_3dof_B.UnitDelay_DSTATE_tmp_g = control_3dof_B.F_trim_c +
    control_3dof_B.q_3_k;
  control_3dof_DW.UnitDelay_DSTATE[0] = control_3dof_B.UnitDelay_DSTATE_tmp_g;
  control_3dof_B.q_3_k = control_3dof_B.q_3[1] /
    control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.q_3[1] = control_3dof_B.q_3_k;
  control_3dof_B.vu[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_3[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[1] = control_3dof_B.e_1[1];
  control_3dof_B.Q_new[4] = control_3dof_B.q_2[1];
  control_3dof_B.Q_new[7] = control_3dof_B.q_3_k;
  control_3dof_B.F_trim_c = (0.0 - control_3dof_DW.obj.Z_3[1]) *
    CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.F_trim[1] = control_3dof_B.F_trim_c;
  control_3dof_B.q_3_k = CONTROL_PARAM.MASS_LOAD * control_3dof_B.u0_p;
  control_3dof_B.F_way[1] = control_3dof_B.q_3_k;
  control_3dof_B.UnitDelay_DSTATE_tmp_m = control_3dof_B.F_trim_c +
    control_3dof_B.q_3_k;
  control_3dof_DW.UnitDelay_DSTATE[1] = control_3dof_B.UnitDelay_DSTATE_tmp_m;
  control_3dof_B.q_3_k = control_3dof_B.q_3[2] /
    control_3dof_B.UnitDelay_DSTATE_tmp;
  control_3dof_B.q_3[2] = control_3dof_B.q_3_k;
  control_3dof_B.vu[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_3[2]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[2] = control_3dof_B.dv1[0];
  control_3dof_B.Q_new[5] = control_3dof_B.q_2[2];
  control_3dof_B.Q_new[8] = control_3dof_B.q_3_k;
  control_3dof_B.F_trim_c = (-9.8 - control_3dof_DW.obj.Z_3[2]) *
    CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.F_trim[2] = control_3dof_B.F_trim_c;
  control_3dof_B.UnitDelay_DSTATE_tmp = CONTROL_PARAM.MASS_LOAD *
    control_3dof_B.e_2_b + control_3dof_B.F_trim_c;
  control_3dof_DW.UnitDelay_DSTATE[2] = control_3dof_B.UnitDelay_DSTATE_tmp;
  _mm_storeu_pd(&control_3dof_B.rtb_control_in3_w[0], _mm_sub_pd(_mm_mul_pd
    (_mm_set_pd(control_3dof_B.vu[0], control_3dof_B.q_3[1]), _mm_set_pd
     (control_3dof_B.q_3_k, control_3dof_B.vu[2])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_3[0], control_3dof_B.vu[1]), _mm_set_pd(control_3dof_B.vu
    [2], control_3dof_B.q_3_k))));
  control_3dof_B.rtb_control_in3_w[2] = control_3dof_B.q_3[0] *
    control_3dof_B.vu[1] - control_3dof_B.vu[0] * control_3dof_B.q_3[1];
  control_3dof_pinv(control_3dof_B.Q_new, control_3dof_B.rtb_vu_tmp);
  control_3dof_B.UnitDelay_DSTATE = control_3dof_DW.UnitDelay_DSTATE[1];
  control_3dof_B.UnitDelay_DSTATE_c = control_3dof_DW.UnitDelay_DSTATE[0];
  control_3dof_B.UnitDelay_DSTATE_f = control_3dof_DW.UnitDelay_DSTATE[2];
  for (control_3dof_B.i = 0; control_3dof_B.i < 3; control_3dof_B.i++) {
    control_3dof_B.rtb_vu_tmp_g = (control_3dof_B.rtb_vu_tmp[control_3dof_B.i +
      3] * control_3dof_B.UnitDelay_DSTATE +
      control_3dof_B.rtb_vu_tmp[control_3dof_B.i] *
      control_3dof_B.UnitDelay_DSTATE_c) +
      control_3dof_B.rtb_vu_tmp[control_3dof_B.i + 6] *
      control_3dof_B.UnitDelay_DSTATE_f;
    control_3dof_B.vu[control_3dof_B.i] = control_3dof_B.rtb_vu_tmp_g;
    control_3dof_B.f_vertical[control_3dof_B.i] = -control_3dof_B.rtb_vu_tmp_g;
    control_3dof_B.x[control_3dof_B.i] = (-control_3dof_B.rtb_vu_tmp_g <
      CONTROL_PARAM.TENSION_MIN - 1.0E-5);
  }

  control_3dof_B.FixPtRelationalOperator = false;
  control_3dof_B.i = 0;
  exitg1 = false;
  while ((!exitg1) && (control_3dof_B.i < 3)) {
    if (control_3dof_B.x[control_3dof_B.i]) {
      control_3dof_B.FixPtRelationalOperator = true;
      exitg1 = true;
    } else {
      control_3dof_B.i++;
    }
  }

  if (!control_3dof_B.FixPtRelationalOperator) {
    control_3dof_B.x[0] = (control_3dof_B.f_vertical[0] >
      CONTROL_PARAM.TENSION_MAX + 1.0E-5);
    control_3dof_B.x[1] = (control_3dof_B.f_vertical[1] >
      CONTROL_PARAM.TENSION_MAX + 1.0E-5);
    control_3dof_B.x[2] = (control_3dof_B.f_vertical[2] >
      CONTROL_PARAM.TENSION_MAX + 1.0E-5);
    control_3dof_B.FixPtRelationalOperator = false;
    control_3dof_B.i = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.i < 3)) {
      if (control_3dof_B.x[control_3dof_B.i]) {
        control_3dof_B.FixPtRelationalOperator = true;
        exitg1 = true;
      } else {
        control_3dof_B.i++;
      }
    }

    control_3dof_B.FixPtRelationalOperator =
      !control_3dof_B.FixPtRelationalOperator;
  } else {
    control_3dof_B.FixPtRelationalOperator = false;
  }

  if (!control_3dof_B.FixPtRelationalOperator) {
    control_3dof_B.UnitDelay_DSTATE = control_3dof_B.F_trim[1];
    control_3dof_B.UnitDelay_DSTATE_c = control_3dof_B.F_trim[0];
    for (control_3dof_B.i = 0; control_3dof_B.i < 3; control_3dof_B.i++) {
      control_3dof_B.rtb_vu_tmp_g = (control_3dof_B.rtb_vu_tmp[control_3dof_B.i
        + 3] * control_3dof_B.UnitDelay_DSTATE +
        control_3dof_B.rtb_vu_tmp[control_3dof_B.i] *
        control_3dof_B.UnitDelay_DSTATE_c) +
        control_3dof_B.rtb_vu_tmp[control_3dof_B.i + 6] *
        control_3dof_B.F_trim_c;
      control_3dof_B.vu[control_3dof_B.i] = control_3dof_B.rtb_vu_tmp_g;
      control_3dof_B.x[control_3dof_B.i] = (-control_3dof_B.rtb_vu_tmp_g <
        CONTROL_PARAM.TENSION_MIN - 1.0E-5);
    }

    control_3dof_B.FixPtRelationalOperator = false;
    control_3dof_B.i = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.i < 3)) {
      if (control_3dof_B.x[control_3dof_B.i]) {
        control_3dof_B.FixPtRelationalOperator = true;
        exitg1 = true;
      } else {
        control_3dof_B.i++;
      }
    }

    if (!control_3dof_B.FixPtRelationalOperator) {
      control_3dof_B.x[0] = (-control_3dof_B.vu[0] > CONTROL_PARAM.TENSION_MAX +
        1.0E-5);
      control_3dof_B.x[1] = (-control_3dof_B.vu[1] > CONTROL_PARAM.TENSION_MAX +
        1.0E-5);
      control_3dof_B.x[2] = (-control_3dof_B.vu[2] > CONTROL_PARAM.TENSION_MAX +
        1.0E-5);
      control_3dof_B.FixPtRelationalOperator = false;
      control_3dof_B.i = 0;
      exitg1 = false;
      while ((!exitg1) && (control_3dof_B.i < 3)) {
        if (control_3dof_B.x[control_3dof_B.i]) {
          control_3dof_B.FixPtRelationalOperator = true;
          exitg1 = true;
        } else {
          control_3dof_B.i++;
        }
      }

      control_3dof_B.FixPtRelationalOperator =
        !control_3dof_B.FixPtRelationalOperator;
    } else {
      control_3dof_B.FixPtRelationalOperator = false;
    }

    if (control_3dof_B.FixPtRelationalOperator) {
      tmp = _mm_add_pd(_mm_loadu_pd(&control_3dof_B.F_trim[0]), _mm_loadu_pd
                       (&control_3dof_B.F_way[0]));
      _mm_storeu_pd(&control_3dof_B.f_vertical[0], tmp);
      control_3dof_B.f_vertical[2] = control_3dof_B.UnitDelay_DSTATE_tmp;
      control_3dof_TCISolver(control_3dof_B.f_vertical, control_3dof_B.F_trim,
        control_3dof_B.Q_new, CONTROL_PARAM.TENSION_MIN,
        CONTROL_PARAM.TENSION_MAX, control_3dof_B.vu,
        control_3dof_DW.UnitDelay_DSTATE);
    } else {
      control_3dof_B.f_vertical[0] = 0.0;
      control_3dof_B.f_vertical[1] = 0.0;
      control_3dof_B.f_vertical[2] = control_3dof_B.F_trim_c;
      control_3dof_TCISolver(control_3dof_B.F_trim, control_3dof_B.f_vertical,
        control_3dof_B.Q_new, CONTROL_PARAM.TENSION_MIN,
        CONTROL_PARAM.TENSION_MAX, control_3dof_B.vu,
        control_3dof_DW.UnitDelay_DSTATE);
    }
  }

  control_3dof_B.F_trim[0] = control_3dof_B.u0 -
    (control_3dof_B.UnitDelay_DSTATE_tmp_g - control_3dof_DW.UnitDelay_DSTATE[0])
    / CONTROL_PARAM.MASS_LOAD;
  control_3dof_DW.UnitDelay_DSTATE[0] /= CONTROL_PARAM.MASS_LOAD;
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_2[0], control_3dof_B.vu[0]), _mm_set_pd(control_3dof_B.vu
    [1], control_3dof_B.e_1[0])));
  control_3dof_B.acc_limit[0] = control_3dof_B.dv1[0];
  control_3dof_B.F_way[0] = control_3dof_B.dv1[1];
  control_3dof_B.F_trim[1] = control_3dof_B.u0_p -
    (control_3dof_B.UnitDelay_DSTATE_tmp_m - control_3dof_DW.UnitDelay_DSTATE[1])
    / CONTROL_PARAM.MASS_LOAD;
  control_3dof_DW.UnitDelay_DSTATE[1] /= CONTROL_PARAM.MASS_LOAD;
  tmp = _mm_mul_pd(_mm_loadu_pd(&control_3dof_B.vu[0]), _mm_set_pd
                   (control_3dof_B.q_2[1], control_3dof_B.e_1[1]));
  _mm_storeu_pd(&control_3dof_B.dv1[0], tmp);

  // MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
  //   Saturate: '<S19>/acc_limit'
  //   UnitDelay: '<S4>/Unit Delay'

  control_3dof_B.acc_limit[1] = control_3dof_B.dv1[0];
  control_3dof_B.F_way[1] = control_3dof_B.dv1[1];
  control_3dof_B.F_trim[2] = (control_3dof_B.e_2_b - 9.8) -
    (control_3dof_B.UnitDelay_DSTATE_tmp - control_3dof_DW.UnitDelay_DSTATE[2]) /
    CONTROL_PARAM.MASS_LOAD;
  control_3dof_DW.UnitDelay_DSTATE[2] /= CONTROL_PARAM.MASS_LOAD;
  tmp = _mm_mul_pd(_mm_loadu_pd(&control_3dof_B.vu[0]), _mm_set_pd
                   (control_3dof_B.q_2[2], control_3dof_B.e_1[2]));
  _mm_storeu_pd(&control_3dof_B.dv1[0], tmp);

  // MATLAB Function: '<S4>/Force Saturation && Disturbution'
  control_3dof_B.acc_limit[2] = control_3dof_B.dv1[0];
  control_3dof_B.F_way[2] = control_3dof_B.dv1[1];
  control_3dof_B.u0 = control_3dof_B.vu[2];
  tmp = _mm_mul_pd(_mm_set1_pd(control_3dof_B.u0), _mm_loadu_pd
                   (&control_3dof_B.q_3[0]));
  _mm_storeu_pd(&control_3dof_B.vu[0], tmp);
  control_3dof_B.vu[2] = control_3dof_B.u0 * control_3dof_B.q_3_k;

  // MATLAB Function: '<S7>/Vertical Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'
  //   Switch: '<S20>/Switch'

  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[0],
    control_3dof_B.e_1, control_3dof_B.e_2, control_3dof_B.F_trim,
    control_3dof_B.f_vertical, control_3dof_B.dot_err_g, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl);

  // MATLAB Function: '<S11>/Vertical Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'
  //   Switch: '<S20>/Switch'

  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[3],
    control_3dof_B.q_2, control_3dof_B.rtb_control_in2_w, control_3dof_B.F_trim,
    control_3dof_Y.force_sp1, control_3dof_B.dot_err_c, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_m);

  // MATLAB Function: '<S15>/Vertical Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'
  //   Switch: '<S20>/Switch'

  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[6],
    control_3dof_B.q_3, control_3dof_B.rtb_control_in3_w, control_3dof_B.F_trim,
    control_3dof_Y.force_sp2, control_3dof_B.dot_err, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_h);

  // MATLAB Function: '<S14>/Parallel Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_ParallelControl(control_3dof_B.q_3,
    control_3dof_B.rtb_control_in3_w, control_3dof_B.vu, control_3dof_B.F_trim,
    control_3dof_Y.force_sp3, CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl_p);

  // Sum: '<S15>/Sum' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S3>/Sum1'

  tmp = _mm_add_pd(_mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp2[0]),
    _mm_loadu_pd(&control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0])),
                   _mm_loadu_pd(&control_3dof_Y.force_sp3[0]));

  // DiscreteIntegrator: '<S11>/Discrete-Time Integrator' incorporates:
  //   Outport: '<Root>/force_sp3'
  //   Sum: '<S3>/Sum1'

  _mm_storeu_pd(&control_3dof_Y.force_sp3[0], tmp);

  // Outport: '<Root>/force_sp3' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S15>/Sum'
  //   Sum: '<S3>/Sum1'

  control_3dof_Y.force_sp3[2] += control_3dof_Y.force_sp2[2] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2];

  // MATLAB Function: '<S10>/Parallel Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_ParallelControl(control_3dof_B.q_2,
    control_3dof_B.rtb_control_in2_w, control_3dof_B.F_way,
    control_3dof_B.F_trim, control_3dof_Y.force_sp2, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_k);

  // Sum: '<S11>/Sum' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S2>/Sum1'

  tmp = _mm_add_pd(_mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp1[0]),
    _mm_loadu_pd(&control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0])),
                   _mm_loadu_pd(&control_3dof_Y.force_sp2[0]));

  // DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
  //   Outport: '<Root>/force_sp2'
  //   Sum: '<S2>/Sum1'

  _mm_storeu_pd(&control_3dof_Y.force_sp2[0], tmp);

  // Outport: '<Root>/force_sp2' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S11>/Sum'
  //   Sum: '<S2>/Sum1'

  control_3dof_Y.force_sp2[2] += control_3dof_Y.force_sp1[2] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2];

  // MATLAB Function: '<S6>/Parallel Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_ParallelControl(control_3dof_B.e_1, control_3dof_B.e_2,
    control_3dof_B.acc_limit, control_3dof_B.F_trim, control_3dof_Y.force_sp1,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl);

  // BusCreator generated from: '<Root>/state' incorporates:
  //   DiscreteIntegrator: '<S22>/Discrete-Time Integrator'
  //   Outport: '<Root>/state'

  control_3dof_Y.state.timestamp =
    control_3dof_DW.DiscreteTimeIntegrator_DSTAT_ps;

  // Update for UnitDelay: '<S21>/Delay Input1' incorporates:
  //   Inport: '<Root>/Dir_sp'
  //
  //  Block description for '<S21>/Delay Input1':
  //
  //   Store in Global RAM

  control_3dof_DW.DelayInput1_DSTATE = control_3dof_U.Dir_sp.timestamp;

  // Update for Delay: '<S20>/Delay'
  control_3dof_DW.icLoad = false;

  // Outport: '<Root>/force_sp1' incorporates:
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S1>/Sum1'
  //   Sum: '<S7>/Sum'

  control_3dof_Y.force_sp1[0] += control_3dof_B.f_vertical[0] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[0];

  // BusCreator generated from: '<Root>/state' incorporates:
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Outport: '<Root>/state'
  //
  control_3dof_Y.state.dL[0] = control_3dof_DW.obj.Z_3[0];

  // Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  tmp = _mm_set1_pd(0.02);

  // Gain: '<S15>/Gain'
  tmp_0 = _mm_set1_pd(CONTROL_PARAM.KQI);

  // Gain: '<S11>/Gain' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  //   Gain: '<S15>/Gain'

  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_0,
    _mm_set_pd(control_3dof_B.dot_err_c[0], control_3dof_B.dot_err[0])), tmp),
    _mm_set_pd(control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0],
               control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0])));

  // Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0] = control_3dof_B.dv1[0];

  // Update for DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0] = control_3dof_B.dv1[1];

  // Update for DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
  //   Gain: '<S7>/Gain'
  //   Sum: '<S1>/Sum1'

  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[0] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[0] * 0.02;

  // Outport: '<Root>/force_sp1' incorporates:
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S1>/Sum1'
  //   Sum: '<S7>/Sum'

  control_3dof_Y.force_sp1[1] += control_3dof_B.f_vertical[1] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[1];

  // BusCreator generated from: '<Root>/state' incorporates:
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Outport: '<Root>/state'
  //
  control_3dof_Y.state.dL[1] = control_3dof_DW.obj.Z_3[1];

  // Gain: '<S11>/Gain' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  //   Gain: '<S15>/Gain'

  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_0,
    _mm_set_pd(control_3dof_B.dot_err_c[1], control_3dof_B.dot_err[1])), tmp),
    _mm_set_pd(control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[1],
               control_3dof_DW.DiscreteTimeIntegrator_DSTATE[1])));

  // Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[1] = control_3dof_B.dv1[0];

  // Update for DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[1] = control_3dof_B.dv1[1];

  // Update for DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
  //   Gain: '<S7>/Gain'
  //   Sum: '<S1>/Sum1'

  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[1] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[1] * 0.02;

  // Outport: '<Root>/force_sp1' incorporates:
  //   DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
  //   Sum: '<S1>/Sum1'
  //   Sum: '<S7>/Sum'

  control_3dof_Y.force_sp1[2] += control_3dof_B.f_vertical[2] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[2];

  // BusCreator generated from: '<Root>/state' incorporates:
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Outport: '<Root>/state'
  //
  control_3dof_Y.state.dL[2] = control_3dof_DW.obj.Z_3[2];

  // Gain: '<S11>/Gain' incorporates:
  //   DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  //   DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  //   Gain: '<S15>/Gain'

  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_0,
    _mm_set_pd(control_3dof_B.dot_err_c[2], control_3dof_B.dot_err[2])), tmp),
    _mm_set_pd(control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2],
               control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2])));

  // Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2] = control_3dof_B.dv1[0];

  // Update for DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2] = control_3dof_B.dv1[1];

  // Update for DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
  //   Gain: '<S7>/Gain'
  //   Sum: '<S1>/Sum1'

  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[2] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[2] * 0.02;
}

// Model initialize function
void control_3dof::initialize()
{
  // InitializeConditions for Delay: '<S20>/Delay'
  control_3dof_DW.icLoad = true;

  // Start for MATLABSystem: '<S4>/Position 2nd ESO'
  //  Constructor
  //  Support name-value pair arguments when constructing object
  control_3dof_DW.obj.isInitialized = 1;

  //         %% Common functions
  //  Perform one-time calculations, such as computing constants
  control_3dof_DW.obj.K_p[0] = CONTROL_PARAM.ESO_PL[0];
  control_3dof_DW.obj.K_v[0] = CONTROL_PARAM.ESO_VL[0];
  control_3dof_DW.obj.Z_1[0] = 0.0;
  control_3dof_DW.obj.Z_2[0] = 0.0;
  control_3dof_DW.obj.Z_3[0] = 0.0;
  control_3dof_DW.obj.K_p[1] = CONTROL_PARAM.ESO_PL[1];
  control_3dof_DW.obj.K_v[1] = CONTROL_PARAM.ESO_VL[1];
  control_3dof_DW.obj.Z_1[1] = 0.0;
  control_3dof_DW.obj.Z_2[1] = 0.0;
  control_3dof_DW.obj.Z_3[1] = 0.0;
  control_3dof_DW.obj.K_p[2] = CONTROL_PARAM.ESO_PL[2];
  control_3dof_DW.obj.K_v[2] = CONTROL_PARAM.ESO_VL[2];
  control_3dof_DW.obj.Z_1[2] = 0.0;
  control_3dof_DW.obj.Z_2[2] = 0.0;
  control_3dof_DW.obj.Z_3[2] = 0.0;

  // InitializeConditions for MATLABSystem: '<S4>/Position 2nd ESO'
  //  Initialize / reset internal or discrete properties
  control_3dof_DW.obj.is_init = false;
}

// Model terminate function
void control_3dof::terminate()
{
  // (no terminate code required)
}

const char_T* control_3dof::RT_MODEL_control_3dof_T::getErrorStatus() const
{
  return (errorStatus);
}

void control_3dof::RT_MODEL_control_3dof_T::setErrorStatus(const char_T* const
  volatile aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
control_3dof::control_3dof() :
  control_3dof_U(),
  control_3dof_Y(),
  control_3dof_B(),
  control_3dof_DW(),
  control_3dof_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
control_3dof::~control_3dof() = default;

// Real-Time Model get method
control_3dof::RT_MODEL_control_3dof_T * control_3dof::getRTM()
{
  return (&control_3dof_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
