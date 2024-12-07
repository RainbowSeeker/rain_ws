//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: control_3dof.cpp
//
// Code generated for Simulink model 'control_3dof'.
//
// Model version                  : 1.442
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Wed Dec  4 20:07:06 2024
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
#include "cmath"
#include "limits"

// Exported block parameters
struct_KOZZ2f41NjtyX9PsCFQIaF CONTROL_PARAM{
  1.0,
  1.5,
  1.0,
  0.2,
  1.0,
  1.0,
  4.0,

  { 0.5, 0.0, -0.1, -0.5, 0.4, -0.1, -0.5, -0.4, -0.1 },
  1.0,
  2.0
} ;                                    // Variable: CONTROL_PARAM
                                          //  Referenced by:
                                          //    '<S4>/Force Disturbution'
                                          //    '<S5>/Kp'
                                          //    '<S6>/Parallel Control'
                                          //    '<S7>/Vertical Control'
                                          //    '<S10>/Parallel Control'
                                          //    '<S11>/Vertical Control'
                                          //    '<S14>/Parallel Control'
                                          //    '<S15>/Vertical Control'
                                          //    '<S19>/Kv'


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
  rty_y[3], real_T rtp_li, real_T rtp_mi, B_ParallelControl_control_3do_T
  *localB)
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
    _mm_storeu_pd(&rty_y[localB->i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd
      (_mm_mul_pd(tmp, _mm_set1_pd(localB->absxk)), _mm_mul_pd(tmp_0,
      _mm_set1_pd(localB->t))), _mm_mul_pd(tmp_1, _mm_set1_pd(localB->a))),
      _mm_mul_pd(_mm_set1_pd(localB->scale), tmp_2)), _mm_add_pd(_mm_add_pd
      (_mm_mul_pd(tmp_3, _mm_set1_pd(localB->rtu_ai_sp)), _mm_mul_pd(tmp_4,
      _mm_set1_pd(localB->rtu_ai_sp_m))), _mm_mul_pd(tmp_5, _mm_set1_pd
      (localB->rtu_ai_sp_c)))));
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    rty_y[localB->i] = (((localB->rtu_q[localB->i + 3] * localB->absxk +
                          localB->rtu_q[localB->i] * localB->t) + localB->
                         rtu_q[localB->i + 6] * localB->a) + localB->scale *
                        rtu_q[localB->i]) + ((localB->rtp_mi[localB->i + 3] *
      localB->rtu_ai_sp + localB->rtp_mi[localB->i] * localB->rtu_ai_sp_m) +
      localB->rtp_mi[localB->i + 6] * localB->rtu_ai_sp_c);
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
  rty_y[3], real_T rty_state[12], real_T rtp_Kq, real_T rtp_Kw, real_T rtp_li,
  real_T rtp_mi, B_VerticalControl_control_3do_T *localB)
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
    _mm_storeu_pd(&localB->w_sp[localB->i], _mm_mul_pd(_mm_set1_pd(-rtp_Kq),
      tmp_5));
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    localB->w_sp[localB->i] = -rtp_Kq * localB->e_q[localB->i];
  }

  for (localB->i_p = 0; localB->i_p < 3; localB->i_p++) {
    for (localB->i = 0; localB->i <= 0; localB->i += 2) {
      tmp_5 = _mm_loadu_pd(&localB->Sqi[localB->i + 3]);
      tmp_3 = _mm_loadu_pd(&localB->Sqi[localB->i]);
      tmp_4 = _mm_loadu_pd(&localB->Sqi[localB->i + 6]);
      _mm_storeu_pd(&localB->Sqi_m[localB->i + 3 * localB->i_p], _mm_add_pd
                    (_mm_add_pd(_mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_p + 1]), tmp_5), _mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_p]), tmp_3)), _mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_p + 2]), tmp_4)));
    }

    for (localB->i = 2; localB->i < 3; localB->i++) {
      localB->Sqi_m[localB->i + 3 * localB->i_p] = (localB->Sqi[3 * localB->i_p
        + 1] * localB->Sqi[localB->i + 3] + localB->Sqi[3 * localB->i_p] *
        localB->Sqi[localB->i]) + localB->Sqi[3 * localB->i_p + 2] * localB->
        Sqi[localB->i + 6];
    }
  }

  localB->a = rtp_mi * rtp_li;
  for (localB->i = 0; localB->i < 3; localB->i++) {
    localB->e_w_tmp = 0.0;
    for (localB->i_p = 0; localB->i_p < 3; localB->i_p++) {
      localB->e_w_tmp_tmp = 3 * localB->i_p + localB->i;
      localB->e_w_tmp += localB->Sqi_m[localB->e_w_tmp_tmp] * localB->
        w_sp[localB->i_p];
      localB->rtp_mi[localB->e_w_tmp_tmp] = (localB->Sqi[localB->i + 3] * rtp_mi
        * localB->Sqi[3 * localB->i_p + 1] + rtp_mi * localB->Sqi[localB->i] *
        localB->Sqi[3 * localB->i_p]) + localB->Sqi[localB->i + 6] * rtp_mi *
        localB->Sqi[3 * localB->i_p + 2];
    }

    localB->e_w_c = rtu_w[localB->i] + localB->e_w_tmp;
    localB->e_w[localB->i] = localB->e_w_c;
    localB->rtp_Kw[localB->i] = -rtp_Kw * localB->e_w_c - localB->e_w_tmp;
  }

  localB->e_w_tmp = localB->rtp_Kw[1];
  localB->e_w_c = localB->rtp_Kw[0];
  localB->rtp_Kw_k = localB->rtp_Kw[2];
  localB->rtu_ai_sp = rtu_ai_sp[1];
  localB->rtu_ai_sp_c = rtu_ai_sp[0];
  localB->rtu_ai_sp_b = rtu_ai_sp[2];
  for (localB->i = 0; localB->i <= 0; localB->i += 2) {
    __m128d tmp;
    __m128d tmp_0;
    __m128d tmp_1;
    __m128d tmp_2;
    tmp_5 = _mm_loadu_pd(&localB->Sqi[localB->i + 3]);
    tmp_3 = _mm_set1_pd(localB->a);
    tmp_4 = _mm_loadu_pd(&localB->Sqi[localB->i]);
    tmp = _mm_loadu_pd(&localB->Sqi[localB->i + 6]);
    tmp_0 = _mm_loadu_pd(&localB->rtp_mi[localB->i + 3]);
    tmp_1 = _mm_loadu_pd(&localB->rtp_mi[localB->i]);
    tmp_2 = _mm_loadu_pd(&localB->rtp_mi[localB->i + 6]);
    _mm_storeu_pd(&rty_y[localB->i], _mm_sub_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (_mm_mul_pd(tmp_5, tmp_3), _mm_set1_pd(localB->e_w_tmp)), _mm_mul_pd
      (_mm_mul_pd(tmp_3, tmp_4), _mm_set1_pd(localB->e_w_c))), _mm_mul_pd
      (_mm_mul_pd(tmp, tmp_3), _mm_set1_pd(localB->rtp_Kw_k))), _mm_add_pd
      (_mm_add_pd(_mm_mul_pd(tmp_0, _mm_set1_pd(localB->rtu_ai_sp)), _mm_mul_pd
                  (tmp_1, _mm_set1_pd(localB->rtu_ai_sp_c))), _mm_mul_pd(tmp_2,
      _mm_set1_pd(localB->rtu_ai_sp_b)))));
    tmp_5 = _mm_loadu_pd(&rtu_q_sp[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i], tmp_5);
    tmp_5 = _mm_loadu_pd(&localB->e_q[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i + 3], tmp_5);
    tmp_5 = _mm_loadu_pd(&localB->w_sp[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i + 6], tmp_5);
    tmp_5 = _mm_loadu_pd(&localB->e_w[localB->i]);
    _mm_storeu_pd(&rty_state[localB->i + 9], tmp_5);
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    rty_y[localB->i] = ((localB->Sqi[localB->i + 3] * localB->a *
                         localB->e_w_tmp + localB->a * localB->Sqi[localB->i] *
                         localB->e_w_c) + localB->Sqi[localB->i + 6] * localB->a
                        * localB->rtp_Kw_k) - ((localB->rtp_mi[localB->i + 3] *
      localB->rtu_ai_sp + localB->rtp_mi[localB->i] * localB->rtu_ai_sp_c) +
      localB->rtp_mi[localB->i + 6] * localB->rtu_ai_sp_b);
    rty_state[localB->i] = rtu_q_sp[localB->i];
    rty_state[localB->i + 3] = localB->e_q[localB->i];
    rty_state[localB->i + 6] = localB->w_sp[localB->i];
    rty_state[localB->i + 9] = localB->e_w[localB->i];
  }
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

// Function for MATLAB Function: '<S4>/Force Disturbution'
void control_3dof::control_3dof_xzlascl(real_T cfrom, real_T cto, int32_T m,
  int32_T n, real_T A[9], int32_T iA0, int32_T lda)
{
  control_3dof_B.cfromc_k = cfrom;
  control_3dof_B.ctoc_c = cto;
  control_3dof_B.notdone_b = true;
  while (control_3dof_B.notdone_b) {
    control_3dof_B.cfrom1_b = control_3dof_B.cfromc_k * 2.0041683600089728E-292;
    control_3dof_B.cto1_p = control_3dof_B.ctoc_c / 4.9896007738368E+291;
    if ((std::abs(control_3dof_B.cfrom1_b) > std::abs(control_3dof_B.ctoc_c)) &&
        (control_3dof_B.ctoc_c != 0.0)) {
      control_3dof_B.mul_c = 2.0041683600089728E-292;
      control_3dof_B.cfromc_k = control_3dof_B.cfrom1_b;
    } else if (std::abs(control_3dof_B.cto1_p) > std::abs
               (control_3dof_B.cfromc_k)) {
      control_3dof_B.mul_c = 4.9896007738368E+291;
      control_3dof_B.ctoc_c = control_3dof_B.cto1_p;
    } else {
      control_3dof_B.mul_c = control_3dof_B.ctoc_c / control_3dof_B.cfromc_k;
      control_3dof_B.notdone_b = false;
    }

    for (control_3dof_B.j_d = 0; control_3dof_B.j_d < n; control_3dof_B.j_d++) {
      control_3dof_B.offset_j = (control_3dof_B.j_d * lda + iA0) - 2;
      control_3dof_B.scalarLB_l = (m / 2) << 1;
      control_3dof_B.vectorUB_d = control_3dof_B.scalarLB_l - 2;
      for (control_3dof_B.b_i_g = 0; control_3dof_B.b_i_g <=
           control_3dof_B.vectorUB_d; control_3dof_B.b_i_g += 2) {
        __m128d tmp;
        control_3dof_B.i1 = (control_3dof_B.b_i_g + control_3dof_B.offset_j) + 1;
        tmp = _mm_loadu_pd(&A[control_3dof_B.i1]);
        _mm_storeu_pd(&A[control_3dof_B.i1], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.mul_c)));
      }

      for (control_3dof_B.b_i_g = control_3dof_B.scalarLB_l;
           control_3dof_B.b_i_g < m; control_3dof_B.b_i_g++) {
        control_3dof_B.i1 = (control_3dof_B.b_i_g + control_3dof_B.offset_j) + 1;
        A[control_3dof_B.i1] *= control_3dof_B.mul_c;
      }
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
real_T control_3dof::control_3dof_xnrm2(int32_T n, const real_T x[9], int32_T
  ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      control_3dof_B.scale_g = 3.3121686421112381E-170;
      control_3dof_B.kend_l = ix0 + n;
      for (control_3dof_B.k_h = ix0; control_3dof_B.k_h < control_3dof_B.kend_l;
           control_3dof_B.k_h++) {
        control_3dof_B.absxk_g = std::abs(x[control_3dof_B.k_h - 1]);
        if (control_3dof_B.absxk_g > control_3dof_B.scale_g) {
          control_3dof_B.t_m = control_3dof_B.scale_g / control_3dof_B.absxk_g;
          y = y * control_3dof_B.t_m * control_3dof_B.t_m + 1.0;
          control_3dof_B.scale_g = control_3dof_B.absxk_g;
        } else {
          control_3dof_B.t_m = control_3dof_B.absxk_g / control_3dof_B.scale_g;
          y += control_3dof_B.t_m * control_3dof_B.t_m;
        }
      }

      y = control_3dof_B.scale_g * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

// Function for MATLAB Function: '<S4>/Force Disturbution'
real_T control_3dof::control_3dof_xnrm2_n(int32_T n, const real_T x[3], int32_T
  ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      control_3dof_B.scale_f = 3.3121686421112381E-170;
      control_3dof_B.kend = ix0 + n;
      for (control_3dof_B.k_b = ix0; control_3dof_B.k_b < control_3dof_B.kend;
           control_3dof_B.k_b++) {
        control_3dof_B.absxk = std::abs(x[control_3dof_B.k_b - 1]);
        if (control_3dof_B.absxk > control_3dof_B.scale_f) {
          control_3dof_B.t = control_3dof_B.scale_f / control_3dof_B.absxk;
          y = y * control_3dof_B.t * control_3dof_B.t + 1.0;
          control_3dof_B.scale_f = control_3dof_B.absxk;
        } else {
          control_3dof_B.t = control_3dof_B.absxk / control_3dof_B.scale_f;
          y += control_3dof_B.t * control_3dof_B.t;
        }
      }

      y = control_3dof_B.scale_f * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
void control_3dof::control_3dof_xaxpy_f(int32_T n, real_T a, const real_T x[9],
  int32_T ix0, real_T y[3], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_b = (n / 2) << 1;
    control_3dof_B.vectorUB_n = control_3dof_B.scalarLB_b - 2;
    for (control_3dof_B.k_o = 0; control_3dof_B.k_o <= control_3dof_B.vectorUB_n;
         control_3dof_B.k_o += 2) {
      __m128d tmp;
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_o) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i3]);
      _mm_storeu_pd(&y[control_3dof_B.i3], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k_o) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k_o = control_3dof_B.scalarLB_b; control_3dof_B.k_o < n;
         control_3dof_B.k_o++) {
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_o) - 1;
      y[control_3dof_B.i3] += x[(ix0 + control_3dof_B.k_o) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
void control_3dof::control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3],
  int32_T ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_d = (n / 2) << 1;
    control_3dof_B.vectorUB_lx = control_3dof_B.scalarLB_d - 2;
    for (control_3dof_B.k = 0; control_3dof_B.k <= control_3dof_B.vectorUB_lx;
         control_3dof_B.k += 2) {
      __m128d tmp;
      control_3dof_B.i2 = (iy0 + control_3dof_B.k) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i2]);
      _mm_storeu_pd(&y[control_3dof_B.i2], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k = control_3dof_B.scalarLB_d; control_3dof_B.k < n;
         control_3dof_B.k++) {
      control_3dof_B.i2 = (iy0 + control_3dof_B.k) - 1;
      y[control_3dof_B.i2] += x[(ix0 + control_3dof_B.k) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

    for (control_3dof_B.j_p = 0; control_3dof_B.j_p < n; control_3dof_B.j_p++) {
      control_3dof_B.offset = (control_3dof_B.j_p * lda + iA0) - 2;
      control_3dof_B.scalarLB = (m / 2) << 1;
      control_3dof_B.vectorUB_l = control_3dof_B.scalarLB - 2;
      for (control_3dof_B.b_i = 0; control_3dof_B.b_i <=
           control_3dof_B.vectorUB_l; control_3dof_B.b_i += 2) {
        __m128d tmp;
        control_3dof_B.i = (control_3dof_B.b_i + control_3dof_B.offset) + 1;
        tmp = _mm_loadu_pd(&A[control_3dof_B.i]);
        _mm_storeu_pd(&A[control_3dof_B.i], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.mul)));
      }

      for (control_3dof_B.b_i = control_3dof_B.scalarLB; control_3dof_B.b_i < m;
           control_3dof_B.b_i++) {
        control_3dof_B.i = (control_3dof_B.b_i + control_3dof_B.offset) + 1;
        A[control_3dof_B.i] *= control_3dof_B.mul;
      }
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

// Function for MATLAB Function: '<S4>/Force Disturbution'
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

// Function for MATLAB Function: '<S4>/Force Disturbution'
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
      _mm_storeu_pd(&control_3dof_B.dv[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.m], control_3dof_B.rt), tmp));
      control_3dof_B.rt = control_3dof_B.dv[0];
      control_3dof_B.smm1 = control_3dof_B.dv[1];
      _mm_storeu_pd(&control_3dof_B.dv[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.qq],
         control_3dof_B.e[control_3dof_B.m]), tmp));
      control_3dof_B.smm1 = ((control_3dof_B.smm1 + control_3dof_B.rt) *
        (control_3dof_B.smm1 - control_3dof_B.rt) + control_3dof_B.dv[0] *
        control_3dof_B.dv[0]) / 2.0;
      control_3dof_B.c = control_3dof_B.rt * control_3dof_B.dv[0];
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

      control_3dof_B.rt = (control_3dof_B.dv[1] + control_3dof_B.rt) *
        (control_3dof_B.dv[1] - control_3dof_B.rt) + control_3dof_B.shift;
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq] / control_3dof_B.r *
        control_3dof_B.dv[1];
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

// Model step function
void control_3dof::step()
{
  __m128d tmp;
  __m128d tmp_0;
  boolean_T exitg1;

  // Gain: '<S5>/Kp' incorporates:
  //   Constant: '<S5>/pL_sp'
  //   Inport: '<Root>/Payload_Out'
  //   Sum: '<S5>/Sum'

  control_3dof_B.absx = (0.0 - control_3dof_U.Payload_Out.pL[0]) *
    CONTROL_PARAM.KP;

  // Saturate: '<S5>/Saturation'
  if (control_3dof_B.absx > 2.0) {
    control_3dof_B.absx = 2.0;
  } else if (control_3dof_B.absx < -2.0) {
    control_3dof_B.absx = -2.0;
  }

  // Sum: '<S19>/Sum' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Saturate: '<S5>/Saturation'

  control_3dof_B.absx -= control_3dof_U.Payload_Out.vL[0];
  control_3dof_Y.force_sp3[0] = control_3dof_B.absx;

  // Bias: '<S19>/g' incorporates:
  //   Gain: '<S19>/Kv'

  control_3dof_B.aL_sp[0] = CONTROL_PARAM.KV * control_3dof_B.absx;

  // MATLAB Function: '<S4>/Q sp'
  control_3dof_B.y[0] = 0.24999999999999994;
  control_3dof_B.y[1] = 0.43301270189221946;
  control_3dof_B.y[2] = 0.8660254037844386;

  // MATLAB Function: '<S4>/Force Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'

  control_3dof_B.A[0] = control_3dof_U.Payload_Out.q_1[0];
  control_3dof_B.A[3] = control_3dof_U.Payload_Out.q_2[0];
  control_3dof_B.A[6] = control_3dof_U.Payload_Out.q_3[0];

  // Gain: '<S5>/Kp' incorporates:
  //   Constant: '<S5>/pL_sp'
  //   Inport: '<Root>/Payload_Out'
  //   Sum: '<S5>/Sum'

  control_3dof_B.absx = (0.0 - control_3dof_U.Payload_Out.pL[1]) *
    CONTROL_PARAM.KP;

  // Saturate: '<S5>/Saturation'
  if (control_3dof_B.absx > 2.0) {
    control_3dof_B.absx = 2.0;
  } else if (control_3dof_B.absx < -2.0) {
    control_3dof_B.absx = -2.0;
  }

  // Sum: '<S19>/Sum' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Saturate: '<S5>/Saturation'

  control_3dof_B.absx -= control_3dof_U.Payload_Out.vL[1];
  control_3dof_Y.force_sp3[1] = control_3dof_B.absx;

  // Bias: '<S19>/g' incorporates:
  //   Gain: '<S19>/Kv'

  control_3dof_B.aL_sp[1] = CONTROL_PARAM.KV * control_3dof_B.absx;

  // MATLAB Function: '<S4>/Q sp'
  control_3dof_B.y[3] = -0.50000000000000011;
  control_3dof_B.y[4] = -0.0;
  control_3dof_B.y[5] = 0.8660254037844386;

  // MATLAB Function: '<S4>/Force Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'

  control_3dof_B.A[1] = control_3dof_U.Payload_Out.q_1[1];
  control_3dof_B.A[4] = control_3dof_U.Payload_Out.q_2[1];
  control_3dof_B.A[7] = control_3dof_U.Payload_Out.q_3[1];

  // Gain: '<S5>/Kp' incorporates:
  //   Constant: '<S5>/pL_sp'
  //   Inport: '<Root>/Payload_Out'
  //   Sum: '<S5>/Sum'

  control_3dof_B.absx = (-10.0 - control_3dof_U.Payload_Out.pL[2]) *
    CONTROL_PARAM.KP;

  // Saturate: '<S5>/Saturation'
  if (control_3dof_B.absx > 2.0) {
    control_3dof_B.absx = 2.0;
  } else if (control_3dof_B.absx < -2.0) {
    control_3dof_B.absx = -2.0;
  }

  // Sum: '<S19>/Sum' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Saturate: '<S5>/Saturation'

  control_3dof_B.absx -= control_3dof_U.Payload_Out.vL[2];
  control_3dof_Y.force_sp3[2] = control_3dof_B.absx;

  // Bias: '<S19>/g' incorporates:
  //   Gain: '<S19>/Kv'

  control_3dof_B.aL_sp[2] = CONTROL_PARAM.KV * control_3dof_B.absx - 9.8;

  // MATLAB Function: '<S4>/Q sp'
  control_3dof_B.y[6] = 0.24999999999999994;
  control_3dof_B.y[7] = -0.43301270189221946;
  control_3dof_B.y[8] = 0.8660254037844386;

  // MATLAB Function: '<S4>/Force Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'

  control_3dof_B.A[2] = control_3dof_U.Payload_Out.q_1[2];
  control_3dof_B.A[5] = control_3dof_U.Payload_Out.q_2[2];
  control_3dof_B.A[8] = control_3dof_U.Payload_Out.q_3[2];
  control_3dof_B.p = true;
  for (control_3dof_B.r_n = 0; control_3dof_B.r_n < 9; control_3dof_B.r_n++) {
    control_3dof_B.a[control_3dof_B.r_n] = 0.0;
    if (control_3dof_B.p) {
      control_3dof_B.absx = control_3dof_B.A[control_3dof_B.r_n];
      if (std::isinf(control_3dof_B.absx) || std::isnan(control_3dof_B.absx)) {
        control_3dof_B.p = false;
      }
    }
  }

  if (!control_3dof_B.p) {
    for (control_3dof_B.r_n = 0; control_3dof_B.r_n < 9; control_3dof_B.r_n++) {
      control_3dof_B.a[control_3dof_B.r_n] = (rtNaN);
    }
  } else {
    control_3dof_svd(control_3dof_B.A, control_3dof_B.U, control_3dof_B.vu,
                     control_3dof_B.V);
    control_3dof_B.absx = std::abs(control_3dof_B.vu[0]);
    if (std::isinf(control_3dof_B.absx) || std::isnan(control_3dof_B.absx)) {
      control_3dof_B.absx = (rtNaN);
    } else if (control_3dof_B.absx < 4.4501477170144028E-308) {
      control_3dof_B.absx = 4.94065645841247E-324;
    } else {
      std::frexp(control_3dof_B.absx, &control_3dof_B.exponent);
      control_3dof_B.absx = std::ldexp(1.0, control_3dof_B.exponent - 53);
    }

    control_3dof_B.absx *= 3.0;
    control_3dof_B.r_n = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.r_n < 3)) {
      if (std::isinf(control_3dof_B.vu[control_3dof_B.r_n]) || std::isnan
          (control_3dof_B.vu[control_3dof_B.r_n])) {
        control_3dof_B.absx = 1.7976931348623157E+308;
        exitg1 = true;
      } else {
        control_3dof_B.r_n++;
      }
    }

    control_3dof_B.r_n = -1;
    control_3dof_B.exponent = 0;
    while ((control_3dof_B.exponent < 3) &&
           (control_3dof_B.vu[control_3dof_B.exponent] > control_3dof_B.absx)) {
      control_3dof_B.r_n++;
      control_3dof_B.exponent++;
    }

    if (control_3dof_B.r_n + 1 > 0) {
      control_3dof_B.vcol = 1;
      for (control_3dof_B.j = 0; control_3dof_B.j <= control_3dof_B.r_n;
           control_3dof_B.j++) {
        control_3dof_B.absx = 1.0 / control_3dof_B.vu[control_3dof_B.j];
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
          control_3dof_B.a[control_3dof_B.vectorUB - 1] = 0.0;
        }
      }

      for (control_3dof_B.j = 0; control_3dof_B.j <= 6; control_3dof_B.j += 3) {
        control_3dof_B.ar = -1;
        control_3dof_B.vcol++;
        control_3dof_B.b = 3 * control_3dof_B.r_n + control_3dof_B.vcol;
        for (control_3dof_B.ib = control_3dof_B.vcol; control_3dof_B.ib <=
             control_3dof_B.b; control_3dof_B.ib += 3) {
          control_3dof_B.exponent = control_3dof_B.j + 3;
          control_3dof_B.vectorUB = control_3dof_B.j + 1;
          for (control_3dof_B.b_ic = control_3dof_B.j + 1; control_3dof_B.b_ic <=
               control_3dof_B.vectorUB; control_3dof_B.b_ic += 2) {
            tmp_0 = _mm_loadu_pd(&control_3dof_B.V[(control_3dof_B.ar +
              control_3dof_B.b_ic) - control_3dof_B.j]);
            tmp = _mm_loadu_pd(&control_3dof_B.a[control_3dof_B.b_ic - 1]);
            _mm_storeu_pd(&control_3dof_B.a[control_3dof_B.b_ic - 1], _mm_add_pd
                          (_mm_mul_pd(tmp_0, _mm_set1_pd
              (control_3dof_B.U[control_3dof_B.ib - 1])), tmp));
          }

          for (control_3dof_B.b_ic = control_3dof_B.exponent;
               control_3dof_B.b_ic <= control_3dof_B.j + 3; control_3dof_B.b_ic
               ++) {
            control_3dof_B.a[control_3dof_B.b_ic - 1] += control_3dof_B.V
              [(control_3dof_B.ar + control_3dof_B.b_ic) - control_3dof_B.j] *
              control_3dof_B.U[control_3dof_B.ib - 1];
          }

          control_3dof_B.ar += 3;
        }
      }
    }
  }

  control_3dof_B.absx = control_3dof_B.aL_sp[0];
  control_3dof_B.rtb_aL_sp_m = control_3dof_B.aL_sp[1];
  control_3dof_B.rtb_aL_sp_c = control_3dof_B.aL_sp[2];
  for (control_3dof_B.r_n = 0; control_3dof_B.r_n < 3; control_3dof_B.r_n++) {
    control_3dof_B.vu[control_3dof_B.r_n] = (CONTROL_PARAM.MASS_LOAD *
      control_3dof_B.absx * control_3dof_B.a[control_3dof_B.r_n] +
      control_3dof_B.a[control_3dof_B.r_n + 3] * (CONTROL_PARAM.MASS_LOAD *
      control_3dof_B.rtb_aL_sp_m)) + control_3dof_B.a[control_3dof_B.r_n + 6] *
      (CONTROL_PARAM.MASS_LOAD * control_3dof_B.rtb_aL_sp_c);
    control_3dof_B.rtb_control_in1_vforce[control_3dof_B.r_n] =
      control_3dof_B.vu[0] * control_3dof_U.Payload_Out.q_1[control_3dof_B.r_n];
  }

  _mm_storeu_pd(&control_3dof_B.rtb_control_in2_vforce[0], _mm_mul_pd(_mm_set_pd
    (control_3dof_B.vu[1], control_3dof_U.Payload_Out.q_2[0]), _mm_set_pd
    (control_3dof_U.Payload_Out.q_2[1], control_3dof_B.vu[1])));
  control_3dof_B.rtb_control_in2_vforce[2] = control_3dof_B.vu[1] *
    control_3dof_U.Payload_Out.q_2[2];
  control_3dof_B.absx = control_3dof_B.vu[2];
  _mm_storeu_pd(&control_3dof_B.vu[0], _mm_mul_pd(_mm_set1_pd
    (control_3dof_B.absx), _mm_loadu_pd(&control_3dof_U.Payload_Out.q_3[0])));
  control_3dof_B.vu[2] = control_3dof_B.absx * control_3dof_U.Payload_Out.q_3[2];

  // MATLAB Function: '<S7>/Vertical Control' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Disturbution'

  control_3dof_VerticalControl(&control_3dof_B.y[0],
    control_3dof_U.Payload_Out.q_1, control_3dof_U.Payload_Out.w_1,
    control_3dof_B.aL_sp, control_3dof_B.y_n, control_3dof_B.state_e,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl);

  // MATLAB Function: '<S11>/Vertical Control' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Disturbution'

  control_3dof_VerticalControl(&control_3dof_B.y[3],
    control_3dof_U.Payload_Out.q_2, control_3dof_U.Payload_Out.w_2,
    control_3dof_B.aL_sp, control_3dof_Y.force_sp1, control_3dof_B.state_e,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_o);

  // MATLAB Function: '<S15>/Vertical Control' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Disturbution'

  control_3dof_VerticalControl(&control_3dof_B.y[6],
    control_3dof_U.Payload_Out.q_3, control_3dof_U.Payload_Out.w_3,
    control_3dof_B.aL_sp, control_3dof_Y.force_sp3, control_3dof_B.state_e,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_i);

  // MATLAB Function: '<S14>/Parallel Control' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Disturbution'

  control_3dof_ParallelControl(control_3dof_U.Payload_Out.q_3,
    control_3dof_U.Payload_Out.w_3, control_3dof_B.vu, control_3dof_B.aL_sp,
    control_3dof_Y.force_sp2, CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl_k);

  // Sum: '<S3>/Sum'
  tmp_0 = _mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp2[0]), _mm_loadu_pd
                     (&control_3dof_Y.force_sp3[0]));

  // Outport: '<Root>/force_sp3' incorporates:
  //   Sum: '<S3>/Sum'

  _mm_storeu_pd(&control_3dof_Y.force_sp3[0], tmp_0);
  control_3dof_Y.force_sp3[2] += control_3dof_Y.force_sp2[2];

  // MATLAB Function: '<S10>/Parallel Control' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Disturbution'

  control_3dof_ParallelControl(control_3dof_U.Payload_Out.q_2,
    control_3dof_U.Payload_Out.w_2, control_3dof_B.rtb_control_in2_vforce,
    control_3dof_B.aL_sp, control_3dof_Y.force_sp2, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_g);

  // Sum: '<S2>/Sum'
  tmp_0 = _mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp2[0]), _mm_loadu_pd
                     (&control_3dof_Y.force_sp1[0]));

  // Outport: '<Root>/force_sp2' incorporates:
  //   Sum: '<S2>/Sum'

  _mm_storeu_pd(&control_3dof_Y.force_sp2[0], tmp_0);
  control_3dof_Y.force_sp2[2] += control_3dof_Y.force_sp1[2];

  // MATLAB Function: '<S6>/Parallel Control' incorporates:
  //   BusCreator generated from: '<S4>/Force Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Disturbution'

  control_3dof_ParallelControl(control_3dof_U.Payload_Out.q_1,
    control_3dof_U.Payload_Out.w_1, control_3dof_B.rtb_control_in1_vforce,
    control_3dof_B.aL_sp, control_3dof_Y.force_sp1, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl);

  // Sum: '<S1>/Sum'
  tmp_0 = _mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp1[0]), _mm_loadu_pd
                     (&control_3dof_B.y_n[0]));

  // Outport: '<Root>/force_sp1' incorporates:
  //   Sum: '<S1>/Sum'

  _mm_storeu_pd(&control_3dof_Y.force_sp1[0], tmp_0);
  control_3dof_Y.force_sp1[2] += control_3dof_B.y_n[2];
}

// Model initialize function
void control_3dof::initialize()
{
  // (no initialization code required)
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
