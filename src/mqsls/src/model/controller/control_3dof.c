/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: control_3dof.c
 *
 * Code generated for Simulink model 'control_3dof'.
 *
 * Model version                  : 1.704
 * Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
 * C/C++ source code generated on : Thu Jan  2 17:06:04 2025
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "control_3dof.h"
#include "rtwtypes.h"
#include <math.h>
#include <emmintrin.h>
#include "math.h"

/* Exported block parameters */
struct_VO89MlSMmM5pi10TnqAIZ CONTROL_PARAM = {
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
} ;                                    /* Variable: CONTROL_PARAM
                                        * Referenced by:
                                        *   '<S4>/Force Saturation && Disturbution'
                                        *   '<S4>/Position 2nd ESO'
                                        *   '<S6>/Parallel Control'
                                        *   '<S7>/Vertical Control'
                                        *   '<S7>/Gain'
                                        *   '<S10>/Parallel Control'
                                        *   '<S11>/Vertical Control'
                                        *   '<S11>/Gain'
                                        *   '<S14>/Parallel Control'
                                        *   '<S15>/Vertical Control'
                                        *   '<S15>/Gain'
                                        *   '<S19>/KP'
                                        *   '<S19>/KV'
                                        */

/* Block signals (default storage) */
B_control_3dof_T control_3dof_B;

/* Block states (default storage) */
DW_control_3dof_T control_3dof_DW;

/* External inputs (root inport signals with default storage) */
ExtU_control_3dof_T control_3dof_U;

/* External outputs (root outports fed by signals with default storage) */
ExtY_control_3dof_T control_3dof_Y;

/* Real-time model */
static RT_MODEL_control_3dof_T control_3dof_M_;
RT_MODEL_control_3dof_T *const control_3dof_M = &control_3dof_M_;
static void control_3dof_ParallelControl(const real_T rtu_q[3], const real_T
  rtu_w[3], const real_T rtu_vforce[3], const real_T rtu_ai_sp[3], real_T
  rty_f_parallel[3], real_T rtp_li, real_T rtp_mi,
  B_ParallelControl_control_3do_T *localB);
static void control_3dof_VerticalControl(const real_T rtu_q_sp[3], const real_T
  rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
  rty_f_vertical[3], real_T rty_dot_err[3], real_T rty_state[12], real_T rtp_Kq,
  real_T rtp_Kw, real_T rtp_li, real_T rtp_mi, B_VerticalControl_control_3do_T
  *localB);

/* Forward declaration for local functions */
static real_T control_3dof_norm(const real_T x[3]);
static real_T control_3dof_xzlangeM(const real_T x[9]);
static void control_3dof_xzlascl(real_T cfrom, real_T cto, int32_T m, int32_T n,
  real_T A[9], int32_T iA0, int32_T lda);
static real_T control_3dof_xnrm2(int32_T n, const real_T x[9], int32_T ix0);
static real_T control_3dof_xdotc(int32_T n, const real_T x[9], int32_T ix0,
  const real_T y[9], int32_T iy0);
static void control_3dof_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[9],
  int32_T iy0);
static real_T control_3dof_xnrm2_n(int32_T n, const real_T x[3], int32_T ix0);
static void control_3dof_xaxpy_f(int32_T n, real_T a, const real_T x[9], int32_T
  ix0, real_T y[3], int32_T iy0);
static void control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3],
  int32_T ix0, real_T y[9], int32_T iy0);
static void control_3dof_xzlascl_l(real_T cfrom, real_T cto, int32_T m, int32_T
  n, real_T A[3], int32_T iA0, int32_T lda);
static void control_3dof_xswap(real_T x[9], int32_T ix0, int32_T iy0);
static void control_3dof_xrotg(real_T *a, real_T *b, real_T *c, real_T *s);
static void control_3dof_xrot(real_T x[9], int32_T ix0, int32_T iy0, real_T c,
  real_T s);
static void control_3dof_svd(const real_T A[9], real_T U[9], real_T s[3], real_T
  V[9]);
static void control_3dof_pinv(const real_T A[9], real_T X[9]);
static real_T rtGetNaN(void);
static real32_T rtGetNaNF(void);
extern real_T rtInf;
extern real_T rtMinusInf;
extern real_T rtNaN;
extern real32_T rtInfF;
extern real32_T rtMinusInfF;
extern real32_T rtNaNF;
static boolean_T rtIsInf(real_T value);
static boolean_T rtIsInfF(real32_T value);
static boolean_T rtIsNaN(real_T value);
static boolean_T rtIsNaNF(real32_T value);
real_T rtNaN = -(real_T)NAN;
real_T rtInf = (real_T)INFINITY;
real_T rtMinusInf = -(real_T)INFINITY;
real32_T rtNaNF = -(real32_T)NAN;
real32_T rtInfF = (real32_T)INFINITY;
real32_T rtMinusInfF = -(real32_T)INFINITY;

/* Return rtNaN needed by the generated code. */
static real_T rtGetNaN(void)
{
  return rtNaN;
}

/* Return rtNaNF needed by the generated code. */
static real32_T rtGetNaNF(void)
{
  return rtNaNF;
}

/* Test if value is infinite */
static boolean_T rtIsInf(real_T value)
{
  return (boolean_T)isinf(value);
}

/* Test if single-precision value is infinite */
static boolean_T rtIsInfF(real32_T value)
{
  return (boolean_T)isinf(value);
}

/* Test if value is not a number */
static boolean_T rtIsNaN(real_T value)
{
  return (boolean_T)(isnan(value) != 0);
}

/* Test if single-precision value is not a number */
static boolean_T rtIsNaNF(real32_T value)
{
  return (boolean_T)(isnan(value) != 0);
}

/*
 * Output and update for atomic system:
 *    '<S6>/Parallel Control'
 *    '<S10>/Parallel Control'
 *    '<S14>/Parallel Control'
 */
static void control_3dof_ParallelControl(const real_T rtu_q[3], const real_T
  rtu_w[3], const real_T rtu_vforce[3], const real_T rtu_ai_sp[3], real_T
  rty_f_parallel[3], real_T rtp_li, real_T rtp_mi,
  B_ParallelControl_control_3do_T *localB)
{
  /* MATLAB Function 'Cable Controller/Parallel Control/Parallel Control': '<S8>:1' */
  /* '<S8>:1:3' f_parallel = q * q' * vforce + mi * li * (norm(w) ^ 2) * q + mi * q * q' * ai_sp; */
  localB->scale = 3.3121686421112381E-170;
  localB->absxk = fabs(rtu_w[0]);
  if (localB->absxk > 3.3121686421112381E-170) {
    localB->a = 1.0;
    localB->scale = localB->absxk;
  } else {
    localB->t = localB->absxk / 3.3121686421112381E-170;
    localB->a = localB->t * localB->t;
  }

  localB->absxk = fabs(rtu_w[1]);
  if (localB->absxk > localB->scale) {
    localB->t = localB->scale / localB->absxk;
    localB->a = localB->a * localB->t * localB->t + 1.0;
    localB->scale = localB->absxk;
  } else {
    localB->t = localB->absxk / localB->scale;
    localB->a += localB->t * localB->t;
  }

  localB->absxk = fabs(rtu_w[2]);
  if (localB->absxk > localB->scale) {
    localB->t = localB->scale / localB->absxk;
    localB->a = localB->a * localB->t * localB->t + 1.0;
    localB->scale = localB->absxk;
  } else {
    localB->t = localB->absxk / localB->scale;
    localB->a += localB->t * localB->t;
  }

  localB->scale *= sqrt(localB->a);
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

/*
 * Output and update for atomic system:
 *    '<S7>/Vertical Control'
 *    '<S11>/Vertical Control'
 *    '<S15>/Vertical Control'
 */
static void control_3dof_VerticalControl(const real_T rtu_q_sp[3], const real_T
  rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
  rty_f_vertical[3], real_T rty_dot_err[3], real_T rty_state[12], real_T rtp_Kq,
  real_T rtp_Kw, real_T rtp_li, real_T rtp_mi, B_VerticalControl_control_3do_T
  *localB)
{
  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;

  /* MATLAB Function 'Cable Controller/Vertical Control/Vertical Control': '<S9>:1' */
  /* '<S9>:1:3' Sqi = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0]; */
  localB->Sqi[0] = 0.0;
  localB->Sqi[3] = -rtu_q[2];
  localB->Sqi[6] = rtu_q[1];
  localB->Sqi[1] = rtu_q[2];
  localB->Sqi[4] = 0.0;
  localB->Sqi[7] = -rtu_q[0];
  localB->Sqi[2] = -rtu_q[1];
  localB->Sqi[5] = rtu_q[0];
  localB->Sqi[8] = 0.0;

  /* '<S9>:1:4' e_q = cross(q_sp, q); */
  _mm_storeu_pd(&localB->e_q[0], _mm_sub_pd(_mm_mul_pd(_mm_set_pd(rtu_q[0],
    rtu_q_sp[1]), _mm_set_pd(rtu_q_sp[2], rtu_q[2])), _mm_mul_pd(_mm_set_pd
    (rtu_q_sp[0], rtu_q[1]), _mm_set_pd(rtu_q[2], rtu_q_sp[2]))));
  localB->e_q[2] = rtu_q_sp[0] * rtu_q[1] - rtu_q[0] * rtu_q_sp[1];

  /* '<S9>:1:7' w_sp = -Kq * e_q; */
  /* '<S9>:1:10' e_w = w + Sqi * Sqi * w_sp; */
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

  /* '<S9>:1:11' dot_q_sp = cross(w_sp, q_sp); */
  /* '<S9>:1:12' f_vertical = mi * li * Sqi * (-Kq * e_q -Kw * e_w - Sqi * Sqi * w_sp - (q' * w_sp) * dot_q_sp) - mi * Sqi * Sqi * ai_sp; */
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

  /* '<S9>:1:13' dot_err = -1 / mi / li * Sqi * (e_q + 2 * e_w); */
  localB->rtu_q = -1.0 / rtp_mi / rtp_li;
  tmp_5 = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(2.0), _mm_loadu_pd(&localB->e_w[0])),
                     _mm_loadu_pd(&localB->e_q[0]));
  _mm_storeu_pd(&localB->dv[0], tmp_5);
  localB->e_w_tmp_tmp = localB->dv[0];
  localB->e_q_idx_1 = localB->dv[1];
  localB->e_q_idx_2 = 2.0 * localB->e_w[2] + localB->e_q[2];

  /* '<S9>:1:15' state = [q_sp; e_q; w_sp; e_w]; */
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

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_norm(const real_T x[3])
{
  real_T absxk;
  real_T scale;
  real_T t;
  real_T y;
  scale = 3.3121686421112381E-170;
  absxk = fabs(x[0]);
  if (absxk > 3.3121686421112381E-170) {
    y = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    y = t * t;
  }

  absxk = fabs(x[1]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  absxk = fabs(x[2]);
  if (absxk > scale) {
    t = scale / absxk;
    y = y * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    y += t * t;
  }

  return scale * sqrt(y);
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_xzlangeM(const real_T x[9])
{
  real_T y;
  int32_T k;
  boolean_T exitg1;
  y = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 9)) {
    real_T absxk;
    absxk = fabs(x[k]);
    if (rtIsNaN(absxk)) {
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

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xzlascl(real_T cfrom, real_T cto, int32_T m, int32_T n,
  real_T A[9], int32_T iA0, int32_T lda)
{
  control_3dof_B.cfromc_g = cfrom;
  control_3dof_B.ctoc_m = cto;
  control_3dof_B.notdone_n = true;
  while (control_3dof_B.notdone_n) {
    control_3dof_B.cfrom1_n = control_3dof_B.cfromc_g * 2.0041683600089728E-292;
    control_3dof_B.cto1_p = control_3dof_B.ctoc_m / 4.9896007738368E+291;
    if ((fabs(control_3dof_B.cfrom1_n) > fabs(control_3dof_B.ctoc_m)) &&
        (control_3dof_B.ctoc_m != 0.0)) {
      control_3dof_B.mul_l = 2.0041683600089728E-292;
      control_3dof_B.cfromc_g = control_3dof_B.cfrom1_n;
    } else if (fabs(control_3dof_B.cto1_p) > fabs(control_3dof_B.cfromc_g)) {
      control_3dof_B.mul_l = 4.9896007738368E+291;
      control_3dof_B.ctoc_m = control_3dof_B.cto1_p;
    } else {
      control_3dof_B.mul_l = control_3dof_B.ctoc_m / control_3dof_B.cfromc_g;
      control_3dof_B.notdone_n = false;
    }

    for (control_3dof_B.j_n = 0; control_3dof_B.j_n < n; control_3dof_B.j_n++) {
      control_3dof_B.offset_b = (control_3dof_B.j_n * lda + iA0) - 2;
      control_3dof_B.scalarLB_l = (m / 2) << 1;
      control_3dof_B.vectorUB_h = control_3dof_B.scalarLB_l - 2;
      for (control_3dof_B.b_i_b = 0; control_3dof_B.b_i_b <=
           control_3dof_B.vectorUB_h; control_3dof_B.b_i_b += 2) {
        __m128d tmp;
        control_3dof_B.i2 = (control_3dof_B.b_i_b + control_3dof_B.offset_b) + 1;
        tmp = _mm_loadu_pd(&A[control_3dof_B.i2]);
        _mm_storeu_pd(&A[control_3dof_B.i2], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.mul_l)));
      }

      for (control_3dof_B.b_i_b = control_3dof_B.scalarLB_l;
           control_3dof_B.b_i_b < m; control_3dof_B.b_i_b++) {
        control_3dof_B.i2 = (control_3dof_B.b_i_b + control_3dof_B.offset_b) + 1;
        A[control_3dof_B.i2] *= control_3dof_B.mul_l;
      }
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_xnrm2(int32_T n, const real_T x[9], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      control_3dof_B.scale_d = 3.3121686421112381E-170;
      control_3dof_B.kend_a = ix0 + n;
      for (control_3dof_B.k_j = ix0; control_3dof_B.k_j < control_3dof_B.kend_a;
           control_3dof_B.k_j++) {
        control_3dof_B.absxk_g = fabs(x[control_3dof_B.k_j - 1]);
        if (control_3dof_B.absxk_g > control_3dof_B.scale_d) {
          control_3dof_B.t_l = control_3dof_B.scale_d / control_3dof_B.absxk_g;
          y = y * control_3dof_B.t_l * control_3dof_B.t_l + 1.0;
          control_3dof_B.scale_d = control_3dof_B.absxk_g;
        } else {
          control_3dof_B.t_l = control_3dof_B.absxk_g / control_3dof_B.scale_d;
          y += control_3dof_B.t_l * control_3dof_B.t_l;
        }
      }

      y = control_3dof_B.scale_d * sqrt(y);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_xdotc(int32_T n, const real_T x[9], int32_T ix0,
  const real_T y[9], int32_T iy0)
{
  real_T d;
  int32_T k;
  d = 0.0;
  if (n >= 1) {
    for (k = 0; k < n; k++) {
      d += x[(ix0 + k) - 1] * y[(iy0 + k) - 1];
    }
  }

  return d;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[9],
  int32_T iy0)
{
  int32_T k;
  if ((n >= 1) && (!(a == 0.0))) {
    for (k = 0; k < n; k++) {
      int32_T tmp;
      tmp = (iy0 + k) - 1;
      y[tmp] += y[(ix0 + k) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_xnrm2_n(int32_T n, const real_T x[3], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      control_3dof_B.scale_j = 3.3121686421112381E-170;
      control_3dof_B.kend = ix0 + n;
      for (control_3dof_B.k_f = ix0; control_3dof_B.k_f < control_3dof_B.kend;
           control_3dof_B.k_f++) {
        control_3dof_B.absxk = fabs(x[control_3dof_B.k_f - 1]);
        if (control_3dof_B.absxk > control_3dof_B.scale_j) {
          control_3dof_B.t = control_3dof_B.scale_j / control_3dof_B.absxk;
          y = y * control_3dof_B.t * control_3dof_B.t + 1.0;
          control_3dof_B.scale_j = control_3dof_B.absxk;
        } else {
          control_3dof_B.t = control_3dof_B.absxk / control_3dof_B.scale_j;
          y += control_3dof_B.t * control_3dof_B.t;
        }
      }

      y = control_3dof_B.scale_j * sqrt(y);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xaxpy_f(int32_T n, real_T a, const real_T x[9], int32_T
  ix0, real_T y[3], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_bj = (n / 2) << 1;
    control_3dof_B.vectorUB_j = control_3dof_B.scalarLB_bj - 2;
    for (control_3dof_B.k_e = 0; control_3dof_B.k_e <= control_3dof_B.vectorUB_j;
         control_3dof_B.k_e += 2) {
      __m128d tmp;
      control_3dof_B.i4 = (iy0 + control_3dof_B.k_e) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i4]);
      _mm_storeu_pd(&y[control_3dof_B.i4], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k_e) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k_e = control_3dof_B.scalarLB_bj; control_3dof_B.k_e < n;
         control_3dof_B.k_e++) {
      control_3dof_B.i4 = (iy0 + control_3dof_B.k_e) - 1;
      y[control_3dof_B.i4] += x[(ix0 + control_3dof_B.k_e) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3],
  int32_T ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_b = (n / 2) << 1;
    control_3dof_B.vectorUB_d = control_3dof_B.scalarLB_b - 2;
    for (control_3dof_B.k = 0; control_3dof_B.k <= control_3dof_B.vectorUB_d;
         control_3dof_B.k += 2) {
      __m128d tmp;
      control_3dof_B.i3 = (iy0 + control_3dof_B.k) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i3]);
      _mm_storeu_pd(&y[control_3dof_B.i3], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k = control_3dof_B.scalarLB_b; control_3dof_B.k < n;
         control_3dof_B.k++) {
      control_3dof_B.i3 = (iy0 + control_3dof_B.k) - 1;
      y[control_3dof_B.i3] += x[(ix0 + control_3dof_B.k) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xzlascl_l(real_T cfrom, real_T cto, int32_T m, int32_T
  n, real_T A[3], int32_T iA0, int32_T lda)
{
  control_3dof_B.cfromc = cfrom;
  control_3dof_B.ctoc = cto;
  control_3dof_B.notdone = true;
  while (control_3dof_B.notdone) {
    control_3dof_B.cfrom1 = control_3dof_B.cfromc * 2.0041683600089728E-292;
    control_3dof_B.cto1 = control_3dof_B.ctoc / 4.9896007738368E+291;
    if ((fabs(control_3dof_B.cfrom1) > fabs(control_3dof_B.ctoc)) &&
        (control_3dof_B.ctoc != 0.0)) {
      control_3dof_B.mul = 2.0041683600089728E-292;
      control_3dof_B.cfromc = control_3dof_B.cfrom1;
    } else if (fabs(control_3dof_B.cto1) > fabs(control_3dof_B.cfromc)) {
      control_3dof_B.mul = 4.9896007738368E+291;
      control_3dof_B.ctoc = control_3dof_B.cto1;
    } else {
      control_3dof_B.mul = control_3dof_B.ctoc / control_3dof_B.cfromc;
      control_3dof_B.notdone = false;
    }

    for (control_3dof_B.j_l = 0; control_3dof_B.j_l < n; control_3dof_B.j_l++) {
      control_3dof_B.offset = (control_3dof_B.j_l * lda + iA0) - 2;
      control_3dof_B.scalarLB = (m / 2) << 1;
      control_3dof_B.vectorUB_o = control_3dof_B.scalarLB - 2;
      for (control_3dof_B.b_i = 0; control_3dof_B.b_i <=
           control_3dof_B.vectorUB_o; control_3dof_B.b_i += 2) {
        __m128d tmp;
        control_3dof_B.i1 = (control_3dof_B.b_i + control_3dof_B.offset) + 1;
        tmp = _mm_loadu_pd(&A[control_3dof_B.i1]);
        _mm_storeu_pd(&A[control_3dof_B.i1], _mm_mul_pd(tmp, _mm_set1_pd
          (control_3dof_B.mul)));
      }

      for (control_3dof_B.b_i = control_3dof_B.scalarLB; control_3dof_B.b_i < m;
           control_3dof_B.b_i++) {
        control_3dof_B.i1 = (control_3dof_B.b_i + control_3dof_B.offset) + 1;
        A[control_3dof_B.i1] *= control_3dof_B.mul;
      }
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xswap(real_T x[9], int32_T ix0, int32_T iy0)
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

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xrotg(real_T *a, real_T *b, real_T *c, real_T *s)
{
  control_3dof_B.roe = *b;
  control_3dof_B.absa = fabs(*a);
  control_3dof_B.absb = fabs(*b);
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
    control_3dof_B.scale *= sqrt(control_3dof_B.ads * control_3dof_B.ads +
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

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xrot(real_T x[9], int32_T ix0, int32_T iy0, real_T c,
  real_T s)
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

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_svd(const real_T A[9], real_T U[9], real_T s[3], real_T
  V[9])
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
      if (fabs(control_3dof_B.nrm) >= 1.0020841800044864E-292) {
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
        if (fabs(control_3dof_B.e[0]) >= 1.0020841800044864E-292) {
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
      control_3dof_B.rt = fabs(control_3dof_B.r);
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
        control_3dof_B.rt = fabs(control_3dof_B.r);
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

    control_3dof_B.nrm = fmax(control_3dof_B.nrm, fmax(fabs
      (control_3dof_B.b_s[control_3dof_B.qq]), fabs
      (control_3dof_B.e[control_3dof_B.qq])));
  }

  while ((control_3dof_B.m + 2 > 0) && (control_3dof_B.qp1 < 75)) {
    control_3dof_B.qq = control_3dof_B.m + 1;
    exitg1 = false;
    while (!(exitg1 || (control_3dof_B.qq == 0))) {
      control_3dof_B.rt = fabs(control_3dof_B.e[control_3dof_B.qq - 1]);
      if (control_3dof_B.rt <= (fabs(control_3dof_B.b_s[control_3dof_B.qq - 1])
           + fabs(control_3dof_B.b_s[control_3dof_B.qq])) *
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
            control_3dof_B.rt = fabs(control_3dof_B.e[control_3dof_B.kase - 1]);
          }

          if (control_3dof_B.kase > control_3dof_B.qq + 1) {
            control_3dof_B.rt += fabs(control_3dof_B.e[control_3dof_B.kase - 2]);
          }

          control_3dof_B.r = fabs(control_3dof_B.b_s[control_3dof_B.kase - 1]);
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
      control_3dof_B.r = fmax(fmax(fmax(fmax(fabs(control_3dof_B.rt), fabs
        (control_3dof_B.b_s[control_3dof_B.m])), fabs
        (control_3dof_B.e[control_3dof_B.m])), fabs
        (control_3dof_B.b_s[control_3dof_B.qq])), fabs
        (control_3dof_B.e[control_3dof_B.qq]));
      tmp = _mm_set1_pd(control_3dof_B.r);
      _mm_storeu_pd(&control_3dof_B.dv2[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.m], control_3dof_B.rt), tmp));
      control_3dof_B.rt = control_3dof_B.dv2[0];
      control_3dof_B.smm1 = control_3dof_B.dv2[1];
      _mm_storeu_pd(&control_3dof_B.dv2[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.qq],
         control_3dof_B.e[control_3dof_B.m]), tmp));
      control_3dof_B.smm1 = ((control_3dof_B.smm1 + control_3dof_B.rt) *
        (control_3dof_B.smm1 - control_3dof_B.rt) + control_3dof_B.dv2[0] *
        control_3dof_B.dv2[0]) / 2.0;
      control_3dof_B.c = control_3dof_B.rt * control_3dof_B.dv2[0];
      control_3dof_B.c *= control_3dof_B.c;
      if ((control_3dof_B.smm1 != 0.0) || (control_3dof_B.c != 0.0)) {
        control_3dof_B.shift = sqrt(control_3dof_B.smm1 * control_3dof_B.smm1 +
          control_3dof_B.c);
        if (control_3dof_B.smm1 < 0.0) {
          control_3dof_B.shift = -control_3dof_B.shift;
        }

        control_3dof_B.shift = control_3dof_B.c / (control_3dof_B.smm1 +
          control_3dof_B.shift);
      } else {
        control_3dof_B.shift = 0.0;
      }

      control_3dof_B.rt = (control_3dof_B.dv2[1] + control_3dof_B.rt) *
        (control_3dof_B.dv2[1] - control_3dof_B.rt) + control_3dof_B.shift;
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq] / control_3dof_B.r *
        control_3dof_B.dv2[1];
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

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_pinv(const real_T A[9], real_T X[9])
{
  __m128d tmp;
  __m128d tmp_0;
  boolean_T exitg1;
  control_3dof_B.p = true;
  for (control_3dof_B.r_d = 0; control_3dof_B.r_d < 9; control_3dof_B.r_d++) {
    X[control_3dof_B.r_d] = 0.0;
    if (control_3dof_B.p) {
      control_3dof_B.absx = A[control_3dof_B.r_d];
      if (rtIsInf(control_3dof_B.absx) || rtIsNaN(control_3dof_B.absx)) {
        control_3dof_B.p = false;
      }
    }
  }

  if (!control_3dof_B.p) {
    for (control_3dof_B.r_d = 0; control_3dof_B.r_d < 9; control_3dof_B.r_d++) {
      X[control_3dof_B.r_d] = (rtNaN);
    }
  } else {
    control_3dof_svd(A, control_3dof_B.U, control_3dof_B.s, control_3dof_B.V);
    control_3dof_B.absx = fabs(control_3dof_B.s[0]);
    if (rtIsInf(control_3dof_B.absx) || rtIsNaN(control_3dof_B.absx)) {
      control_3dof_B.absx = (rtNaN);
    } else if (control_3dof_B.absx < 4.4501477170144028E-308) {
      control_3dof_B.absx = 4.94065645841247E-324;
    } else {
      frexp(control_3dof_B.absx, &control_3dof_B.exponent);
      control_3dof_B.absx = ldexp(1.0, control_3dof_B.exponent - 53);
    }

    control_3dof_B.absx *= 3.0;
    control_3dof_B.r_d = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.r_d < 3)) {
      if (rtIsInf(control_3dof_B.s[control_3dof_B.r_d]) || rtIsNaN
          (control_3dof_B.s[control_3dof_B.r_d])) {
        control_3dof_B.absx = 1.7976931348623157E+308;
        exitg1 = true;
      } else {
        control_3dof_B.r_d++;
      }
    }

    control_3dof_B.r_d = -1;
    control_3dof_B.exponent = 0;
    while ((control_3dof_B.exponent < 3) &&
           (control_3dof_B.s[control_3dof_B.exponent] > control_3dof_B.absx)) {
      control_3dof_B.r_d++;
      control_3dof_B.exponent++;
    }

    if (control_3dof_B.r_d + 1 > 0) {
      control_3dof_B.vcol = 1;
      for (control_3dof_B.j = 0; control_3dof_B.j <= control_3dof_B.r_d;
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
        control_3dof_B.b = 3 * control_3dof_B.r_d + control_3dof_B.vcol;
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

/* Model step function */
void control_3dof_step(void)
{
  static const real_T b[3] = { 0.0, 0.0, -9.8 };

  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  static const int8_T a[18] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0,
    0, -1 };

  static const real_T b_0[3] = { 0.0, 0.0, -9.8 };

  boolean_T exitg1;

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[0] = control_3dof_U.Dir_sp.q_sp1[0];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[3] = control_3dof_U.Dir_sp.q_sp2[0];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[6] = control_3dof_U.Dir_sp.q_sp3[0];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[1] = control_3dof_U.Dir_sp.q_sp1[1];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[4] = control_3dof_U.Dir_sp.q_sp2[1];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[7] = control_3dof_U.Dir_sp.q_sp3[1];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[2] = control_3dof_U.Dir_sp.q_sp1[2];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[5] = control_3dof_U.Dir_sp.q_sp2[2];

  /* SignalConversion generated from: '<S20>/Vector Concatenate' incorporates:
   *  Inport: '<Root>/Dir_sp'
   */
  control_3dof_B.Q_new[8] = control_3dof_U.Dir_sp.q_sp3[2];

  /* RelationalOperator: '<S21>/FixPt Relational Operator' incorporates:
   *  Inport: '<Root>/Dir_sp'
   *  UnitDelay: '<S21>/Delay Input1'
   *
   * Block description for '<S21>/Delay Input1':
   *
   *  Store in Global RAM
   */
  control_3dof_B.FixPtRelationalOperator = (control_3dof_U.Dir_sp.timestamp !=
    control_3dof_DW.DelayInput1_DSTATE);
  for (control_3dof_B.i = 0; control_3dof_B.i < 9; control_3dof_B.i++) {
    /* Delay: '<S20>/Delay' incorporates:
     *  Concatenate: '<S20>/Vector Concatenate'
     *  Switch: '<S20>/Switch'
     */
    if (control_3dof_DW.icLoad) {
      control_3dof_DW.Delay_DSTATE[control_3dof_B.i] =
        control_3dof_B.Q_new[control_3dof_B.i];
    }

    /* End of Delay: '<S20>/Delay' */

    /* Switch: '<S20>/Switch' incorporates:
     *  Concatenate: '<S20>/Vector Concatenate'
     */
    if (control_3dof_B.FixPtRelationalOperator) {
      control_3dof_DW.Delay_DSTATE[control_3dof_B.i] =
        control_3dof_B.Q_new[control_3dof_B.i];
    }

    /* End of Switch: '<S20>/Switch' */
  }

  /* MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  Inport: '<Root>/Payload_Out'
   */
  /* 'pos_2nd_eso:13' (3, 1) */
  /* 'pos_2nd_eso:13' K_p */
  /* 'pos_2nd_eso:15' (3, 1) */
  /* 'pos_2nd_eso:15' K_v */
  control_3dof_DW.obj.K_p[0] = CONTROL_PARAM.ESO_PL[0];
  control_3dof_DW.obj.K_v[0] = CONTROL_PARAM.ESO_VL[0];
  control_3dof_DW.obj.K_p[1] = CONTROL_PARAM.ESO_PL[1];
  control_3dof_DW.obj.K_v[1] = CONTROL_PARAM.ESO_VL[1];
  control_3dof_DW.obj.K_p[2] = CONTROL_PARAM.ESO_PL[2];
  control_3dof_DW.obj.K_v[2] = CONTROL_PARAM.ESO_VL[2];

  /* 'pos_2nd_eso:56' if ~obj.is_init */
  if (!control_3dof_DW.obj.is_init) {
    /* 'pos_2nd_eso:57' obj.is_init = true; */
    /* 'pos_2nd_eso:34' logical */
    /* 'pos_2nd_eso:34' (1, 1) */
    /* 'pos_2nd_eso:34' is_init */
    control_3dof_DW.obj.is_init = true;

    /* 'pos_2nd_eso:58' obj.Z_1 = pos; */
    /* 'pos_2nd_eso:59' obj.Z_2 = vel; */
    /* 'pos_2nd_eso:60' obj.Z_3 = zeros(size(pos)); */
    control_3dof_DW.obj.Z_1[0] = control_3dof_U.Payload_Out.pL[0];
    control_3dof_DW.obj.Z_2[0] = control_3dof_U.Payload_Out.vL[0];
    control_3dof_DW.obj.Z_3[0] = 0.0;
    control_3dof_DW.obj.Z_1[1] = control_3dof_U.Payload_Out.pL[1];
    control_3dof_DW.obj.Z_2[1] = control_3dof_U.Payload_Out.vL[1];
    control_3dof_DW.obj.Z_3[1] = 0.0;
    control_3dof_DW.obj.Z_1[2] = control_3dof_U.Payload_Out.pL[2];
    control_3dof_DW.obj.Z_2[2] = control_3dof_U.Payload_Out.vL[2];
    control_3dof_DW.obj.Z_3[2] = 0.0;
  } else {
    /* 'pos_2nd_eso:61' else */
    /* 'pos_2nd_eso:62' obj.Z_1 = obj.Z_1; */
    /* 'pos_2nd_eso:63' obj.Z_2 = obj.Z_2; */
    /* 'pos_2nd_eso:64' obj.Z_3 = obj.Z_3; */
  }

  /* Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  Inport: '<Root>/Payload_Out'
   */
  /* 'pos_2nd_eso:66' sts = getSampleTime(obj); */
  /* 'pos_2nd_eso:67' dt = sts.SampleTime; */
  /* 'pos_2nd_eso:68' g = [0; 0; 9.8]; */
  /*  Extended state observer */
  /* 'pos_2nd_eso:71' e_1 = pos - obj.Z_1; */
  /* 'pos_2nd_eso:72' e_2 = vel - obj.Z_2; */
  /* 'pos_2nd_eso:73' dot_Z_1 = obj.Z_2 + obj.K_p(1) * e_1; */
  /* 'pos_2nd_eso:74' dot_Z_2 = acc + obj.Z_3 + g + ... */
  /* 'pos_2nd_eso:75'                       obj.K_p(2) * e_1 + obj.K_v(2) * e_2; */
  /* 'pos_2nd_eso:76' dot_Z_3 = obj.K_p(3) * e_1 + obj.K_v(3) * e_2; */
  /*  Update states */
  /* 'pos_2nd_eso:79' obj.Z_1 = obj.Z_1 + dt * dot_Z_1; */
  /* 'pos_2nd_eso:80' obj.Z_2 = obj.Z_2 + dt * dot_Z_2; */
  /* 'pos_2nd_eso:81' obj.Z_3 = obj.Z_3 + dt * dot_Z_3; */
  /* 'pos_2nd_eso:83' z1 = obj.Z_1; */
  /* 'pos_2nd_eso:84' z2 = obj.Z_2; */
  /* 'pos_2nd_eso:85' z3 = obj.Z_3; */
  /* MATLAB Function 'Payload Controller/Force Saturation && Disturbution': '<S18>:1' */
  /* '<S18>:1:3' q_1 = Payload_Out.pL - Payload_Out.p_1; */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_sub_pd(_mm_set_pd
    (control_3dof_U.Payload_Out.vL[0], control_3dof_U.Payload_Out.pL[0]),
    _mm_set_pd(control_3dof_DW.obj.Z_2[0], control_3dof_DW.obj.Z_1[0])));

  /* MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  UnitDelay: '<S4>/Unit Delay'
   * */
  control_3dof_DW.obj.Z_1[0] += (control_3dof_DW.obj.K_p[0] *
    control_3dof_B.dv1[0] + control_3dof_DW.obj.Z_2[0]) * 0.02;
  control_3dof_DW.obj.Z_2[0] += (((control_3dof_DW.UnitDelay_DSTATE[0] +
    control_3dof_DW.obj.Z_3[0]) + control_3dof_B.dv1[0] *
    control_3dof_DW.obj.K_p[1]) + control_3dof_DW.obj.K_v[1] *
    control_3dof_B.dv1[1]) * 0.02;
  control_3dof_DW.obj.Z_3[0] += (control_3dof_B.dv1[0] *
    control_3dof_DW.obj.K_p[2] + control_3dof_B.dv1[1] *
    control_3dof_DW.obj.K_v[2]) * 0.02;

  /* Gain: '<S19>/KP' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Sum: '<S19>/Sum'
   * */
  control_3dof_B.upper = (control_3dof_U.Traj_sp.pos_sp[0] -
    control_3dof_U.Payload_Out.pL[0]) * CONTROL_PARAM.KP;

  /* Saturate: '<S19>/vel_limit' */
  if (control_3dof_B.upper > 5.0) {
    control_3dof_B.upper = 5.0;
  } else if (control_3dof_B.upper < -5.0) {
    control_3dof_B.upper = -5.0;
  }

  /* Sum: '<S19>/Sum1' incorporates:
   *  Gain: '<S19>/KV'
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Saturate: '<S19>/vel_limit'
   *  Sum: '<S19>/Sum2'
   *  Sum: '<S19>/Sum3'
   * */
  control_3dof_B.upper = ((control_3dof_B.upper + control_3dof_U.Traj_sp.vel_sp
    [0]) - control_3dof_U.Payload_Out.vL[0]) * CONTROL_PARAM.KV +
    control_3dof_U.Traj_sp.acc_ff[0];

  /* Saturate: '<S19>/acc_limit' */
  if (control_3dof_B.upper > 3.0) {
    control_3dof_B.acc_limit[0] = 3.0;
  } else if (control_3dof_B.upper < -3.0) {
    control_3dof_B.acc_limit[0] = -3.0;
  } else {
    control_3dof_B.acc_limit[0] = control_3dof_B.upper;
  }

  /* Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   * */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_sub_pd(_mm_loadu_pd
    (&control_3dof_U.Payload_Out.pL[0]), _mm_set_pd(control_3dof_DW.obj.Z_1[1],
    control_3dof_U.Payload_Out.p_1[0])));

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' */
  control_3dof_B.e_1[0] = control_3dof_B.dv1[0];

  /* MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  UnitDelay: '<S4>/Unit Delay'
   * */
  control_3dof_B.e_2_c = control_3dof_U.Payload_Out.vL[1] -
    control_3dof_DW.obj.Z_2[1];
  control_3dof_DW.obj.Z_1[1] += (control_3dof_DW.obj.K_p[0] *
    control_3dof_B.dv1[1] + control_3dof_DW.obj.Z_2[1]) * 0.02;
  control_3dof_DW.obj.Z_2[1] += (((control_3dof_DW.UnitDelay_DSTATE[1] +
    control_3dof_DW.obj.Z_3[1]) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.dv1[1]) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_c) *
    0.02;
  control_3dof_DW.obj.Z_3[1] += (control_3dof_B.dv1[1] *
    control_3dof_DW.obj.K_p[2] + control_3dof_DW.obj.K_v[2] *
    control_3dof_B.e_2_c) * 0.02;

  /* Gain: '<S19>/KP' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Sum: '<S19>/Sum'
   * */
  control_3dof_B.lower = (control_3dof_U.Traj_sp.pos_sp[1] -
    control_3dof_U.Payload_Out.pL[1]) * CONTROL_PARAM.KP;

  /* Saturate: '<S19>/vel_limit' */
  if (control_3dof_B.lower > 5.0) {
    control_3dof_B.lower = 5.0;
  } else if (control_3dof_B.lower < -5.0) {
    control_3dof_B.lower = -5.0;
  }

  /* Sum: '<S19>/Sum1' incorporates:
   *  Gain: '<S19>/KV'
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Saturate: '<S19>/vel_limit'
   *  Sum: '<S19>/Sum2'
   *  Sum: '<S19>/Sum3'
   * */
  control_3dof_B.lower = ((control_3dof_B.lower + control_3dof_U.Traj_sp.vel_sp
    [1]) - control_3dof_U.Payload_Out.vL[1]) * CONTROL_PARAM.KV +
    control_3dof_U.Traj_sp.acc_ff[1];

  /* Saturate: '<S19>/acc_limit' */
  if (control_3dof_B.lower > 3.0) {
    control_3dof_B.acc_limit[1] = 3.0;
  } else if (control_3dof_B.lower < -3.0) {
    control_3dof_B.acc_limit[1] = -3.0;
  } else {
    control_3dof_B.acc_limit[1] = control_3dof_B.lower;
  }

  /* Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   * */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_sub_pd(_mm_loadu_pd
    (&control_3dof_U.Payload_Out.pL[1]), _mm_set_pd(control_3dof_DW.obj.Z_1[2],
    control_3dof_U.Payload_Out.p_1[1])));

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' */
  control_3dof_B.e_1[1] = control_3dof_B.dv1[0];

  /* MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  UnitDelay: '<S4>/Unit Delay'
   * */
  control_3dof_B.e_2_c = control_3dof_U.Payload_Out.vL[2] -
    control_3dof_DW.obj.Z_2[2];
  control_3dof_DW.obj.Z_1[2] += (control_3dof_DW.obj.K_p[0] *
    control_3dof_B.dv1[1] + control_3dof_DW.obj.Z_2[2]) * 0.02;
  control_3dof_DW.obj.Z_2[2] += ((((control_3dof_DW.UnitDelay_DSTATE[2] +
    control_3dof_DW.obj.Z_3[2]) + 9.8) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.dv1[1]) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_c) *
    0.02;
  control_3dof_DW.obj.Z_3[2] += (control_3dof_B.dv1[1] *
    control_3dof_DW.obj.K_p[2] + control_3dof_DW.obj.K_v[2] *
    control_3dof_B.e_2_c) * 0.02;

  /* Gain: '<S19>/KP' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Sum: '<S19>/Sum'
   * */
  control_3dof_B.u0 = (control_3dof_U.Traj_sp.pos_sp[2] -
                       control_3dof_U.Payload_Out.pL[2]) * CONTROL_PARAM.KP;

  /* Saturate: '<S19>/vel_limit' */
  if (control_3dof_B.u0 > 5.0) {
    control_3dof_B.u0 = 5.0;
  } else if (control_3dof_B.u0 < -5.0) {
    control_3dof_B.u0 = -5.0;
  }

  /* Sum: '<S19>/Sum1' incorporates:
   *  Gain: '<S19>/KV'
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Saturate: '<S19>/vel_limit'
   *  Sum: '<S19>/Sum2'
   *  Sum: '<S19>/Sum3'
   * */
  control_3dof_B.u0 = ((control_3dof_B.u0 + control_3dof_U.Traj_sp.vel_sp[2]) -
                       control_3dof_U.Payload_Out.vL[2]) * CONTROL_PARAM.KV +
    control_3dof_U.Traj_sp.acc_ff[2];

  /* Saturate: '<S19>/acc_limit' */
  if (control_3dof_B.u0 > 3.0) {
    control_3dof_B.acc_limit[2] = 3.0;
  } else if (control_3dof_B.u0 < -3.0) {
    control_3dof_B.acc_limit[2] = -3.0;
  } else {
    control_3dof_B.acc_limit[2] = control_3dof_B.u0;
  }

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  control_3dof_B.e_1[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_1[2];

  /* '<S18>:1:4' q_1 = q_1 / norm(q_1); */
  control_3dof_B.F_trim_g = control_3dof_norm(control_3dof_B.e_1);

  /* '<S18>:1:5' w_1 = cross(q_1, (Payload_Out.vL - Payload_Out.v_1) / li); */
  /* '<S18>:1:8' q_2 = Payload_Out.pL - Payload_Out.p_2; */
  control_3dof_B.e_1[0] /= control_3dof_B.F_trim_g;
  control_3dof_B.F_var[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_1[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_2[0];
  control_3dof_B.e_1[1] = control_3dof_B.dv1[0] / control_3dof_B.F_trim_g;
  control_3dof_B.F_var[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_1[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_2[1];

  /* BusCreator generated from: '<S4>/Force Saturation && Disturbution' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_div_pd(_mm_sub_pd(_mm_set_pd
    (control_3dof_U.Payload_Out.vL[2], control_3dof_U.Payload_Out.pL[2]),
    _mm_set_pd(control_3dof_U.Payload_Out.v_1[2],
               control_3dof_U.Payload_Out.p_1[2])), _mm_set_pd
    (CONTROL_PARAM.CABLE_LEN, control_3dof_B.F_trim_g)));

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Constant: '<S4>/-g'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Sum: '<S4>/Sum'
   * */
  control_3dof_B.e_1[2] = control_3dof_B.dv1[0];
  control_3dof_B.q_2[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_2[2];
  _mm_storeu_pd(&control_3dof_B.e_2[0], _mm_sub_pd(_mm_mul_pd(_mm_set_pd
    (control_3dof_B.F_var[0], control_3dof_B.e_1[1]), _mm_set_pd
    (control_3dof_B.dv1[0], control_3dof_B.dv1[1])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.e_1[0], control_3dof_B.F_var[1]), _mm_set_pd
    (control_3dof_B.dv1[1], control_3dof_B.dv1[0]))));
  control_3dof_B.e_2[2] = control_3dof_B.e_1[0] * control_3dof_B.F_var[1] -
    control_3dof_B.F_var[0] * control_3dof_B.e_1[1];

  /* '<S18>:1:9' q_2 = q_2 / norm(q_2); */
  control_3dof_B.F_trim_g = control_3dof_norm(control_3dof_B.q_2);

  /* '<S18>:1:10' w_2 = cross(q_2, (Payload_Out.vL - Payload_Out.v_2) / li); */
  /* '<S18>:1:13' q_3 = Payload_Out.pL - Payload_Out.p_3; */
  control_3dof_B.q_2[0] /= control_3dof_B.F_trim_g;
  control_3dof_B.F_var[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_2[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_3[0];
  control_3dof_B.q_2[1] /= control_3dof_B.F_trim_g;
  control_3dof_B.F_var[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_2[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_3[1];
  control_3dof_B.q_2[2] /= control_3dof_B.F_trim_g;
  control_3dof_B.F_var[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_2[2]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_3[2];
  _mm_storeu_pd(&control_3dof_B.rtb_control_in2_w[0], _mm_sub_pd(_mm_mul_pd
    (_mm_set_pd(control_3dof_B.F_var[0], control_3dof_B.q_2[1]), _mm_set_pd
     (control_3dof_B.q_2[2], control_3dof_B.F_var[2])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_2[0], control_3dof_B.F_var[1]), _mm_set_pd
    (control_3dof_B.F_var[2], control_3dof_B.q_2[2]))));
  control_3dof_B.rtb_control_in2_w[2] = control_3dof_B.q_2[0] *
    control_3dof_B.F_var[1] - control_3dof_B.F_var[0] * control_3dof_B.q_2[1];

  /* '<S18>:1:14' q_3 = q_3 / norm(q_3); */
  control_3dof_B.F_trim_g = control_3dof_norm(control_3dof_B.q_3);

  /* '<S18>:1:15' w_3 = cross(q_3, (Payload_Out.vL - Payload_Out.v_3) / li); */
  /* '<S18>:1:17' Q = [q_1 q_2 q_3]; */
  /* '<S18>:1:20' F_trim = mL * a_trim; */
  /* '<S18>:1:21' F_way  = mL * a_way; */
  control_3dof_B.e_2_c = control_3dof_B.q_3[0] / control_3dof_B.F_trim_g;
  control_3dof_B.q_3[0] = control_3dof_B.e_2_c;
  control_3dof_B.F_var[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_3[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[0] = control_3dof_B.e_1[0];
  control_3dof_B.Q_new[3] = control_3dof_B.q_2[0];
  control_3dof_B.Q_new[6] = control_3dof_B.e_2_c;
  control_3dof_B.F_trim[0] = (0.0 - control_3dof_DW.obj.Z_3[0]) *
    CONTROL_PARAM.MASS_LOAD;

  /* Saturate: '<S19>/acc_limit' */
  if (control_3dof_B.upper > 3.0) {
    control_3dof_B.upper = 3.0;
  } else if (control_3dof_B.upper < -3.0) {
    control_3dof_B.upper = -3.0;
  }

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Constant: '<S4>/-g'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Saturate: '<S19>/acc_limit'
   *  Sum: '<S4>/Sum'
   * */
  control_3dof_B.F_way[0] = CONTROL_PARAM.MASS_LOAD * control_3dof_B.upper;
  control_3dof_B.e_2_c = control_3dof_B.q_3[1] / control_3dof_B.F_trim_g;
  control_3dof_B.q_3[1] = control_3dof_B.e_2_c;
  control_3dof_B.F_var[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_3[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[1] = control_3dof_B.e_1[1];
  control_3dof_B.Q_new[4] = control_3dof_B.q_2[1];
  control_3dof_B.Q_new[7] = control_3dof_B.e_2_c;
  control_3dof_B.F_trim[1] = (0.0 - control_3dof_DW.obj.Z_3[1]) *
    CONTROL_PARAM.MASS_LOAD;

  /* Saturate: '<S19>/acc_limit' */
  if (control_3dof_B.lower > 3.0) {
    control_3dof_B.lower = 3.0;
  } else if (control_3dof_B.lower < -3.0) {
    control_3dof_B.lower = -3.0;
  }

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Constant: '<S4>/-g'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Saturate: '<S19>/acc_limit'
   *  Sum: '<S4>/Sum'
   * */
  control_3dof_B.F_way[1] = CONTROL_PARAM.MASS_LOAD * control_3dof_B.lower;
  control_3dof_B.e_2_c = control_3dof_B.q_3[2] / control_3dof_B.F_trim_g;
  control_3dof_B.q_3[2] = control_3dof_B.e_2_c;
  control_3dof_B.F_var[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_3[2]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[2] = control_3dof_B.dv1[0];
  control_3dof_B.Q_new[5] = control_3dof_B.q_2[2];
  control_3dof_B.Q_new[8] = control_3dof_B.e_2_c;
  control_3dof_B.F_trim[2] = (-9.8 - control_3dof_DW.obj.Z_3[2]) *
    CONTROL_PARAM.MASS_LOAD;

  /* Saturate: '<S19>/acc_limit' */
  if (control_3dof_B.u0 > 3.0) {
    control_3dof_B.u0 = 3.0;
  } else if (control_3dof_B.u0 < -3.0) {
    control_3dof_B.u0 = -3.0;
  }

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  Saturate: '<S19>/acc_limit'
   *  UnitDelay: '<S4>/Unit Delay'
   */
  control_3dof_B.F_way[2] = CONTROL_PARAM.MASS_LOAD * control_3dof_B.u0;
  _mm_storeu_pd(&control_3dof_B.rtb_control_in3_w[0], _mm_sub_pd(_mm_mul_pd
    (_mm_set_pd(control_3dof_B.F_var[0], control_3dof_B.q_3[1]), _mm_set_pd
     (control_3dof_B.e_2_c, control_3dof_B.F_var[2])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_3[0], control_3dof_B.F_var[1]), _mm_set_pd
    (control_3dof_B.F_var[2], control_3dof_B.e_2_c))));
  control_3dof_B.rtb_control_in3_w[2] = control_3dof_B.q_3[0] *
    control_3dof_B.F_var[1] - control_3dof_B.F_var[0] * control_3dof_B.q_3[1];

  /* '<S18>:1:22' T_min = min_tension * ones(3, 1); */
  /* '<S18>:1:23' T_max = max_tension * ones(3, 1); */
  /* '<S18>:1:24' [vu, F_actual] = PWAS(F_trim, F_way, Q, T_min, T_max); */
  /* 'PWAS:3' if ~isvector(T_max) || ~isvector(T_min) */
  /* 'PWAS:7' N = size(Q, 2); */
  /* 'PWAS:8' W = -Q; */
  /* 'PWAS:9' A_wrench = [eye(N); -eye(N)] * pinv(W); */
  for (control_3dof_B.i = 0; control_3dof_B.i <= 6; control_3dof_B.i += 2) {
    tmp_2 = _mm_loadu_pd(&control_3dof_B.Q_new[control_3dof_B.i]);
    _mm_storeu_pd(&control_3dof_B.rtb_Q_new_m[control_3dof_B.i], _mm_mul_pd
                  (tmp_2, _mm_set1_pd(-1.0)));
  }

  for (control_3dof_B.i = 8; control_3dof_B.i < 9; control_3dof_B.i++) {
    control_3dof_B.rtb_Q_new_m[control_3dof_B.i] =
      -control_3dof_B.Q_new[control_3dof_B.i];
  }

  control_3dof_pinv(control_3dof_B.rtb_Q_new_m, control_3dof_B.dv);
  for (control_3dof_B.i = 0; control_3dof_B.i < 6; control_3dof_B.i++) {
    control_3dof_B.a = a[control_3dof_B.i + 6];
    control_3dof_B.a_j = a[control_3dof_B.i];
    control_3dof_B.a_o = a[control_3dof_B.i + 12];
    for (control_3dof_B.i_d = 0; control_3dof_B.i_d < 3; control_3dof_B.i_d++) {
      control_3dof_B.A_wrench[control_3dof_B.i + 6 * control_3dof_B.i_d] =
        (control_3dof_B.dv[3 * control_3dof_B.i_d + 1] * (real_T)
         control_3dof_B.a + control_3dof_B.dv[3 * control_3dof_B.i_d] * (real_T)
         control_3dof_B.a_j) + control_3dof_B.dv[3 * control_3dof_B.i_d + 2] *
        (real_T)control_3dof_B.a_o;
    }
  }

  /* 'PWAS:10' b_wrench = [T_max; -T_min]; */
  control_3dof_B.b_wrench[0] = CONTROL_PARAM.TENSION_MAX;
  control_3dof_B.b_wrench[3] = -CONTROL_PARAM.TENSION_MIN;
  control_3dof_B.b_wrench[1] = CONTROL_PARAM.TENSION_MAX;
  control_3dof_B.b_wrench[4] = -CONTROL_PARAM.TENSION_MIN;
  control_3dof_B.b_wrench[2] = CONTROL_PARAM.TENSION_MAX;
  control_3dof_B.b_wrench[5] = -CONTROL_PARAM.TENSION_MIN;

  /* 'PWAS:11' if WrenchIsValid(F_trim, A_wrench, b_wrench) */
  /* 'PWAS:54' retval = all(A_w * F_sp <= b_w); */
  control_3dof_B.F_trim_g = control_3dof_B.F_trim[1];
  control_3dof_B.upper = control_3dof_B.F_trim[0];
  control_3dof_B.lower = control_3dof_B.F_trim[2];
  for (control_3dof_B.i = 0; control_3dof_B.i <= 4; control_3dof_B.i += 2) {
    tmp_2 = _mm_loadu_pd(&control_3dof_B.A_wrench[control_3dof_B.i + 6]);
    tmp_0 = _mm_loadu_pd(&control_3dof_B.A_wrench[control_3dof_B.i]);
    tmp_1 = _mm_loadu_pd(&control_3dof_B.A_wrench[control_3dof_B.i + 12]);
    _mm_storeu_pd(&control_3dof_B.A[control_3dof_B.i], _mm_add_pd(_mm_add_pd
      (_mm_mul_pd(tmp_2, _mm_set1_pd(control_3dof_B.F_trim_g)), _mm_mul_pd(tmp_0,
      _mm_set1_pd(control_3dof_B.upper))), _mm_mul_pd(tmp_1, _mm_set1_pd
      (control_3dof_B.lower))));
  }

  control_3dof_B.FixPtRelationalOperator = true;
  control_3dof_B.i = 0;
  exitg1 = false;
  while ((!exitg1) && (control_3dof_B.i < 6)) {
    if (!(control_3dof_B.A[control_3dof_B.i] <=
          control_3dof_B.b_wrench[control_3dof_B.i])) {
      control_3dof_B.FixPtRelationalOperator = false;
      exitg1 = true;
    } else {
      control_3dof_B.i++;
    }
  }

  if (control_3dof_B.FixPtRelationalOperator) {
    /* 'PWAS:13' F_fixed = F_trim; */
    /* 'PWAS:14' F_var = F_way; */
    control_3dof_DW.UnitDelay_DSTATE[0] = control_3dof_B.F_trim[0];
    control_3dof_B.F_var[0] = control_3dof_B.F_way[0];
    control_3dof_DW.UnitDelay_DSTATE[1] = control_3dof_B.F_trim[1];
    control_3dof_B.F_var[1] = control_3dof_B.F_way[1];
    control_3dof_DW.UnitDelay_DSTATE[2] = control_3dof_B.F_trim[2];
    control_3dof_B.F_var[2] = control_3dof_B.F_way[2];
  } else {
    control_3dof_B.F_var[0] = 0.0;
    control_3dof_B.F_var[1] = 0.0;
    control_3dof_B.F_var[2] = control_3dof_B.F_trim[2];

    /* 'PWAS:54' retval = all(A_w * F_sp <= b_w); */
    for (control_3dof_B.i = 0; control_3dof_B.i <= 4; control_3dof_B.i += 2) {
      tmp_2 = _mm_loadu_pd(&control_3dof_B.A_wrench[control_3dof_B.i + 6]);
      tmp_0 = _mm_set1_pd(0.0);
      tmp_1 = _mm_loadu_pd(&control_3dof_B.A_wrench[control_3dof_B.i]);
      tmp = _mm_loadu_pd(&control_3dof_B.A_wrench[control_3dof_B.i + 12]);
      _mm_storeu_pd(&control_3dof_B.A[control_3dof_B.i], _mm_add_pd(_mm_add_pd
        (_mm_mul_pd(tmp_2, tmp_0), _mm_mul_pd(tmp_1, tmp_0)), _mm_mul_pd(tmp,
        _mm_set1_pd(control_3dof_B.F_var[2]))));
    }

    control_3dof_B.FixPtRelationalOperator = true;
    control_3dof_B.i = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.i < 6)) {
      if (!(control_3dof_B.A[control_3dof_B.i] <=
            control_3dof_B.b_wrench[control_3dof_B.i])) {
        control_3dof_B.FixPtRelationalOperator = false;
        exitg1 = true;
      } else {
        control_3dof_B.i++;
      }
    }

    if (control_3dof_B.FixPtRelationalOperator) {
      /* 'PWAS:15' elseif WrenchIsValid([0;0;F_trim(3)], A_wrench, b_wrench) */
      /* 'PWAS:17' F_fixed = [0;0;F_trim(3)]; */
      control_3dof_DW.UnitDelay_DSTATE[0] = 0.0;
      control_3dof_DW.UnitDelay_DSTATE[1] = 0.0;
      control_3dof_DW.UnitDelay_DSTATE[2] = control_3dof_B.F_trim[2];

      /* 'PWAS:18' F_var = [F_trim(1);F_trim(2);0]; */
      control_3dof_B.F_var[0] = control_3dof_B.F_trim[0];
      control_3dof_B.F_var[1] = control_3dof_B.F_trim[1];
      control_3dof_B.F_var[2] = 0.0;
    } else {
      /* 'PWAS:19' else */
      /* 'PWAS:21' F_fixed = [0;0;0]; */
      control_3dof_DW.UnitDelay_DSTATE[0] = 0.0;
      control_3dof_DW.UnitDelay_DSTATE[1] = 0.0;
      control_3dof_DW.UnitDelay_DSTATE[2] = 0.0;

      /* 'PWAS:22' F_var = [0;0;F_trim(3)]; */
    }
  }

  /* 'PWAS:25' alpha = PWAS_LOPT(F_fixed, F_var, A_wrench, b_wrench); */
  /* 'PWAS:32' A = A_w * F_var; */
  /* 'PWAS:33' B = b_w - A_w * F_fixed; */
  /* 'PWAS:35' upper = 1; */
  control_3dof_B.upper = 1.0;

  /* 'PWAS:36' lower = 0; */
  control_3dof_B.lower = 0.0;

  /* 'PWAS:37' for i = 1:length(A) */
  control_3dof_B.u0 = control_3dof_B.F_var[0];
  control_3dof_B.F_var_b = control_3dof_B.F_var[1];
  control_3dof_B.F_var_p = control_3dof_B.F_var[2];
  control_3dof_B.UnitDelay_DSTATE = control_3dof_DW.UnitDelay_DSTATE[0];
  control_3dof_B.UnitDelay_DSTATE_c = control_3dof_DW.UnitDelay_DSTATE[1];
  control_3dof_B.UnitDelay_DSTATE_f = control_3dof_DW.UnitDelay_DSTATE[2];
  for (control_3dof_B.i = 0; control_3dof_B.i < 6; control_3dof_B.i++) {
    control_3dof_B.A_wrench_k = control_3dof_B.A_wrench[control_3dof_B.i];
    control_3dof_B.A_c = control_3dof_B.A_wrench_k * control_3dof_B.u0;
    control_3dof_B.F_trim_g = control_3dof_B.A_wrench_k *
      control_3dof_B.UnitDelay_DSTATE;
    control_3dof_B.A_wrench_k = control_3dof_B.A_wrench[control_3dof_B.i + 6];
    control_3dof_B.A_c += control_3dof_B.A_wrench_k * control_3dof_B.F_var_b;
    control_3dof_B.F_trim_g += control_3dof_B.A_wrench_k *
      control_3dof_B.UnitDelay_DSTATE_c;
    control_3dof_B.A_wrench_k = control_3dof_B.A_wrench[control_3dof_B.i + 12];
    control_3dof_B.A_c += control_3dof_B.A_wrench_k * control_3dof_B.F_var_p;
    control_3dof_B.F_trim_g = control_3dof_B.b_wrench[control_3dof_B.i] -
      (control_3dof_B.A_wrench_k * control_3dof_B.UnitDelay_DSTATE_f +
       control_3dof_B.F_trim_g);

    /* 'PWAS:38' if A(i) > 0 */
    if (control_3dof_B.A_c > 0.0) {
      /* 'PWAS:39' upper = min(upper, B(i) / A(i)); */
      control_3dof_B.upper = fmin(control_3dof_B.upper, control_3dof_B.F_trim_g /
        control_3dof_B.A_c);
    } else if (control_3dof_B.A_c < 0.0) {
      /* 'PWAS:40' elseif A(i) < 0 */
      /* 'PWAS:41' lower = max(lower, B(i) / A(i)); */
      control_3dof_B.lower = fmax(control_3dof_B.lower, control_3dof_B.F_trim_g /
        control_3dof_B.A_c);
    }
  }

  /* 'PWAS:45' if lower > upper */
  if (!(control_3dof_B.lower > control_3dof_B.upper)) {
    /* 'PWAS:47' else */
    /* 'PWAS:48' alpha = upper; */
    control_3dof_B.alpha = control_3dof_B.upper;
  } else {
    /* 'PWAS:46' error('No Way'); */
  }

  /* 'PWAS:27' F_actual = F_fixed + alpha * F_var; */
  tmp_2 = _mm_add_pd(_mm_mul_pd(_mm_set1_pd(control_3dof_B.alpha), _mm_loadu_pd(
    &control_3dof_B.F_var[0])), _mm_loadu_pd(&control_3dof_DW.UnitDelay_DSTATE[0]));
  _mm_storeu_pd(&control_3dof_DW.UnitDelay_DSTATE[0], tmp_2);
  control_3dof_DW.UnitDelay_DSTATE[2] += control_3dof_B.alpha *
    control_3dof_B.F_var[2];

  /* 'PWAS:28' result = pinv(Q) * F_actual; */
  control_3dof_pinv(control_3dof_B.Q_new, control_3dof_B.dv);

  /* '<S18>:1:27' aL_sp = [0;0;-9.8] + a_way - (F_trim + F_way - F_actual) / mL; */
  for (control_3dof_B.i = 0; control_3dof_B.i <= 0; control_3dof_B.i += 2) {
    tmp_2 = _mm_loadu_pd(&control_3dof_B.dv[control_3dof_B.i + 3]);
    tmp_0 = _mm_loadu_pd(&control_3dof_B.dv[control_3dof_B.i]);
    tmp_1 = _mm_loadu_pd(&control_3dof_B.dv[control_3dof_B.i + 6]);
    _mm_storeu_pd(&control_3dof_B.F_var[control_3dof_B.i], _mm_add_pd(_mm_add_pd
      (_mm_mul_pd(tmp_2, _mm_set1_pd(control_3dof_DW.UnitDelay_DSTATE[1])),
       _mm_mul_pd(tmp_0, _mm_set1_pd(control_3dof_DW.UnitDelay_DSTATE[0]))),
      _mm_mul_pd(tmp_1, _mm_set1_pd(control_3dof_DW.UnitDelay_DSTATE[2]))));
    tmp_2 = _mm_loadu_pd(&control_3dof_B.acc_limit[control_3dof_B.i]);
    tmp_0 = _mm_loadu_pd(&control_3dof_B.F_trim[control_3dof_B.i]);
    tmp_1 = _mm_loadu_pd(&control_3dof_B.F_way[control_3dof_B.i]);
    tmp = _mm_loadu_pd(&control_3dof_DW.UnitDelay_DSTATE[control_3dof_B.i]);
    _mm_storeu_pd(&control_3dof_B.F_trim[control_3dof_B.i], _mm_sub_pd
                  (_mm_add_pd(_mm_loadu_pd(&b[control_3dof_B.i]), tmp_2),
                   _mm_div_pd(_mm_sub_pd(_mm_add_pd(tmp_0, tmp_1), tmp),
      _mm_set1_pd(CONTROL_PARAM.MASS_LOAD))));
  }

  for (control_3dof_B.i = 2; control_3dof_B.i < 3; control_3dof_B.i++) {
    control_3dof_B.F_var[control_3dof_B.i] = (control_3dof_B.dv[control_3dof_B.i
      + 3] * control_3dof_DW.UnitDelay_DSTATE[1] +
      control_3dof_B.dv[control_3dof_B.i] * control_3dof_DW.UnitDelay_DSTATE[0])
      + control_3dof_B.dv[control_3dof_B.i + 6] *
      control_3dof_DW.UnitDelay_DSTATE[2];
    control_3dof_B.F_trim[control_3dof_B.i] = (b_0[control_3dof_B.i] +
      control_3dof_B.acc_limit[control_3dof_B.i]) -
      ((control_3dof_B.F_trim[control_3dof_B.i] +
        control_3dof_B.F_way[control_3dof_B.i]) -
       control_3dof_DW.UnitDelay_DSTATE[control_3dof_B.i]) /
      CONTROL_PARAM.MASS_LOAD;
  }

  /* '<S18>:1:29' F_actual_div_m = F_actual / mL; */
  /* '<S18>:1:31' control_in1.timestamp = uint64(1); */
  /* '<S18>:1:32' control_in1.vforce = vu(1) * q_1; */
  /* '<S18>:1:33' control_in1.q_sp = Q_sp(:, 1); */
  /* '<S18>:1:34' control_in1.q = q_1; */
  /* '<S18>:1:35' control_in1.w = w_1; */
  /* '<S18>:1:36' control_in1.ai_sp = aL_sp; */
  /* '<S18>:1:37' control_in1.pi = Payload_Out.p_1; */
  /* '<S18>:1:38' control_in1.vi = Payload_Out.v_1; */
  /* '<S18>:1:40' control_in2.timestamp = uint64(1); */
  /* '<S18>:1:41' control_in2.vforce = vu(2) * q_2; */
  control_3dof_DW.UnitDelay_DSTATE[0] /= CONTROL_PARAM.MASS_LOAD;
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_2[0], control_3dof_B.F_var[0]), _mm_set_pd
    (control_3dof_B.F_var[1], control_3dof_B.e_1[0])));
  control_3dof_B.acc_limit[0] = control_3dof_B.dv1[0];
  control_3dof_B.F_way[0] = control_3dof_B.dv1[1];
  control_3dof_DW.UnitDelay_DSTATE[1] /= CONTROL_PARAM.MASS_LOAD;
  tmp_2 = _mm_mul_pd(_mm_loadu_pd(&control_3dof_B.F_var[0]), _mm_set_pd
                     (control_3dof_B.q_2[1], control_3dof_B.e_1[1]));
  _mm_storeu_pd(&control_3dof_B.dv1[0], tmp_2);

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  UnitDelay: '<S4>/Unit Delay'
   */
  control_3dof_B.acc_limit[1] = control_3dof_B.dv1[0];
  control_3dof_B.F_way[1] = control_3dof_B.dv1[1];
  control_3dof_DW.UnitDelay_DSTATE[2] /= CONTROL_PARAM.MASS_LOAD;
  tmp_2 = _mm_mul_pd(_mm_loadu_pd(&control_3dof_B.F_var[0]), _mm_set_pd
                     (control_3dof_B.q_2[2], control_3dof_B.e_1[2]));
  _mm_storeu_pd(&control_3dof_B.dv1[0], tmp_2);

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' */
  control_3dof_B.acc_limit[2] = control_3dof_B.dv1[0];
  control_3dof_B.F_way[2] = control_3dof_B.dv1[1];

  /* '<S18>:1:42' control_in2.q_sp = Q_sp(:, 2); */
  /* '<S18>:1:43' control_in2.q = q_2; */
  /* '<S18>:1:44' control_in2.w = w_2; */
  /* '<S18>:1:45' control_in2.ai_sp = aL_sp; */
  /* '<S18>:1:46' control_in2.pi = Payload_Out.p_2; */
  /* '<S18>:1:47' control_in2.vi = Payload_Out.v_2; */
  /* '<S18>:1:49' control_in3.timestamp = uint64(1); */
  /* '<S18>:1:50' control_in3.vforce = vu(3) * q_3; */
  control_3dof_B.u0 = control_3dof_B.F_var[2];
  tmp_2 = _mm_mul_pd(_mm_set1_pd(control_3dof_B.u0), _mm_loadu_pd
                     (&control_3dof_B.q_3[0]));
  _mm_storeu_pd(&control_3dof_B.F_var[0], tmp_2);
  control_3dof_B.F_var[2] = control_3dof_B.u0 * control_3dof_B.e_2_c;

  /* MATLAB Function: '<S7>/Vertical Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  Switch: '<S20>/Switch'
   */
  /* '<S18>:1:51' control_in3.q_sp = Q_sp(:, 3); */
  /* '<S18>:1:52' control_in3.q = q_3; */
  /* '<S18>:1:53' control_in3.w = w_3; */
  /* '<S18>:1:54' control_in3.ai_sp = aL_sp; */
  /* '<S18>:1:55' control_in3.pi = Payload_Out.p_3; */
  /* '<S18>:1:56' control_in3.vi = Payload_Out.v_3; */
  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[0],
    control_3dof_B.e_1, control_3dof_B.e_2, control_3dof_B.F_trim,
    control_3dof_B.f_vertical, control_3dof_B.dot_err_g, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl);

  /* MATLAB Function: '<S11>/Vertical Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  Switch: '<S20>/Switch'
   */
  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[3],
    control_3dof_B.q_2, control_3dof_B.rtb_control_in2_w, control_3dof_B.F_trim,
    control_3dof_Y.force_sp1, control_3dof_B.dot_err_c, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_m);

  /* MATLAB Function: '<S15>/Vertical Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  Switch: '<S20>/Switch'
   */
  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[6],
    control_3dof_B.q_3, control_3dof_B.rtb_control_in3_w, control_3dof_B.F_trim,
    control_3dof_Y.force_sp2, control_3dof_B.dot_err, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_h);

  /* MATLAB Function: '<S14>/Parallel Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   */
  control_3dof_ParallelControl(control_3dof_B.q_3,
    control_3dof_B.rtb_control_in3_w, control_3dof_B.F_var,
    control_3dof_B.F_trim, control_3dof_Y.force_sp3, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_p);

  /* Sum: '<S15>/Sum' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S3>/Sum1'
   */
  tmp_2 = _mm_add_pd(_mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp2[0]),
    _mm_loadu_pd(&control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0])),
                     _mm_loadu_pd(&control_3dof_Y.force_sp3[0]));

  /* DiscreteIntegrator: '<S11>/Discrete-Time Integrator' incorporates:
   *  Outport: '<Root>/force_sp3'
   *  Sum: '<S3>/Sum1'
   */
  _mm_storeu_pd(&control_3dof_Y.force_sp3[0], tmp_2);

  /* Outport: '<Root>/force_sp3' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S15>/Sum'
   *  Sum: '<S3>/Sum1'
   */
  control_3dof_Y.force_sp3[2] += control_3dof_Y.force_sp2[2] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2];

  /* MATLAB Function: '<S10>/Parallel Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   */
  control_3dof_ParallelControl(control_3dof_B.q_2,
    control_3dof_B.rtb_control_in2_w, control_3dof_B.F_way,
    control_3dof_B.F_trim, control_3dof_Y.force_sp2, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_k);

  /* Sum: '<S11>/Sum' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S2>/Sum1'
   */
  tmp_2 = _mm_add_pd(_mm_add_pd(_mm_loadu_pd(&control_3dof_Y.force_sp1[0]),
    _mm_loadu_pd(&control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0])),
                     _mm_loadu_pd(&control_3dof_Y.force_sp2[0]));

  /* DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
   *  Outport: '<Root>/force_sp2'
   *  Sum: '<S2>/Sum1'
   */
  _mm_storeu_pd(&control_3dof_Y.force_sp2[0], tmp_2);

  /* Outport: '<Root>/force_sp2' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S11>/Sum'
   *  Sum: '<S2>/Sum1'
   */
  control_3dof_Y.force_sp2[2] += control_3dof_Y.force_sp1[2] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2];

  /* MATLAB Function: '<S6>/Parallel Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   */
  control_3dof_ParallelControl(control_3dof_B.e_1, control_3dof_B.e_2,
    control_3dof_B.acc_limit, control_3dof_B.F_trim, control_3dof_Y.force_sp1,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl);

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  DiscreteIntegrator: '<S22>/Discrete-Time Integrator'
   *  Outport: '<Root>/state'
   */
  control_3dof_Y.state.timestamp =
    control_3dof_DW.DiscreteTimeIntegrator_DSTAT_ps;

  /* Update for UnitDelay: '<S21>/Delay Input1' incorporates:
   *  Inport: '<Root>/Dir_sp'
   *
   * Block description for '<S21>/Delay Input1':
   *
   *  Store in Global RAM
   */
  control_3dof_DW.DelayInput1_DSTATE = control_3dof_U.Dir_sp.timestamp;

  /* Update for Delay: '<S20>/Delay' */
  control_3dof_DW.icLoad = false;

  /* Outport: '<Root>/force_sp1' incorporates:
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S1>/Sum1'
   *  Sum: '<S7>/Sum'
   */
  control_3dof_Y.force_sp1[0] += control_3dof_B.f_vertical[0] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[0];

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Outport: '<Root>/state'
   * */
  control_3dof_Y.state.dL[0] = control_3dof_DW.obj.Z_3[0];

  /* Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator' */
  tmp_2 = _mm_set1_pd(0.02);

  /* Gain: '<S15>/Gain' */
  tmp_0 = _mm_set1_pd(CONTROL_PARAM.KQI);

  /* Gain: '<S11>/Gain' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
   *  Gain: '<S15>/Gain'
   */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_0,
    _mm_set_pd(control_3dof_B.dot_err_c[0], control_3dof_B.dot_err[0])), tmp_2),
    _mm_set_pd(control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0],
               control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0])));

  /* Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator' */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0] = control_3dof_B.dv1[0];

  /* Update for DiscreteIntegrator: '<S11>/Discrete-Time Integrator' */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0] = control_3dof_B.dv1[1];

  /* Update for DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S7>/Gain'
   *  Sum: '<S1>/Sum1'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[0] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[0] * 0.02;

  /* Outport: '<Root>/force_sp1' incorporates:
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S1>/Sum1'
   *  Sum: '<S7>/Sum'
   */
  control_3dof_Y.force_sp1[1] += control_3dof_B.f_vertical[1] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[1];

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Outport: '<Root>/state'
   * */
  control_3dof_Y.state.dL[1] = control_3dof_DW.obj.Z_3[1];

  /* Gain: '<S11>/Gain' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
   *  Gain: '<S15>/Gain'
   */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_0,
    _mm_set_pd(control_3dof_B.dot_err_c[1], control_3dof_B.dot_err[1])), tmp_2),
    _mm_set_pd(control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[1],
               control_3dof_DW.DiscreteTimeIntegrator_DSTATE[1])));

  /* Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator' */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[1] = control_3dof_B.dv1[0];

  /* Update for DiscreteIntegrator: '<S11>/Discrete-Time Integrator' */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[1] = control_3dof_B.dv1[1];

  /* Update for DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S7>/Gain'
   *  Sum: '<S1>/Sum1'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[1] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[1] * 0.02;

  /* Outport: '<Root>/force_sp1' incorporates:
   *  DiscreteIntegrator: '<S7>/Discrete-Time Integrator'
   *  Sum: '<S1>/Sum1'
   *  Sum: '<S7>/Sum'
   */
  control_3dof_Y.force_sp1[2] += control_3dof_B.f_vertical[2] +
    control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[2];

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Outport: '<Root>/state'
   * */
  control_3dof_Y.state.dL[2] = control_3dof_DW.obj.Z_3[2];

  /* Gain: '<S11>/Gain' incorporates:
   *  DiscreteIntegrator: '<S11>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S15>/Discrete-Time Integrator'
   *  Gain: '<S15>/Gain'
   */
  _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_add_pd(_mm_mul_pd(_mm_mul_pd(tmp_0,
    _mm_set_pd(control_3dof_B.dot_err_c[2], control_3dof_B.dot_err[2])), tmp_2),
    _mm_set_pd(control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2],
               control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2])));

  /* Update for DiscreteIntegrator: '<S15>/Discrete-Time Integrator' */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2] = control_3dof_B.dv1[0];

  /* Update for DiscreteIntegrator: '<S11>/Discrete-Time Integrator' */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2] = control_3dof_B.dv1[1];

  /* Update for DiscreteIntegrator: '<S7>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S7>/Gain'
   *  Sum: '<S1>/Sum1'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[2] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[2] * 0.02;
}

/* Model initialize function */
void control_3dof_initialize(void)
{
  /* InitializeConditions for Delay: '<S20>/Delay' */
  control_3dof_DW.icLoad = true;

  /* Start for MATLABSystem: '<S4>/Position 2nd ESO' */
  /*  Constructor */
  /* 'pos_2nd_eso:1' matlab.System */
  /*  Support name-value pair arguments when constructing object */
  /* 'pos_2nd_eso:41' setProperties(obj,nargin,varargin{:}) */
  /* 'pos_2nd_eso:13' (3, 1) */
  /* 'pos_2nd_eso:13' K_p */
  /* 'pos_2nd_eso:15' (3, 1) */
  /* 'pos_2nd_eso:15' K_v */
  control_3dof_DW.obj.isInitialized = 1;

  /*         %% Common functions */
  /*  Perform one-time calculations, such as computing constants */
  /* 'pos_2nd_eso:49' obj.Z_1 = pos; */
  /* 'pos_2nd_eso:50' obj.Z_2 = vel; */
  /* 'pos_2nd_eso:51' obj.Z_3 = zeros(size(pos)); */
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

  /* InitializeConditions for MATLABSystem: '<S4>/Position 2nd ESO' */
  /* 'pos_2nd_eso:52' obj.is_init = false; */
  /* 'pos_2nd_eso:34' logical */
  /* 'pos_2nd_eso:34' (1, 1) */
  /* 'pos_2nd_eso:34' is_init */
  /*  Initialize / reset internal or discrete properties */
  /* 'pos_2nd_eso:90' obj.is_init = false; */
  /* 'pos_2nd_eso:34' logical */
  /* 'pos_2nd_eso:34' (1, 1) */
  /* 'pos_2nd_eso:34' is_init */
  control_3dof_DW.obj.is_init = false;
}

/* Model terminate function */
void control_3dof_terminate(void)
{
  /* (no terminate code required) */
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
