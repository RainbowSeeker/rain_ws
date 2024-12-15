//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: control_3dof.cpp
//
// Code generated for Simulink model 'control_3dof'.
//
// Model version                  : 1.656
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Sun Dec 15 15:40:18 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "control_3dof.h"
#include "rtwtypes.h"
#include <emmintrin.h>
#include <cmath>
#include <cstring>
#include "cmath"
#include "limits"

// Exported block parameters
struct_n4M4IuMQiYjb6H8OxPPEC CONTROL_PARAM{
  1.0,
  1.5,
  1.0,
  1.0,
  20.0,

  { 3.0, 5.0, 10.0 },

  { 0.0, 20.0, 20.0 },

  { 3.0, 5.0, 10.0 },

  { 0.0, 20.0, 20.0 },
  0.2,
  1.0,
  1.0,
  2.0
} ;                                    // Variable: CONTROL_PARAM
                                          //  Referenced by:
                                          //    '<S1>/Gain'
                                          //    '<S1>/Gain1'
                                          //    '<S1>/MATLAB System'
                                          //    '<S2>/Gain'
                                          //    '<S2>/Gain1'
                                          //    '<S2>/MATLAB System'
                                          //    '<S3>/Gain'
                                          //    '<S3>/Gain1'
                                          //    '<S3>/MATLAB System'
                                          //    '<S4>/Force Saturation && Disturbution'
                                          //    '<S4>/Position 2nd ESO'
                                          //    '<S6>/Parallel Control'
                                          //    '<S7>/Vertical Control'
                                          //    '<S11>/Parallel Control'
                                          //    '<S12>/Vertical Control'
                                          //    '<S16>/Parallel Control'
                                          //    '<S17>/Vertical Control'
                                          //    '<S22>/KP'
                                          //    '<S22>/KV'


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
//    '<S1>/MATLAB Function'
//    '<S2>/MATLAB Function'
//    '<S3>/MATLAB Function'
//
void control_3dof::control_3dof_MATLABFunction(const real_T rtu_z3[3], const
  real_T rtu_q[3], real_T rty_y[3])
{
  real_T tmp[2];
  real_T b_idx_2;
  _mm_storeu_pd(&tmp[0], _mm_sub_pd(_mm_mul_pd(_mm_set_pd(rtu_z3[0], rtu_q[1]),
    _mm_set_pd(rtu_q[2], rtu_z3[2])), _mm_mul_pd(_mm_set_pd(rtu_q[0], rtu_z3[1]),
    _mm_set_pd(rtu_z3[2], rtu_q[2]))));
  b_idx_2 = rtu_q[0] * rtu_z3[1] - rtu_z3[0] * rtu_q[1];
  rty_y[0] = -(rtu_q[1] * b_idx_2 - tmp[1] * rtu_q[2]);
  rty_y[1] = -(tmp[0] * rtu_q[2] - rtu_q[0] * b_idx_2);
  rty_y[2] = -(rtu_q[0] * tmp[1] - tmp[0] * rtu_q[1]);
}

// System initialize for atomic system:
void control_3dof::control_3dof_MATLABSystem_Init(DW_MATLABSystem_control_3dof_T
  *localDW)
{
  // Start for MATLABSystem: '<S1>/MATLAB System'
  //  Constructor
  //  Support name-value pair arguments when constructing object
  localDW->objisempty = true;
  localDW->obj.isInitialized = 1;

  //         %% Common functions
  //  Perform one-time calculations, such as computing constants
  localDW->obj.K_p[0] = CONTROL_PARAM.ESO_PI[0];
  localDW->obj.K_v[0] = CONTROL_PARAM.ESO_VI[0];
  localDW->obj.Z_1[0] = 0.0;
  localDW->obj.Z_2[0] = 0.0;
  localDW->obj.Z_3[0] = 0.0;
  localDW->obj.K_p[1] = CONTROL_PARAM.ESO_PI[1];
  localDW->obj.K_v[1] = CONTROL_PARAM.ESO_VI[1];
  localDW->obj.Z_1[1] = 0.0;
  localDW->obj.Z_2[1] = 0.0;
  localDW->obj.Z_3[1] = 0.0;
  localDW->obj.K_p[2] = CONTROL_PARAM.ESO_PI[2];
  localDW->obj.K_v[2] = CONTROL_PARAM.ESO_VI[2];
  localDW->obj.Z_1[2] = 0.0;
  localDW->obj.Z_2[2] = 0.0;
  localDW->obj.Z_3[2] = 0.0;

  // InitializeConditions for MATLABSystem: '<S1>/MATLAB System'
  //  Initialize / reset internal or discrete properties
  localDW->obj.is_init = false;
}

// Output and update for atomic system:
void control_3dof::control_3dof_MATLABSystem(const real_T rtu_0[3], const real_T
  rtu_1[3], const real_T rtu_2[3], B_MATLABSystem_control_3dof_T *localB,
  DW_MATLABSystem_control_3dof_T *localDW)
{
  real_T tmp[2];
  real_T obj_Z_3;

  // MATLABSystem: '<S1>/MATLAB System'
  localDW->obj.K_p[0] = CONTROL_PARAM.ESO_PI[0];
  localDW->obj.K_v[0] = CONTROL_PARAM.ESO_VI[0];
  localDW->obj.K_p[1] = CONTROL_PARAM.ESO_PI[1];
  localDW->obj.K_v[1] = CONTROL_PARAM.ESO_VI[1];
  localDW->obj.K_p[2] = CONTROL_PARAM.ESO_PI[2];
  localDW->obj.K_v[2] = CONTROL_PARAM.ESO_VI[2];
  if (!localDW->obj.is_init) {
    localDW->obj.is_init = true;
    localDW->obj.Z_1[0] = rtu_1[0];
    localDW->obj.Z_2[0] = rtu_2[0];
    localDW->obj.Z_3[0] = 0.0;
    localDW->obj.Z_1[1] = rtu_1[1];
    localDW->obj.Z_2[1] = rtu_2[1];
    localDW->obj.Z_3[1] = 0.0;
    localDW->obj.Z_1[2] = rtu_1[2];
    localDW->obj.Z_2[2] = rtu_2[2];
    localDW->obj.Z_3[2] = 0.0;
  }

  // Start for MATLABSystem: '<S1>/MATLAB System'
  //  Extended state observer
  //  Update states
  _mm_storeu_pd(&tmp[0], _mm_sub_pd(_mm_set_pd(rtu_2[0], rtu_1[0]), _mm_set_pd
    (localDW->obj.Z_2[0], localDW->obj.Z_1[0])));

  // MATLABSystem: '<S1>/MATLAB System'
  localDW->obj.Z_1[0] += (localDW->obj.K_p[0] * tmp[0] + localDW->obj.Z_2[0]) *
    0.02;
  localDW->obj.Z_2[0] += (((rtu_0[0] + localDW->obj.Z_3[0]) + tmp[0] *
    localDW->obj.K_p[1]) + localDW->obj.K_v[1] * tmp[1]) * 0.02;
  obj_Z_3 = (tmp[0] * localDW->obj.K_p[2] + tmp[1] * localDW->obj.K_v[2]) * 0.02
    + localDW->obj.Z_3[0];
  localDW->obj.Z_3[0] = obj_Z_3;

  // MATLABSystem: '<S1>/MATLAB System'
  localB->MATLABSystem_o3[0] = obj_Z_3;

  // Start for MATLABSystem: '<S1>/MATLAB System'
  _mm_storeu_pd(&tmp[0], _mm_sub_pd(_mm_set_pd(rtu_2[1], rtu_1[1]), _mm_set_pd
    (localDW->obj.Z_2[1], localDW->obj.Z_1[1])));

  // MATLABSystem: '<S1>/MATLAB System'
  localDW->obj.Z_1[1] += (localDW->obj.K_p[0] * tmp[0] + localDW->obj.Z_2[1]) *
    0.02;
  localDW->obj.Z_2[1] += (((rtu_0[1] + localDW->obj.Z_3[1]) + tmp[0] *
    localDW->obj.K_p[1]) + localDW->obj.K_v[1] * tmp[1]) * 0.02;
  obj_Z_3 = (tmp[0] * localDW->obj.K_p[2] + tmp[1] * localDW->obj.K_v[2]) * 0.02
    + localDW->obj.Z_3[1];
  localDW->obj.Z_3[1] = obj_Z_3;

  // MATLABSystem: '<S1>/MATLAB System'
  localB->MATLABSystem_o3[1] = obj_Z_3;

  // Start for MATLABSystem: '<S1>/MATLAB System'
  _mm_storeu_pd(&tmp[0], _mm_sub_pd(_mm_set_pd(rtu_2[2], rtu_1[2]), _mm_set_pd
    (localDW->obj.Z_2[2], localDW->obj.Z_1[2])));

  // MATLABSystem: '<S1>/MATLAB System'
  localDW->obj.Z_1[2] += (localDW->obj.K_p[0] * tmp[0] + localDW->obj.Z_2[2]) *
    0.02;
  localDW->obj.Z_2[2] += ((((rtu_0[2] + localDW->obj.Z_3[2]) + 9.8) + tmp[0] *
    localDW->obj.K_p[1]) + localDW->obj.K_v[1] * tmp[1]) * 0.02;
  obj_Z_3 = (tmp[0] * localDW->obj.K_p[2] + tmp[1] * localDW->obj.K_v[2]) * 0.02
    + localDW->obj.Z_3[2];
  localDW->obj.Z_3[2] = obj_Z_3;

  // MATLABSystem: '<S1>/MATLAB System'
  localB->MATLABSystem_o3[2] = obj_Z_3;
}

//
// Output and update for atomic system:
//    '<S6>/Parallel Control'
//    '<S11>/Parallel Control'
//    '<S16>/Parallel Control'
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
//    '<S12>/Vertical Control'
//    '<S17>/Vertical Control'
//
void control_3dof::control_3dof_VerticalControl(const real_T rtu_q_sp[3], const
  real_T rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
  rty_f_vertical[3], real_T rty_state[12], real_T rtp_Kq, real_T rtp_Kw, real_T
  rtp_li, real_T rtp_mi, B_VerticalControl_control_3do_T *localB)
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

  for (localB->i_b = 0; localB->i_b < 3; localB->i_b++) {
    for (localB->i = 0; localB->i <= 0; localB->i += 2) {
      tmp_5 = _mm_loadu_pd(&localB->Sqi[localB->i + 3]);
      tmp_3 = _mm_loadu_pd(&localB->Sqi[localB->i]);
      tmp_4 = _mm_loadu_pd(&localB->Sqi[localB->i + 6]);
      _mm_storeu_pd(&localB->Sqi_m[localB->i + 3 * localB->i_b], _mm_add_pd
                    (_mm_add_pd(_mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_b + 1]), tmp_5), _mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_b]), tmp_3)), _mm_mul_pd(_mm_set1_pd(localB->Sqi[3 *
        localB->i_b + 2]), tmp_4)));
    }

    for (localB->i = 2; localB->i < 3; localB->i++) {
      localB->Sqi_m[localB->i + 3 * localB->i_b] = (localB->Sqi[3 * localB->i_b
        + 1] * localB->Sqi[localB->i + 3] + localB->Sqi[3 * localB->i_b] *
        localB->Sqi[localB->i]) + localB->Sqi[3 * localB->i_b + 2] * localB->
        Sqi[localB->i + 6];
    }
  }

  localB->a = rtp_mi * rtp_li;
  localB->rtu_q = 0.0;
  for (localB->i = 0; localB->i < 3; localB->i++) {
    localB->e_w_tmp_tmp = (localB->Sqi_m[localB->i + 3] * localB->w_sp[1] +
      localB->Sqi_m[localB->i] * localB->w_sp[0]) + localB->Sqi_m[localB->i + 6]
      * localB->w_sp[2];
    localB->e_w_tmp[localB->i] = localB->e_w_tmp_tmp;
    localB->e_w[localB->i] = rtu_w[localB->i] + localB->e_w_tmp_tmp;
    localB->rtu_q += rtu_q[localB->i] * localB->w_sp[localB->i];
  }

  _mm_storeu_pd(&localB->b_a[0], _mm_mul_pd(_mm_sub_pd(_mm_mul_pd(_mm_set_pd
    (rtu_q_sp[0], localB->w_sp[1]), _mm_set_pd(localB->w_sp[2], rtu_q_sp[2])),
    _mm_mul_pd(_mm_set_pd(localB->w_sp[0], rtu_q_sp[1]), _mm_set_pd(rtu_q_sp[2],
    localB->w_sp[2]))), _mm_set1_pd(localB->rtu_q)));
  localB->b_a[2] = (localB->w_sp[0] * rtu_q_sp[1] - rtu_q_sp[0] * localB->w_sp[1])
    * localB->rtu_q;
  for (localB->i_b = 0; localB->i_b < 3; localB->i_b++) {
    localB->rtp_Kw[localB->i_b] = (-rtp_Kw * localB->e_w[localB->i_b] -
      localB->e_w_tmp[localB->i_b]) - localB->b_a[localB->i_b];
    for (localB->i = 0; localB->i < 3; localB->i++) {
      localB->Sqi_m[localB->i_b + 3 * localB->i] = (localB->Sqi[localB->i_b + 3]
        * rtp_mi * localB->Sqi[3 * localB->i + 1] + rtp_mi * localB->Sqi
        [localB->i_b] * localB->Sqi[3 * localB->i]) + localB->Sqi[localB->i_b +
        6] * rtp_mi * localB->Sqi[3 * localB->i + 2];
    }
  }

  localB->rtu_q = localB->rtp_Kw[1];
  localB->e_w_tmp_tmp = localB->rtp_Kw[0];
  localB->rtp_Kw_c = localB->rtp_Kw[2];
  localB->rtu_ai_sp = rtu_ai_sp[1];
  localB->rtu_ai_sp_k = rtu_ai_sp[0];
  localB->rtu_ai_sp_c = rtu_ai_sp[2];
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
      (_mm_mul_pd(_mm_mul_pd(tmp_5, tmp_3), _mm_set1_pd(localB->rtu_q)),
       _mm_mul_pd(_mm_mul_pd(tmp_3, tmp_4), _mm_set1_pd(localB->e_w_tmp_tmp))),
      _mm_mul_pd(_mm_mul_pd(tmp, tmp_3), _mm_set1_pd(localB->rtp_Kw_c))),
      _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_0, _mm_set1_pd(localB->rtu_ai_sp)),
      _mm_mul_pd(tmp_1, _mm_set1_pd(localB->rtu_ai_sp_k))), _mm_mul_pd(tmp_2,
      _mm_set1_pd(localB->rtu_ai_sp_c)))));
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
    rty_f_vertical[localB->i] = ((localB->Sqi[localB->i + 3] * localB->a *
      localB->rtu_q + localB->a * localB->Sqi[localB->i] * localB->e_w_tmp_tmp)
      + localB->Sqi[localB->i + 6] * localB->a * localB->rtp_Kw_c) -
      ((localB->Sqi_m[localB->i + 3] * localB->rtu_ai_sp + localB->Sqi_m
        [localB->i] * localB->rtu_ai_sp_k) + localB->Sqi_m[localB->i + 6] *
       localB->rtu_ai_sp_c);
    rty_state[localB->i] = rtu_q_sp[localB->i];
    rty_state[localB->i + 3] = localB->e_q[localB->i];
    rty_state[localB->i + 6] = localB->w_sp[localB->i];
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
  control_3dof_B.cfromc_c = cfrom;
  control_3dof_B.ctoc_b = cto;
  control_3dof_B.notdone_d = true;
  while (control_3dof_B.notdone_d) {
    control_3dof_B.cfrom1_p = control_3dof_B.cfromc_c * 2.0041683600089728E-292;
    control_3dof_B.cto1_c = control_3dof_B.ctoc_b / 4.9896007738368E+291;
    if ((std::abs(control_3dof_B.cfrom1_p) > std::abs(control_3dof_B.ctoc_b)) &&
        (control_3dof_B.ctoc_b != 0.0)) {
      control_3dof_B.mul_f = 2.0041683600089728E-292;
      control_3dof_B.cfromc_c = control_3dof_B.cfrom1_p;
    } else if (std::abs(control_3dof_B.cto1_c) > std::abs
               (control_3dof_B.cfromc_c)) {
      control_3dof_B.mul_f = 4.9896007738368E+291;
      control_3dof_B.ctoc_b = control_3dof_B.cto1_c;
    } else {
      control_3dof_B.mul_f = control_3dof_B.ctoc_b / control_3dof_B.cfromc_c;
      control_3dof_B.notdone_d = false;
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
          (control_3dof_B.mul_f)));
      }

      for (control_3dof_B.b_i_g = control_3dof_B.scalarLB_l;
           control_3dof_B.b_i_g < m; control_3dof_B.b_i_g++) {
        control_3dof_B.i1 = (control_3dof_B.b_i_g + control_3dof_B.offset_j) + 1;
        A[control_3dof_B.i1] *= control_3dof_B.mul_f;
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
      control_3dof_B.scale_g1 = 3.3121686421112381E-170;
      control_3dof_B.kend_h = ix0 + n;
      for (control_3dof_B.k_bn = ix0; control_3dof_B.k_bn <
           control_3dof_B.kend_h; control_3dof_B.k_bn++) {
        control_3dof_B.absxk_m = std::abs(x[control_3dof_B.k_bn - 1]);
        if (control_3dof_B.absxk_m > control_3dof_B.scale_g1) {
          control_3dof_B.t_n = control_3dof_B.scale_g1 / control_3dof_B.absxk_m;
          y = y * control_3dof_B.t_n * control_3dof_B.t_n + 1.0;
          control_3dof_B.scale_g1 = control_3dof_B.absxk_m;
        } else {
          control_3dof_B.t_n = control_3dof_B.absxk_m / control_3dof_B.scale_g1;
          y += control_3dof_B.t_n * control_3dof_B.t_n;
        }
      }

      y = control_3dof_B.scale_g1 * std::sqrt(y);
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
      for (control_3dof_B.k_l = ix0; control_3dof_B.k_l < control_3dof_B.kend;
           control_3dof_B.k_l++) {
        control_3dof_B.absxk = std::abs(x[control_3dof_B.k_l - 1]);
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
    control_3dof_B.scalarLB_n = (n / 2) << 1;
    control_3dof_B.vectorUB_b = control_3dof_B.scalarLB_n - 2;
    for (control_3dof_B.k_b = 0; control_3dof_B.k_b <= control_3dof_B.vectorUB_b;
         control_3dof_B.k_b += 2) {
      __m128d tmp;
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_b) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i3]);
      _mm_storeu_pd(&y[control_3dof_B.i3], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k_b) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k_b = control_3dof_B.scalarLB_n; control_3dof_B.k_b < n;
         control_3dof_B.k_b++) {
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_b) - 1;
      y[control_3dof_B.i3] += x[(ix0 + control_3dof_B.k_b) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S4>/Force Saturation && Disturbution'
void control_3dof::control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3],
  int32_T ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    control_3dof_B.scalarLB_lx = (n / 2) << 1;
    control_3dof_B.vectorUB_o = control_3dof_B.scalarLB_lx - 2;
    for (control_3dof_B.k_d = 0; control_3dof_B.k_d <= control_3dof_B.vectorUB_o;
         control_3dof_B.k_d += 2) {
      __m128d tmp;
      control_3dof_B.i2 = (iy0 + control_3dof_B.k_d) - 1;
      tmp = _mm_loadu_pd(&y[control_3dof_B.i2]);
      _mm_storeu_pd(&y[control_3dof_B.i2], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
        (&x[(ix0 + control_3dof_B.k_d) - 1]), _mm_set1_pd(a)), tmp));
    }

    for (control_3dof_B.k_d = control_3dof_B.scalarLB_lx; control_3dof_B.k_d < n;
         control_3dof_B.k_d++) {
      control_3dof_B.i2 = (iy0 + control_3dof_B.k_d) - 1;
      y[control_3dof_B.i2] += x[(ix0 + control_3dof_B.k_d) - 1] * a;
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

    for (control_3dof_B.j = 0; control_3dof_B.j < n; control_3dof_B.j++) {
      control_3dof_B.offset = (control_3dof_B.j * lda + iA0) - 2;
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
      _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.m], control_3dof_B.rt), tmp));
      control_3dof_B.rt = control_3dof_B.dv1[0];
      control_3dof_B.smm1 = control_3dof_B.dv1[1];
      _mm_storeu_pd(&control_3dof_B.dv1[0], _mm_div_pd(_mm_set_pd
        (control_3dof_B.b_s[control_3dof_B.qq],
         control_3dof_B.e[control_3dof_B.m]), tmp));
      control_3dof_B.smm1 = ((control_3dof_B.smm1 + control_3dof_B.rt) *
        (control_3dof_B.smm1 - control_3dof_B.rt) + control_3dof_B.dv1[0] *
        control_3dof_B.dv1[0]) / 2.0;
      control_3dof_B.c = control_3dof_B.rt * control_3dof_B.dv1[0];
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

      control_3dof_B.rt = (control_3dof_B.dv1[1] + control_3dof_B.rt) *
        (control_3dof_B.dv1[1] - control_3dof_B.rt) + control_3dof_B.shift;
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq] / control_3dof_B.r *
        control_3dof_B.dv1[1];
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
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  int32_T exitg1;
  boolean_T exitg2;

  // MATLAB Function: '<S4>/Q sp'
  control_3dof_B.y[0] = 0.24999999999999994;
  control_3dof_B.y[1] = 0.43301270189221946;
  control_3dof_B.y[2] = 0.8660254037844386;

  // Gain: '<S22>/KP' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   Sum: '<S22>/Sum'

  control_3dof_B.u0 = (control_3dof_U.Traj_sp.pos_sp[0] -
                       control_3dof_U.Payload_Out.pL[0]) * CONTROL_PARAM.KP;

  // Saturate: '<S22>/vel_limit'
  if (control_3dof_B.u0 > 2.0) {
    control_3dof_B.u0 = 2.0;
  } else if (control_3dof_B.u0 < -2.0) {
    control_3dof_B.u0 = -2.0;
  }

  // Gain: '<S22>/KV' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   Saturate: '<S22>/vel_limit'
  //   Sum: '<S22>/Sum2'
  //   Sum: '<S22>/Sum3'

  control_3dof_B.u0 = ((control_3dof_B.u0 + control_3dof_U.Traj_sp.vel_sp[0]) -
                       control_3dof_U.Payload_Out.vL[0]) * CONTROL_PARAM.KV;

  // Saturate: '<S22>/acc_limit'
  if (control_3dof_B.u0 > 5.0) {
    control_3dof_B.u0 = 5.0;
  } else if (control_3dof_B.u0 < -5.0) {
    control_3dof_B.u0 = -5.0;
  }

  // Bias: '<S22>/-g' incorporates:
  //   Inport: '<Root>/Traj_sp'
  //   Saturate: '<S22>/acc_limit'
  //   Sum: '<S1>/Sum'
  //   Sum: '<S22>/Sum1'

  control_3dof_Y.force_sp1[0] = control_3dof_B.u0 +
    control_3dof_U.Traj_sp.acc_ff[0];

  // MATLABSystem: '<S4>/Position 2nd ESO'
  control_3dof_DW.obj.K_p[0] = CONTROL_PARAM.ESO_PL[0];
  control_3dof_DW.obj.K_v[0] = CONTROL_PARAM.ESO_VL[0];

  // MATLAB Function: '<S4>/Q sp'
  control_3dof_B.y[3] = -0.50000000000000011;
  control_3dof_B.y[4] = -0.0;
  control_3dof_B.y[5] = 0.8660254037844386;

  // Gain: '<S22>/KP' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   Sum: '<S22>/Sum'

  control_3dof_B.u0 = (control_3dof_U.Traj_sp.pos_sp[1] -
                       control_3dof_U.Payload_Out.pL[1]) * CONTROL_PARAM.KP;

  // Saturate: '<S22>/vel_limit'
  if (control_3dof_B.u0 > 2.0) {
    control_3dof_B.u0 = 2.0;
  } else if (control_3dof_B.u0 < -2.0) {
    control_3dof_B.u0 = -2.0;
  }

  // Gain: '<S22>/KV' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   Saturate: '<S22>/vel_limit'
  //   Sum: '<S22>/Sum2'
  //   Sum: '<S22>/Sum3'

  control_3dof_B.u0 = ((control_3dof_B.u0 + control_3dof_U.Traj_sp.vel_sp[1]) -
                       control_3dof_U.Payload_Out.vL[1]) * CONTROL_PARAM.KV;

  // Saturate: '<S22>/acc_limit'
  if (control_3dof_B.u0 > 5.0) {
    control_3dof_B.u0 = 5.0;
  } else if (control_3dof_B.u0 < -5.0) {
    control_3dof_B.u0 = -5.0;
  }

  // Bias: '<S22>/-g' incorporates:
  //   Inport: '<Root>/Traj_sp'
  //   Saturate: '<S22>/acc_limit'
  //   Sum: '<S1>/Sum'
  //   Sum: '<S22>/Sum1'

  control_3dof_Y.force_sp1[1] = control_3dof_B.u0 +
    control_3dof_U.Traj_sp.acc_ff[1];

  // MATLABSystem: '<S4>/Position 2nd ESO'
  control_3dof_DW.obj.K_p[1] = CONTROL_PARAM.ESO_PL[1];
  control_3dof_DW.obj.K_v[1] = CONTROL_PARAM.ESO_VL[1];

  // MATLAB Function: '<S4>/Q sp'
  control_3dof_B.y[6] = 0.24999999999999994;
  control_3dof_B.y[7] = -0.43301270189221946;
  control_3dof_B.y[8] = 0.8660254037844386;

  // Gain: '<S22>/KP' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   Sum: '<S22>/Sum'

  control_3dof_B.u0 = (control_3dof_U.Traj_sp.pos_sp[2] -
                       control_3dof_U.Payload_Out.pL[2]) * CONTROL_PARAM.KP;

  // Saturate: '<S22>/vel_limit'
  if (control_3dof_B.u0 > 2.0) {
    control_3dof_B.u0 = 2.0;
  } else if (control_3dof_B.u0 < -2.0) {
    control_3dof_B.u0 = -2.0;
  }

  // Gain: '<S22>/KV' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   Inport: '<Root>/Traj_sp'
  //   Saturate: '<S22>/vel_limit'
  //   Sum: '<S22>/Sum2'
  //   Sum: '<S22>/Sum3'

  control_3dof_B.u0 = ((control_3dof_B.u0 + control_3dof_U.Traj_sp.vel_sp[2]) -
                       control_3dof_U.Payload_Out.vL[2]) * CONTROL_PARAM.KV;

  // Saturate: '<S22>/acc_limit'
  if (control_3dof_B.u0 > 5.0) {
    control_3dof_B.u0 = 5.0;
  } else if (control_3dof_B.u0 < -5.0) {
    control_3dof_B.u0 = -5.0;
  }

  // Bias: '<S22>/-g' incorporates:
  //   Inport: '<Root>/Traj_sp'
  //   Saturate: '<S22>/acc_limit'
  //   Sum: '<S1>/Sum'
  //   Sum: '<S22>/Sum1'

  control_3dof_Y.force_sp1[2] = (control_3dof_B.u0 +
    control_3dof_U.Traj_sp.acc_ff[2]) - 9.8;

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'

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
  _mm_storeu_pd(&control_3dof_B.dv[0], _mm_sub_pd(_mm_set_pd
    (control_3dof_U.Payload_Out.vL[0], control_3dof_U.Payload_Out.pL[0]),
    _mm_set_pd(control_3dof_DW.obj.Z_2[0], control_3dof_DW.obj.Z_1[0])));

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_DW.obj.Z_1[0] += (control_3dof_DW.obj.K_p[0] * control_3dof_B.dv
    [0] + control_3dof_DW.obj.Z_2[0]) * 0.02;
  control_3dof_DW.obj.Z_2[0] += (((control_3dof_DW.UnitDelay_DSTATE[0] +
    control_3dof_DW.obj.Z_3[0]) + control_3dof_B.dv[0] *
    control_3dof_DW.obj.K_p[1]) + control_3dof_DW.obj.K_v[1] *
    control_3dof_B.dv[1]) * 0.02;
  control_3dof_DW.obj.Z_3[0] += (control_3dof_B.dv[0] * control_3dof_DW.obj.K_p
    [2] + control_3dof_B.dv[1] * control_3dof_DW.obj.K_v[2]) * 0.02;

  // Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //
  _mm_storeu_pd(&control_3dof_B.dv[0], _mm_sub_pd(_mm_loadu_pd
    (&control_3dof_U.Payload_Out.pL[0]), _mm_set_pd(control_3dof_DW.obj.Z_1[1],
    control_3dof_U.Payload_Out.p_1[0])));

  // MATLAB Function: '<S4>/Force Saturation && Disturbution'
  control_3dof_B.e_1[0] = control_3dof_B.dv[0];

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_B.u0 = control_3dof_U.Payload_Out.vL[1] -
    control_3dof_DW.obj.Z_2[1];
  control_3dof_DW.obj.Z_1[1] += (control_3dof_DW.obj.K_p[0] * control_3dof_B.dv
    [1] + control_3dof_DW.obj.Z_2[1]) * 0.02;
  control_3dof_DW.obj.Z_2[1] += (((control_3dof_DW.UnitDelay_DSTATE[1] +
    control_3dof_DW.obj.Z_3[1]) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.dv[1]) + control_3dof_DW.obj.K_v[1] * control_3dof_B.u0) *
    0.02;
  control_3dof_DW.obj.Z_3[1] += (control_3dof_B.dv[1] * control_3dof_DW.obj.K_p
    [2] + control_3dof_DW.obj.K_v[2] * control_3dof_B.u0) * 0.02;

  // Start for MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //
  _mm_storeu_pd(&control_3dof_B.dv[0], _mm_sub_pd(_mm_loadu_pd
    (&control_3dof_U.Payload_Out.pL[1]), _mm_set_pd(control_3dof_DW.obj.Z_1[2],
    control_3dof_U.Payload_Out.p_1[1])));

  // MATLAB Function: '<S4>/Force Saturation && Disturbution'
  control_3dof_B.e_1[1] = control_3dof_B.dv[0];

  // MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_B.u0 = control_3dof_U.Payload_Out.vL[2] -
    control_3dof_DW.obj.Z_2[2];
  control_3dof_DW.obj.Z_1[2] += (control_3dof_DW.obj.K_p[0] * control_3dof_B.dv
    [1] + control_3dof_DW.obj.Z_2[2]) * 0.02;
  control_3dof_DW.obj.Z_2[2] += ((((control_3dof_DW.UnitDelay_DSTATE[2] +
    control_3dof_DW.obj.Z_3[2]) + 9.8) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.dv[1]) + control_3dof_DW.obj.K_v[1] * control_3dof_B.u0) *
    0.02;
  control_3dof_DW.obj.Z_3[2] += (control_3dof_B.dv[1] * control_3dof_DW.obj.K_p
    [2] + control_3dof_DW.obj.K_v[2] * control_3dof_B.u0) * 0.02;

  // MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //
  control_3dof_B.e_1[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_1[2];
  control_3dof_B.u0 = control_3dof_norm(control_3dof_B.e_1);
  control_3dof_B.e_1[0] /= control_3dof_B.u0;
  control_3dof_B.vu[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_1[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_2[0];
  control_3dof_B.e_1[1] = control_3dof_B.dv[0] / control_3dof_B.u0;
  control_3dof_B.vu[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_1[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_2[1];

  // BusCreator generated from: '<S4>/Force Saturation && Disturbution' incorporates:
  //   Inport: '<Root>/Payload_Out'
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //
  _mm_storeu_pd(&control_3dof_B.dv[0], _mm_div_pd(_mm_sub_pd(_mm_set_pd
    (control_3dof_U.Payload_Out.vL[2], control_3dof_U.Payload_Out.pL[2]),
    _mm_set_pd(control_3dof_U.Payload_Out.v_1[2],
               control_3dof_U.Payload_Out.p_1[2])), _mm_set_pd
    (CONTROL_PARAM.CABLE_LEN, control_3dof_B.u0)));

  // MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'
  //   MATLABSystem: '<S4>/Position 2nd ESO'
  //   Sum: '<S4>/Sum'
  //   UnitDelay: '<S4>/Unit Delay'
  //
  control_3dof_B.e_1[2] = control_3dof_B.dv[0];
  control_3dof_B.q_2[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_2[2];
  _mm_storeu_pd(&control_3dof_B.e_2[0], _mm_sub_pd(_mm_mul_pd(_mm_set_pd
    (control_3dof_B.vu[0], control_3dof_B.e_1[1]), _mm_set_pd(control_3dof_B.dv
    [0], control_3dof_B.dv[1])), _mm_mul_pd(_mm_set_pd(control_3dof_B.e_1[0],
    control_3dof_B.vu[1]), _mm_set_pd(control_3dof_B.dv[1], control_3dof_B.dv[0]))));
  control_3dof_B.e_2[2] = control_3dof_B.e_1[0] * control_3dof_B.vu[1] -
    control_3dof_B.vu[0] * control_3dof_B.e_1[1];
  control_3dof_B.u0 = control_3dof_norm(control_3dof_B.q_2);
  control_3dof_B.q_2[0] /= control_3dof_B.u0;
  control_3dof_B.vu[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_2[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_3[0];
  control_3dof_B.q_2[1] /= control_3dof_B.u0;
  control_3dof_B.vu[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_2[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_3[1];
  control_3dof_B.q_2[2] /= control_3dof_B.u0;
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
  control_3dof_B.u0 = control_3dof_norm(control_3dof_B.q_3);
  control_3dof_B.q_3[0] /= control_3dof_B.u0;
  control_3dof_B.vu[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_3[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[1] /= control_3dof_B.u0;
  control_3dof_B.vu[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_3[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_3[2] /= control_3dof_B.u0;
  control_3dof_B.vu[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_3[2]) / CONTROL_PARAM.CABLE_LEN;
  _mm_storeu_pd(&control_3dof_B.rtb_control_in3_w[0], _mm_sub_pd(_mm_mul_pd
    (_mm_set_pd(control_3dof_B.vu[0], control_3dof_B.q_3[1]), _mm_set_pd
     (control_3dof_B.q_3[2], control_3dof_B.vu[2])), _mm_mul_pd(_mm_set_pd
    (control_3dof_B.q_3[0], control_3dof_B.vu[1]), _mm_set_pd(control_3dof_B.vu
    [2], control_3dof_B.q_3[2]))));
  control_3dof_B.rtb_control_in3_w[2] = control_3dof_B.q_3[0] *
    control_3dof_B.vu[1] - control_3dof_B.vu[0] * control_3dof_B.q_3[1];
  control_3dof_B.Q[0] = control_3dof_B.e_1[0];
  control_3dof_B.Q[3] = control_3dof_B.q_2[0];
  control_3dof_B.Q[6] = control_3dof_B.q_3[0];
  control_3dof_B.u0 = (control_3dof_Y.force_sp1[0] - control_3dof_DW.obj.Z_3[0])
    * CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.F_sp[0] = control_3dof_B.u0;
  control_3dof_B.found[0] = false;
  control_3dof_B.vu_m[0] = 0.0;
  control_3dof_B.F_d[0] = control_3dof_B.u0;
  control_3dof_B.Q[1] = control_3dof_B.e_1[1];
  control_3dof_B.Q[4] = control_3dof_B.q_2[1];
  control_3dof_B.Q[7] = control_3dof_B.q_3[1];
  control_3dof_B.u0 = (control_3dof_Y.force_sp1[1] - control_3dof_DW.obj.Z_3[1])
    * CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.F_sp[1] = control_3dof_B.u0;
  control_3dof_B.found[1] = false;
  control_3dof_B.vu_m[1] = 0.0;
  control_3dof_B.F_d[1] = control_3dof_B.u0;
  control_3dof_B.Q[2] = control_3dof_B.dv[0];
  control_3dof_B.Q[5] = control_3dof_B.q_2[2];
  control_3dof_B.Q[8] = control_3dof_B.q_3[2];
  control_3dof_B.u0 = (control_3dof_Y.force_sp1[2] - control_3dof_DW.obj.Z_3[2])
    * CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.F_sp[2] = control_3dof_B.u0;
  control_3dof_B.found[2] = false;
  control_3dof_B.vu_m[2] = 0.0;
  control_3dof_B.F_d[2] = control_3dof_B.u0;
  std::memcpy(&control_3dof_B.Q_new[0], &control_3dof_B.Q[0], 9U * sizeof(real_T));
  do {
    exitg1 = 0;
    control_3dof_B.is_valid = true;
    for (control_3dof_B.r_p = 0; control_3dof_B.r_p < 9; control_3dof_B.r_p++) {
      control_3dof_B.a[control_3dof_B.r_p] = 0.0;
      if (control_3dof_B.is_valid) {
        control_3dof_B.absx = control_3dof_B.Q_new[control_3dof_B.r_p];
        if (std::isinf(control_3dof_B.absx) || std::isnan(control_3dof_B.absx))
        {
          control_3dof_B.is_valid = false;
        }
      }
    }

    if (!control_3dof_B.is_valid) {
      for (control_3dof_B.r_p = 0; control_3dof_B.r_p < 9; control_3dof_B.r_p++)
      {
        control_3dof_B.a[control_3dof_B.r_p] = (rtNaN);
      }
    } else {
      control_3dof_svd(control_3dof_B.Q_new, control_3dof_B.U, control_3dof_B.vu,
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
      control_3dof_B.r_p = 0;
      exitg2 = false;
      while ((!exitg2) && (control_3dof_B.r_p < 3)) {
        if (std::isinf(control_3dof_B.vu[control_3dof_B.r_p]) || std::isnan
            (control_3dof_B.vu[control_3dof_B.r_p])) {
          control_3dof_B.absx = 1.7976931348623157E+308;
          exitg2 = true;
        } else {
          control_3dof_B.r_p++;
        }
      }

      control_3dof_B.r_p = -1;
      control_3dof_B.k = 0;
      while ((control_3dof_B.k < 3) && (control_3dof_B.vu[control_3dof_B.k] >
              control_3dof_B.absx)) {
        control_3dof_B.r_p++;
        control_3dof_B.k++;
      }

      if (control_3dof_B.r_p + 1 > 0) {
        control_3dof_B.vcol = 1;
        for (control_3dof_B.c_j = 0; control_3dof_B.c_j <= control_3dof_B.r_p;
             control_3dof_B.c_j++) {
          control_3dof_B.absx = 1.0 / control_3dof_B.vu[control_3dof_B.c_j];
          control_3dof_B.k = 2 + control_3dof_B.vcol;
          control_3dof_B.vectorUB = control_3dof_B.vcol;
          for (control_3dof_B.ar = control_3dof_B.vcol; control_3dof_B.ar <=
               control_3dof_B.vectorUB; control_3dof_B.ar += 2) {
            tmp_2 = _mm_loadu_pd(&control_3dof_B.V[control_3dof_B.ar - 1]);
            _mm_storeu_pd(&control_3dof_B.V[control_3dof_B.ar - 1], _mm_mul_pd
                          (tmp_2, _mm_set1_pd(control_3dof_B.absx)));
          }

          for (control_3dof_B.ar = control_3dof_B.k; control_3dof_B.ar <=
               control_3dof_B.vcol + 2; control_3dof_B.ar++) {
            control_3dof_B.V[control_3dof_B.ar - 1] *= control_3dof_B.absx;
          }

          control_3dof_B.vcol += 3;
        }

        control_3dof_B.vcol = 0;
        for (control_3dof_B.k = 0; control_3dof_B.k <= 6; control_3dof_B.k += 3)
        {
          for (control_3dof_B.vectorUB = control_3dof_B.k + 1;
               control_3dof_B.vectorUB <= control_3dof_B.k + 3;
               control_3dof_B.vectorUB++) {
            control_3dof_B.a[control_3dof_B.vectorUB - 1] = 0.0;
          }
        }

        for (control_3dof_B.c_j = 0; control_3dof_B.c_j <= 6; control_3dof_B.c_j
             += 3) {
          control_3dof_B.ar = -1;
          control_3dof_B.vcol++;
          control_3dof_B.b = 3 * control_3dof_B.r_p + control_3dof_B.vcol;
          for (control_3dof_B.ib = control_3dof_B.vcol; control_3dof_B.ib <=
               control_3dof_B.b; control_3dof_B.ib += 3) {
            control_3dof_B.k = control_3dof_B.c_j + 3;
            control_3dof_B.vectorUB = control_3dof_B.c_j + 1;
            for (control_3dof_B.b_ic = control_3dof_B.c_j + 1;
                 control_3dof_B.b_ic <= control_3dof_B.vectorUB;
                 control_3dof_B.b_ic += 2) {
              tmp_2 = _mm_loadu_pd(&control_3dof_B.V[(control_3dof_B.ar +
                control_3dof_B.b_ic) - control_3dof_B.c_j]);
              tmp_1 = _mm_loadu_pd(&control_3dof_B.a[control_3dof_B.b_ic - 1]);
              _mm_storeu_pd(&control_3dof_B.a[control_3dof_B.b_ic - 1],
                            _mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
                (control_3dof_B.U[control_3dof_B.ib - 1])), tmp_1));
            }

            for (control_3dof_B.b_ic = control_3dof_B.k; control_3dof_B.b_ic <=
                 control_3dof_B.c_j + 3; control_3dof_B.b_ic++) {
              control_3dof_B.a[control_3dof_B.b_ic - 1] += control_3dof_B.V
                [(control_3dof_B.ar + control_3dof_B.b_ic) - control_3dof_B.c_j]
                * control_3dof_B.U[control_3dof_B.ib - 1];
            }

            control_3dof_B.ar += 3;
          }
        }
      }
    }

    control_3dof_B.absx = control_3dof_B.F_d[1];
    control_3dof_B.F_d_c = control_3dof_B.F_d[0];
    control_3dof_B.F_d_k = control_3dof_B.F_d[2];
    for (control_3dof_B.r_p = 0; control_3dof_B.r_p <= 0; control_3dof_B.r_p +=
         2) {
      tmp_2 = _mm_loadu_pd(&control_3dof_B.a[control_3dof_B.r_p + 3]);
      tmp_1 = _mm_loadu_pd(&control_3dof_B.a[control_3dof_B.r_p]);
      tmp_0 = _mm_loadu_pd(&control_3dof_B.a[control_3dof_B.r_p + 6]);
      _mm_storeu_pd(&control_3dof_B.vu[control_3dof_B.r_p], _mm_add_pd
                    (_mm_add_pd(_mm_mul_pd(tmp_2, _mm_set1_pd
        (control_3dof_B.absx)), _mm_mul_pd(tmp_1, _mm_set1_pd
        (control_3dof_B.F_d_c))), _mm_mul_pd(tmp_0, _mm_set1_pd
        (control_3dof_B.F_d_k))));
    }

    for (control_3dof_B.r_p = 2; control_3dof_B.r_p < 3; control_3dof_B.r_p++) {
      control_3dof_B.vu[control_3dof_B.r_p] =
        (control_3dof_B.a[control_3dof_B.r_p + 3] * control_3dof_B.absx +
         control_3dof_B.a[control_3dof_B.r_p] * control_3dof_B.F_d_c) +
        control_3dof_B.a[control_3dof_B.r_p + 6] * control_3dof_B.F_d_k;
    }

    control_3dof_B.is_valid = true;
    control_3dof_B.r_p = 0;
    control_3dof_B.k = 0;
    exitg2 = false;
    while ((!exitg2) && (control_3dof_B.k < 3)) {
      control_3dof_B.r_p = control_3dof_B.k;
      if (control_3dof_B.found[control_3dof_B.k]) {
        control_3dof_B.k++;
      } else if (-control_3dof_B.vu[control_3dof_B.k] <
                 CONTROL_PARAM.TENSION_MIN - 1.0E-5) {
        control_3dof_B.vu[control_3dof_B.k] = -CONTROL_PARAM.TENSION_MIN;
        control_3dof_B.is_valid = false;
        exitg2 = true;
      } else if (-control_3dof_B.vu[control_3dof_B.k] >
                 CONTROL_PARAM.TENSION_MAX + 1.0E-5) {
        control_3dof_B.vu[control_3dof_B.k] = -CONTROL_PARAM.TENSION_MAX;
        control_3dof_B.is_valid = false;
        exitg2 = true;
      } else {
        control_3dof_B.k++;
      }
    }

    if (control_3dof_B.is_valid) {
      exitg1 = 1;
    } else {
      control_3dof_B.found[control_3dof_B.r_p] = true;
      control_3dof_B.vu_m[control_3dof_B.r_p] =
        control_3dof_B.vu[control_3dof_B.r_p];
      control_3dof_B.F_d[0] = 0.0;
      control_3dof_B.F_d[1] = 0.0;
      control_3dof_B.F_d[2] = control_3dof_B.u0;
      for (control_3dof_B.k = 0; control_3dof_B.k < 3; control_3dof_B.k++) {
        if (control_3dof_B.found[control_3dof_B.k]) {
          control_3dof_B.absx = control_3dof_B.vu_m[control_3dof_B.k];
          tmp_2 = _mm_sub_pd(_mm_loadu_pd(&control_3dof_B.F_d[0]), _mm_mul_pd
                             (_mm_loadu_pd(&control_3dof_B.Q[3 *
            control_3dof_B.k]), _mm_set1_pd(control_3dof_B.absx)));
          _mm_storeu_pd(&control_3dof_B.F_d[0], tmp_2);
          control_3dof_B.F_d[2] -= control_3dof_B.Q[3 * control_3dof_B.k + 2] *
            control_3dof_B.absx;
        }
      }

      control_3dof_B.Q_new[3 * control_3dof_B.r_p] = -control_3dof_B.F_sp[0];
      control_3dof_B.Q_new[3 * control_3dof_B.r_p + 1] = -control_3dof_B.F_sp[1];
      control_3dof_B.Q_new[3 * control_3dof_B.r_p + 2] = -0.0;
    }
  } while (exitg1 == 0);

  if (!control_3dof_B.found[0]) {
    control_3dof_B.vu_m[0] = control_3dof_B.vu[0];
  }

  if (!control_3dof_B.found[1]) {
    control_3dof_B.vu_m[1] = control_3dof_B.vu[1];
  }

  if (!control_3dof_B.found[2]) {
    control_3dof_B.vu_m[2] = control_3dof_B.vu[2];
  }

  control_3dof_B.absx = control_3dof_B.vu_m[1];
  control_3dof_B.u0 = control_3dof_B.vu_m[0];
  control_3dof_B.F_d_c = control_3dof_B.vu_m[2];
  for (control_3dof_B.r_p = 0; control_3dof_B.r_p <= 0; control_3dof_B.r_p += 2)
  {
    tmp_2 = _mm_loadu_pd(&control_3dof_B.Q[control_3dof_B.r_p + 3]);
    tmp_1 = _mm_set1_pd(control_3dof_B.absx);
    tmp_0 = _mm_loadu_pd(&control_3dof_B.Q[control_3dof_B.r_p]);
    tmp_3 = _mm_set1_pd(control_3dof_B.u0);
    tmp = _mm_loadu_pd(&control_3dof_B.Q[control_3dof_B.r_p + 6]);

    // UnitDelay: '<S4>/Unit Delay'
    tmp_2 = _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_2, tmp_1), _mm_mul_pd(tmp_0,
      tmp_3)), _mm_mul_pd(tmp, _mm_set1_pd(control_3dof_B.F_d_c)));
    tmp_0 = _mm_loadu_pd(&control_3dof_B.F_sp[control_3dof_B.r_p]);
    tmp = _mm_loadu_pd(&control_3dof_Y.force_sp1[control_3dof_B.r_p]);
    tmp_4 = _mm_set1_pd(CONTROL_PARAM.MASS_LOAD);
    _mm_storeu_pd(&control_3dof_B.F_sp[control_3dof_B.r_p], _mm_sub_pd(tmp,
      _mm_div_pd(_mm_sub_pd(tmp_0, tmp_2), tmp_4)));

    // UnitDelay: '<S4>/Unit Delay'
    _mm_storeu_pd(&control_3dof_DW.UnitDelay_DSTATE[control_3dof_B.r_p],
                  _mm_div_pd(tmp_2, tmp_4));
    tmp_2 = _mm_loadu_pd(&control_3dof_B.e_1[control_3dof_B.r_p]);
    _mm_storeu_pd(&control_3dof_B.vu[control_3dof_B.r_p], _mm_mul_pd(tmp_3,
      tmp_2));
    tmp_2 = _mm_loadu_pd(&control_3dof_B.q_2[control_3dof_B.r_p]);
    _mm_storeu_pd(&control_3dof_B.rtb_control_in2_vforce[control_3dof_B.r_p],
                  _mm_mul_pd(tmp_1, tmp_2));
  }

  for (control_3dof_B.r_p = 2; control_3dof_B.r_p < 3; control_3dof_B.r_p++) {
    control_3dof_B.F_d_k = (control_3dof_B.Q[control_3dof_B.r_p + 3] *
      control_3dof_B.absx + control_3dof_B.Q[control_3dof_B.r_p] *
      control_3dof_B.u0) + control_3dof_B.Q[control_3dof_B.r_p + 6] *
      control_3dof_B.F_d_c;
    control_3dof_B.F_sp[control_3dof_B.r_p] =
      control_3dof_Y.force_sp1[control_3dof_B.r_p] -
      (control_3dof_B.F_sp[control_3dof_B.r_p] - control_3dof_B.F_d_k) /
      CONTROL_PARAM.MASS_LOAD;

    // UnitDelay: '<S4>/Unit Delay'
    control_3dof_DW.UnitDelay_DSTATE[control_3dof_B.r_p] = control_3dof_B.F_d_k /
      CONTROL_PARAM.MASS_LOAD;
    _mm_storeu_pd(&control_3dof_B.dv[0], _mm_mul_pd(_mm_set_pd
      (control_3dof_B.absx, control_3dof_B.u0), _mm_set_pd
      (control_3dof_B.q_2[control_3dof_B.r_p],
       control_3dof_B.e_1[control_3dof_B.r_p])));
    control_3dof_B.vu[control_3dof_B.r_p] = control_3dof_B.dv[0];
    control_3dof_B.rtb_control_in2_vforce[control_3dof_B.r_p] =
      control_3dof_B.dv[1];
  }

  control_3dof_B.absx = control_3dof_B.vu_m[2];
  tmp_2 = _mm_mul_pd(_mm_set1_pd(control_3dof_B.absx), _mm_loadu_pd
                     (&control_3dof_B.q_3[0]));
  _mm_storeu_pd(&control_3dof_B.vu_m[0], tmp_2);
  control_3dof_B.vu_m[2] = control_3dof_B.absx * control_3dof_B.q_3[2];

  // MATLAB Function: '<S7>/Vertical Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_VerticalControl(&control_3dof_B.y[0], control_3dof_B.e_1,
    control_3dof_B.e_2, control_3dof_B.F_sp, control_3dof_B.F_d,
    control_3dof_B.state_e, CONTROL_PARAM.KQ, CONTROL_PARAM.KW,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_VerticalControl);

  // MATLAB Function: '<S12>/Vertical Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_VerticalControl(&control_3dof_B.y[3], control_3dof_B.q_2,
    control_3dof_B.rtb_control_in2_w, control_3dof_B.F_sp,
    control_3dof_B.f_parallel, control_3dof_B.state_e, CONTROL_PARAM.KQ,
    CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_VerticalControl_h);

  // MATLAB Function: '<S17>/Vertical Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_VerticalControl(&control_3dof_B.y[6], control_3dof_B.q_3,
    control_3dof_B.rtb_control_in3_w, control_3dof_B.F_sp,
    control_3dof_Y.force_sp1, control_3dof_B.state_e, CONTROL_PARAM.KQ,
    CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_VerticalControl_f);

  // UnitDelay: '<S3>/Unit Delay' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'

  control_3dof_MATLABSystem(control_3dof_DW.UnitDelay_DSTATE_h,
    control_3dof_U.Payload_Out.p_3, control_3dof_U.Payload_Out.v_3,
    &control_3dof_B.MATLABSystem_a, &control_3dof_DW.MATLABSystem_a);

  // MATLAB Function: '<S3>/MATLAB Function' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_MATLABFunction(control_3dof_B.MATLABSystem_a.MATLABSystem_o3,
    control_3dof_B.q_3, control_3dof_B.y_a);

  // MATLAB Function: '<S16>/Parallel Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_ParallelControl(control_3dof_B.q_3,
    control_3dof_B.rtb_control_in3_w, control_3dof_B.vu_m, control_3dof_B.F_sp,
    control_3dof_B.y_b, CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl_a);

  // Gain: '<S3>/Gain1' incorporates:
  //   Gain: '<S1>/Gain1'
  //   Gain: '<S2>/Gain1'

  control_3dof_B.u0 = 1.0 / CONTROL_PARAM.MASS_UAV;

  // Sum: '<S3>/Sum1' incorporates:
  //   Gain: '<S3>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.y_b[0] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_a[0]) + control_3dof_Y.force_sp1[0];

  // Outport: '<Root>/force_sp3'
  control_3dof_Y.force_sp3[0] = control_3dof_B.absx;

  // Gain: '<S3>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   Sum: '<S3>/Sum'
  //   UnitDelay: '<S3>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_h[0] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Sum: '<S3>/Sum1' incorporates:
  //   Gain: '<S3>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.y_b[1] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_a[1]) + control_3dof_Y.force_sp1[1];

  // Outport: '<Root>/force_sp3'
  control_3dof_Y.force_sp3[1] = control_3dof_B.absx;

  // Gain: '<S3>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   Sum: '<S3>/Sum'
  //   UnitDelay: '<S3>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_h[1] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Sum: '<S3>/Sum1' incorporates:
  //   Gain: '<S3>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.y_b[2] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_a[2]) + control_3dof_Y.force_sp1[2];

  // Outport: '<Root>/force_sp3'
  control_3dof_Y.force_sp3[2] = control_3dof_B.absx;

  // Gain: '<S3>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   Sum: '<S3>/Sum'
  //   UnitDelay: '<S3>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_h[2] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // UnitDelay: '<S2>/Unit Delay' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'

  control_3dof_MATLABSystem(control_3dof_DW.UnitDelay_DSTATE_e,
    control_3dof_U.Payload_Out.p_2, control_3dof_U.Payload_Out.v_2,
    &control_3dof_B.MATLABSystem_j, &control_3dof_DW.MATLABSystem_j);

  // MATLAB Function: '<S2>/MATLAB Function' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_MATLABFunction(control_3dof_B.MATLABSystem_j.MATLABSystem_o3,
    control_3dof_B.q_2, control_3dof_B.y_b);

  // MATLAB Function: '<S11>/Parallel Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_ParallelControl(control_3dof_B.q_2,
    control_3dof_B.rtb_control_in2_w, control_3dof_B.rtb_control_in2_vforce,
    control_3dof_B.F_sp, control_3dof_B.y_a, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_h);

  // Sum: '<S2>/Sum1' incorporates:
  //   Gain: '<S2>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.y_a[0] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_b[0]) + control_3dof_B.f_parallel[0];

  // Outport: '<Root>/force_sp2'
  control_3dof_Y.force_sp2[0] = control_3dof_B.absx;

  // Gain: '<S2>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   Sum: '<S2>/Sum'
  //   UnitDelay: '<S2>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_e[0] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Sum: '<S2>/Sum1' incorporates:
  //   Gain: '<S2>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.y_a[1] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_b[1]) + control_3dof_B.f_parallel[1];

  // Outport: '<Root>/force_sp2'
  control_3dof_Y.force_sp2[1] = control_3dof_B.absx;

  // Gain: '<S2>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   Sum: '<S2>/Sum'
  //   UnitDelay: '<S2>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_e[1] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Sum: '<S2>/Sum1' incorporates:
  //   Gain: '<S2>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.y_a[2] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_b[2]) + control_3dof_B.f_parallel[2];

  // Outport: '<Root>/force_sp2'
  control_3dof_Y.force_sp2[2] = control_3dof_B.absx;

  // Gain: '<S2>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   Sum: '<S2>/Sum'
  //   UnitDelay: '<S2>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_e[2] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // UnitDelay: '<S1>/Unit Delay' incorporates:
  //   BusCreator generated from: '<S4>/Force Saturation && Disturbution'
  //   Inport: '<Root>/Payload_Out'

  control_3dof_MATLABSystem(control_3dof_DW.UnitDelay_DSTATE_f,
    control_3dof_U.Payload_Out.p_1, control_3dof_U.Payload_Out.v_1,
    &control_3dof_B.MATLABSystem, &control_3dof_DW.MATLABSystem);

  // MATLAB Function: '<S1>/MATLAB Function' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_MATLABFunction(control_3dof_B.MATLABSystem.MATLABSystem_o3,
    control_3dof_B.e_1, control_3dof_B.y_a);

  // MATLAB Function: '<S6>/Parallel Control' incorporates:
  //   MATLAB Function: '<S4>/Force Saturation && Disturbution'

  control_3dof_ParallelControl(control_3dof_B.e_1, control_3dof_B.e_2,
    control_3dof_B.vu, control_3dof_B.F_sp, control_3dof_B.f_parallel,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl);

  // Sum: '<S1>/Sum1' incorporates:
  //   Gain: '<S1>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.f_parallel[0] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_a[0]) + control_3dof_B.F_d[0];
  control_3dof_Y.force_sp1[0] = control_3dof_B.absx;

  // Gain: '<S1>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   UnitDelay: '<S1>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_f[0] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Sum: '<S1>/Sum1' incorporates:
  //   Gain: '<S1>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.f_parallel[1] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_a[1]) + control_3dof_B.F_d[1];
  control_3dof_Y.force_sp1[1] = control_3dof_B.absx;

  // Gain: '<S1>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   UnitDelay: '<S1>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_f[1] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Sum: '<S1>/Sum1' incorporates:
  //   Gain: '<S1>/Gain'
  //   Sum: '<S1>/Sum'

  control_3dof_B.absx = (control_3dof_B.f_parallel[2] - CONTROL_PARAM.MASS_UAV *
    control_3dof_B.y_a[2]) + control_3dof_B.F_d[2];
  control_3dof_Y.force_sp1[2] = control_3dof_B.absx;

  // Gain: '<S1>/Gain1' incorporates:
  //   Sum: '<S1>/Sum'
  //   UnitDelay: '<S1>/Unit Delay'

  control_3dof_DW.UnitDelay_DSTATE_f[2] = control_3dof_B.u0 *
    control_3dof_B.absx;

  // Update for DiscreteIntegrator: '<S20>/Discrete-Time Integrator' incorporates:
  //   Constant: '<S20>/Constant'

  control_3dof_DW.DiscreteTimeIntegrator_DSTATE += 0.02;
}

// Model initialize function
void control_3dof::initialize()
{
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
  control_3dof_MATLABSystem_Init(&control_3dof_DW.MATLABSystem_a);
  control_3dof_MATLABSystem_Init(&control_3dof_DW.MATLABSystem_j);
  control_3dof_MATLABSystem_Init(&control_3dof_DW.MATLABSystem);
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
