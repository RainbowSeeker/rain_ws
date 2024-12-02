//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: payload_controller.cpp
//
// Code generated for Simulink model 'payload_controller'.
//
// Model version                  : 1.266
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Wed Nov 27 19:34:30 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "payload_controller.h"
#include "rtwtypes.h"
#include <cmath>
#include <emmintrin.h>
#include "cmath"
#include "limits"

// Exported block parameters
struct_yHOZp7GlS2mevFwObUORyE CONTROL_PARAM{
  1.0,
  1.5,
  1.0,
  0.2,
  4.0,
  1.0,
  5.0
} ;                                    // Variable: CONTROL_PARAM
                                          //  Referenced by:
                                          //    '<S7>/Kp'
                                          //    '<S11>/Parallel Control'
                                          //    '<S12>/Vertical Control'
                                          //    '<S15>/Parallel Control'
                                          //    '<S16>/Vertical Control'
                                          //    '<S19>/Parallel Control'
                                          //    '<S20>/Vertical Control'
                                          //    '<S23>/mL'
                                          //    '<S24>/Kv'


// Invariant block signals (default storage)
const payload_controller::ConstB_payload_controller_T payload_controller_ConstB{
  {
    0.0,
    0.0,
    -9.8
  }
  // '<S24>/-g'
};

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
//    '<S1>/q1'
//    '<S1>/q2'
//    '<S1>/q3'
//
void payload_controller::payload_controller_q1(real_T rtu_theta, real_T rtu_psi,
  real_T rty_q[3])
{
  real_T tmp;
  tmp = std::cos(rtu_theta);
  rty_q[0] = -std::cos(rtu_psi) * tmp;
  rty_q[1] = -std::sin(rtu_psi) * tmp;
  rty_q[2] = std::sin(rtu_theta);
}

//
// Output and update for atomic system:
//    '<S11>/Parallel Control'
//    '<S15>/Parallel Control'
//    '<S19>/Parallel Control'
//
void payload_controller::payload_control_ParallelControl(const real_T rtu_q[3],
  const real_T rtu_w[3], real_T rtu_vforce, const real_T rtu_aL[3], real_T
  rty_y[3], real_T rtp_li, real_T rtp_mi, B_ParallelControl_payload_con_T
  *localB)
{
  __m128d tmp_2;
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
    tmp_2 = _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(rtp_mi), _mm_loadu_pd(&rtu_q[0])),
                       _mm_set1_pd(rtu_q[localB->i]));
    _mm_storeu_pd(&localB->rtp_mi[3 * localB->i], tmp_2);
    localB->rtp_mi[3 * localB->i + 2] = rtp_mi * rtu_q[2] * rtu_q[localB->i];
  }

  localB->absxk = rtu_aL[1];
  localB->t = rtu_aL[0];
  localB->a = rtu_aL[2];
  for (localB->i = 0; localB->i <= 0; localB->i += 2) {
    __m128d tmp;
    __m128d tmp_0;
    __m128d tmp_1;
    tmp_2 = _mm_loadu_pd(&localB->rtp_mi[localB->i + 3]);
    tmp = _mm_loadu_pd(&localB->rtp_mi[localB->i]);
    tmp_0 = _mm_loadu_pd(&localB->rtp_mi[localB->i + 6]);
    tmp_1 = _mm_loadu_pd(&rtu_q[localB->i]);
    _mm_storeu_pd(&rty_y[localB->i], _mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd
      (tmp_2, _mm_set1_pd(localB->absxk)), _mm_mul_pd(tmp, _mm_set1_pd(localB->t))),
      _mm_mul_pd(tmp_0, _mm_set1_pd(localB->a))), _mm_add_pd(_mm_mul_pd
      (_mm_set1_pd(rtu_vforce), tmp_1), _mm_mul_pd(_mm_set1_pd(localB->scale),
      tmp_1))));
  }

  for (localB->i = 2; localB->i < 3; localB->i++) {
    localB->rtu_q = rtu_q[localB->i];
    rty_y[localB->i] = ((localB->rtp_mi[localB->i + 3] * localB->absxk +
                         localB->rtp_mi[localB->i] * localB->t) + localB->
                        rtp_mi[localB->i + 6] * localB->a) + (rtu_vforce *
      localB->rtu_q + localB->scale * localB->rtu_q);
  }
}

//
// Output and update for atomic system:
//    '<S12>/Vertical Control'
//    '<S16>/Vertical Control'
//    '<S20>/Vertical Control'
//
void payload_controller::payload_control_VerticalControl(const real_T rtu_q_sp[3],
  const real_T rtu_q[3], const real_T rtu_w[3], const real_T rtu_aL_sp[3],
  real_T rty_y[3], real_T rty_state[6], real_T rtp_Kq, real_T rtp_Kw, real_T
  rtp_li, real_T rtp_mi, B_VerticalControl_payload_con_T *localB)
{
  __m128d tmp;
  real_T a;
  real_T e_w;
  real_T e_w_tmp;
  localB->Sqi[0] = 0.0;
  localB->Sqi[3] = -rtu_q[2];
  localB->Sqi[6] = rtu_q[1];
  localB->Sqi[1] = rtu_q[2];
  localB->Sqi[4] = 0.0;
  localB->Sqi[7] = -rtu_q[0];
  localB->Sqi[2] = -rtu_q[1];
  localB->Sqi[5] = rtu_q[0];
  localB->Sqi[8] = 0.0;
  tmp = _mm_sub_pd(_mm_mul_pd(_mm_set_pd(rtu_q[0], rtu_q_sp[1]), _mm_set_pd
    (rtu_q_sp[2], rtu_q[2])), _mm_mul_pd(_mm_set_pd(rtu_q_sp[0], rtu_q[1]),
    _mm_set_pd(rtu_q[2], rtu_q_sp[2])));
  _mm_storeu_pd(&localB->e_q[0], tmp);
  localB->e_q[2] = rtu_q_sp[0] * rtu_q[1] - rtu_q[0] * rtu_q_sp[1];
  a = rtp_mi * rtp_li;
  for (int32_T i{0}; i < 3; i++) {
    e_w_tmp = 0.0;
    for (int32_T i_0{0}; i_0 < 3; i_0++) {
      e_w_tmp += ((localB->Sqi[3 * i_0 + 1] * localB->Sqi[i + 3] + localB->Sqi[3
                   * i_0] * localB->Sqi[i]) + localB->Sqi[3 * i_0 + 2] *
                  localB->Sqi[i + 6]) * (-rtp_Kq * localB->e_q[i_0]);
    }

    e_w = rtu_w[i] + e_w_tmp;
    localB->e_w[i] = e_w;
    localB->rtp_Kw[i] = -rtp_Kw * e_w - e_w_tmp;
  }

  for (int32_T i{0}; i < 3; i++) {
    e_w_tmp = 0.0;
    e_w = 0.0;
    for (int32_T i_0{0}; i_0 < 3; i_0++) {
      e_w_tmp += localB->Sqi[3 * i_0 + i] * a * localB->rtp_Kw[i_0];
      e_w += ((localB->Sqi[i + 3] * rtp_mi * localB->Sqi[3 * i_0 + 1] + rtp_mi *
               localB->Sqi[i] * localB->Sqi[3 * i_0]) + localB->Sqi[i + 6] *
              rtp_mi * localB->Sqi[3 * i_0 + 2]) * rtu_aL_sp[i_0];
    }

    rty_y[i] = e_w_tmp - e_w;
    rty_state[i] = localB->e_q[i];
    rty_state[i + 3] = localB->e_w[i];
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
real_T payload_controller::payload_controller_xzlangeM(const real_T x[9])
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

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xzlascl(real_T cfrom, real_T cto,
  int32_T m, int32_T n, real_T A[9], int32_T iA0, int32_T lda)
{
  payload_controller_B.cfromc_c = cfrom;
  payload_controller_B.ctoc_k = cto;
  payload_controller_B.notdone_l = true;
  while (payload_controller_B.notdone_l) {
    payload_controller_B.cfrom1_c = payload_controller_B.cfromc_c *
      2.0041683600089728E-292;
    payload_controller_B.cto1_b = payload_controller_B.ctoc_k /
      4.9896007738368E+291;
    if ((std::abs(payload_controller_B.cfrom1_c) > std::abs
         (payload_controller_B.ctoc_k)) && (payload_controller_B.ctoc_k != 0.0))
    {
      payload_controller_B.mul_p = 2.0041683600089728E-292;
      payload_controller_B.cfromc_c = payload_controller_B.cfrom1_c;
    } else if (std::abs(payload_controller_B.cto1_b) > std::abs
               (payload_controller_B.cfromc_c)) {
      payload_controller_B.mul_p = 4.9896007738368E+291;
      payload_controller_B.ctoc_k = payload_controller_B.cto1_b;
    } else {
      payload_controller_B.mul_p = payload_controller_B.ctoc_k /
        payload_controller_B.cfromc_c;
      payload_controller_B.notdone_l = false;
    }

    for (payload_controller_B.j_l = 0; payload_controller_B.j_l < n;
         payload_controller_B.j_l++) {
      payload_controller_B.offset_p = (payload_controller_B.j_l * lda + iA0) - 2;
      payload_controller_B.scalarLB_d = (m / 2) << 1;
      payload_controller_B.vectorUB_g = payload_controller_B.scalarLB_d - 2;
      for (payload_controller_B.b_i_j = 0; payload_controller_B.b_i_j <=
           payload_controller_B.vectorUB_g; payload_controller_B.b_i_j += 2) {
        __m128d tmp;
        payload_controller_B.i1 = (payload_controller_B.b_i_j +
          payload_controller_B.offset_p) + 1;
        tmp = _mm_loadu_pd(&A[payload_controller_B.i1]);
        _mm_storeu_pd(&A[payload_controller_B.i1], _mm_mul_pd(tmp, _mm_set1_pd
          (payload_controller_B.mul_p)));
      }

      for (payload_controller_B.b_i_j = payload_controller_B.scalarLB_d;
           payload_controller_B.b_i_j < m; payload_controller_B.b_i_j++) {
        payload_controller_B.i1 = (payload_controller_B.b_i_j +
          payload_controller_B.offset_p) + 1;
        A[payload_controller_B.i1] *= payload_controller_B.mul_p;
      }
    }
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
real_T payload_controller::payload_controller_xnrm2(int32_T n, const real_T x[9],
  int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      payload_controller_B.scale_f = 3.3121686421112381E-170;
      payload_controller_B.kend_n = ix0 + n;
      for (payload_controller_B.k_bs = ix0; payload_controller_B.k_bs <
           payload_controller_B.kend_n; payload_controller_B.k_bs++) {
        payload_controller_B.absxk_g = std::abs(x[payload_controller_B.k_bs - 1]);
        if (payload_controller_B.absxk_g > payload_controller_B.scale_f) {
          payload_controller_B.t_g = payload_controller_B.scale_f /
            payload_controller_B.absxk_g;
          y = y * payload_controller_B.t_g * payload_controller_B.t_g + 1.0;
          payload_controller_B.scale_f = payload_controller_B.absxk_g;
        } else {
          payload_controller_B.t_g = payload_controller_B.absxk_g /
            payload_controller_B.scale_f;
          y += payload_controller_B.t_g * payload_controller_B.t_g;
        }
      }

      y = payload_controller_B.scale_f * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
real_T payload_controller::payload_controller_xdotc(int32_T n, const real_T x[9],
  int32_T ix0, const real_T y[9], int32_T iy0)
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

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xaxpy(int32_T n, real_T a, int32_T
  ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (int32_T k{0}; k < n; k++) {
      int32_T tmp;
      tmp = (iy0 + k) - 1;
      y[tmp] += y[(ix0 + k) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
real_T payload_controller::payload_controller_xnrm2_d(int32_T n, const real_T x
  [3], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = std::abs(x[ix0 - 1]);
    } else {
      payload_controller_B.scale_c = 3.3121686421112381E-170;
      payload_controller_B.kend = ix0 + n;
      for (payload_controller_B.k_b = ix0; payload_controller_B.k_b <
           payload_controller_B.kend; payload_controller_B.k_b++) {
        payload_controller_B.absxk = std::abs(x[payload_controller_B.k_b - 1]);
        if (payload_controller_B.absxk > payload_controller_B.scale_c) {
          payload_controller_B.t = payload_controller_B.scale_c /
            payload_controller_B.absxk;
          y = y * payload_controller_B.t * payload_controller_B.t + 1.0;
          payload_controller_B.scale_c = payload_controller_B.absxk;
        } else {
          payload_controller_B.t = payload_controller_B.absxk /
            payload_controller_B.scale_c;
          y += payload_controller_B.t * payload_controller_B.t;
        }
      }

      y = payload_controller_B.scale_c * std::sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xaxpy_g(int32_T n, real_T a, const
  real_T x[9], int32_T ix0, real_T y[3], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    payload_controller_B.scalarLB_lx = (n / 2) << 1;
    payload_controller_B.vectorUB_o = payload_controller_B.scalarLB_lx - 2;
    for (payload_controller_B.k_d = 0; payload_controller_B.k_d <=
         payload_controller_B.vectorUB_o; payload_controller_B.k_d += 2) {
      __m128d tmp;
      payload_controller_B.i3 = (iy0 + payload_controller_B.k_d) - 1;
      tmp = _mm_loadu_pd(&y[payload_controller_B.i3]);
      _mm_storeu_pd(&y[payload_controller_B.i3], _mm_add_pd(_mm_mul_pd
        (_mm_loadu_pd(&x[(ix0 + payload_controller_B.k_d) - 1]), _mm_set1_pd(a)),
        tmp));
    }

    for (payload_controller_B.k_d = payload_controller_B.scalarLB_lx;
         payload_controller_B.k_d < n; payload_controller_B.k_d++) {
      payload_controller_B.i3 = (iy0 + payload_controller_B.k_d) - 1;
      y[payload_controller_B.i3] += x[(ix0 + payload_controller_B.k_d) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xaxpy_g3(int32_T n, real_T a, const
  real_T x[3], int32_T ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    payload_controller_B.scalarLB_l = (n / 2) << 1;
    payload_controller_B.vectorUB_d = payload_controller_B.scalarLB_l - 2;
    for (payload_controller_B.k = 0; payload_controller_B.k <=
         payload_controller_B.vectorUB_d; payload_controller_B.k += 2) {
      __m128d tmp;
      payload_controller_B.i2 = (iy0 + payload_controller_B.k) - 1;
      tmp = _mm_loadu_pd(&y[payload_controller_B.i2]);
      _mm_storeu_pd(&y[payload_controller_B.i2], _mm_add_pd(_mm_mul_pd
        (_mm_loadu_pd(&x[(ix0 + payload_controller_B.k) - 1]), _mm_set1_pd(a)),
        tmp));
    }

    for (payload_controller_B.k = payload_controller_B.scalarLB_l;
         payload_controller_B.k < n; payload_controller_B.k++) {
      payload_controller_B.i2 = (iy0 + payload_controller_B.k) - 1;
      y[payload_controller_B.i2] += x[(ix0 + payload_controller_B.k) - 1] * a;
    }
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xzlascl_a(real_T cfrom, real_T cto,
  int32_T m, int32_T n, real_T A[3], int32_T iA0, int32_T lda)
{
  payload_controller_B.cfromc = cfrom;
  payload_controller_B.ctoc = cto;
  payload_controller_B.notdone = true;
  while (payload_controller_B.notdone) {
    payload_controller_B.cfrom1 = payload_controller_B.cfromc *
      2.0041683600089728E-292;
    payload_controller_B.cto1 = payload_controller_B.ctoc / 4.9896007738368E+291;
    if ((std::abs(payload_controller_B.cfrom1) > std::abs
         (payload_controller_B.ctoc)) && (payload_controller_B.ctoc != 0.0)) {
      payload_controller_B.mul = 2.0041683600089728E-292;
      payload_controller_B.cfromc = payload_controller_B.cfrom1;
    } else if (std::abs(payload_controller_B.cto1) > std::abs
               (payload_controller_B.cfromc)) {
      payload_controller_B.mul = 4.9896007738368E+291;
      payload_controller_B.ctoc = payload_controller_B.cto1;
    } else {
      payload_controller_B.mul = payload_controller_B.ctoc /
        payload_controller_B.cfromc;
      payload_controller_B.notdone = false;
    }

    for (payload_controller_B.j_m = 0; payload_controller_B.j_m < n;
         payload_controller_B.j_m++) {
      payload_controller_B.offset = (payload_controller_B.j_m * lda + iA0) - 2;
      payload_controller_B.scalarLB = (m / 2) << 1;
      payload_controller_B.vectorUB_n = payload_controller_B.scalarLB - 2;
      for (payload_controller_B.b_i = 0; payload_controller_B.b_i <=
           payload_controller_B.vectorUB_n; payload_controller_B.b_i += 2) {
        __m128d tmp;
        payload_controller_B.i = (payload_controller_B.b_i +
          payload_controller_B.offset) + 1;
        tmp = _mm_loadu_pd(&A[payload_controller_B.i]);
        _mm_storeu_pd(&A[payload_controller_B.i], _mm_mul_pd(tmp, _mm_set1_pd
          (payload_controller_B.mul)));
      }

      for (payload_controller_B.b_i = payload_controller_B.scalarLB;
           payload_controller_B.b_i < m; payload_controller_B.b_i++) {
        payload_controller_B.i = (payload_controller_B.b_i +
          payload_controller_B.offset) + 1;
        A[payload_controller_B.i] *= payload_controller_B.mul;
      }
    }
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xswap(real_T x[9], int32_T ix0,
  int32_T iy0)
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

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xrotg(real_T *a, real_T *b, real_T
  *c, real_T *s)
{
  payload_controller_B.roe = *b;
  payload_controller_B.absa = std::abs(*a);
  payload_controller_B.absb = std::abs(*b);
  if (payload_controller_B.absa > payload_controller_B.absb) {
    payload_controller_B.roe = *a;
  }

  payload_controller_B.scale = payload_controller_B.absa +
    payload_controller_B.absb;
  if (payload_controller_B.scale == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    payload_controller_B.ads = payload_controller_B.absa /
      payload_controller_B.scale;
    payload_controller_B.bds = payload_controller_B.absb /
      payload_controller_B.scale;
    payload_controller_B.scale *= std::sqrt(payload_controller_B.ads *
      payload_controller_B.ads + payload_controller_B.bds *
      payload_controller_B.bds);
    if (payload_controller_B.roe < 0.0) {
      payload_controller_B.scale = -payload_controller_B.scale;
    }

    *c = *a / payload_controller_B.scale;
    *s = *b / payload_controller_B.scale;
    if (payload_controller_B.absa > payload_controller_B.absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = payload_controller_B.scale;
  }
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_xrot(real_T x[9], int32_T ix0,
  int32_T iy0, real_T c, real_T s)
{
  payload_controller_B.temp = x[iy0 - 1];
  payload_controller_B.temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = payload_controller_B.temp * c - payload_controller_B.temp_tmp * s;
  x[ix0 - 1] = payload_controller_B.temp_tmp * c + payload_controller_B.temp * s;
  payload_controller_B.temp = x[ix0] * c + x[iy0] * s;
  x[iy0] = x[iy0] * c - x[ix0] * s;
  x[ix0] = payload_controller_B.temp;
  payload_controller_B.temp = x[iy0 + 1];
  payload_controller_B.temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = payload_controller_B.temp * c - payload_controller_B.temp_tmp * s;
  x[ix0 + 1] = payload_controller_B.temp_tmp * c + payload_controller_B.temp * s;
}

// Function for MATLAB Function: '<S23>/Moore–Penrose inverse'
void payload_controller::payload_controller_svd(const real_T A[9], real_T U[9],
  real_T s[3], real_T V[9])
{
  __m128d tmp;
  boolean_T exitg1;
  payload_controller_B.b_s[0] = 0.0;
  payload_controller_B.e[0] = 0.0;
  payload_controller_B.work[0] = 0.0;
  payload_controller_B.b_s[1] = 0.0;
  payload_controller_B.e[1] = 0.0;
  payload_controller_B.work[1] = 0.0;
  payload_controller_B.b_s[2] = 0.0;
  payload_controller_B.e[2] = 0.0;
  payload_controller_B.work[2] = 0.0;
  for (payload_controller_B.qjj = 0; payload_controller_B.qjj < 9;
       payload_controller_B.qjj++) {
    payload_controller_B.b_A[payload_controller_B.qjj] =
      A[payload_controller_B.qjj];
    U[payload_controller_B.qjj] = 0.0;
    payload_controller_B.Vf[payload_controller_B.qjj] = 0.0;
  }

  payload_controller_B.doscale = false;
  payload_controller_B.anrm = payload_controller_xzlangeM(A);
  payload_controller_B.cscale = payload_controller_B.anrm;
  if ((payload_controller_B.anrm > 0.0) && (payload_controller_B.anrm <
       6.7178761075670888E-139)) {
    payload_controller_B.doscale = true;
    payload_controller_B.cscale = 6.7178761075670888E-139;
    payload_controller_xzlascl(payload_controller_B.anrm,
      payload_controller_B.cscale, 3, 3, payload_controller_B.b_A, 1, 3);
  } else if (payload_controller_B.anrm > 1.4885657073574029E+138) {
    payload_controller_B.doscale = true;
    payload_controller_B.cscale = 1.4885657073574029E+138;
    payload_controller_xzlascl(payload_controller_B.anrm,
      payload_controller_B.cscale, 3, 3, payload_controller_B.b_A, 1, 3);
  }

  for (payload_controller_B.m = 0; payload_controller_B.m < 2;
       payload_controller_B.m++) {
    payload_controller_B.qp1 = payload_controller_B.m + 2;
    payload_controller_B.qq_tmp = 3 * payload_controller_B.m +
      payload_controller_B.m;
    payload_controller_B.qq = payload_controller_B.qq_tmp + 1;
    payload_controller_B.apply_transform = false;
    payload_controller_B.nrm = payload_controller_xnrm2(3 -
      payload_controller_B.m, payload_controller_B.b_A,
      payload_controller_B.qq_tmp + 1);
    if (payload_controller_B.nrm > 0.0) {
      payload_controller_B.apply_transform = true;
      if (payload_controller_B.b_A[payload_controller_B.qq_tmp] < 0.0) {
        payload_controller_B.nrm = -payload_controller_B.nrm;
      }

      payload_controller_B.b_s[payload_controller_B.m] =
        payload_controller_B.nrm;
      if (std::abs(payload_controller_B.nrm) >= 1.0020841800044864E-292) {
        payload_controller_B.nrm = 1.0 / payload_controller_B.nrm;
        payload_controller_B.d = (payload_controller_B.qq_tmp -
          payload_controller_B.m) + 3;
        payload_controller_B.qjj = ((((payload_controller_B.d -
          payload_controller_B.qq_tmp) / 2) << 1) + payload_controller_B.qq_tmp)
          + 1;
        payload_controller_B.kase = payload_controller_B.qjj - 2;
        for (payload_controller_B.e_k = payload_controller_B.qq;
             payload_controller_B.e_k <= payload_controller_B.kase;
             payload_controller_B.e_k += 2) {
          tmp = _mm_loadu_pd(&payload_controller_B.b_A[payload_controller_B.e_k
                             - 1]);
          _mm_storeu_pd(&payload_controller_B.b_A[payload_controller_B.e_k - 1],
                        _mm_mul_pd(tmp, _mm_set1_pd(payload_controller_B.nrm)));
        }

        for (payload_controller_B.e_k = payload_controller_B.qjj;
             payload_controller_B.e_k <= payload_controller_B.d;
             payload_controller_B.e_k++) {
          payload_controller_B.b_A[payload_controller_B.e_k - 1] *=
            payload_controller_B.nrm;
        }
      } else {
        payload_controller_B.d = (payload_controller_B.qq_tmp -
          payload_controller_B.m) + 3;
        payload_controller_B.qjj = ((((payload_controller_B.d -
          payload_controller_B.qq_tmp) / 2) << 1) + payload_controller_B.qq_tmp)
          + 1;
        payload_controller_B.kase = payload_controller_B.qjj - 2;
        for (payload_controller_B.e_k = payload_controller_B.qq;
             payload_controller_B.e_k <= payload_controller_B.kase;
             payload_controller_B.e_k += 2) {
          tmp = _mm_loadu_pd(&payload_controller_B.b_A[payload_controller_B.e_k
                             - 1]);
          _mm_storeu_pd(&payload_controller_B.b_A[payload_controller_B.e_k - 1],
                        _mm_div_pd(tmp, _mm_set1_pd
            (payload_controller_B.b_s[payload_controller_B.m])));
        }

        for (payload_controller_B.e_k = payload_controller_B.qjj;
             payload_controller_B.e_k <= payload_controller_B.d;
             payload_controller_B.e_k++) {
          payload_controller_B.b_A[payload_controller_B.e_k - 1] /=
            payload_controller_B.b_s[payload_controller_B.m];
        }
      }

      payload_controller_B.b_A[payload_controller_B.qq_tmp]++;
      payload_controller_B.b_s[payload_controller_B.m] =
        -payload_controller_B.b_s[payload_controller_B.m];
    } else {
      payload_controller_B.b_s[payload_controller_B.m] = 0.0;
    }

    for (payload_controller_B.kase = payload_controller_B.qp1;
         payload_controller_B.kase < 4; payload_controller_B.kase++) {
      payload_controller_B.qjj = (payload_controller_B.kase - 1) * 3 +
        payload_controller_B.m;
      if (payload_controller_B.apply_transform) {
        payload_controller_xaxpy(3 - payload_controller_B.m,
          -(payload_controller_xdotc(3 - payload_controller_B.m,
          payload_controller_B.b_A, payload_controller_B.qq_tmp + 1,
          payload_controller_B.b_A, payload_controller_B.qjj + 1) /
            payload_controller_B.b_A[payload_controller_B.qq_tmp]),
          payload_controller_B.qq_tmp + 1, payload_controller_B.b_A,
          payload_controller_B.qjj + 1);
      }

      payload_controller_B.e[payload_controller_B.kase - 1] =
        payload_controller_B.b_A[payload_controller_B.qjj];
    }

    for (payload_controller_B.kase = payload_controller_B.m + 1;
         payload_controller_B.kase < 4; payload_controller_B.kase++) {
      payload_controller_B.qjj = (3 * payload_controller_B.m +
        payload_controller_B.kase) - 1;
      U[payload_controller_B.qjj] =
        payload_controller_B.b_A[payload_controller_B.qjj];
    }

    if (payload_controller_B.m + 1 <= 1) {
      payload_controller_B.nrm = payload_controller_xnrm2_d(2,
        payload_controller_B.e, 2);
      if (payload_controller_B.nrm == 0.0) {
        payload_controller_B.e[0] = 0.0;
      } else {
        if (payload_controller_B.e[1] < 0.0) {
          payload_controller_B.e[0] = -payload_controller_B.nrm;
        } else {
          payload_controller_B.e[0] = payload_controller_B.nrm;
        }

        payload_controller_B.nrm = payload_controller_B.e[0];
        if (std::abs(payload_controller_B.e[0]) >= 1.0020841800044864E-292) {
          payload_controller_B.nrm = 1.0 / payload_controller_B.e[0];
          payload_controller_B.qjj = ((((2 - payload_controller_B.m) / 2) << 1)
            + payload_controller_B.m) + 2;
          payload_controller_B.kase = payload_controller_B.qjj - 2;
          for (payload_controller_B.qq = payload_controller_B.qp1;
               payload_controller_B.qq <= payload_controller_B.kase;
               payload_controller_B.qq += 2) {
            tmp = _mm_loadu_pd(&payload_controller_B.e[payload_controller_B.qq -
                               1]);
            _mm_storeu_pd(&payload_controller_B.e[payload_controller_B.qq - 1],
                          _mm_mul_pd(tmp, _mm_set1_pd(payload_controller_B.nrm)));
          }

          for (payload_controller_B.qq = payload_controller_B.qjj;
               payload_controller_B.qq < 4; payload_controller_B.qq++) {
            payload_controller_B.e[payload_controller_B.qq - 1] *=
              payload_controller_B.nrm;
          }
        } else {
          payload_controller_B.qjj = ((((2 - payload_controller_B.m) / 2) << 1)
            + payload_controller_B.m) + 2;
          payload_controller_B.kase = payload_controller_B.qjj - 2;
          for (payload_controller_B.qq = payload_controller_B.qp1;
               payload_controller_B.qq <= payload_controller_B.kase;
               payload_controller_B.qq += 2) {
            tmp = _mm_loadu_pd(&payload_controller_B.e[payload_controller_B.qq -
                               1]);
            _mm_storeu_pd(&payload_controller_B.e[payload_controller_B.qq - 1],
                          _mm_div_pd(tmp, _mm_set1_pd(payload_controller_B.nrm)));
          }

          for (payload_controller_B.qq = payload_controller_B.qjj;
               payload_controller_B.qq < 4; payload_controller_B.qq++) {
            payload_controller_B.e[payload_controller_B.qq - 1] /=
              payload_controller_B.nrm;
          }
        }

        payload_controller_B.e[1]++;
        payload_controller_B.e[0] = -payload_controller_B.e[0];
        for (payload_controller_B.qjj = payload_controller_B.qp1;
             payload_controller_B.qjj < 4; payload_controller_B.qjj++) {
          payload_controller_B.work[payload_controller_B.qjj - 1] = 0.0;
        }

        for (payload_controller_B.qjj = payload_controller_B.qp1;
             payload_controller_B.qjj < 4; payload_controller_B.qjj++) {
          payload_controller_xaxpy_g(2,
            payload_controller_B.e[payload_controller_B.qjj - 1],
            payload_controller_B.b_A, 3 * (payload_controller_B.qjj - 1) + 2,
            payload_controller_B.work, 2);
        }

        for (payload_controller_B.qjj = payload_controller_B.qp1;
             payload_controller_B.qjj < 4; payload_controller_B.qjj++) {
          payload_controller_xaxpy_g3(2,
            -payload_controller_B.e[payload_controller_B.qjj - 1] /
            payload_controller_B.e[1], payload_controller_B.work, 2,
            payload_controller_B.b_A, 3 * (payload_controller_B.qjj - 1) + 2);
        }
      }

      for (payload_controller_B.qjj = payload_controller_B.qp1;
           payload_controller_B.qjj < 4; payload_controller_B.qjj++) {
        payload_controller_B.Vf[payload_controller_B.qjj - 1] =
          payload_controller_B.e[payload_controller_B.qjj - 1];
      }
    }
  }

  payload_controller_B.m = 1;
  payload_controller_B.b_s[2] = payload_controller_B.b_A[8];
  payload_controller_B.e[1] = payload_controller_B.b_A[7];
  payload_controller_B.e[2] = 0.0;
  U[6] = 0.0;
  U[7] = 0.0;
  U[8] = 1.0;
  for (payload_controller_B.qp1 = 1; payload_controller_B.qp1 >= 0;
       payload_controller_B.qp1--) {
    payload_controller_B.qq = 3 * payload_controller_B.qp1 +
      payload_controller_B.qp1;
    if (payload_controller_B.b_s[payload_controller_B.qp1] != 0.0) {
      for (payload_controller_B.kase = payload_controller_B.qp1 + 2;
           payload_controller_B.kase < 4; payload_controller_B.kase++) {
        payload_controller_B.qjj = ((payload_controller_B.kase - 1) * 3 +
          payload_controller_B.qp1) + 1;
        payload_controller_xaxpy(3 - payload_controller_B.qp1,
          -(payload_controller_xdotc(3 - payload_controller_B.qp1, U,
          payload_controller_B.qq + 1, U, payload_controller_B.qjj) /
            U[payload_controller_B.qq]), payload_controller_B.qq + 1, U,
          payload_controller_B.qjj);
      }

      for (payload_controller_B.kase = payload_controller_B.qp1 + 1;
           payload_controller_B.kase < 4; payload_controller_B.kase++) {
        payload_controller_B.qjj = (3 * payload_controller_B.qp1 +
          payload_controller_B.kase) - 1;
        U[payload_controller_B.qjj] = -U[payload_controller_B.qjj];
      }

      U[payload_controller_B.qq]++;
      if (payload_controller_B.qp1 - 1 >= 0) {
        U[3 * payload_controller_B.qp1] = 0.0;
      }
    } else {
      U[3 * payload_controller_B.qp1] = 0.0;
      U[3 * payload_controller_B.qp1 + 1] = 0.0;
      U[3 * payload_controller_B.qp1 + 2] = 0.0;
      U[payload_controller_B.qq] = 1.0;
    }
  }

  for (payload_controller_B.qp1 = 2; payload_controller_B.qp1 >= 0;
       payload_controller_B.qp1--) {
    if ((payload_controller_B.qp1 + 1 <= 1) && (payload_controller_B.e[0] != 0.0))
    {
      payload_controller_xaxpy(2, -(payload_controller_xdotc(2,
        payload_controller_B.Vf, 2, payload_controller_B.Vf, 5) /
        payload_controller_B.Vf[1]), 2, payload_controller_B.Vf, 5);
      payload_controller_xaxpy(2, -(payload_controller_xdotc(2,
        payload_controller_B.Vf, 2, payload_controller_B.Vf, 8) /
        payload_controller_B.Vf[1]), 2, payload_controller_B.Vf, 8);
    }

    payload_controller_B.Vf[3 * payload_controller_B.qp1] = 0.0;
    payload_controller_B.Vf[3 * payload_controller_B.qp1 + 1] = 0.0;
    payload_controller_B.Vf[3 * payload_controller_B.qp1 + 2] = 0.0;
    payload_controller_B.Vf[payload_controller_B.qp1 + 3 *
      payload_controller_B.qp1] = 1.0;
  }

  payload_controller_B.qp1 = 0;
  payload_controller_B.nrm = 0.0;
  for (payload_controller_B.qq = 0; payload_controller_B.qq < 3;
       payload_controller_B.qq++) {
    payload_controller_B.r = payload_controller_B.b_s[payload_controller_B.qq];
    if (payload_controller_B.r != 0.0) {
      payload_controller_B.rt = std::abs(payload_controller_B.r);
      payload_controller_B.r /= payload_controller_B.rt;
      payload_controller_B.b_s[payload_controller_B.qq] =
        payload_controller_B.rt;
      if (payload_controller_B.qq + 1 < 3) {
        payload_controller_B.e[payload_controller_B.qq] /=
          payload_controller_B.r;
      }

      payload_controller_B.qq_tmp = 3 * payload_controller_B.qq + 1;
      payload_controller_B.qjj = 2 + payload_controller_B.qq_tmp;
      payload_controller_B.kase = payload_controller_B.qq_tmp;
      for (payload_controller_B.d = payload_controller_B.qq_tmp;
           payload_controller_B.d <= payload_controller_B.kase;
           payload_controller_B.d += 2) {
        tmp = _mm_loadu_pd(&U[payload_controller_B.d - 1]);
        _mm_storeu_pd(&U[payload_controller_B.d - 1], _mm_mul_pd(tmp,
          _mm_set1_pd(payload_controller_B.r)));
      }

      for (payload_controller_B.d = payload_controller_B.qjj;
           payload_controller_B.d <= payload_controller_B.qq_tmp + 2;
           payload_controller_B.d++) {
        U[payload_controller_B.d - 1] *= payload_controller_B.r;
      }
    }

    if (payload_controller_B.qq + 1 < 3) {
      payload_controller_B.r = payload_controller_B.e[payload_controller_B.qq];
      if (payload_controller_B.r != 0.0) {
        payload_controller_B.rt = std::abs(payload_controller_B.r);
        payload_controller_B.r = payload_controller_B.rt /
          payload_controller_B.r;
        payload_controller_B.e[payload_controller_B.qq] =
          payload_controller_B.rt;
        payload_controller_B.b_s[payload_controller_B.qq + 1] *=
          payload_controller_B.r;
        payload_controller_B.qq_tmp = (payload_controller_B.qq + 1) * 3 + 1;
        payload_controller_B.qjj = 2 + payload_controller_B.qq_tmp;
        payload_controller_B.kase = payload_controller_B.qq_tmp;
        for (payload_controller_B.d = payload_controller_B.qq_tmp;
             payload_controller_B.d <= payload_controller_B.kase;
             payload_controller_B.d += 2) {
          tmp = _mm_loadu_pd(&payload_controller_B.Vf[payload_controller_B.d - 1]);
          _mm_storeu_pd(&payload_controller_B.Vf[payload_controller_B.d - 1],
                        _mm_mul_pd(tmp, _mm_set1_pd(payload_controller_B.r)));
        }

        for (payload_controller_B.d = payload_controller_B.qjj;
             payload_controller_B.d <= payload_controller_B.qq_tmp + 2;
             payload_controller_B.d++) {
          payload_controller_B.Vf[payload_controller_B.d - 1] *=
            payload_controller_B.r;
        }
      }
    }

    payload_controller_B.nrm = std::fmax(payload_controller_B.nrm, std::fmax(std::
      abs(payload_controller_B.b_s[payload_controller_B.qq]), std::abs
      (payload_controller_B.e[payload_controller_B.qq])));
  }

  while ((payload_controller_B.m + 2 > 0) && (payload_controller_B.qp1 < 75)) {
    payload_controller_B.qq = payload_controller_B.m + 1;
    exitg1 = false;
    while (!(exitg1 || (payload_controller_B.qq == 0))) {
      payload_controller_B.rt = std::abs
        (payload_controller_B.e[payload_controller_B.qq - 1]);
      if (payload_controller_B.rt <= (std::abs
           (payload_controller_B.b_s[payload_controller_B.qq - 1]) + std::abs
           (payload_controller_B.b_s[payload_controller_B.qq])) *
          2.2204460492503131E-16) {
        payload_controller_B.e[payload_controller_B.qq - 1] = 0.0;
        exitg1 = true;
      } else if ((payload_controller_B.rt <= 1.0020841800044864E-292) ||
                 ((payload_controller_B.qp1 > 20) && (payload_controller_B.rt <=
        2.2204460492503131E-16 * payload_controller_B.nrm))) {
        payload_controller_B.e[payload_controller_B.qq - 1] = 0.0;
        exitg1 = true;
      } else {
        payload_controller_B.qq--;
      }
    }

    if (payload_controller_B.m + 1 == payload_controller_B.qq) {
      payload_controller_B.kase = 4;
    } else {
      payload_controller_B.qjj = payload_controller_B.m + 2;
      payload_controller_B.kase = payload_controller_B.m + 2;
      exitg1 = false;
      while ((!exitg1) && (payload_controller_B.kase >= payload_controller_B.qq))
      {
        payload_controller_B.qjj = payload_controller_B.kase;
        if (payload_controller_B.kase == payload_controller_B.qq) {
          exitg1 = true;
        } else {
          payload_controller_B.rt = 0.0;
          if (payload_controller_B.kase < payload_controller_B.m + 2) {
            payload_controller_B.rt = std::abs
              (payload_controller_B.e[payload_controller_B.kase - 1]);
          }

          if (payload_controller_B.kase > payload_controller_B.qq + 1) {
            payload_controller_B.rt += std::abs
              (payload_controller_B.e[payload_controller_B.kase - 2]);
          }

          payload_controller_B.r = std::abs
            (payload_controller_B.b_s[payload_controller_B.kase - 1]);
          if ((payload_controller_B.r <= 2.2204460492503131E-16 *
               payload_controller_B.rt) || (payload_controller_B.r <=
               1.0020841800044864E-292)) {
            payload_controller_B.b_s[payload_controller_B.kase - 1] = 0.0;
            exitg1 = true;
          } else {
            payload_controller_B.kase--;
          }
        }
      }

      if (payload_controller_B.qjj == payload_controller_B.qq) {
        payload_controller_B.kase = 3;
      } else if (payload_controller_B.m + 2 == payload_controller_B.qjj) {
        payload_controller_B.kase = 1;
      } else {
        payload_controller_B.kase = 2;
        payload_controller_B.qq = payload_controller_B.qjj;
      }
    }

    switch (payload_controller_B.kase) {
     case 1:
      payload_controller_B.rt = payload_controller_B.e[payload_controller_B.m];
      payload_controller_B.e[payload_controller_B.m] = 0.0;
      for (payload_controller_B.qjj = payload_controller_B.m + 1;
           payload_controller_B.qjj >= payload_controller_B.qq + 1;
           payload_controller_B.qjj--) {
        payload_controller_xrotg
          (&payload_controller_B.b_s[payload_controller_B.qjj - 1],
           &payload_controller_B.rt, &payload_controller_B.r,
           &payload_controller_B.smm1);
        if (payload_controller_B.qjj > payload_controller_B.qq + 1) {
          payload_controller_B.rt = -payload_controller_B.smm1 *
            payload_controller_B.e[0];
          payload_controller_B.e[0] *= payload_controller_B.r;
        }

        payload_controller_xrot(payload_controller_B.Vf, 3 *
          (payload_controller_B.qjj - 1) + 1, 3 * (payload_controller_B.m + 1) +
          1, payload_controller_B.r, payload_controller_B.smm1);
      }
      break;

     case 2:
      payload_controller_B.rt = payload_controller_B.e[payload_controller_B.qq -
        1];
      payload_controller_B.e[payload_controller_B.qq - 1] = 0.0;
      for (payload_controller_B.qjj = payload_controller_B.qq + 1;
           payload_controller_B.qjj <= payload_controller_B.m + 2;
           payload_controller_B.qjj++) {
        payload_controller_xrotg
          (&payload_controller_B.b_s[payload_controller_B.qjj - 1],
           &payload_controller_B.rt, &payload_controller_B.smm1,
           &payload_controller_B.c);
        payload_controller_B.r = payload_controller_B.e[payload_controller_B.qjj
          - 1];
        payload_controller_B.rt = -payload_controller_B.c *
          payload_controller_B.r;
        payload_controller_B.e[payload_controller_B.qjj - 1] =
          payload_controller_B.r * payload_controller_B.smm1;
        payload_controller_xrot(U, 3 * (payload_controller_B.qjj - 1) + 1, 3 *
          (payload_controller_B.qq - 1) + 1, payload_controller_B.smm1,
          payload_controller_B.c);
      }
      break;

     case 3:
      payload_controller_B.rt = payload_controller_B.b_s[payload_controller_B.m
        + 1];
      payload_controller_B.r = std::fmax(std::fmax(std::fmax(std::fmax(std::abs
        (payload_controller_B.rt), std::abs
        (payload_controller_B.b_s[payload_controller_B.m])), std::abs
        (payload_controller_B.e[payload_controller_B.m])), std::abs
        (payload_controller_B.b_s[payload_controller_B.qq])), std::abs
        (payload_controller_B.e[payload_controller_B.qq]));
      tmp = _mm_set1_pd(payload_controller_B.r);
      _mm_storeu_pd(&payload_controller_B.dv[0], _mm_div_pd(_mm_set_pd
        (payload_controller_B.b_s[payload_controller_B.m],
         payload_controller_B.rt), tmp));
      payload_controller_B.rt = payload_controller_B.dv[0];
      payload_controller_B.smm1 = payload_controller_B.dv[1];
      _mm_storeu_pd(&payload_controller_B.dv[0], _mm_div_pd(_mm_set_pd
        (payload_controller_B.b_s[payload_controller_B.qq],
         payload_controller_B.e[payload_controller_B.m]), tmp));
      payload_controller_B.smm1 = ((payload_controller_B.smm1 +
        payload_controller_B.rt) * (payload_controller_B.smm1 -
        payload_controller_B.rt) + payload_controller_B.dv[0] *
        payload_controller_B.dv[0]) / 2.0;
      payload_controller_B.c = payload_controller_B.rt *
        payload_controller_B.dv[0];
      payload_controller_B.c *= payload_controller_B.c;
      if ((payload_controller_B.smm1 != 0.0) || (payload_controller_B.c != 0.0))
      {
        payload_controller_B.shift = std::sqrt(payload_controller_B.smm1 *
          payload_controller_B.smm1 + payload_controller_B.c);
        if (payload_controller_B.smm1 < 0.0) {
          payload_controller_B.shift = -payload_controller_B.shift;
        }

        payload_controller_B.shift = payload_controller_B.c /
          (payload_controller_B.smm1 + payload_controller_B.shift);
      } else {
        payload_controller_B.shift = 0.0;
      }

      payload_controller_B.rt = (payload_controller_B.dv[1] +
        payload_controller_B.rt) * (payload_controller_B.dv[1] -
        payload_controller_B.rt) + payload_controller_B.shift;
      payload_controller_B.r = payload_controller_B.e[payload_controller_B.qq] /
        payload_controller_B.r * payload_controller_B.dv[1];
      for (payload_controller_B.qq_tmp = payload_controller_B.qq + 1;
           payload_controller_B.qq_tmp <= payload_controller_B.m + 1;
           payload_controller_B.qq_tmp++) {
        payload_controller_xrotg(&payload_controller_B.rt,
          &payload_controller_B.r, &payload_controller_B.smm1,
          &payload_controller_B.c);
        if (payload_controller_B.qq_tmp > payload_controller_B.qq + 1) {
          payload_controller_B.e[0] = payload_controller_B.rt;
        }

        payload_controller_B.r =
          payload_controller_B.e[payload_controller_B.qq_tmp - 1];
        payload_controller_B.shift =
          payload_controller_B.b_s[payload_controller_B.qq_tmp - 1];
        payload_controller_B.e[payload_controller_B.qq_tmp - 1] =
          payload_controller_B.r * payload_controller_B.smm1 -
          payload_controller_B.shift * payload_controller_B.c;
        payload_controller_B.rt = payload_controller_B.c *
          payload_controller_B.b_s[payload_controller_B.qq_tmp];
        payload_controller_B.b_s[payload_controller_B.qq_tmp] *=
          payload_controller_B.smm1;
        payload_controller_B.qjj = (payload_controller_B.qq_tmp - 1) * 3 + 1;
        payload_controller_B.kase = 3 * payload_controller_B.qq_tmp + 1;
        payload_controller_xrot(payload_controller_B.Vf,
          payload_controller_B.qjj, payload_controller_B.kase,
          payload_controller_B.smm1, payload_controller_B.c);
        payload_controller_B.b_s[payload_controller_B.qq_tmp - 1] =
          payload_controller_B.shift * payload_controller_B.smm1 +
          payload_controller_B.r * payload_controller_B.c;
        payload_controller_xrotg
          (&payload_controller_B.b_s[payload_controller_B.qq_tmp - 1],
           &payload_controller_B.rt, &payload_controller_B.smm1,
           &payload_controller_B.c);
        payload_controller_B.shift =
          payload_controller_B.e[payload_controller_B.qq_tmp - 1];
        payload_controller_B.rt = payload_controller_B.shift *
          payload_controller_B.smm1 + payload_controller_B.c *
          payload_controller_B.b_s[payload_controller_B.qq_tmp];
        payload_controller_B.b_s[payload_controller_B.qq_tmp] =
          payload_controller_B.shift * -payload_controller_B.c +
          payload_controller_B.smm1 *
          payload_controller_B.b_s[payload_controller_B.qq_tmp];
        payload_controller_B.r = payload_controller_B.c *
          payload_controller_B.e[payload_controller_B.qq_tmp];
        payload_controller_B.e[payload_controller_B.qq_tmp] *=
          payload_controller_B.smm1;
        payload_controller_xrot(U, payload_controller_B.qjj,
          payload_controller_B.kase, payload_controller_B.smm1,
          payload_controller_B.c);
      }

      payload_controller_B.e[payload_controller_B.m] = payload_controller_B.rt;
      payload_controller_B.qp1++;
      break;

     default:
      if (payload_controller_B.b_s[payload_controller_B.qq] < 0.0) {
        payload_controller_B.b_s[payload_controller_B.qq] =
          -payload_controller_B.b_s[payload_controller_B.qq];
        payload_controller_B.qp1 = 3 * payload_controller_B.qq + 1;
        payload_controller_B.qjj = 2 + payload_controller_B.qp1;
        payload_controller_B.kase = payload_controller_B.qp1;
        for (payload_controller_B.qq_tmp = payload_controller_B.qp1;
             payload_controller_B.qq_tmp <= payload_controller_B.kase;
             payload_controller_B.qq_tmp += 2) {
          tmp = _mm_loadu_pd
            (&payload_controller_B.Vf[payload_controller_B.qq_tmp - 1]);
          _mm_storeu_pd(&payload_controller_B.Vf[payload_controller_B.qq_tmp - 1],
                        _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
        }

        for (payload_controller_B.qq_tmp = payload_controller_B.qjj;
             payload_controller_B.qq_tmp <= payload_controller_B.qp1 + 2;
             payload_controller_B.qq_tmp++) {
          payload_controller_B.Vf[payload_controller_B.qq_tmp - 1] =
            -payload_controller_B.Vf[payload_controller_B.qq_tmp - 1];
        }
      }

      payload_controller_B.qp1 = payload_controller_B.qq + 1;
      while ((payload_controller_B.qq + 1 < 3) &&
             (payload_controller_B.b_s[payload_controller_B.qq] <
              payload_controller_B.b_s[payload_controller_B.qp1])) {
        payload_controller_B.rt =
          payload_controller_B.b_s[payload_controller_B.qq];
        payload_controller_B.b_s[payload_controller_B.qq] =
          payload_controller_B.b_s[payload_controller_B.qp1];
        payload_controller_B.b_s[payload_controller_B.qp1] =
          payload_controller_B.rt;
        payload_controller_B.qjj = 3 * payload_controller_B.qq + 1;
        payload_controller_B.kase = (payload_controller_B.qq + 1) * 3 + 1;
        payload_controller_xswap(payload_controller_B.Vf,
          payload_controller_B.qjj, payload_controller_B.kase);
        payload_controller_xswap(U, payload_controller_B.qjj,
          payload_controller_B.kase);
        payload_controller_B.qq = payload_controller_B.qp1;
        payload_controller_B.qp1++;
      }

      payload_controller_B.qp1 = 0;
      payload_controller_B.m--;
      break;
    }
  }

  s[0] = payload_controller_B.b_s[0];
  s[1] = payload_controller_B.b_s[1];
  s[2] = payload_controller_B.b_s[2];
  if (payload_controller_B.doscale) {
    payload_controller_xzlascl_a(payload_controller_B.cscale,
      payload_controller_B.anrm, 3, 1, s, 1, 3);
  }

  for (payload_controller_B.m = 0; payload_controller_B.m < 3;
       payload_controller_B.m++) {
    V[3 * payload_controller_B.m] = payload_controller_B.Vf[3 *
      payload_controller_B.m];
    payload_controller_B.qp1 = 3 * payload_controller_B.m + 1;
    V[payload_controller_B.qp1] =
      payload_controller_B.Vf[payload_controller_B.qp1];
    payload_controller_B.qp1 = 3 * payload_controller_B.m + 2;
    V[payload_controller_B.qp1] =
      payload_controller_B.Vf[payload_controller_B.qp1];
  }
}

// Model step function
void payload_controller::step()
{
  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  boolean_T exitg1;

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[0] = payload_controller_U.Payload_Out1.q_1[0];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[3] = payload_controller_U.Payload_Out1.q_2[0];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[6] = payload_controller_U.Payload_Out1.q_3[0];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[1] = payload_controller_U.Payload_Out1.q_1[1];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[4] = payload_controller_U.Payload_Out1.q_2[1];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[7] = payload_controller_U.Payload_Out1.q_3[1];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[2] = payload_controller_U.Payload_Out1.q_1[2];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[5] = payload_controller_U.Payload_Out1.q_2[2];

  // SignalConversion generated from: '<S23>/Vector Concatenate' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_controller_B.Q[8] = payload_controller_U.Payload_Out1.q_3[2];

  // MATLAB Function: '<S23>/Moore–Penrose inverse' incorporates:
  //   Concatenate: '<S23>/Vector Concatenate'

  payload_controller_B.p = true;
  for (payload_controller_B.c_k = 0; payload_controller_B.c_k < 9;
       payload_controller_B.c_k++) {
    payload_controller_B.y[payload_controller_B.c_k] = 0.0;
    if (payload_controller_B.p) {
      // Concatenate: '<S23>/Vector Concatenate'
      payload_controller_B.absx =
        payload_controller_B.Q[payload_controller_B.c_k];
      if (std::isinf(payload_controller_B.absx) || std::isnan
          (payload_controller_B.absx)) {
        payload_controller_B.p = false;
      }
    }
  }

  if (!payload_controller_B.p) {
    for (payload_controller_B.c_k = 0; payload_controller_B.c_k < 9;
         payload_controller_B.c_k++) {
      payload_controller_B.y[payload_controller_B.c_k] = (rtNaN);
    }
  } else {
    payload_controller_svd(payload_controller_B.Q, payload_controller_B.U,
      payload_controller_B.aL_sp, payload_controller_B.V);
    payload_controller_B.absx = std::abs(payload_controller_B.aL_sp[0]);
    if (std::isinf(payload_controller_B.absx) || std::isnan
        (payload_controller_B.absx)) {
      payload_controller_B.absx = (rtNaN);
    } else if (payload_controller_B.absx < 4.4501477170144028E-308) {
      payload_controller_B.absx = 4.94065645841247E-324;
    } else {
      std::frexp(payload_controller_B.absx, &payload_controller_B.exponent);
      payload_controller_B.absx = std::ldexp(1.0, payload_controller_B.exponent
        - 53);
    }

    payload_controller_B.absx *= 3.0;
    payload_controller_B.c_k = 0;
    exitg1 = false;
    while ((!exitg1) && (payload_controller_B.c_k < 3)) {
      if (std::isinf(payload_controller_B.aL_sp[payload_controller_B.c_k]) ||
          std::isnan(payload_controller_B.aL_sp[payload_controller_B.c_k])) {
        payload_controller_B.absx = 1.7976931348623157E+308;
        exitg1 = true;
      } else {
        payload_controller_B.c_k++;
      }
    }

    payload_controller_B.c_k = -1;
    payload_controller_B.exponent = 0;
    while ((payload_controller_B.exponent < 3) &&
           (payload_controller_B.aL_sp[payload_controller_B.exponent] >
            payload_controller_B.absx)) {
      payload_controller_B.c_k++;
      payload_controller_B.exponent++;
    }

    if (payload_controller_B.c_k + 1 > 0) {
      payload_controller_B.vcol = 1;
      for (payload_controller_B.j = 0; payload_controller_B.j <=
           payload_controller_B.c_k; payload_controller_B.j++) {
        payload_controller_B.absx = 1.0 /
          payload_controller_B.aL_sp[payload_controller_B.j];
        payload_controller_B.exponent = 2 + payload_controller_B.vcol;
        payload_controller_B.vectorUB = payload_controller_B.vcol;
        for (payload_controller_B.ar = payload_controller_B.vcol;
             payload_controller_B.ar <= payload_controller_B.vectorUB;
             payload_controller_B.ar += 2) {
          tmp_1 = _mm_loadu_pd(&payload_controller_B.V[payload_controller_B.ar -
                               1]);
          _mm_storeu_pd(&payload_controller_B.V[payload_controller_B.ar - 1],
                        _mm_mul_pd(tmp_1, _mm_set1_pd(payload_controller_B.absx)));
        }

        for (payload_controller_B.ar = payload_controller_B.exponent;
             payload_controller_B.ar <= payload_controller_B.vcol + 2;
             payload_controller_B.ar++) {
          payload_controller_B.V[payload_controller_B.ar - 1] *=
            payload_controller_B.absx;
        }

        payload_controller_B.vcol += 3;
      }

      payload_controller_B.vcol = 0;
      for (payload_controller_B.exponent = 0; payload_controller_B.exponent <= 6;
           payload_controller_B.exponent += 3) {
        for (payload_controller_B.vectorUB = payload_controller_B.exponent + 1;
             payload_controller_B.vectorUB <= payload_controller_B.exponent + 3;
             payload_controller_B.vectorUB++) {
          payload_controller_B.y[payload_controller_B.vectorUB - 1] = 0.0;
        }
      }

      for (payload_controller_B.j = 0; payload_controller_B.j <= 6;
           payload_controller_B.j += 3) {
        payload_controller_B.ar = -1;
        payload_controller_B.vcol++;
        payload_controller_B.b = 3 * payload_controller_B.c_k +
          payload_controller_B.vcol;
        for (payload_controller_B.ib = payload_controller_B.vcol;
             payload_controller_B.ib <= payload_controller_B.b;
             payload_controller_B.ib += 3) {
          payload_controller_B.exponent = payload_controller_B.j + 3;
          payload_controller_B.vectorUB = payload_controller_B.j + 1;
          for (payload_controller_B.b_ic = payload_controller_B.j + 1;
               payload_controller_B.b_ic <= payload_controller_B.vectorUB;
               payload_controller_B.b_ic += 2) {
            tmp_1 = _mm_loadu_pd(&payload_controller_B.V
                                 [(payload_controller_B.ar +
              payload_controller_B.b_ic) - payload_controller_B.j]);
            tmp_0 = _mm_loadu_pd
              (&payload_controller_B.y[payload_controller_B.b_ic - 1]);
            _mm_storeu_pd(&payload_controller_B.y[payload_controller_B.b_ic - 1],
                          _mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd
              (payload_controller_B.U[payload_controller_B.ib - 1])), tmp_0));
          }

          for (payload_controller_B.b_ic = payload_controller_B.exponent;
               payload_controller_B.b_ic <= payload_controller_B.j + 3;
               payload_controller_B.b_ic++) {
            payload_controller_B.y[payload_controller_B.b_ic - 1] +=
              payload_controller_B.V[(payload_controller_B.ar +
              payload_controller_B.b_ic) - payload_controller_B.j] *
              payload_controller_B.U[payload_controller_B.ib - 1];
          }

          payload_controller_B.ar += 3;
        }
      }
    }
  }

  // End of MATLAB Function: '<S23>/Moore–Penrose inverse'

  // Gain: '<S7>/Kp' incorporates:
  //   Constant: '<S7>/pL_sp'
  //   Inport: '<Root>/Payload_Out1'
  //   Sum: '<S7>/Sum'

  payload_controller_B.absx = (0.0 - payload_controller_U.Payload_Out1.pL[0]) *
    CONTROL_PARAM.KP;

  // Saturate: '<S7>/Saturation'
  if (payload_controller_B.absx > 2.0) {
    payload_controller_B.absx = 2.0;
  } else if (payload_controller_B.absx < -2.0) {
    payload_controller_B.absx = -2.0;
  }

  // Gain: '<S24>/Kv' incorporates:
  //   Inport: '<Root>/Payload_Out1'
  //   Saturate: '<S7>/Saturation'
  //   Sum: '<S24>/Sum2'

  payload_controller_B.absx = (payload_controller_B.absx -
    payload_controller_U.Payload_Out1.vL[0]) * CONTROL_PARAM.KV;

  // Saturate: '<S24>/Saturation1'
  if (payload_controller_B.absx > 1.0) {
    payload_controller_B.absx = 1.0;
  } else if (payload_controller_B.absx < -1.0) {
    payload_controller_B.absx = -1.0;
  }

  payload_controller_Y.force_sp3[0] = payload_controller_B.absx;

  // Sum: '<S24>/Add'
  payload_controller_B.aL_sp[0] = payload_controller_B.absx +
    payload_controller_ConstB.g[0];

  // Gain: '<S7>/Kp' incorporates:
  //   Constant: '<S7>/pL_sp'
  //   Inport: '<Root>/Payload_Out1'
  //   Sum: '<S7>/Sum'

  payload_controller_B.absx = (0.0 - payload_controller_U.Payload_Out1.pL[1]) *
    CONTROL_PARAM.KP;

  // Saturate: '<S7>/Saturation'
  if (payload_controller_B.absx > 2.0) {
    payload_controller_B.absx = 2.0;
  } else if (payload_controller_B.absx < -2.0) {
    payload_controller_B.absx = -2.0;
  }

  // Gain: '<S24>/Kv' incorporates:
  //   Inport: '<Root>/Payload_Out1'
  //   Saturate: '<S7>/Saturation'
  //   Sum: '<S24>/Sum2'

  payload_controller_B.absx = (payload_controller_B.absx -
    payload_controller_U.Payload_Out1.vL[1]) * CONTROL_PARAM.KV;

  // Saturate: '<S24>/Saturation1'
  if (payload_controller_B.absx > 1.0) {
    payload_controller_B.absx = 1.0;
  } else if (payload_controller_B.absx < -1.0) {
    payload_controller_B.absx = -1.0;
  }

  payload_controller_Y.force_sp3[1] = payload_controller_B.absx;

  // Sum: '<S24>/Add'
  payload_controller_B.aL_sp[1] = payload_controller_B.absx +
    payload_controller_ConstB.g[1];

  // Gain: '<S7>/Kp' incorporates:
  //   Constant: '<S7>/pL_sp'
  //   Inport: '<Root>/Payload_Out1'
  //   Sum: '<S7>/Sum'

  payload_controller_B.absx = (-10.0 - payload_controller_U.Payload_Out1.pL[2]) *
    CONTROL_PARAM.KP;

  // Saturate: '<S7>/Saturation'
  if (payload_controller_B.absx > 2.0) {
    payload_controller_B.absx = 2.0;
  } else if (payload_controller_B.absx < -2.0) {
    payload_controller_B.absx = -2.0;
  }

  // Gain: '<S24>/Kv' incorporates:
  //   Inport: '<Root>/Payload_Out1'
  //   Saturate: '<S7>/Saturation'
  //   Sum: '<S24>/Sum2'

  payload_controller_B.absx = (payload_controller_B.absx -
    payload_controller_U.Payload_Out1.vL[2]) * CONTROL_PARAM.KV;

  // Saturate: '<S24>/Saturation1'
  if (payload_controller_B.absx > 1.0) {
    payload_controller_B.absx = 1.0;
  } else if (payload_controller_B.absx < -1.0) {
    payload_controller_B.absx = -1.0;
  }

  payload_controller_Y.force_sp3[2] = payload_controller_B.absx;

  // Sum: '<S24>/Add'
  payload_controller_B.aL_sp[2] = payload_controller_B.absx +
    payload_controller_ConstB.g[2];

  // MATLAB Function: '<S1>/q1' incorporates:
  //   Constant: '<S1>/Constant'
  //   Constant: '<S1>/Constant1'

  payload_controller_q1(1.0471975511965976, -2.0943951023931953,
                        payload_controller_B.y_f);

  // MATLAB Function: '<S12>/Vertical Control' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_control_VerticalControl(payload_controller_B.y_f,
    payload_controller_U.Payload_Out1.q_1, payload_controller_U.Payload_Out1.w_1,
    payload_controller_B.aL_sp, payload_controller_Y.force_sp1,
    payload_controller_B.state_e, CONTROL_PARAM.KQ, CONTROL_PARAM.KW,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &payload_controller_B.sf_VerticalControl);

  // MATLAB Function: '<S1>/q2' incorporates:
  //   Constant: '<S1>/Constant'
  //   Constant: '<S1>/Constant2'

  payload_controller_q1(1.0471975511965976, 0.0, payload_controller_Y.force_sp2);

  // MATLAB Function: '<S16>/Vertical Control' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_control_VerticalControl(payload_controller_Y.force_sp2,
    payload_controller_U.Payload_Out1.q_2, payload_controller_U.Payload_Out1.w_2,
    payload_controller_B.aL_sp, payload_controller_B.y_f,
    payload_controller_B.state_e, CONTROL_PARAM.KQ, CONTROL_PARAM.KW,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &payload_controller_B.sf_VerticalControl_o);

  // MATLAB Function: '<S1>/q3' incorporates:
  //   Constant: '<S1>/Constant'
  //   Constant: '<S1>/Constant3'

  payload_controller_q1(1.0471975511965976, 2.0943951023931953,
                        payload_controller_Y.force_sp3);

  // MATLAB Function: '<S20>/Vertical Control' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_control_VerticalControl(payload_controller_Y.force_sp3,
    payload_controller_U.Payload_Out1.q_3, payload_controller_U.Payload_Out1.w_3,
    payload_controller_B.aL_sp, payload_controller_Y.force_sp2,
    payload_controller_B.state_e, CONTROL_PARAM.KQ, CONTROL_PARAM.KW,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &payload_controller_B.sf_VerticalControl_i);

  // Gain: '<S23>/mL'
  tmp_1 = _mm_mul_pd(_mm_set1_pd(CONTROL_PARAM.MASS_LOAD), _mm_loadu_pd
                     (&payload_controller_B.aL_sp[0]));
  _mm_storeu_pd(&payload_controller_Y.force_sp3[0], tmp_1);
  payload_controller_Y.force_sp3[2] = CONTROL_PARAM.MASS_LOAD *
    payload_controller_B.aL_sp[2];

  // Product: '<S23>/Matrix Multiply'
  payload_controller_B.absx = payload_controller_Y.force_sp3[1];
  payload_controller_B.force_sp3 = payload_controller_Y.force_sp3[0];
  payload_controller_B.force_sp3_m = payload_controller_Y.force_sp3[2];
  for (payload_controller_B.c_k = 0; payload_controller_B.c_k <= 0;
       payload_controller_B.c_k += 2) {
    tmp_1 = _mm_loadu_pd(&payload_controller_B.y[payload_controller_B.c_k + 3]);
    tmp_0 = _mm_loadu_pd(&payload_controller_B.y[payload_controller_B.c_k]);
    tmp = _mm_loadu_pd(&payload_controller_B.y[payload_controller_B.c_k + 6]);
    _mm_storeu_pd(&payload_controller_B.MatrixMultiply[payload_controller_B.c_k],
                  _mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp_1, _mm_set1_pd
      (payload_controller_B.absx)), _mm_mul_pd(tmp_0, _mm_set1_pd
      (payload_controller_B.force_sp3))), _mm_mul_pd(tmp, _mm_set1_pd
      (payload_controller_B.force_sp3_m))));
  }

  for (payload_controller_B.c_k = 2; payload_controller_B.c_k < 3;
       payload_controller_B.c_k++) {
    payload_controller_B.MatrixMultiply[payload_controller_B.c_k] =
      (payload_controller_B.y[payload_controller_B.c_k + 3] *
       payload_controller_B.absx +
       payload_controller_B.y[payload_controller_B.c_k] *
       payload_controller_B.force_sp3) +
      payload_controller_B.y[payload_controller_B.c_k + 6] *
      payload_controller_B.force_sp3_m;
  }

  // End of Product: '<S23>/Matrix Multiply'

  // MATLAB Function: '<S19>/Parallel Control' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_control_ParallelControl(payload_controller_U.Payload_Out1.q_3,
    payload_controller_U.Payload_Out1.w_3, payload_controller_B.MatrixMultiply[2],
    payload_controller_B.aL_sp, payload_controller_Y.force_sp3,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &payload_controller_B.sf_ParallelControl_k);

  // Sum: '<S4>/Sum'
  tmp_1 = _mm_add_pd(_mm_loadu_pd(&payload_controller_Y.force_sp3[0]),
                     _mm_loadu_pd(&payload_controller_Y.force_sp2[0]));

  // Outport: '<Root>/force_sp3' incorporates:
  //   Sum: '<S4>/Sum'

  _mm_storeu_pd(&payload_controller_Y.force_sp3[0], tmp_1);
  payload_controller_Y.force_sp3[2] += payload_controller_Y.force_sp2[2];

  // MATLAB Function: '<S11>/Parallel Control' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_control_ParallelControl(payload_controller_U.Payload_Out1.q_1,
    payload_controller_U.Payload_Out1.w_1, payload_controller_B.MatrixMultiply[0],
    payload_controller_B.aL_sp, payload_controller_Y.force_sp2,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &payload_controller_B.sf_ParallelControl);

  // Sum: '<S2>/Sum'
  tmp_1 = _mm_add_pd(_mm_loadu_pd(&payload_controller_Y.force_sp2[0]),
                     _mm_loadu_pd(&payload_controller_Y.force_sp1[0]));

  // Outport: '<Root>/force_sp1' incorporates:
  //   Sum: '<S2>/Sum'

  _mm_storeu_pd(&payload_controller_Y.force_sp1[0], tmp_1);
  payload_controller_Y.force_sp1[2] += payload_controller_Y.force_sp2[2];

  // MATLAB Function: '<S15>/Parallel Control' incorporates:
  //   Inport: '<Root>/Payload_Out1'

  payload_control_ParallelControl(payload_controller_U.Payload_Out1.q_2,
    payload_controller_U.Payload_Out1.w_2, payload_controller_B.MatrixMultiply[1],
    payload_controller_B.aL_sp, payload_controller_Y.force_sp2,
    CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &payload_controller_B.sf_ParallelControl_g);

  // Sum: '<S3>/Sum'
  tmp_1 = _mm_add_pd(_mm_loadu_pd(&payload_controller_Y.force_sp2[0]),
                     _mm_loadu_pd(&payload_controller_B.y_f[0]));

  // Outport: '<Root>/force_sp2' incorporates:
  //   Sum: '<S3>/Sum'

  _mm_storeu_pd(&payload_controller_Y.force_sp2[0], tmp_1);
  payload_controller_Y.force_sp2[2] += payload_controller_B.y_f[2];
}

// Model initialize function
void payload_controller::initialize()
{
  // (no initialization code required)
}

// Model terminate function
void payload_controller::terminate()
{
  // (no terminate code required)
}

const char_T* payload_controller::RT_MODEL_payload_controller_T::getErrorStatus()
  const
{
  return (errorStatus);
}

void payload_controller::RT_MODEL_payload_controller_T::setErrorStatus(const
  char_T* const volatile aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

// Constructor
payload_controller::payload_controller() :
  payload_controller_U(),
  payload_controller_Y(),
  payload_controller_B(),
  payload_controller_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
// Currently there is no destructor body generated.
payload_controller::~payload_controller() = default;

// Real-Time Model get method
payload_controller::RT_MODEL_payload_controller_T * payload_controller::getRTM()
{
  return (&payload_controller_M);
}

//
// File trailer for generated code.
//
// [EOF]
//
