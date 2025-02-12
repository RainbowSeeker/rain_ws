/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: control_3dof.c
 *
 * Code generated for Simulink model 'control_3dof'.
 *
 * Model version                  : 1.756
 * Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
 * C/C++ source code generated on : Wed Feb 12 10:21:57 2025
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Custom Processor->Custom Processor
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "control_3dof.h"
#include "rtwtypes.h"
#include <math.h>
#include "math.h"
#ifndef UCHAR_MAX
#include <limits.h>
#endif

#if ( UCHAR_MAX != (0xFFU) ) || ( SCHAR_MAX != (0x7F) )
#error Code was generated for compiler with different sized uchar/char. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( USHRT_MAX != (0xFFFFU) ) || ( SHRT_MAX != (0x7FFF) )
#error Code was generated for compiler with different sized ushort/short. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

#if ( UINT_MAX != (0xFFFFFFFFU) ) || ( INT_MAX != (0x7FFFFFFF) )
#error Code was generated for compiler with different sized uint/int. \
Consider adjusting Test hardware word size settings on the \
Hardware Implementation pane to match your compiler word sizes as \
defined in limits.h of the compiler. Alternatively, you can \
select the Test hardware is the same as production hardware option and \
select the Enable portable word sizes option on the Code Generation > \
Verification pane for ERT based targets, which will disable the \
preprocessor word size checks.
#endif

/* Skipping ulong/long check: insufficient preprocessor integer range. */

/* Exported block parameters */
struct_e7GeJ2KaE57MoQuFkrFyOD CONTROL_PARAM = {
  1.0,
  1.5,
  1.0,
  0.0,
  20.0,

  { 3.0, 5.0, 10.0 },

  { 0.0, 0.0, 0.0 },

  { 3.0, 5.0, 10.0 },

  { 0.0, 20.0, 30.0 },
  0.2,
  0.05,
  1.0,
  0.5,
  3.0,
  0.0
} ;                                    /* Variable: CONTROL_PARAM
                                        * Referenced by:
                                        *   '<S1>/Gain'
                                        *   '<S1>/Gain1'
                                        *   '<S1>/MATLAB System'
                                        *   '<S2>/Gain'
                                        *   '<S2>/Gain1'
                                        *   '<S2>/MATLAB System'
                                        *   '<S3>/Gain'
                                        *   '<S3>/Gain1'
                                        *   '<S3>/MATLAB System'
                                        *   '<S4>/Force Saturation && Disturbution'
                                        *   '<S4>/Position 2nd ESO'
                                        *   '<S5>/Parallel Control'
                                        *   '<S6>/Vertical Control'
                                        *   '<S6>/Gain'
                                        *   '<S9>/Parallel Control'
                                        *   '<S10>/Vertical Control'
                                        *   '<S10>/Gain'
                                        *   '<S13>/Parallel Control'
                                        *   '<S14>/Vertical Control'
                                        *   '<S14>/Gain'
                                        *   '<S18>/KV'
                                        *   '<S21>/KP'
                                        *   '<S21>/KPI'
                                        */

struct_s1X7MNNw1KKuxbD9HCjUEH CONTROL_EXPORT = {
  50
} ;                                    /* Variable: CONTROL_EXPORT
                                        * Referenced by: '<S22>/Period [ms]'
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
static void control_3dof_MATLABSystem_Init(DW_MATLABSystem_control_3dof_T
  *localDW);
static void control_3dof_MATLABSystem(const real_T rtu_0[3], const real_T rtu_1
  [3], const real_T rtu_2[3], B_MATLABSystem_control_3dof_T *localB,
  DW_MATLABSystem_control_3dof_T *localDW);
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

/* System initialize for atomic system: */
static void control_3dof_MATLABSystem_Init(DW_MATLABSystem_control_3dof_T
  *localDW)
{
  /* Start for MATLABSystem: '<S1>/MATLAB System' */
  /*  Constructor */
  /* 'pos_2nd_eso:1' matlab.System */
  /*  Support name-value pair arguments when constructing object */
  /* 'pos_2nd_eso:43' setProperties(obj,nargin,varargin{:}) */
  localDW->objisempty = true;

  /* 'pos_2nd_eso:13' (3, 1) */
  /* 'pos_2nd_eso:13' K_p */
  /* 'pos_2nd_eso:15' (3, 1) */
  /* 'pos_2nd_eso:15' K_v */
  localDW->obj.isInitialized = 1;

  /*         %% Common functions */
  /*  Perform one-time calculations, such as computing constants */
  /* 'pos_2nd_eso:51' obj.Z_1 = pos; */
  /* 'pos_2nd_eso:52' obj.Z_2 = vel; */
  /* 'pos_2nd_eso:53' obj.Z_3 = zeros(size(pos)); */
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

  /* InitializeConditions for MATLABSystem: '<S1>/MATLAB System' */
  /* 'pos_2nd_eso:54' obj.is_init = false; */
  /* 'pos_2nd_eso:36' logical */
  /* 'pos_2nd_eso:36' (1, 1) */
  /* 'pos_2nd_eso:36' is_init */
  /*  Initialize / reset internal or discrete properties */
  /* 'pos_2nd_eso:92' obj.is_init = false; */
  /* 'pos_2nd_eso:36' logical */
  /* 'pos_2nd_eso:36' (1, 1) */
  /* 'pos_2nd_eso:36' is_init */
  localDW->obj.is_init = false;
}

/* Output and update for atomic system: */
static void control_3dof_MATLABSystem(const real_T rtu_0[3], const real_T rtu_1
  [3], const real_T rtu_2[3], B_MATLABSystem_control_3dof_T *localB,
  DW_MATLABSystem_control_3dof_T *localDW)
{
  real_T e_1;
  real_T e_2;

  /* MATLABSystem: '<S1>/MATLAB System' */
  /* 'pos_2nd_eso:13' (3, 1) */
  /* 'pos_2nd_eso:13' K_p */
  /* 'pos_2nd_eso:15' (3, 1) */
  /* 'pos_2nd_eso:15' K_v */
  localDW->obj.K_p[0] = CONTROL_PARAM.ESO_PI[0];
  localDW->obj.K_v[0] = CONTROL_PARAM.ESO_VI[0];
  localDW->obj.K_p[1] = CONTROL_PARAM.ESO_PI[1];
  localDW->obj.K_v[1] = CONTROL_PARAM.ESO_VI[1];
  localDW->obj.K_p[2] = CONTROL_PARAM.ESO_PI[2];
  localDW->obj.K_v[2] = CONTROL_PARAM.ESO_VI[2];

  /* 'pos_2nd_eso:58' if ~obj.is_init */
  if (!localDW->obj.is_init) {
    /* 'pos_2nd_eso:59' obj.is_init = true; */
    /* 'pos_2nd_eso:36' logical */
    /* 'pos_2nd_eso:36' (1, 1) */
    /* 'pos_2nd_eso:36' is_init */
    localDW->obj.is_init = true;

    /* 'pos_2nd_eso:60' obj.Z_1 = pos; */
    /* 'pos_2nd_eso:61' obj.Z_2 = vel; */
    /* 'pos_2nd_eso:62' obj.Z_3 = zeros(size(pos)); */
    localDW->obj.Z_1[0] = rtu_1[0];
    localDW->obj.Z_2[0] = rtu_2[0];
    localDW->obj.Z_3[0] = 0.0;
    localDW->obj.Z_1[1] = rtu_1[1];
    localDW->obj.Z_2[1] = rtu_2[1];
    localDW->obj.Z_3[1] = 0.0;
    localDW->obj.Z_1[2] = rtu_1[2];
    localDW->obj.Z_2[2] = rtu_2[2];
    localDW->obj.Z_3[2] = 0.0;
  } else {
    /* 'pos_2nd_eso:63' else */
    /* 'pos_2nd_eso:64' obj.Z_1 = obj.Z_1; */
    /* 'pos_2nd_eso:65' obj.Z_2 = obj.Z_2; */
    /* 'pos_2nd_eso:66' obj.Z_3 = obj.Z_3; */
  }

  /* 'pos_2nd_eso:68' sts = getSampleTime(obj); */
  /* 'pos_2nd_eso:69' dt = sts.SampleTime; */
  /* 'pos_2nd_eso:70' g = [0; 0; 9.8]; */
  /*  Extended state observer */
  /* 'pos_2nd_eso:73' e_1 = pos - obj.Z_1; */
  /* 'pos_2nd_eso:74' e_2 = vel - obj.Z_2; */
  /* 'pos_2nd_eso:75' dot_Z_1 = obj.Z_2 + obj.K_p(1) * e_1; */
  /* 'pos_2nd_eso:76' dot_Z_2 = acc + obj.Z_3 + g + ... */
  /* 'pos_2nd_eso:77'                       obj.K_p(2) * e_1 + obj.K_v(2) * e_2; */
  /* 'pos_2nd_eso:78' dot_Z_3 = obj.K_p(3) * e_1 + obj.K_v(3) * e_2; */
  /*  Update states */
  /* 'pos_2nd_eso:81' obj.Z_1 = obj.Z_1 + dt * dot_Z_1; */
  /* 'pos_2nd_eso:82' obj.Z_2 = obj.Z_2 + dt * dot_Z_2; */
  /* 'pos_2nd_eso:83' obj.Z_3 = obj.Z_3 + dt * dot_Z_3; */
  /* 'pos_2nd_eso:85' z1 = obj.Z_1; */
  /* 'pos_2nd_eso:86' z2 = obj.Z_2; */
  /* 'pos_2nd_eso:87' z3 = obj.Z_3; */
  e_1 = rtu_1[0] - localDW->obj.Z_1[0];
  e_2 = rtu_2[0] - localDW->obj.Z_2[0];
  localDW->obj.Z_1[0] += (localDW->obj.K_p[0] * e_1 + localDW->obj.Z_2[0]) *
    0.05;
  localDW->obj.Z_2[0] += (((rtu_0[0] + localDW->obj.Z_3[0]) + localDW->obj.K_p[1]
    * e_1) + localDW->obj.K_v[1] * e_2) * 0.05;
  e_1 = (localDW->obj.K_p[2] * e_1 + localDW->obj.K_v[2] * e_2) * 0.05 +
    localDW->obj.Z_3[0];
  localDW->obj.Z_3[0] = e_1;

  /* MATLABSystem: '<S1>/MATLAB System' */
  localB->MATLABSystem_o3[0] = e_1;

  /* MATLABSystem: '<S1>/MATLAB System' */
  e_1 = rtu_1[1] - localDW->obj.Z_1[1];
  e_2 = rtu_2[1] - localDW->obj.Z_2[1];
  localDW->obj.Z_1[1] += (localDW->obj.K_p[0] * e_1 + localDW->obj.Z_2[1]) *
    0.05;
  localDW->obj.Z_2[1] += (((rtu_0[1] + localDW->obj.Z_3[1]) + localDW->obj.K_p[1]
    * e_1) + localDW->obj.K_v[1] * e_2) * 0.05;
  e_1 = (localDW->obj.K_p[2] * e_1 + localDW->obj.K_v[2] * e_2) * 0.05 +
    localDW->obj.Z_3[1];
  localDW->obj.Z_3[1] = e_1;

  /* MATLABSystem: '<S1>/MATLAB System' */
  localB->MATLABSystem_o3[1] = e_1;

  /* MATLABSystem: '<S1>/MATLAB System' */
  e_1 = rtu_1[2] - localDW->obj.Z_1[2];
  e_2 = rtu_2[2] - localDW->obj.Z_2[2];
  localDW->obj.Z_1[2] += (localDW->obj.K_p[0] * e_1 + localDW->obj.Z_2[2]) *
    0.05;
  localDW->obj.Z_2[2] += ((((rtu_0[2] + localDW->obj.Z_3[2]) + 9.8) +
    localDW->obj.K_p[1] * e_1) + localDW->obj.K_v[1] * e_2) * 0.05;
  e_1 = (localDW->obj.K_p[2] * e_1 + localDW->obj.K_v[2] * e_2) * 0.05 +
    localDW->obj.Z_3[2];
  localDW->obj.Z_3[2] = e_1;

  /* MATLABSystem: '<S1>/MATLAB System' */
  localB->MATLABSystem_o3[2] = e_1;
}

/*
 * Output and update for atomic system:
 *    '<S5>/Parallel Control'
 *    '<S9>/Parallel Control'
 *    '<S13>/Parallel Control'
 */
static void control_3dof_ParallelControl(const real_T rtu_q[3], const real_T
  rtu_w[3], const real_T rtu_vforce[3], const real_T rtu_ai_sp[3], real_T
  rty_f_parallel[3], real_T rtp_li, real_T rtp_mi,
  B_ParallelControl_control_3do_T *localB)
{
  real_T a;
  real_T absxk;
  real_T rtu_ai_sp_0;
  real_T rtu_ai_sp_1;
  real_T rtu_ai_sp_2;
  real_T scale;
  real_T t;
  int32_T i;

  /* MATLAB Function 'Cable Controller/Parallel Control/Parallel Control': '<S7>:1' */
  /* '<S7>:1:3' f_parallel = q * q' * vforce + mi * li * (norm(w) ^ 2) * q + mi * q * q' * ai_sp; */
  scale = 3.3121686421112381E-170;
  absxk = fabs(rtu_w[0]);
  if (absxk > 3.3121686421112381E-170) {
    a = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    a = t * t;
  }

  absxk = fabs(rtu_w[1]);
  if (absxk > scale) {
    t = scale / absxk;
    a = a * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    a += t * t;
  }

  absxk = fabs(rtu_w[2]);
  if (absxk > scale) {
    t = scale / absxk;
    a = a * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    a += t * t;
  }

  scale *= sqrt(a);
  scale = rtp_mi * rtp_li * (scale * scale);
  for (i = 0; i < 3; i++) {
    int32_T rtu_q_tmp;
    localB->rtu_q[3 * i] = rtu_q[0] * rtu_q[i];
    localB->rtp_mi[3 * i] = rtp_mi * rtu_q[0] * rtu_q[i];
    rtu_q_tmp = 3 * i + 1;
    localB->rtu_q[rtu_q_tmp] = rtu_q[1] * rtu_q[i];
    localB->rtp_mi[rtu_q_tmp] = rtp_mi * rtu_q[1] * rtu_q[i];
    rtu_q_tmp = 3 * i + 2;
    localB->rtu_q[rtu_q_tmp] = rtu_q[2] * rtu_q[i];
    localB->rtp_mi[rtu_q_tmp] = rtp_mi * rtu_q[2] * rtu_q[i];
  }

  absxk = rtu_vforce[1];
  t = rtu_vforce[0];
  a = rtu_vforce[2];
  rtu_ai_sp_0 = rtu_ai_sp[1];
  rtu_ai_sp_1 = rtu_ai_sp[0];
  rtu_ai_sp_2 = rtu_ai_sp[2];
  for (i = 0; i < 3; i++) {
    rty_f_parallel[i] = (((localB->rtu_q[i + 3] * absxk + localB->rtu_q[i] * t)
                          + localB->rtu_q[i + 6] * a) + scale * rtu_q[i]) +
      ((localB->rtp_mi[i + 3] * rtu_ai_sp_0 + localB->rtp_mi[i] * rtu_ai_sp_1) +
       localB->rtp_mi[i + 6] * rtu_ai_sp_2);
  }
}

/*
 * Output and update for atomic system:
 *    '<S6>/Vertical Control'
 *    '<S10>/Vertical Control'
 *    '<S14>/Vertical Control'
 */
static void control_3dof_VerticalControl(const real_T rtu_q_sp[3], const real_T
  rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
  rty_f_vertical[3], real_T rty_dot_err[3], real_T rty_state[12], real_T rtp_Kq,
  real_T rtp_Kw, real_T rtp_li, real_T rtp_mi, B_VerticalControl_control_3do_T
  *localB)
{
  real_T rtu_ai_sp_0;
  real_T rtu_ai_sp_1;
  real_T rtu_ai_sp_2;
  real_T w_sp_tmp;
  real_T w_sp_tmp_0;
  real_T w_sp_tmp_1;
  int32_T i;
  int32_T i_0;

  /* MATLAB Function 'Cable Controller/Vertical Control/Vertical Control': '<S8>:1' */
  /* '<S8>:1:3' Sqi = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0]; */
  localB->Sqi[0] = 0.0;
  localB->Sqi[3] = -rtu_q[2];
  localB->Sqi[6] = rtu_q[1];
  localB->Sqi[1] = rtu_q[2];
  localB->Sqi[4] = 0.0;
  localB->Sqi[7] = -rtu_q[0];
  localB->Sqi[2] = -rtu_q[1];
  localB->Sqi[5] = rtu_q[0];
  localB->Sqi[8] = 0.0;

  /* '<S8>:1:4' e_q = cross(q_sp, q); */
  localB->e_q[0] = rtu_q_sp[1] * rtu_q[2] - rtu_q[1] * rtu_q_sp[2];
  localB->e_q[1] = rtu_q[0] * rtu_q_sp[2] - rtu_q_sp[0] * rtu_q[2];
  localB->e_q[2] = rtu_q_sp[0] * rtu_q[1] - rtu_q[0] * rtu_q_sp[1];

  /* '<S8>:1:7' w_sp = -Kq * e_q; */
  /* '<S8>:1:10' e_w = w + Sqi * Sqi * w_sp; */
  for (i = 0; i < 3; i++) {
    localB->w_sp_tmp[i] = -rtp_Kq * localB->e_q[i];
    for (i_0 = 0; i_0 < 3; i_0++) {
      localB->Sqi_m[i + 3 * i_0] = (localB->Sqi[3 * i_0 + 1] * localB->Sqi[i + 3]
        + localB->Sqi[3 * i_0] * localB->Sqi[i]) + localB->Sqi[3 * i_0 + 2] *
        localB->Sqi[i + 6];
    }
  }

  /* '<S8>:1:11' dot_q_sp = cross(w_sp, q_sp); */
  /* '<S8>:1:12' f_vertical = mi * li * Sqi * (-Kq * e_q -Kw * e_w - Sqi * Sqi * w_sp - (q' * w_sp) * dot_q_sp) - mi * Sqi * Sqi * ai_sp; */
  localB->a = rtp_mi * rtp_li;
  localB->rtu_q = 0.0;
  for (i = 0; i < 3; i++) {
    localB->e_w_tmp_k = (localB->Sqi_m[i + 3] * localB->w_sp_tmp[1] +
                         localB->Sqi_m[i] * localB->w_sp_tmp[0]) + localB->
      Sqi_m[i + 6] * localB->w_sp_tmp[2];
    localB->e_w_tmp[i] = localB->e_w_tmp_k;
    localB->e_w[i] = rtu_w[i] + localB->e_w_tmp_k;
    localB->rtu_q += rtu_q[i] * localB->w_sp_tmp[i];
  }

  localB->b_a[0] = (localB->w_sp_tmp[1] * rtu_q_sp[2] - rtu_q_sp[1] *
                    localB->w_sp_tmp[2]) * localB->rtu_q;
  localB->b_a[1] = (rtu_q_sp[0] * localB->w_sp_tmp[2] - localB->w_sp_tmp[0] *
                    rtu_q_sp[2]) * localB->rtu_q;
  localB->b_a[2] = (localB->w_sp_tmp[0] * rtu_q_sp[1] - rtu_q_sp[0] *
                    localB->w_sp_tmp[1]) * localB->rtu_q;
  for (i_0 = 0; i_0 < 3; i_0++) {
    localB->w_sp_tmp_c[i_0] = ((localB->w_sp_tmp[i_0] - rtp_Kw * localB->e_w[i_0])
      - localB->e_w_tmp[i_0]) - localB->b_a[i_0];
    for (i = 0; i < 3; i++) {
      localB->Sqi_m[i_0 + 3 * i] = (localB->Sqi[i_0 + 3] * rtp_mi * localB->Sqi
        [3 * i + 1] + rtp_mi * localB->Sqi[i_0] * localB->Sqi[3 * i]) +
        localB->Sqi[i_0 + 6] * rtp_mi * localB->Sqi[3 * i + 2];
    }
  }

  /* '<S8>:1:13' dot_err = -1 / mi / li * Sqi * (e_q + 2 * e_w); */
  localB->rtu_q = -1.0 / rtp_mi / rtp_li;
  localB->e_w_tmp_k = 2.0 * localB->e_w[0] + localB->e_q[0];
  localB->e_q_idx_1 = 2.0 * localB->e_w[1] + localB->e_q[1];
  localB->e_q_idx_2 = 2.0 * localB->e_w[2] + localB->e_q[2];

  /* '<S8>:1:15' state = [q_sp; e_q; w_sp; e_w]; */
  w_sp_tmp = localB->w_sp_tmp_c[1];
  w_sp_tmp_0 = localB->w_sp_tmp_c[0];
  w_sp_tmp_1 = localB->w_sp_tmp_c[2];
  rtu_ai_sp_0 = rtu_ai_sp[1];
  rtu_ai_sp_1 = rtu_ai_sp[0];
  rtu_ai_sp_2 = rtu_ai_sp[2];
  for (i = 0; i < 3; i++) {
    real_T Sqi;
    localB->Sqi_c = localB->Sqi[i + 3];
    localB->Sqi_b = localB->Sqi[i];
    Sqi = localB->Sqi[i + 6];
    rty_f_vertical[i] = ((localB->Sqi_c * localB->a * w_sp_tmp + localB->a *
                          localB->Sqi_b * w_sp_tmp_0) + Sqi * localB->a *
                         w_sp_tmp_1) - ((localB->Sqi_m[i + 3] * rtu_ai_sp_0 +
      localB->Sqi_m[i] * rtu_ai_sp_1) + localB->Sqi_m[i + 6] * rtu_ai_sp_2);
    rty_dot_err[i] = (localB->Sqi_c * localB->rtu_q * localB->e_q_idx_1 +
                      localB->rtu_q * localB->Sqi_b * localB->e_w_tmp_k) + Sqi *
      localB->rtu_q * localB->e_q_idx_2;
    rty_state[i] = rtu_q_sp[i];
    rty_state[i + 3] = localB->e_q[i];
    rty_state[i + 6] = localB->w_sp_tmp[i];
    rty_state[i + 9] = localB->e_w[i];
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_norm(const real_T x[3])
{
  real_T y;
  control_3dof_B.scale_dy = 3.3121686421112381E-170;
  control_3dof_B.absxk_lx = fabs(x[0]);
  if (control_3dof_B.absxk_lx > 3.3121686421112381E-170) {
    y = 1.0;
    control_3dof_B.scale_dy = control_3dof_B.absxk_lx;
  } else {
    control_3dof_B.t_o = control_3dof_B.absxk_lx / 3.3121686421112381E-170;
    y = control_3dof_B.t_o * control_3dof_B.t_o;
  }

  control_3dof_B.absxk_lx = fabs(x[1]);
  if (control_3dof_B.absxk_lx > control_3dof_B.scale_dy) {
    control_3dof_B.t_o = control_3dof_B.scale_dy / control_3dof_B.absxk_lx;
    y = y * control_3dof_B.t_o * control_3dof_B.t_o + 1.0;
    control_3dof_B.scale_dy = control_3dof_B.absxk_lx;
  } else {
    control_3dof_B.t_o = control_3dof_B.absxk_lx / control_3dof_B.scale_dy;
    y += control_3dof_B.t_o * control_3dof_B.t_o;
  }

  control_3dof_B.absxk_lx = fabs(x[2]);
  if (control_3dof_B.absxk_lx > control_3dof_B.scale_dy) {
    control_3dof_B.t_o = control_3dof_B.scale_dy / control_3dof_B.absxk_lx;
    y = y * control_3dof_B.t_o * control_3dof_B.t_o + 1.0;
    control_3dof_B.scale_dy = control_3dof_B.absxk_lx;
  } else {
    control_3dof_B.t_o = control_3dof_B.absxk_lx / control_3dof_B.scale_dy;
    y += control_3dof_B.t_o * control_3dof_B.t_o;
  }

  return control_3dof_B.scale_dy * sqrt(y);
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_xzlangeM(const real_T x[9])
{
  real_T y;
  boolean_T exitg1;
  y = 0.0;
  control_3dof_B.k_f = 0;
  exitg1 = false;
  while ((!exitg1) && (control_3dof_B.k_f < 9)) {
    control_3dof_B.absxk_b = fabs(x[control_3dof_B.k_f]);
    if (rtIsNaN(control_3dof_B.absxk_b)) {
      y = (rtNaN);
      exitg1 = true;
    } else {
      if (control_3dof_B.absxk_b > y) {
        y = control_3dof_B.absxk_b;
      }

      control_3dof_B.k_f++;
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xzlascl(real_T cfrom, real_T cto, int32_T m, int32_T n,
  real_T A[9], int32_T iA0, int32_T lda)
{
  control_3dof_B.cfromc_m = cfrom;
  control_3dof_B.ctoc_n = cto;
  control_3dof_B.notdone_n = true;
  while (control_3dof_B.notdone_n) {
    control_3dof_B.cfrom1_p = control_3dof_B.cfromc_m * 2.0041683600089728E-292;
    control_3dof_B.cto1_l = control_3dof_B.ctoc_n / 4.9896007738368E+291;
    if ((fabs(control_3dof_B.cfrom1_p) > fabs(control_3dof_B.ctoc_n)) &&
        (control_3dof_B.ctoc_n != 0.0)) {
      control_3dof_B.mul_j = 2.0041683600089728E-292;
      control_3dof_B.cfromc_m = control_3dof_B.cfrom1_p;
    } else if (fabs(control_3dof_B.cto1_l) > fabs(control_3dof_B.cfromc_m)) {
      control_3dof_B.mul_j = 4.9896007738368E+291;
      control_3dof_B.ctoc_n = control_3dof_B.cto1_l;
    } else {
      control_3dof_B.mul_j = control_3dof_B.ctoc_n / control_3dof_B.cfromc_m;
      control_3dof_B.notdone_n = false;
    }

    for (control_3dof_B.j_d = 0; control_3dof_B.j_d < n; control_3dof_B.j_d++) {
      control_3dof_B.offset_b = (control_3dof_B.j_d * lda + iA0) - 2;
      for (control_3dof_B.b_i_e = 0; control_3dof_B.b_i_e < m;
           control_3dof_B.b_i_e++) {
        control_3dof_B.i1 = (control_3dof_B.b_i_e + control_3dof_B.offset_b) + 1;
        A[control_3dof_B.i1] *= control_3dof_B.mul_j;
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
      control_3dof_B.scale_g = 3.3121686421112381E-170;
      control_3dof_B.kend_b = ix0 + n;
      for (control_3dof_B.k_j = ix0; control_3dof_B.k_j < control_3dof_B.kend_b;
           control_3dof_B.k_j++) {
        control_3dof_B.absxk_l = fabs(x[control_3dof_B.k_j - 1]);
        if (control_3dof_B.absxk_l > control_3dof_B.scale_g) {
          control_3dof_B.t_d = control_3dof_B.scale_g / control_3dof_B.absxk_l;
          y = y * control_3dof_B.t_d * control_3dof_B.t_d + 1.0;
          control_3dof_B.scale_g = control_3dof_B.absxk_l;
        } else {
          control_3dof_B.t_d = control_3dof_B.absxk_l / control_3dof_B.scale_g;
          y += control_3dof_B.t_d * control_3dof_B.t_d;
        }
      }

      y = control_3dof_B.scale_g * sqrt(y);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static real_T control_3dof_xdotc(int32_T n, const real_T x[9], int32_T ix0,
  const real_T y[9], int32_T iy0)
{
  real_T d;
  d = 0.0;
  if (n >= 1) {
    for (control_3dof_B.k_o = 0; control_3dof_B.k_o < n; control_3dof_B.k_o++) {
      d += x[(ix0 + control_3dof_B.k_o) - 1] * y[(iy0 + control_3dof_B.k_o) - 1];
    }
  }

  return d;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[9],
  int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (control_3dof_B.k_a = 0; control_3dof_B.k_a < n; control_3dof_B.k_a++) {
      control_3dof_B.i2 = (iy0 + control_3dof_B.k_a) - 1;
      y[control_3dof_B.i2] += y[(ix0 + control_3dof_B.k_a) - 1] * a;
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
      control_3dof_B.scale_d = 3.3121686421112381E-170;
      control_3dof_B.kend = ix0 + n;
      for (control_3dof_B.k = ix0; control_3dof_B.k < control_3dof_B.kend;
           control_3dof_B.k++) {
        control_3dof_B.absxk = fabs(x[control_3dof_B.k - 1]);
        if (control_3dof_B.absxk > control_3dof_B.scale_d) {
          control_3dof_B.t = control_3dof_B.scale_d / control_3dof_B.absxk;
          y = y * control_3dof_B.t * control_3dof_B.t + 1.0;
          control_3dof_B.scale_d = control_3dof_B.absxk;
        } else {
          control_3dof_B.t = control_3dof_B.absxk / control_3dof_B.scale_d;
          y += control_3dof_B.t * control_3dof_B.t;
        }
      }

      y = control_3dof_B.scale_d * sqrt(y);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xaxpy_f(int32_T n, real_T a, const real_T x[9], int32_T
  ix0, real_T y[3], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (control_3dof_B.k_jz = 0; control_3dof_B.k_jz < n; control_3dof_B.k_jz++)
    {
      control_3dof_B.i4 = (iy0 + control_3dof_B.k_jz) - 1;
      y[control_3dof_B.i4] += x[(ix0 + control_3dof_B.k_jz) - 1] * a;
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3],
  int32_T ix0, real_T y[9], int32_T iy0)
{
  if ((n >= 1) && (!(a == 0.0))) {
    for (control_3dof_B.k_ju = 0; control_3dof_B.k_ju < n; control_3dof_B.k_ju++)
    {
      control_3dof_B.i3 = (iy0 + control_3dof_B.k_ju) - 1;
      y[control_3dof_B.i3] += x[(ix0 + control_3dof_B.k_ju) - 1] * a;
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
      for (control_3dof_B.b_i = 0; control_3dof_B.b_i < m; control_3dof_B.b_i++)
      {
        control_3dof_B.i_h = (control_3dof_B.b_i + control_3dof_B.offset) + 1;
        A[control_3dof_B.i_h] *= control_3dof_B.mul;
      }
    }
  }
}

/* Function for MATLAB Function: '<S4>/Force Saturation && Disturbution' */
static void control_3dof_xswap(real_T x[9], int32_T ix0, int32_T iy0)
{
  control_3dof_B.temp_n = x[ix0 - 1];
  x[ix0 - 1] = x[iy0 - 1];
  x[iy0 - 1] = control_3dof_B.temp_n;
  control_3dof_B.temp_n = x[ix0];
  x[ix0] = x[iy0];
  x[iy0] = control_3dof_B.temp_n;
  control_3dof_B.temp_n = x[ix0 + 1];
  x[ix0 + 1] = x[iy0 + 1];
  x[iy0 + 1] = control_3dof_B.temp_n;
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
  for (control_3dof_B.qs = 0; control_3dof_B.qs < 9; control_3dof_B.qs++) {
    control_3dof_B.b_A[control_3dof_B.qs] = A[control_3dof_B.qs];
    U[control_3dof_B.qs] = 0.0;
    control_3dof_B.Vf[control_3dof_B.qs] = 0.0;
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
    control_3dof_B.qs = 3 * control_3dof_B.m + control_3dof_B.m;
    control_3dof_B.qq = control_3dof_B.qs + 1;
    control_3dof_B.apply_transform = false;
    control_3dof_B.nrm = control_3dof_xnrm2(3 - control_3dof_B.m,
      control_3dof_B.b_A, control_3dof_B.qs + 1);
    if (control_3dof_B.nrm > 0.0) {
      control_3dof_B.apply_transform = true;
      if (control_3dof_B.b_A[control_3dof_B.qs] < 0.0) {
        control_3dof_B.nrm = -control_3dof_B.nrm;
      }

      control_3dof_B.b_s[control_3dof_B.m] = control_3dof_B.nrm;
      if (fabs(control_3dof_B.nrm) >= 1.0020841800044864E-292) {
        control_3dof_B.nrm = 1.0 / control_3dof_B.nrm;
        control_3dof_B.qjj = (control_3dof_B.qs - control_3dof_B.m) + 3;
        for (control_3dof_B.e_k = control_3dof_B.qq; control_3dof_B.e_k <=
             control_3dof_B.qjj; control_3dof_B.e_k++) {
          control_3dof_B.b_A[control_3dof_B.e_k - 1] *= control_3dof_B.nrm;
        }
      } else {
        control_3dof_B.qjj = (control_3dof_B.qs - control_3dof_B.m) + 3;
        for (control_3dof_B.e_k = control_3dof_B.qq; control_3dof_B.e_k <=
             control_3dof_B.qjj; control_3dof_B.e_k++) {
          control_3dof_B.b_A[control_3dof_B.e_k - 1] /=
            control_3dof_B.b_s[control_3dof_B.m];
        }
      }

      control_3dof_B.b_A[control_3dof_B.qs]++;
      control_3dof_B.b_s[control_3dof_B.m] =
        -control_3dof_B.b_s[control_3dof_B.m];
    } else {
      control_3dof_B.b_s[control_3dof_B.m] = 0.0;
    }

    for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
         control_3dof_B.qq++) {
      control_3dof_B.qjj = (control_3dof_B.qq - 1) * 3 + control_3dof_B.m;
      if (control_3dof_B.apply_transform) {
        control_3dof_xaxpy(3 - control_3dof_B.m, -(control_3dof_xdotc(3 -
          control_3dof_B.m, control_3dof_B.b_A, control_3dof_B.qs + 1,
          control_3dof_B.b_A, control_3dof_B.qjj + 1) /
          control_3dof_B.b_A[control_3dof_B.qs]), control_3dof_B.qs + 1,
                           control_3dof_B.b_A, control_3dof_B.qjj + 1);
      }

      control_3dof_B.e[control_3dof_B.qq - 1] =
        control_3dof_B.b_A[control_3dof_B.qjj];
    }

    for (control_3dof_B.qq = control_3dof_B.m + 1; control_3dof_B.qq < 4;
         control_3dof_B.qq++) {
      control_3dof_B.qs = (3 * control_3dof_B.m + control_3dof_B.qq) - 1;
      U[control_3dof_B.qs] = control_3dof_B.b_A[control_3dof_B.qs];
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
          for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
               control_3dof_B.qq++) {
            control_3dof_B.e[control_3dof_B.qq - 1] *= control_3dof_B.nrm;
          }
        } else {
          for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
               control_3dof_B.qq++) {
            control_3dof_B.e[control_3dof_B.qq - 1] /= control_3dof_B.nrm;
          }
        }

        control_3dof_B.e[1]++;
        control_3dof_B.e[0] = -control_3dof_B.e[0];
        for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
             control_3dof_B.qq++) {
          control_3dof_B.work[control_3dof_B.qq - 1] = 0.0;
        }

        for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
             control_3dof_B.qq++) {
          control_3dof_xaxpy_f(2, control_3dof_B.e[control_3dof_B.qq - 1],
                               control_3dof_B.b_A, 3 * (control_3dof_B.qq - 1) +
                               2, control_3dof_B.work, 2);
        }

        for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
             control_3dof_B.qq++) {
          control_3dof_xaxpy_fh(2, -control_3dof_B.e[control_3dof_B.qq - 1] /
                                control_3dof_B.e[1], control_3dof_B.work, 2,
                                control_3dof_B.b_A, 3 * (control_3dof_B.qq - 1)
                                + 2);
        }
      }

      for (control_3dof_B.qq = control_3dof_B.qp1; control_3dof_B.qq < 4;
           control_3dof_B.qq++) {
        control_3dof_B.Vf[control_3dof_B.qq - 1] =
          control_3dof_B.e[control_3dof_B.qq - 1];
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
      for (control_3dof_B.qs = control_3dof_B.qp1 + 2; control_3dof_B.qs < 4;
           control_3dof_B.qs++) {
        control_3dof_B.qjj = ((control_3dof_B.qs - 1) * 3 + control_3dof_B.qp1)
          + 1;
        control_3dof_xaxpy(3 - control_3dof_B.qp1, -(control_3dof_xdotc(3 -
          control_3dof_B.qp1, U, control_3dof_B.qq + 1, U, control_3dof_B.qjj) /
          U[control_3dof_B.qq]), control_3dof_B.qq + 1, U, control_3dof_B.qjj);
      }

      for (control_3dof_B.qjj = control_3dof_B.qp1 + 1; control_3dof_B.qjj < 4;
           control_3dof_B.qjj++) {
        control_3dof_B.qs = (3 * control_3dof_B.qp1 + control_3dof_B.qjj) - 1;
        U[control_3dof_B.qs] = -U[control_3dof_B.qs];
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

      control_3dof_B.qs = 3 * control_3dof_B.qq + 1;
      for (control_3dof_B.qjj = control_3dof_B.qs; control_3dof_B.qjj <=
           control_3dof_B.qs + 2; control_3dof_B.qjj++) {
        U[control_3dof_B.qjj - 1] *= control_3dof_B.r;
      }
    }

    if (control_3dof_B.qq + 1 < 3) {
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq];
      if (control_3dof_B.r != 0.0) {
        control_3dof_B.rt = fabs(control_3dof_B.r);
        control_3dof_B.r = control_3dof_B.rt / control_3dof_B.r;
        control_3dof_B.e[control_3dof_B.qq] = control_3dof_B.rt;
        control_3dof_B.b_s[control_3dof_B.qq + 1] *= control_3dof_B.r;
        control_3dof_B.qs = (control_3dof_B.qq + 1) * 3 + 1;
        for (control_3dof_B.qjj = control_3dof_B.qs; control_3dof_B.qjj <=
             control_3dof_B.qs + 2; control_3dof_B.qjj++) {
          control_3dof_B.Vf[control_3dof_B.qjj - 1] *= control_3dof_B.r;
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
      control_3dof_B.qjj = 4;
    } else {
      control_3dof_B.qs = control_3dof_B.m + 2;
      control_3dof_B.qjj = control_3dof_B.m + 2;
      exitg1 = false;
      while ((!exitg1) && (control_3dof_B.qjj >= control_3dof_B.qq)) {
        control_3dof_B.qs = control_3dof_B.qjj;
        if (control_3dof_B.qjj == control_3dof_B.qq) {
          exitg1 = true;
        } else {
          control_3dof_B.rt = 0.0;
          if (control_3dof_B.qjj < control_3dof_B.m + 2) {
            control_3dof_B.rt = fabs(control_3dof_B.e[control_3dof_B.qjj - 1]);
          }

          if (control_3dof_B.qjj > control_3dof_B.qq + 1) {
            control_3dof_B.rt += fabs(control_3dof_B.e[control_3dof_B.qjj - 2]);
          }

          control_3dof_B.r = fabs(control_3dof_B.b_s[control_3dof_B.qjj - 1]);
          if ((control_3dof_B.r <= 2.2204460492503131E-16 * control_3dof_B.rt) ||
              (control_3dof_B.r <= 1.0020841800044864E-292)) {
            control_3dof_B.b_s[control_3dof_B.qjj - 1] = 0.0;
            exitg1 = true;
          } else {
            control_3dof_B.qjj--;
          }
        }
      }

      if (control_3dof_B.qs == control_3dof_B.qq) {
        control_3dof_B.qjj = 3;
      } else if (control_3dof_B.m + 2 == control_3dof_B.qs) {
        control_3dof_B.qjj = 1;
      } else {
        control_3dof_B.qjj = 2;
        control_3dof_B.qq = control_3dof_B.qs;
      }
    }

    switch (control_3dof_B.qjj) {
     case 1:
      control_3dof_B.rt = control_3dof_B.e[control_3dof_B.m];
      control_3dof_B.e[control_3dof_B.m] = 0.0;
      for (control_3dof_B.qs = control_3dof_B.m + 1; control_3dof_B.qs >=
           control_3dof_B.qq + 1; control_3dof_B.qs--) {
        control_3dof_xrotg(&control_3dof_B.b_s[control_3dof_B.qs - 1],
                           &control_3dof_B.rt, &control_3dof_B.r,
                           &control_3dof_B.sqds);
        if (control_3dof_B.qs > control_3dof_B.qq + 1) {
          control_3dof_B.rt = -control_3dof_B.sqds * control_3dof_B.e[0];
          control_3dof_B.e[0] *= control_3dof_B.r;
        }

        control_3dof_xrot(control_3dof_B.Vf, 3 * (control_3dof_B.qs - 1) + 1, 3 *
                          (control_3dof_B.m + 1) + 1, control_3dof_B.r,
                          control_3dof_B.sqds);
      }
      break;

     case 2:
      control_3dof_B.rt = control_3dof_B.e[control_3dof_B.qq - 1];
      control_3dof_B.e[control_3dof_B.qq - 1] = 0.0;
      for (control_3dof_B.qs = control_3dof_B.qq + 1; control_3dof_B.qs <=
           control_3dof_B.m + 2; control_3dof_B.qs++) {
        control_3dof_xrotg(&control_3dof_B.b_s[control_3dof_B.qs - 1],
                           &control_3dof_B.rt, &control_3dof_B.sqds,
                           &control_3dof_B.smm1);
        control_3dof_B.r = control_3dof_B.e[control_3dof_B.qs - 1];
        control_3dof_B.rt = -control_3dof_B.smm1 * control_3dof_B.r;
        control_3dof_B.e[control_3dof_B.qs - 1] = control_3dof_B.r *
          control_3dof_B.sqds;
        control_3dof_xrot(U, 3 * (control_3dof_B.qs - 1) + 1, 3 *
                          (control_3dof_B.qq - 1) + 1, control_3dof_B.sqds,
                          control_3dof_B.smm1);
      }
      break;

     case 3:
      control_3dof_B.rt = control_3dof_B.b_s[control_3dof_B.m + 1];
      control_3dof_B.r = fmax(fmax(fmax(fmax(fabs(control_3dof_B.rt), fabs
        (control_3dof_B.b_s[control_3dof_B.m])), fabs
        (control_3dof_B.e[control_3dof_B.m])), fabs
        (control_3dof_B.b_s[control_3dof_B.qq])), fabs
        (control_3dof_B.e[control_3dof_B.qq]));
      control_3dof_B.rt /= control_3dof_B.r;
      control_3dof_B.smm1 = control_3dof_B.b_s[control_3dof_B.m] /
        control_3dof_B.r;
      control_3dof_B.emm1 = control_3dof_B.e[control_3dof_B.m] /
        control_3dof_B.r;
      control_3dof_B.sqds = control_3dof_B.b_s[control_3dof_B.qq] /
        control_3dof_B.r;
      control_3dof_B.smm1 = ((control_3dof_B.smm1 + control_3dof_B.rt) *
        (control_3dof_B.smm1 - control_3dof_B.rt) + control_3dof_B.emm1 *
        control_3dof_B.emm1) / 2.0;
      control_3dof_B.emm1 *= control_3dof_B.rt;
      control_3dof_B.emm1 *= control_3dof_B.emm1;
      if ((control_3dof_B.smm1 != 0.0) || (control_3dof_B.emm1 != 0.0)) {
        control_3dof_B.shift = sqrt(control_3dof_B.smm1 * control_3dof_B.smm1 +
          control_3dof_B.emm1);
        if (control_3dof_B.smm1 < 0.0) {
          control_3dof_B.shift = -control_3dof_B.shift;
        }

        control_3dof_B.shift = control_3dof_B.emm1 / (control_3dof_B.smm1 +
          control_3dof_B.shift);
      } else {
        control_3dof_B.shift = 0.0;
      }

      control_3dof_B.rt = (control_3dof_B.sqds + control_3dof_B.rt) *
        (control_3dof_B.sqds - control_3dof_B.rt) + control_3dof_B.shift;
      control_3dof_B.r = control_3dof_B.e[control_3dof_B.qq] / control_3dof_B.r *
        control_3dof_B.sqds;
      for (control_3dof_B.e_k = control_3dof_B.qq + 1; control_3dof_B.e_k <=
           control_3dof_B.m + 1; control_3dof_B.e_k++) {
        control_3dof_xrotg(&control_3dof_B.rt, &control_3dof_B.r,
                           &control_3dof_B.sqds, &control_3dof_B.smm1);
        if (control_3dof_B.e_k > control_3dof_B.qq + 1) {
          control_3dof_B.e[0] = control_3dof_B.rt;
        }

        control_3dof_B.r = control_3dof_B.e[control_3dof_B.e_k - 1];
        control_3dof_B.emm1 = control_3dof_B.b_s[control_3dof_B.e_k - 1];
        control_3dof_B.e[control_3dof_B.e_k - 1] = control_3dof_B.r *
          control_3dof_B.sqds - control_3dof_B.emm1 * control_3dof_B.smm1;
        control_3dof_B.rt = control_3dof_B.smm1 *
          control_3dof_B.b_s[control_3dof_B.e_k];
        control_3dof_B.b_s[control_3dof_B.e_k] *= control_3dof_B.sqds;
        control_3dof_B.qs = (control_3dof_B.e_k - 1) * 3 + 1;
        control_3dof_B.qjj = 3 * control_3dof_B.e_k + 1;
        control_3dof_xrot(control_3dof_B.Vf, control_3dof_B.qs,
                          control_3dof_B.qjj, control_3dof_B.sqds,
                          control_3dof_B.smm1);
        control_3dof_B.b_s[control_3dof_B.e_k - 1] = control_3dof_B.emm1 *
          control_3dof_B.sqds + control_3dof_B.r * control_3dof_B.smm1;
        control_3dof_xrotg(&control_3dof_B.b_s[control_3dof_B.e_k - 1],
                           &control_3dof_B.rt, &control_3dof_B.sqds,
                           &control_3dof_B.smm1);
        control_3dof_B.emm1 = control_3dof_B.e[control_3dof_B.e_k - 1];
        control_3dof_B.rt = control_3dof_B.emm1 * control_3dof_B.sqds +
          control_3dof_B.smm1 * control_3dof_B.b_s[control_3dof_B.e_k];
        control_3dof_B.b_s[control_3dof_B.e_k] = control_3dof_B.emm1 *
          -control_3dof_B.smm1 + control_3dof_B.sqds *
          control_3dof_B.b_s[control_3dof_B.e_k];
        control_3dof_B.r = control_3dof_B.smm1 *
          control_3dof_B.e[control_3dof_B.e_k];
        control_3dof_B.e[control_3dof_B.e_k] *= control_3dof_B.sqds;
        control_3dof_xrot(U, control_3dof_B.qs, control_3dof_B.qjj,
                          control_3dof_B.sqds, control_3dof_B.smm1);
      }

      control_3dof_B.e[control_3dof_B.m] = control_3dof_B.rt;
      control_3dof_B.qp1++;
      break;

     default:
      if (control_3dof_B.b_s[control_3dof_B.qq] < 0.0) {
        control_3dof_B.b_s[control_3dof_B.qq] =
          -control_3dof_B.b_s[control_3dof_B.qq];
        control_3dof_B.qp1 = 3 * control_3dof_B.qq + 1;
        for (control_3dof_B.qs = control_3dof_B.qp1; control_3dof_B.qs <=
             control_3dof_B.qp1 + 2; control_3dof_B.qs++) {
          control_3dof_B.Vf[control_3dof_B.qs - 1] =
            -control_3dof_B.Vf[control_3dof_B.qs - 1];
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
        control_3dof_B.qs = 3 * control_3dof_B.qq + 1;
        control_3dof_B.qjj = (control_3dof_B.qq + 1) * 3 + 1;
        control_3dof_xswap(control_3dof_B.Vf, control_3dof_B.qs,
                           control_3dof_B.qjj);
        control_3dof_xswap(U, control_3dof_B.qs, control_3dof_B.qjj);
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
  boolean_T exitg1;
  control_3dof_B.p = true;
  for (control_3dof_B.r_b = 0; control_3dof_B.r_b < 9; control_3dof_B.r_b++) {
    X[control_3dof_B.r_b] = 0.0;
    if (control_3dof_B.p) {
      control_3dof_B.absx = A[control_3dof_B.r_b];
      if (rtIsInf(control_3dof_B.absx) || rtIsNaN(control_3dof_B.absx)) {
        control_3dof_B.p = false;
      }
    }
  }

  if (!control_3dof_B.p) {
    for (control_3dof_B.r_b = 0; control_3dof_B.r_b < 9; control_3dof_B.r_b++) {
      X[control_3dof_B.r_b] = (rtNaN);
    }
  } else {
    control_3dof_svd(A, control_3dof_B.U, control_3dof_B.s, control_3dof_B.V);
    control_3dof_B.absx = fabs(control_3dof_B.s[0]);
    if (rtIsInf(control_3dof_B.absx) || rtIsNaN(control_3dof_B.absx)) {
      control_3dof_B.absx = (rtNaN);
    } else if (control_3dof_B.absx < 4.4501477170144028E-308) {
      control_3dof_B.absx = 4.94065645841247E-324;
    } else {
      frexp(control_3dof_B.absx, &control_3dof_B.vcol);
      control_3dof_B.absx = ldexp(1.0, control_3dof_B.vcol - 53);
    }

    control_3dof_B.absx *= 3.0;
    control_3dof_B.r_b = 0;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.r_b < 3)) {
      if (rtIsInf(control_3dof_B.s[control_3dof_B.r_b]) || rtIsNaN
          (control_3dof_B.s[control_3dof_B.r_b])) {
        control_3dof_B.absx = 1.7976931348623157E+308;
        exitg1 = true;
      } else {
        control_3dof_B.r_b++;
      }
    }

    control_3dof_B.r_b = -1;
    control_3dof_B.vcol = 0;
    while ((control_3dof_B.vcol < 3) && (control_3dof_B.s[control_3dof_B.vcol] >
            control_3dof_B.absx)) {
      control_3dof_B.r_b++;
      control_3dof_B.vcol++;
    }

    if (control_3dof_B.r_b + 1 > 0) {
      control_3dof_B.vcol = 1;
      for (control_3dof_B.j = 0; control_3dof_B.j <= control_3dof_B.r_b;
           control_3dof_B.j++) {
        control_3dof_B.absx = 1.0 / control_3dof_B.s[control_3dof_B.j];
        for (control_3dof_B.ar = control_3dof_B.vcol; control_3dof_B.ar <=
             control_3dof_B.vcol + 2; control_3dof_B.ar++) {
          control_3dof_B.V[control_3dof_B.ar - 1] *= control_3dof_B.absx;
        }

        control_3dof_B.vcol += 3;
      }

      control_3dof_B.vcol = 0;
      for (control_3dof_B.j = 0; control_3dof_B.j <= 6; control_3dof_B.j += 3) {
        for (control_3dof_B.ar = control_3dof_B.j + 1; control_3dof_B.ar <=
             control_3dof_B.j + 3; control_3dof_B.ar++) {
          X[control_3dof_B.ar - 1] = 0.0;
        }
      }

      for (control_3dof_B.j = 0; control_3dof_B.j <= 6; control_3dof_B.j += 3) {
        control_3dof_B.ar = -1;
        control_3dof_B.vcol++;
        control_3dof_B.b = 3 * control_3dof_B.r_b + control_3dof_B.vcol;
        for (control_3dof_B.ib = control_3dof_B.vcol; control_3dof_B.ib <=
             control_3dof_B.b; control_3dof_B.ib += 3) {
          for (control_3dof_B.b_ic = control_3dof_B.j + 1; control_3dof_B.b_ic <=
               control_3dof_B.j + 3; control_3dof_B.b_ic++) {
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
  static const int8_T a[18] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0, 1, 0,
    0, -1 };

  static const real_T b[3] = { 0.0, 0.0, -9.8 };

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

  /* RelationalOperator: '<S23>/FixPt Relational Operator' incorporates:
   *  Inport: '<Root>/Dir_sp'
   *  UnitDelay: '<S23>/Delay Input1'
   *
   * Block description for '<S23>/Delay Input1':
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
   *  UnitDelay: '<S4>/Unit Delay'
   * */
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

  /* 'pos_2nd_eso:58' if ~obj.is_init */
  if (!control_3dof_DW.obj.is_init) {
    /* 'pos_2nd_eso:59' obj.is_init = true; */
    /* 'pos_2nd_eso:36' logical */
    /* 'pos_2nd_eso:36' (1, 1) */
    /* 'pos_2nd_eso:36' is_init */
    control_3dof_DW.obj.is_init = true;

    /* 'pos_2nd_eso:60' obj.Z_1 = pos; */
    /* 'pos_2nd_eso:61' obj.Z_2 = vel; */
    /* 'pos_2nd_eso:62' obj.Z_3 = zeros(size(pos)); */
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
    /* 'pos_2nd_eso:63' else */
    /* 'pos_2nd_eso:64' obj.Z_1 = obj.Z_1; */
    /* 'pos_2nd_eso:65' obj.Z_2 = obj.Z_2; */
    /* 'pos_2nd_eso:66' obj.Z_3 = obj.Z_3; */
  }

  /* 'pos_2nd_eso:68' sts = getSampleTime(obj); */
  /* 'pos_2nd_eso:69' dt = sts.SampleTime; */
  /* 'pos_2nd_eso:70' g = [0; 0; 9.8]; */
  /*  Extended state observer */
  /* 'pos_2nd_eso:73' e_1 = pos - obj.Z_1; */
  /* 'pos_2nd_eso:74' e_2 = vel - obj.Z_2; */
  /* 'pos_2nd_eso:75' dot_Z_1 = obj.Z_2 + obj.K_p(1) * e_1; */
  /* 'pos_2nd_eso:76' dot_Z_2 = acc + obj.Z_3 + g + ... */
  /* 'pos_2nd_eso:77'                       obj.K_p(2) * e_1 + obj.K_v(2) * e_2; */
  /* 'pos_2nd_eso:78' dot_Z_3 = obj.K_p(3) * e_1 + obj.K_v(3) * e_2; */
  /*  Update states */
  /* 'pos_2nd_eso:81' obj.Z_1 = obj.Z_1 + dt * dot_Z_1; */
  /* 'pos_2nd_eso:82' obj.Z_2 = obj.Z_2 + dt * dot_Z_2; */
  /* 'pos_2nd_eso:83' obj.Z_3 = obj.Z_3 + dt * dot_Z_3; */
  /* 'pos_2nd_eso:85' z1 = obj.Z_1; */
  /* 'pos_2nd_eso:86' z2 = obj.Z_2; */
  /* 'pos_2nd_eso:87' z3 = obj.Z_3; */
  /* MATLAB Function 'Payload Controller/Force Saturation && Disturbution': '<S17>:1' */
  /* '<S17>:1:3' q_1 = Payload_Out.pL - Payload_Out.p_1; */
  control_3dof_B.e_1 = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_DW.obj.Z_1[0];
  control_3dof_B.e_2_k = control_3dof_U.Payload_Out.vL[0] -
    control_3dof_DW.obj.Z_2[0];
  control_3dof_DW.obj.Z_1[0] += (control_3dof_DW.obj.K_p[0] * control_3dof_B.e_1
    + control_3dof_DW.obj.Z_2[0]) * 0.05;
  control_3dof_DW.obj.Z_2[0] += (((control_3dof_DW.UnitDelay_DSTATE[0] +
    control_3dof_DW.obj.Z_3[0]) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.e_1) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_k) *
    0.05;
  control_3dof_DW.obj.Z_3[0] += (control_3dof_DW.obj.K_p[2] * control_3dof_B.e_1
    + control_3dof_DW.obj.K_v[2] * control_3dof_B.e_2_k) * 0.05;

  /* Sum: '<S18>/Sum' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  control_3dof_B.e_1_idx_0 = control_3dof_U.Traj_sp.pos_sp[0] -
    control_3dof_U.Payload_Out.pL[0];

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  control_3dof_B.e_2[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_1[0];

  /* MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  UnitDelay: '<S4>/Unit Delay'
   * */
  control_3dof_B.e_1 = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_DW.obj.Z_1[1];
  control_3dof_B.e_2_k = control_3dof_U.Payload_Out.vL[1] -
    control_3dof_DW.obj.Z_2[1];
  control_3dof_DW.obj.Z_1[1] += (control_3dof_DW.obj.K_p[0] * control_3dof_B.e_1
    + control_3dof_DW.obj.Z_2[1]) * 0.05;
  control_3dof_DW.obj.Z_2[1] += (((control_3dof_DW.UnitDelay_DSTATE[1] +
    control_3dof_DW.obj.Z_3[1]) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.e_1) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_k) *
    0.05;
  control_3dof_DW.obj.Z_3[1] += (control_3dof_DW.obj.K_p[2] * control_3dof_B.e_1
    + control_3dof_DW.obj.K_v[2] * control_3dof_B.e_2_k) * 0.05;

  /* Sum: '<S18>/Sum' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  control_3dof_B.e_1_idx_1 = control_3dof_U.Traj_sp.pos_sp[1] -
    control_3dof_U.Payload_Out.pL[1];

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  control_3dof_B.e_2[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_1[1];

  /* MATLABSystem: '<S4>/Position 2nd ESO' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  UnitDelay: '<S4>/Unit Delay'
   * */
  control_3dof_B.e_1 = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_DW.obj.Z_1[2];
  control_3dof_B.e_2_k = control_3dof_U.Payload_Out.vL[2] -
    control_3dof_DW.obj.Z_2[2];
  control_3dof_DW.obj.Z_1[2] += (control_3dof_DW.obj.K_p[0] * control_3dof_B.e_1
    + control_3dof_DW.obj.Z_2[2]) * 0.05;
  control_3dof_DW.obj.Z_2[2] += ((((control_3dof_DW.UnitDelay_DSTATE[2] +
    control_3dof_DW.obj.Z_3[2]) + 9.8) + control_3dof_DW.obj.K_p[1] *
    control_3dof_B.e_1) + control_3dof_DW.obj.K_v[1] * control_3dof_B.e_2_k) *
    0.05;
  control_3dof_DW.obj.Z_3[2] += (control_3dof_DW.obj.K_p[2] * control_3dof_B.e_1
    + control_3dof_DW.obj.K_v[2] * control_3dof_B.e_2_k) * 0.05;

  /* Sum: '<S18>/Sum' incorporates:
   *  Inport: '<Root>/Payload_Out'
   *  Inport: '<Root>/Traj_sp'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   * */
  control_3dof_B.e_1 = control_3dof_U.Traj_sp.pos_sp[2] -
    control_3dof_U.Payload_Out.pL[2];

  /* MATLAB Function: '<S4>/Force Saturation && Disturbution' incorporates:
   *  Bias: '<S4>/-g'
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Gain: '<S4>/Gain'
   *  Inport: '<Root>/Payload_Out'
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Saturate: '<S18>/acc_limit'
   *  UnitDelay: '<S4>/Unit Delay'
   * */
  control_3dof_B.e_2_k = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_1[2];
  control_3dof_B.e_2[2] = control_3dof_B.e_2_k;

  /* '<S17>:1:4' q_1 = q_1 / norm(q_1); */
  control_3dof_B.F_trim_g = control_3dof_norm(control_3dof_B.e_2);

  /* '<S17>:1:5' w_1 = cross(q_1, (Payload_Out.vL - Payload_Out.v_1) / li); */
  /* '<S17>:1:8' q_2 = Payload_Out.pL - Payload_Out.p_2; */
  control_3dof_B.e_2[0] /= control_3dof_B.F_trim_g;
  control_3dof_B.F_var[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_1[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[0] = control_3dof_U.Payload_Out.pL[0] -
    control_3dof_U.Payload_Out.p_2[0];
  control_3dof_B.e_2[1] /= control_3dof_B.F_trim_g;
  control_3dof_B.F_var[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_1[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[1] = control_3dof_U.Payload_Out.pL[1] -
    control_3dof_U.Payload_Out.p_2[1];
  control_3dof_B.e_2[2] = control_3dof_B.e_2_k / control_3dof_B.F_trim_g;
  control_3dof_B.F_var[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_1[2]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.q_2[2] = control_3dof_U.Payload_Out.pL[2] -
    control_3dof_U.Payload_Out.p_2[2];
  control_3dof_B.rtb_control_in1_w[0] = control_3dof_B.e_2[1] *
    control_3dof_B.F_var[2] - control_3dof_B.F_var[1] * control_3dof_B.e_2[2];
  control_3dof_B.rtb_control_in1_w[1] = control_3dof_B.F_var[0] *
    control_3dof_B.e_2[2] - control_3dof_B.e_2[0] * control_3dof_B.F_var[2];
  control_3dof_B.rtb_control_in1_w[2] = control_3dof_B.e_2[0] *
    control_3dof_B.F_var[1] - control_3dof_B.F_var[0] * control_3dof_B.e_2[1];

  /* '<S17>:1:9' q_2 = q_2 / norm(q_2); */
  control_3dof_B.F_trim_g = control_3dof_norm(control_3dof_B.q_2);

  /* '<S17>:1:10' w_2 = cross(q_2, (Payload_Out.vL - Payload_Out.v_2) / li); */
  /* '<S17>:1:13' q_3 = Payload_Out.pL - Payload_Out.p_3; */
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
  control_3dof_B.rtb_control_in2_w[0] = control_3dof_B.q_2[1] *
    control_3dof_B.F_var[2] - control_3dof_B.F_var[1] * control_3dof_B.q_2[2];
  control_3dof_B.rtb_control_in2_w[1] = control_3dof_B.F_var[0] *
    control_3dof_B.q_2[2] - control_3dof_B.q_2[0] * control_3dof_B.F_var[2];
  control_3dof_B.rtb_control_in2_w[2] = control_3dof_B.q_2[0] *
    control_3dof_B.F_var[1] - control_3dof_B.F_var[0] * control_3dof_B.q_2[1];

  /* '<S17>:1:14' q_3 = q_3 / norm(q_3); */
  control_3dof_B.F_trim_g = control_3dof_norm(control_3dof_B.q_3);

  /* '<S17>:1:15' w_3 = cross(q_3, (Payload_Out.vL - Payload_Out.v_3) / li); */
  /* '<S17>:1:17' Q = [q_1 q_2 q_3]; */
  /* '<S17>:1:20' F_trim = mL * a_trim; */
  control_3dof_B.e_2_k = control_3dof_B.q_3[0] / control_3dof_B.F_trim_g;
  control_3dof_B.q_3[0] = control_3dof_B.e_2_k;
  control_3dof_B.F_var[0] = (control_3dof_U.Payload_Out.vL[0] -
    control_3dof_U.Payload_Out.v_3[0]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[0] = control_3dof_B.e_2[0];
  control_3dof_B.Q_new[3] = control_3dof_B.q_2[0];
  control_3dof_B.Q_new[6] = control_3dof_B.e_2_k;
  control_3dof_B.F_trim[0] = CONTROL_PARAM.MASS_LOAD * -control_3dof_DW.obj.Z_3
    [0];
  control_3dof_B.e_2_k = control_3dof_B.q_3[1] / control_3dof_B.F_trim_g;
  control_3dof_B.q_3[1] = control_3dof_B.e_2_k;
  control_3dof_B.F_var[1] = (control_3dof_U.Payload_Out.vL[1] -
    control_3dof_U.Payload_Out.v_3[1]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[1] = control_3dof_B.e_2[1];
  control_3dof_B.Q_new[4] = control_3dof_B.q_2[1];
  control_3dof_B.Q_new[7] = control_3dof_B.e_2_k;
  control_3dof_B.F_trim[1] = CONTROL_PARAM.MASS_LOAD * -control_3dof_DW.obj.Z_3
    [1];
  control_3dof_B.e_2_k = control_3dof_B.q_3[2] / control_3dof_B.F_trim_g;
  control_3dof_B.q_3[2] = control_3dof_B.e_2_k;
  control_3dof_B.F_var[2] = (control_3dof_U.Payload_Out.vL[2] -
    control_3dof_U.Payload_Out.v_3[2]) / CONTROL_PARAM.CABLE_LEN;
  control_3dof_B.Q_new[2] = control_3dof_B.e_2[2];
  control_3dof_B.Q_new[5] = control_3dof_B.q_2[2];
  control_3dof_B.Q_new[8] = control_3dof_B.e_2_k;
  control_3dof_B.F_trim[2] = (-control_3dof_DW.obj.Z_3[2] - 9.8) *
    CONTROL_PARAM.MASS_LOAD;
  control_3dof_B.rtb_control_in3_w[0] = control_3dof_B.q_3[1] *
    control_3dof_B.F_var[2] - control_3dof_B.F_var[1] * control_3dof_B.e_2_k;
  control_3dof_B.rtb_control_in3_w[1] = control_3dof_B.F_var[0] *
    control_3dof_B.e_2_k - control_3dof_B.q_3[0] * control_3dof_B.F_var[2];
  control_3dof_B.rtb_control_in3_w[2] = control_3dof_B.q_3[0] *
    control_3dof_B.F_var[1] - control_3dof_B.F_var[0] * control_3dof_B.q_3[1];

  /* '<S17>:1:21' F_way  = mL * a_way; */
  /* '<S17>:1:22' T_min = min_tension * ones(3, 1); */
  /* '<S17>:1:23' T_max = max_tension * ones(3, 1); */
  /* '<S17>:1:24' [vu, F_actual] = PWAS(F_trim, F_way, Q, T_min, T_max); */
  /* 'PWAS:3' if ~isvector(T_max) || ~isvector(T_min) */
  /* 'PWAS:7' N = size(Q, 2); */
  /* 'PWAS:8' W = -Q; */
  /* 'PWAS:9' A_wrench = [eye(N); -eye(N)] * pinv(W); */
  for (control_3dof_B.i = 0; control_3dof_B.i < 9; control_3dof_B.i++) {
    control_3dof_B.rtb_Q_new_m[control_3dof_B.i] =
      -control_3dof_B.Q_new[control_3dof_B.i];
  }

  control_3dof_pinv(control_3dof_B.rtb_Q_new_m, control_3dof_B.dv);
  for (control_3dof_B.i = 0; control_3dof_B.i < 3; control_3dof_B.i++) {
    control_3dof_B.F_trim_g = control_3dof_B.dv[3 * control_3dof_B.i + 1];
    control_3dof_B.F_trim_g1 = control_3dof_B.dv[3 * control_3dof_B.i];
    control_3dof_B.upper = control_3dof_B.dv[3 * control_3dof_B.i + 2];
    for (control_3dof_B.c_k = 0; control_3dof_B.c_k < 6; control_3dof_B.c_k++) {
      control_3dof_B.A_wrench_tmp[control_3dof_B.c_k + 6 * control_3dof_B.i] =
        ((real_T)a[control_3dof_B.c_k + 6] * control_3dof_B.F_trim_g +
         control_3dof_B.F_trim_g1 * (real_T)a[control_3dof_B.c_k]) + (real_T)
        a[control_3dof_B.c_k + 12] * control_3dof_B.upper;
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
  /* 'PWAS:50' retval = all(A_w * F_sp <= b_w); */
  control_3dof_B.F_trim_g = control_3dof_B.F_trim[1];
  control_3dof_B.F_trim_g1 = control_3dof_B.F_trim[0];
  control_3dof_B.upper = control_3dof_B.F_trim[2];
  for (control_3dof_B.i = 0; control_3dof_B.i < 6; control_3dof_B.i++) {
    control_3dof_B.B_tmp[control_3dof_B.i] =
      (control_3dof_B.A_wrench_tmp[control_3dof_B.i + 6] *
       control_3dof_B.F_trim_g + control_3dof_B.A_wrench_tmp[control_3dof_B.i] *
       control_3dof_B.F_trim_g1) + control_3dof_B.A_wrench_tmp[control_3dof_B.i
      + 12] * control_3dof_B.upper;
  }

  control_3dof_B.FixPtRelationalOperator = true;
  control_3dof_B.i = 0;
  exitg1 = false;
  while ((!exitg1) && (control_3dof_B.i < 6)) {
    if (!(control_3dof_B.B_tmp[control_3dof_B.i] <=
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

    /* Sum: '<S21>/Add' incorporates:
     *  DiscreteIntegrator: '<S21>/Discrete-Time Integrator'
     *  Gain: '<S21>/KP'
     */
    control_3dof_B.F_trim_g = CONTROL_PARAM.KP * control_3dof_B.e_1_idx_0 +
      control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0];

    /* Saturate: '<S18>/vel_limit' */
    if (control_3dof_B.F_trim_g > 5.0) {
      control_3dof_B.F_trim_g = 5.0;
    } else if (control_3dof_B.F_trim_g < -5.0) {
      control_3dof_B.F_trim_g = -5.0;
    }

    /* Sum: '<S18>/Sum1' incorporates:
     *  Gain: '<S18>/KV'
     *  Inport: '<Root>/Traj_sp'
     *  Saturate: '<S18>/vel_limit'
     *  Sum: '<S18>/Sum2'
     *  Sum: '<S18>/Sum3'
     */
    control_3dof_B.F_trim_g = ((control_3dof_B.F_trim_g +
      control_3dof_U.Traj_sp.vel_sp[0]) - control_3dof_U.Payload_Out.vL[0]) *
      CONTROL_PARAM.KV + control_3dof_U.Traj_sp.acc_ff[0];

    /* Saturate: '<S18>/acc_limit' */
    if (control_3dof_B.F_trim_g > 3.0) {
      control_3dof_B.F_trim_g = 3.0;
    } else if (control_3dof_B.F_trim_g < -3.0) {
      control_3dof_B.F_trim_g = -3.0;
    }

    control_3dof_B.F_var[0] = CONTROL_PARAM.MASS_LOAD * control_3dof_B.F_trim_g;
    control_3dof_DW.UnitDelay_DSTATE[1] = control_3dof_B.F_trim[1];

    /* Sum: '<S21>/Add' incorporates:
     *  DiscreteIntegrator: '<S21>/Discrete-Time Integrator'
     *  Gain: '<S21>/KP'
     *  Saturate: '<S18>/acc_limit'
     */
    control_3dof_B.F_trim_g = CONTROL_PARAM.KP * control_3dof_B.e_1_idx_1 +
      control_3dof_DW.DiscreteTimeIntegrator_DSTATE[1];

    /* Saturate: '<S18>/vel_limit' */
    if (control_3dof_B.F_trim_g > 5.0) {
      control_3dof_B.F_trim_g = 5.0;
    } else if (control_3dof_B.F_trim_g < -5.0) {
      control_3dof_B.F_trim_g = -5.0;
    }

    /* Sum: '<S18>/Sum1' incorporates:
     *  Gain: '<S18>/KV'
     *  Inport: '<Root>/Traj_sp'
     *  Saturate: '<S18>/vel_limit'
     *  Sum: '<S18>/Sum2'
     *  Sum: '<S18>/Sum3'
     */
    control_3dof_B.F_trim_g = ((control_3dof_B.F_trim_g +
      control_3dof_U.Traj_sp.vel_sp[1]) - control_3dof_U.Payload_Out.vL[1]) *
      CONTROL_PARAM.KV + control_3dof_U.Traj_sp.acc_ff[1];

    /* Saturate: '<S18>/acc_limit' */
    if (control_3dof_B.F_trim_g > 3.0) {
      control_3dof_B.F_trim_g = 3.0;
    } else if (control_3dof_B.F_trim_g < -3.0) {
      control_3dof_B.F_trim_g = -3.0;
    }

    control_3dof_B.F_var[1] = CONTROL_PARAM.MASS_LOAD * control_3dof_B.F_trim_g;
    control_3dof_DW.UnitDelay_DSTATE[2] = control_3dof_B.F_trim[2];

    /* Sum: '<S21>/Add' incorporates:
     *  DiscreteIntegrator: '<S21>/Discrete-Time Integrator'
     *  Gain: '<S21>/KP'
     *  Saturate: '<S18>/acc_limit'
     *  Sum: '<S18>/Sum'
     */
    control_3dof_B.F_trim_g = CONTROL_PARAM.KP * control_3dof_B.e_1 +
      control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2];

    /* Saturate: '<S18>/vel_limit' */
    if (control_3dof_B.F_trim_g > 5.0) {
      control_3dof_B.F_trim_g = 5.0;
    } else if (control_3dof_B.F_trim_g < -5.0) {
      control_3dof_B.F_trim_g = -5.0;
    }

    /* Sum: '<S18>/Sum1' incorporates:
     *  Gain: '<S18>/KV'
     *  Inport: '<Root>/Traj_sp'
     *  Saturate: '<S18>/vel_limit'
     *  Sum: '<S18>/Sum2'
     *  Sum: '<S18>/Sum3'
     */
    control_3dof_B.F_trim_g = ((control_3dof_B.F_trim_g +
      control_3dof_U.Traj_sp.vel_sp[2]) - control_3dof_U.Payload_Out.vL[2]) *
      CONTROL_PARAM.KV + control_3dof_U.Traj_sp.acc_ff[2];

    /* Saturate: '<S18>/acc_limit' */
    if (control_3dof_B.F_trim_g > 3.0) {
      control_3dof_B.F_trim_g = 3.0;
    } else if (control_3dof_B.F_trim_g < -3.0) {
      control_3dof_B.F_trim_g = -3.0;
    }

    control_3dof_B.F_var[2] = CONTROL_PARAM.MASS_LOAD * control_3dof_B.F_trim_g;
  } else {
    /* 'PWAS:15' else */
    /* 'PWAS:17' F_fixed = [0;0;0]; */
    control_3dof_DW.UnitDelay_DSTATE[0] = 0.0;
    control_3dof_DW.UnitDelay_DSTATE[1] = 0.0;
    control_3dof_DW.UnitDelay_DSTATE[2] = 0.0;

    /* 'PWAS:18' F_var = [F_trim(1);F_trim(2);F_trim(3)]; */
    control_3dof_B.F_var[0] = control_3dof_B.F_trim[0];
    control_3dof_B.F_var[1] = control_3dof_B.F_trim[1];
    control_3dof_B.F_var[2] = control_3dof_B.F_trim[2];
  }

  /* 'PWAS:21' alpha = PWAS_LOPT(F_fixed, F_var, A_wrench, b_wrench); */
  /* 'PWAS:28' A = A_w * F_var; */
  /* 'PWAS:29' B = b_w - A_w * F_fixed; */
  /* 'PWAS:31' upper = 1; */
  control_3dof_B.upper = 1.0;

  /* 'PWAS:32' lower = 0; */
  control_3dof_B.lower = 0.0;

  /* 'PWAS:33' for i = 1:length(A) */
  control_3dof_B.F_var_c = control_3dof_B.F_var[0];
  control_3dof_B.F_var_b = control_3dof_B.F_var[1];
  control_3dof_B.F_var_p = control_3dof_B.F_var[2];
  control_3dof_B.UnitDelay_DSTATE = control_3dof_DW.UnitDelay_DSTATE[0];
  control_3dof_B.UnitDelay_DSTATE_c = control_3dof_DW.UnitDelay_DSTATE[1];
  control_3dof_B.UnitDelay_DSTATE_f = control_3dof_DW.UnitDelay_DSTATE[2];
  for (control_3dof_B.i = 0; control_3dof_B.i < 6; control_3dof_B.i++) {
    control_3dof_B.A_wrench_tmp_c = control_3dof_B.A_wrench_tmp[control_3dof_B.i];
    control_3dof_B.F_trim_g1 = control_3dof_B.A_wrench_tmp_c *
      control_3dof_B.F_var_c;
    control_3dof_B.F_trim_g = control_3dof_B.A_wrench_tmp_c *
      control_3dof_B.UnitDelay_DSTATE;
    control_3dof_B.A_wrench_tmp_c = control_3dof_B.A_wrench_tmp[control_3dof_B.i
      + 6];
    control_3dof_B.F_trim_g1 += control_3dof_B.A_wrench_tmp_c *
      control_3dof_B.F_var_b;
    control_3dof_B.F_trim_g += control_3dof_B.A_wrench_tmp_c *
      control_3dof_B.UnitDelay_DSTATE_c;
    control_3dof_B.A_wrench_tmp_c = control_3dof_B.A_wrench_tmp[control_3dof_B.i
      + 12];
    control_3dof_B.F_trim_g1 += control_3dof_B.A_wrench_tmp_c *
      control_3dof_B.F_var_p;
    control_3dof_B.F_trim_g = control_3dof_B.b_wrench[control_3dof_B.i] -
      (control_3dof_B.A_wrench_tmp_c * control_3dof_B.UnitDelay_DSTATE_f +
       control_3dof_B.F_trim_g);

    /* 'PWAS:34' if A(i) > 0 */
    if (control_3dof_B.F_trim_g1 > 0.0) {
      /* 'PWAS:35' upper = min(upper, B(i) / A(i)); */
      control_3dof_B.upper = fmin(control_3dof_B.upper, control_3dof_B.F_trim_g /
        control_3dof_B.F_trim_g1);
    } else if (control_3dof_B.F_trim_g1 < 0.0) {
      /* 'PWAS:36' elseif A(i) < 0 */
      /* 'PWAS:37' lower = max(lower, B(i) / A(i)); */
      control_3dof_B.lower = fmax(control_3dof_B.lower, control_3dof_B.F_trim_g /
        control_3dof_B.F_trim_g1);
    }
  }

  /* 'PWAS:41' if lower > upper */
  if (!(control_3dof_B.lower > control_3dof_B.upper)) {
    /* 'PWAS:43' else */
    /* 'PWAS:44' alpha = upper; */
    control_3dof_B.alpha = control_3dof_B.upper;
  } else {
    /* 'PWAS:42' error('No Way'); */
  }

  /* 'PWAS:23' F_actual = F_fixed + alpha * F_var; */
  control_3dof_DW.UnitDelay_DSTATE[0] += control_3dof_B.alpha *
    control_3dof_B.F_var[0];
  control_3dof_DW.UnitDelay_DSTATE[1] += control_3dof_B.alpha *
    control_3dof_B.F_var[1];
  control_3dof_DW.UnitDelay_DSTATE[2] += control_3dof_B.alpha *
    control_3dof_B.F_var[2];

  /* 'PWAS:24' result = pinv(Q) * F_actual; */
  control_3dof_pinv(control_3dof_B.Q_new, control_3dof_B.dv);

  /* '<S17>:1:27' aL_sp = [0;0;-9.8] + (F_actual - F_trim) / mL; */
  for (control_3dof_B.i = 0; control_3dof_B.i < 3; control_3dof_B.i++) {
    control_3dof_B.F_var[control_3dof_B.i] = (control_3dof_B.dv[control_3dof_B.i
      + 3] * control_3dof_DW.UnitDelay_DSTATE[1] +
      control_3dof_B.dv[control_3dof_B.i] * control_3dof_DW.UnitDelay_DSTATE[0])
      + control_3dof_B.dv[control_3dof_B.i + 6] *
      control_3dof_DW.UnitDelay_DSTATE[2];
    control_3dof_B.F_trim[control_3dof_B.i] =
      (control_3dof_DW.UnitDelay_DSTATE[control_3dof_B.i] -
       control_3dof_B.F_trim[control_3dof_B.i]) / CONTROL_PARAM.MASS_LOAD +
      b[control_3dof_B.i];
  }

  /* '<S17>:1:29' F_actual_div_m = F_actual / mL; */
  control_3dof_DW.UnitDelay_DSTATE[0] /= CONTROL_PARAM.MASS_LOAD;
  control_3dof_DW.UnitDelay_DSTATE[1] /= CONTROL_PARAM.MASS_LOAD;
  control_3dof_DW.UnitDelay_DSTATE[2] /= CONTROL_PARAM.MASS_LOAD;

  /* '<S17>:1:31' margin = margin_lp(F_trim, -Q, T_min, T_max); */
  /* 'margin_lp:3' N = size(W, 2); */
  /* 'margin_lp:5' A_wrench = [eye(N); -eye(N)] * pinv(W); */
  /* 'margin_lp:6' b_wrench = [T_max; -T_min]; */
  /* 'margin_lp:9' margin = min((b_wrench - A_wrench * center) ./ vecnorm(A_wrench, 2, 2)); */
  for (control_3dof_B.i = 0; control_3dof_B.i < 6; control_3dof_B.i++) {
    control_3dof_B.xv[0] = control_3dof_B.A_wrench_tmp[control_3dof_B.i];
    control_3dof_B.xv[1] = control_3dof_B.A_wrench_tmp[control_3dof_B.i + 6];
    control_3dof_B.xv[2] = control_3dof_B.A_wrench_tmp[control_3dof_B.i + 12];
    control_3dof_B.B[control_3dof_B.i] = control_3dof_norm(control_3dof_B.xv);
  }

  control_3dof_B.b_wrench[0] = (CONTROL_PARAM.TENSION_MAX -
    control_3dof_B.B_tmp[0]) / control_3dof_B.B[0];
  control_3dof_B.b_wrench[3] = (-CONTROL_PARAM.TENSION_MIN -
    control_3dof_B.B_tmp[3]) / control_3dof_B.B[3];
  control_3dof_B.b_wrench[1] = (CONTROL_PARAM.TENSION_MAX -
    control_3dof_B.B_tmp[1]) / control_3dof_B.B[1];
  control_3dof_B.b_wrench[4] = (-CONTROL_PARAM.TENSION_MIN -
    control_3dof_B.B_tmp[4]) / control_3dof_B.B[4];
  control_3dof_B.b_wrench[2] = (CONTROL_PARAM.TENSION_MAX -
    control_3dof_B.B_tmp[2]) / control_3dof_B.B[2];
  control_3dof_B.b_wrench[5] = (-CONTROL_PARAM.TENSION_MIN -
    control_3dof_B.B_tmp[5]) / control_3dof_B.B[5];
  if (!rtIsNaN(control_3dof_B.b_wrench[0])) {
    control_3dof_B.i = 1;
  } else {
    control_3dof_B.i = 0;
    control_3dof_B.c_k = 2;
    exitg1 = false;
    while ((!exitg1) && (control_3dof_B.c_k < 7)) {
      if (!rtIsNaN(control_3dof_B.b_wrench[control_3dof_B.c_k - 1])) {
        control_3dof_B.i = control_3dof_B.c_k;
        exitg1 = true;
      } else {
        control_3dof_B.c_k++;
      }
    }
  }

  if (control_3dof_B.i == 0) {
    control_3dof_Y.state.margin = control_3dof_B.b_wrench[0];
  } else {
    control_3dof_B.alpha = control_3dof_B.b_wrench[control_3dof_B.i - 1];
    for (control_3dof_B.c_k = control_3dof_B.i + 1; control_3dof_B.c_k < 7;
         control_3dof_B.c_k++) {
      control_3dof_B.F_trim_g1 = control_3dof_B.b_wrench[control_3dof_B.c_k - 1];
      if (control_3dof_B.alpha > control_3dof_B.F_trim_g1) {
        control_3dof_B.alpha = control_3dof_B.F_trim_g1;
      }
    }

    control_3dof_Y.state.margin = control_3dof_B.alpha;
  }

  /* '<S17>:1:34' control_in1.timestamp = uint64(1); */
  /* '<S17>:1:35' control_in1.vforce = vu(1) * q_1; */
  /* '<S17>:1:36' control_in1.q_sp = Q_sp(:, 1); */
  /* '<S17>:1:37' control_in1.q = q_1; */
  /* '<S17>:1:38' control_in1.w = w_1; */
  /* '<S17>:1:39' control_in1.ai_sp = aL_sp; */
  /* '<S17>:1:40' control_in1.pi = Payload_Out.p_1; */
  /* '<S17>:1:41' control_in1.vi = Payload_Out.v_1; */
  /* '<S17>:1:43' control_in2.timestamp = uint64(1); */
  /* '<S17>:1:44' control_in2.vforce = vu(2) * q_2; */
  control_3dof_B.xv[0] = control_3dof_B.F_var[0] * control_3dof_B.e_2[0];
  control_3dof_B.rtb_control_in2_vforce[0] = control_3dof_B.q_2[0] *
    control_3dof_B.F_var[1];
  control_3dof_B.xv[1] = control_3dof_B.F_var[0] * control_3dof_B.e_2[1];
  control_3dof_B.rtb_control_in2_vforce[1] = control_3dof_B.F_var[1] *
    control_3dof_B.q_2[1];
  control_3dof_B.xv[2] = control_3dof_B.F_var[0] * control_3dof_B.e_2[2];
  control_3dof_B.rtb_control_in2_vforce[2] = control_3dof_B.F_var[1] *
    control_3dof_B.q_2[2];

  /* '<S17>:1:45' control_in2.q_sp = Q_sp(:, 2); */
  /* '<S17>:1:46' control_in2.q = q_2; */
  /* '<S17>:1:47' control_in2.w = w_2; */
  /* '<S17>:1:48' control_in2.ai_sp = aL_sp; */
  /* '<S17>:1:49' control_in2.pi = Payload_Out.p_2; */
  /* '<S17>:1:50' control_in2.vi = Payload_Out.v_2; */
  /* '<S17>:1:52' control_in3.timestamp = uint64(1); */
  /* '<S17>:1:53' control_in3.vforce = vu(3) * q_3; */
  control_3dof_B.F_var_c = control_3dof_B.F_var[2];
  control_3dof_B.F_var[0] = control_3dof_B.F_var_c * control_3dof_B.q_3[0];
  control_3dof_B.F_var[1] = control_3dof_B.F_var_c * control_3dof_B.q_3[1];
  control_3dof_B.F_var[2] = control_3dof_B.F_var_c * control_3dof_B.e_2_k;

  /* MATLAB Function: '<S6>/Vertical Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  Switch: '<S20>/Switch'
   */
  /* '<S17>:1:54' control_in3.q_sp = Q_sp(:, 3); */
  /* '<S17>:1:55' control_in3.q = q_3; */
  /* '<S17>:1:56' control_in3.w = w_3; */
  /* '<S17>:1:57' control_in3.ai_sp = aL_sp; */
  /* '<S17>:1:58' control_in3.pi = Payload_Out.p_3; */
  /* '<S17>:1:59' control_in3.vi = Payload_Out.v_3; */
  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[0],
    control_3dof_B.e_2, control_3dof_B.rtb_control_in1_w, control_3dof_B.F_trim,
    control_3dof_B.f_vertical, control_3dof_B.dot_err_g, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl);

  /* MATLAB Function: '<S10>/Vertical Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  Switch: '<S20>/Switch'
   */
  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[3],
    control_3dof_B.q_2, control_3dof_B.rtb_control_in2_w, control_3dof_B.F_trim,
    control_3dof_B.f_parallel, control_3dof_B.dot_err_c, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_m);

  /* MATLAB Function: '<S14>/Vertical Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   *  Switch: '<S20>/Switch'
   */
  control_3dof_VerticalControl(&control_3dof_DW.Delay_DSTATE[6],
    control_3dof_B.q_3, control_3dof_B.rtb_control_in3_w, control_3dof_B.F_trim,
    control_3dof_Y.force_sp1, control_3dof_B.dot_err, control_3dof_B.state_f,
    CONTROL_PARAM.KQ, CONTROL_PARAM.KW, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_VerticalControl_h);

  /* UnitDelay: '<S3>/Unit Delay' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   */
  control_3dof_MATLABSystem(control_3dof_DW.UnitDelay_DSTATE_l,
    control_3dof_U.Payload_Out.p_3, control_3dof_U.Payload_Out.v_3,
    &control_3dof_B.MATLABSystem_a, &control_3dof_DW.MATLABSystem_a);

  /* MATLAB Function: '<S13>/Parallel Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   */
  control_3dof_ParallelControl(control_3dof_B.q_3,
    control_3dof_B.rtb_control_in3_w, control_3dof_B.F_var,
    control_3dof_B.F_trim, control_3dof_Y.force_sp2, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_p);

  /* Gain: '<S3>/Gain1' incorporates:
   *  Gain: '<S1>/Gain1'
   *  Gain: '<S2>/Gain1'
   */
  control_3dof_B.F_trim_g = 1.0 / CONTROL_PARAM.MASS_UAV;

  /* Sum: '<S3>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S14>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S3>/Gain'
   *  MATLABSystem: '<S3>/MATLAB System'
   *  Sum: '<S14>/Sum'
   *  Sum: '<S1>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem_a.MATLABSystem_o3[0] + control_3dof_Y.force_sp2
    [0]) + (control_3dof_Y.force_sp1[0] +
            control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[0]);
  control_3dof_Y.force_sp3[0] = control_3dof_B.e_2_k;

  /* Gain: '<S3>/Gain1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  Sum: '<S3>/Sum'
   *  UnitDelay: '<S3>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_l[0] = (control_3dof_B.e_2_k -
    control_3dof_B.F_var[0]) * control_3dof_B.F_trim_g;

  /* Sum: '<S3>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S14>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S3>/Gain'
   *  MATLABSystem: '<S3>/MATLAB System'
   *  Sum: '<S14>/Sum'
   *  Sum: '<S1>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem_a.MATLABSystem_o3[1] + control_3dof_Y.force_sp2
    [1]) + (control_3dof_Y.force_sp1[1] +
            control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[1]);
  control_3dof_Y.force_sp3[1] = control_3dof_B.e_2_k;

  /* Gain: '<S3>/Gain1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  Sum: '<S3>/Sum'
   *  UnitDelay: '<S3>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_l[1] = (control_3dof_B.e_2_k -
    control_3dof_B.F_var[1]) * control_3dof_B.F_trim_g;

  /* Sum: '<S3>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S14>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S3>/Gain'
   *  MATLABSystem: '<S3>/MATLAB System'
   *  Sum: '<S14>/Sum'
   *  Sum: '<S1>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem_a.MATLABSystem_o3[2] + control_3dof_Y.force_sp2
    [2]) + (control_3dof_Y.force_sp1[2] +
            control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[2]);
  control_3dof_Y.force_sp3[2] = control_3dof_B.e_2_k;

  /* Gain: '<S3>/Gain1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  Sum: '<S3>/Sum'
   *  UnitDelay: '<S3>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_l[2] = (control_3dof_B.e_2_k -
    control_3dof_B.F_var[2]) * control_3dof_B.F_trim_g;

  /* UnitDelay: '<S2>/Unit Delay' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   */
  control_3dof_MATLABSystem(control_3dof_DW.UnitDelay_DSTATE_a,
    control_3dof_U.Payload_Out.p_2, control_3dof_U.Payload_Out.v_2,
    &control_3dof_B.MATLABSystem_g, &control_3dof_DW.MATLABSystem_g);

  /* MATLAB Function: '<S9>/Parallel Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   */
  control_3dof_ParallelControl(control_3dof_B.q_2,
    control_3dof_B.rtb_control_in2_w, control_3dof_B.rtb_control_in2_vforce,
    control_3dof_B.F_trim, control_3dof_Y.force_sp1, CONTROL_PARAM.CABLE_LEN,
    CONTROL_PARAM.MASS_UAV, &control_3dof_B.sf_ParallelControl_k);

  /* Sum: '<S2>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S2>/Gain'
   *  MATLABSystem: '<S2>/MATLAB System'
   *  Sum: '<S10>/Sum'
   *  Sum: '<S1>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem_g.MATLABSystem_o3[0] + control_3dof_Y.force_sp1
    [0]) + (control_3dof_B.f_parallel[0] +
            control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0]);
  control_3dof_Y.force_sp2[0] = control_3dof_B.e_2_k;

  /* Gain: '<S2>/Gain1' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Sum: '<S2>/Sum'
   *  UnitDelay: '<S2>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_a[0] = (control_3dof_B.e_2_k -
    control_3dof_B.rtb_control_in2_vforce[0]) * control_3dof_B.F_trim_g;

  /* Sum: '<S2>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S2>/Gain'
   *  MATLABSystem: '<S2>/MATLAB System'
   *  Sum: '<S10>/Sum'
   *  Sum: '<S1>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem_g.MATLABSystem_o3[1] + control_3dof_Y.force_sp1
    [1]) + (control_3dof_B.f_parallel[1] +
            control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[1]);
  control_3dof_Y.force_sp2[1] = control_3dof_B.e_2_k;

  /* Gain: '<S2>/Gain1' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Sum: '<S2>/Sum'
   *  UnitDelay: '<S2>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_a[1] = (control_3dof_B.e_2_k -
    control_3dof_B.rtb_control_in2_vforce[1]) * control_3dof_B.F_trim_g;

  /* Sum: '<S2>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S10>/Discrete-Time Integrator'
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S2>/Gain'
   *  MATLABSystem: '<S2>/MATLAB System'
   *  Sum: '<S10>/Sum'
   *  Sum: '<S1>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem_g.MATLABSystem_o3[2] + control_3dof_Y.force_sp1
    [2]) + (control_3dof_B.f_parallel[2] +
            control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2]);
  control_3dof_Y.force_sp2[2] = control_3dof_B.e_2_k;

  /* Gain: '<S2>/Gain1' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Sum: '<S2>/Sum'
   *  UnitDelay: '<S2>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_a[2] = (control_3dof_B.e_2_k -
    control_3dof_B.rtb_control_in2_vforce[2]) * control_3dof_B.F_trim_g;

  /* UnitDelay: '<S1>/Unit Delay' incorporates:
   *  BusCreator generated from: '<S4>/Force Saturation && Disturbution'
   *  Inport: '<Root>/Payload_Out'
   */
  control_3dof_MATLABSystem(control_3dof_DW.UnitDelay_DSTATE_f,
    control_3dof_U.Payload_Out.p_1, control_3dof_U.Payload_Out.v_1,
    &control_3dof_B.MATLABSystem, &control_3dof_DW.MATLABSystem);

  /* MATLAB Function: '<S5>/Parallel Control' incorporates:
   *  MATLAB Function: '<S4>/Force Saturation && Disturbution'
   */
  control_3dof_ParallelControl(control_3dof_B.e_2,
    control_3dof_B.rtb_control_in1_w, control_3dof_B.xv, control_3dof_B.F_trim,
    control_3dof_B.f_parallel, CONTROL_PARAM.CABLE_LEN, CONTROL_PARAM.MASS_UAV,
    &control_3dof_B.sf_ParallelControl);

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  DiscreteIntegrator: '<S22>/Accumulator'
   *  Outport: '<Root>/state'
   */
  control_3dof_Y.state.timestamp = control_3dof_DW.Accumulator_DSTATE;

  /* Sum: '<S1>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S1>/Gain'
   *  MATLABSystem: '<S1>/MATLAB System'
   *  Sum: '<S1>/Sum'
   *  Sum: '<S6>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem.MATLABSystem_o3[0] + control_3dof_B.f_parallel[0])
    + (control_3dof_B.f_vertical[0] +
       control_3dof_DW.DiscreteTimeIntegrator_DSTAT_pb[0]);
  control_3dof_Y.force_sp1[0] = control_3dof_B.e_2_k;

  /* Gain: '<S1>/Gain1' incorporates:
   *  Sum: '<S1>/Sum'
   *  UnitDelay: '<S1>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_f[0] = (control_3dof_B.e_2_k -
    control_3dof_B.xv[0]) * control_3dof_B.F_trim_g;

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Outport: '<Root>/state'
   * */
  control_3dof_Y.state.dL[0] = control_3dof_DW.obj.Z_3[0];

  /* Sum: '<S1>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S1>/Gain'
   *  MATLABSystem: '<S1>/MATLAB System'
   *  Sum: '<S1>/Sum'
   *  Sum: '<S6>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem.MATLABSystem_o3[1] + control_3dof_B.f_parallel[1])
    + (control_3dof_B.f_vertical[1] +
       control_3dof_DW.DiscreteTimeIntegrator_DSTAT_pb[1]);
  control_3dof_Y.force_sp1[1] = control_3dof_B.e_2_k;

  /* Gain: '<S1>/Gain1' incorporates:
   *  Sum: '<S1>/Sum'
   *  UnitDelay: '<S1>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_f[1] = (control_3dof_B.e_2_k -
    control_3dof_B.xv[1]) * control_3dof_B.F_trim_g;

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Outport: '<Root>/state'
   * */
  control_3dof_Y.state.dL[1] = control_3dof_DW.obj.Z_3[1];

  /* Sum: '<S1>/Sum1' incorporates:
   *  DiscreteIntegrator: '<S6>/Discrete-Time Integrator'
   *  Gain: '<S1>/Gain'
   *  MATLABSystem: '<S1>/MATLAB System'
   *  Sum: '<S1>/Sum'
   *  Sum: '<S6>/Sum'
   */
  control_3dof_B.e_2_k = (-CONTROL_PARAM.MASS_UAV *
    control_3dof_B.MATLABSystem.MATLABSystem_o3[2] + control_3dof_B.f_parallel[2])
    + (control_3dof_B.f_vertical[2] +
       control_3dof_DW.DiscreteTimeIntegrator_DSTAT_pb[2]);
  control_3dof_Y.force_sp1[2] = control_3dof_B.e_2_k;

  /* Gain: '<S1>/Gain1' incorporates:
   *  Sum: '<S1>/Sum'
   *  UnitDelay: '<S1>/Unit Delay'
   */
  control_3dof_DW.UnitDelay_DSTATE_f[2] = (control_3dof_B.e_2_k -
    control_3dof_B.xv[2]) * control_3dof_B.F_trim_g;

  /* BusCreator generated from: '<Root>/state' incorporates:
   *  MATLABSystem: '<S4>/Position 2nd ESO'
   *  Outport: '<Root>/state'
   * */
  control_3dof_Y.state.dL[2] = control_3dof_DW.obj.Z_3[2];

  /* Update for UnitDelay: '<S23>/Delay Input1' incorporates:
   *  Inport: '<Root>/Dir_sp'
   *
   * Block description for '<S23>/Delay Input1':
   *
   *  Store in Global RAM
   */
  control_3dof_DW.DelayInput1_DSTATE = control_3dof_U.Dir_sp.timestamp;

  /* Update for Delay: '<S20>/Delay' */
  control_3dof_DW.icLoad = false;

  /* Update for DiscreteIntegrator: '<S21>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S21>/KPI'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[0] += CONTROL_PARAM.KPI *
    control_3dof_B.e_1_idx_0 * 0.05;

  /* Update for DiscreteIntegrator: '<S14>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S14>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[0] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err[0] * 0.05;

  /* Update for DiscreteIntegrator: '<S10>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S10>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[0] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_c[0] * 0.05;

  /* Update for DiscreteIntegrator: '<S6>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S6>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTAT_pb[0] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[0] * 0.05;

  /* Update for DiscreteIntegrator: '<S21>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S21>/KPI'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[1] += CONTROL_PARAM.KPI *
    control_3dof_B.e_1_idx_1 * 0.05;

  /* Update for DiscreteIntegrator: '<S14>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S14>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[1] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err[1] * 0.05;

  /* Update for DiscreteIntegrator: '<S10>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S10>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[1] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_c[1] * 0.05;

  /* Update for DiscreteIntegrator: '<S6>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S6>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTAT_pb[1] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[1] * 0.05;

  /* Update for DiscreteIntegrator: '<S21>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S21>/KPI'
   *  Sum: '<S18>/Sum'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE[2] += CONTROL_PARAM.KPI *
    control_3dof_B.e_1 * 0.05;

  /* Update for DiscreteIntegrator: '<S14>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S14>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_p[2] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err[2] * 0.05;

  /* Update for DiscreteIntegrator: '<S10>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S10>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTATE_l[2] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_c[2] * 0.05;

  /* Update for DiscreteIntegrator: '<S6>/Discrete-Time Integrator' incorporates:
   *  Gain: '<S6>/Gain'
   */
  control_3dof_DW.DiscreteTimeIntegrator_DSTAT_pb[2] += CONTROL_PARAM.KQI *
    control_3dof_B.dot_err_g[2] * 0.05;

  /* Switch: '<S22>/Switch' incorporates:
   *  Constant: '<S22>/Period [ms]'
   *  SampleTimeMath: '<S22>/Ts [ms]'
   *
   * About '<S22>/Ts [ms]':
   *  y = K where K = ( w * Ts )
   *   */
  if (CONTROL_EXPORT.period > 0) {
    control_3dof_B.i = CONTROL_EXPORT.period;
  } else {
    control_3dof_B.i = 50;
  }

  /* Update for DiscreteIntegrator: '<S22>/Accumulator' incorporates:
   *  DataTypeConversion: '<S22>/Data Type Conversion'
   *  Switch: '<S22>/Switch'
   */
  control_3dof_DW.Accumulator_DSTATE += (uint64_T)control_3dof_B.i;
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
  /* 'pos_2nd_eso:43' setProperties(obj,nargin,varargin{:}) */
  /* 'pos_2nd_eso:13' (3, 1) */
  /* 'pos_2nd_eso:13' K_p */
  /* 'pos_2nd_eso:15' (3, 1) */
  /* 'pos_2nd_eso:15' K_v */
  control_3dof_DW.obj.isInitialized = 1;

  /*         %% Common functions */
  /*  Perform one-time calculations, such as computing constants */
  /* 'pos_2nd_eso:51' obj.Z_1 = pos; */
  /* 'pos_2nd_eso:52' obj.Z_2 = vel; */
  /* 'pos_2nd_eso:53' obj.Z_3 = zeros(size(pos)); */
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
  /* 'pos_2nd_eso:54' obj.is_init = false; */
  /* 'pos_2nd_eso:36' logical */
  /* 'pos_2nd_eso:36' (1, 1) */
  /* 'pos_2nd_eso:36' is_init */
  /*  Initialize / reset internal or discrete properties */
  /* 'pos_2nd_eso:92' obj.is_init = false; */
  /* 'pos_2nd_eso:36' logical */
  /* 'pos_2nd_eso:36' (1, 1) */
  /* 'pos_2nd_eso:36' is_init */
  control_3dof_DW.obj.is_init = false;
  control_3dof_MATLABSystem_Init(&control_3dof_DW.MATLABSystem_a);
  control_3dof_MATLABSystem_Init(&control_3dof_DW.MATLABSystem_g);
  control_3dof_MATLABSystem_Init(&control_3dof_DW.MATLABSystem);
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
