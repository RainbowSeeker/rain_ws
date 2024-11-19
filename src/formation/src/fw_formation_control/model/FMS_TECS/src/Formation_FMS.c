/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Formation_FMS.c
 *
 * Code generated for Simulink model 'Formation_FMS'.
 *
 * Model version                  : 1.259
 * Simulink Coder version         : 9.8 (R2022b) 13-May-2022
 * C/C++ source code generated on : Fri Aug 16 11:40:06 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Formation_FMS.h"
#include "rtwtypes.h"
#include <math.h>
#include <string.h>
#include "norm_FqL2z8vo.h"
#include "rt_nonfinite.h"
#include "rt_atan2f_snf.h"
#include "zero_crossing_types.h"

/* Named constants for Chart: '<Root>/FMS State Machine' */
#define Formation_FMS_IN_Follower      ((uint8_T)1U)
#define Formation_FMS_IN_FormAssemble  ((uint8_T)1U)
#define Formation_FMS_IN_FormDisband   ((uint8_T)2U)
#define Formation_FMS_IN_FormMission   ((uint8_T)3U)
#define Formation_FMS_IN_Formation     ((uint8_T)1U)
#define Formation_FMS_IN_Hold          ((uint8_T)1U)
#define Formation_FMS_IN_InvalidMode   ((uint8_T)4U)
#define Formation_FMS_IN_NextWP        ((uint8_T)2U)
#define Formation_FMS_IN_Standby       ((uint8_T)2U)
#define Formation_FMS_IN_WaitForUpdate ((uint8_T)3U)
#define Formation_FMS_IN_Waypoint      ((uint8_T)4U)
#define Formation_FMS_IN_Waypoint_g    ((uint8_T)3U)
#define Formation_FM_IN_NO_ACTIVE_CHILD ((uint8_T)0U)

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         (*((rtm)->errorStatus))
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    (*((rtm)->errorStatus) = (val))
#endif

#ifndef rtmGetErrorStatusPointer
#define rtmGetErrorStatusPointer(rtm)  (rtm)->errorStatus
#endif

#ifndef rtmSetErrorStatusPointer
#define rtmSetErrorStatusPointer(rtm, val) ((rtm)->errorStatus = (val))
#endif

/* Forward declaration for local functions */
static void Formation_exit_internal_Vehicle(DW_Formation_FMS_f_T *localDW);
static real_T Formation_FMS_atan3_9FIHpp9F(real_T x, real_T y, real_T x0, real_T
  b_y0, B_Formation_FMS_c_T *localB);
static void Formation_FMS_Dubins(const captured_var_Formation_FMS_T *PhiMaximum,
  const captured_var_Formation_FMS_T *rad2deg, real_T xs, real_T ys, real_T
  psi_s, real_T xf, real_T yf, real_T psi_f, real_T v, real_T xts[4], real_T
  yts[4], real_T xtf[4], real_T ytf[4], real_T cs[4], real_T cf[4], real_T lt[4],
  real_T l[4], real_T pos[4], real_T pof[4], real_T *r, real_T m_data[], int32_T
  *m_size, real_T n_data[], int32_T *n_size, B_Formation_FMS_c_T *localB);
static void Format_enter_internal_Formation(const INS_Out_Bus
  *BusConversion_InsertedFor_FMS_c, const Formation_Cross_Bus
  *BusConversion_InsertedFor_FMS_p, B_Formation_FMS_c_T *localB,
  DW_Formation_FMS_f_T *localDW);

/*
 * Output and update for atomic system:
 *    '<S30>/NearbyRefWP'
 *    '<S163>/NearbyRefWP'
 */
void Formation_FMS_NearbyRefWP(const real32_T rtu_P2[2], const real32_T rtu_P3[2],
  real32_T rtu_L1, real32_T rty_P[2], real32_T *rty_d)
{
  real32_T dis;
  real32_T y_idx_0;
  dis = rtu_P2[0] - rtu_P3[0];
  y_idx_0 = dis * dis;
  dis = rtu_P2[1] - rtu_P3[1];
  dis = (real32_T)sqrt(dis * dis + y_idx_0);
  if (dis <= rtu_L1) {
    *rty_d = dis;
    rty_P[0] = rtu_P2[0];
    rty_P[1] = rtu_P2[1];
  } else {
    *rty_d = -1.0F;
    rty_P[0] = 0.0F;
    rty_P[1] = 0.0F;
  }
}

/*
 * Output and update for atomic system:
 *    '<S30>/OutRegionRegWP'
 *    '<S163>/OutRegionRegWP'
 */
void Formation_FMS_OutRegionRegWP(const real32_T rtu_P1[2], const real32_T
  rtu_P2[2], const real32_T rtu_P3[2], real32_T rty_P[2])
{
  real32_T a;
  real32_T u;
  real32_T y_idx_0;
  rty_P[0] = rtu_P2[0] - rtu_P1[0];
  u = (rtu_P3[0] - rtu_P1[0]) * rty_P[0];
  y_idx_0 = rty_P[0] * rty_P[0];
  rty_P[1] = rtu_P2[1] - rtu_P1[1];
  u = ((rtu_P3[1] - rtu_P1[1]) * rty_P[1] + u) / (rty_P[1] * rty_P[1] + y_idx_0);
  if (u <= 0.0F) {
    rty_P[0] = rtu_P1[0];
    rty_P[1] = rtu_P1[1];
  } else if (u >= 1.0F) {
    rty_P[0] = rtu_P2[0];
    rty_P[1] = rtu_P2[1];
  } else {
    a = (u * rty_P[0] + rtu_P1[0]) - rtu_P3[0];
    y_idx_0 = a * a;
    a = (u * rty_P[1] + rtu_P1[1]) - rtu_P3[1];
    u += (real32_T)sqrt(a * a + y_idx_0) * 0.5774F / (real32_T)sqrt(rty_P[0] *
      rty_P[0] + rty_P[1] * rty_P[1]);
    if (u > 1.0F) {
      u = 1.0F;
    }

    rty_P[0] = u * rty_P[0] + rtu_P1[0];
    rty_P[1] = u * rty_P[1] + rtu_P1[1];
  }
}

/*
 * Output and update for atomic system:
 *    '<S30>/SearchL1RefWP'
 *    '<S163>/SearchL1RefWP'
 */
void Formation_FMS_SearchL1RefWP(const real32_T rtu_P1[2], const real32_T
  rtu_P2[2], const real32_T rtu_P3[2], real32_T rtu_L1, real32_T rty_P[2],
  real32_T *rty_u)
{
  real32_T A;
  real32_T B;
  real32_T D;
  real32_T a_tmp;
  real32_T b_a_tmp;
  real32_T u;
  real32_T u1_tmp;
  boolean_T guard1;
  a_tmp = rtu_P2[0] - rtu_P1[0];
  b_a_tmp = rtu_P2[1] - rtu_P1[1];
  A = a_tmp * a_tmp + b_a_tmp * b_a_tmp;
  B = ((rtu_P1[0] - rtu_P3[0]) * a_tmp + (rtu_P1[1] - rtu_P3[1]) * b_a_tmp) *
    2.0F;
  D = B * B - (((((rtu_P3[0] * rtu_P3[0] + rtu_P3[1] * rtu_P3[1]) + rtu_P1[0] *
                  rtu_P1[0]) + rtu_P1[1] * rtu_P1[1]) - (rtu_P3[0] * rtu_P1[0] +
    rtu_P3[1] * rtu_P1[1]) * 2.0F) - rtu_L1 * rtu_L1) * (4.0F * A);
  u = -1.0F;
  rty_P[0] = 0.0F;
  rty_P[1] = 0.0F;
  guard1 = false;
  if (D > 0.0F) {
    u1_tmp = (real32_T)sqrt(D);
    D = (-B + u1_tmp) / (2.0F * A);
    A = (-B - u1_tmp) / (2.0F * A);
    if ((D >= 0.0F) && (D <= 1.0F) && (A >= 0.0F) && (A <= 1.0F)) {
      if (D > A) {
        u = D;
      } else {
        u = A;
      }

      guard1 = true;
    } else if ((D >= 0.0F) && (D <= 1.0F)) {
      u = D;
      guard1 = true;
    } else if ((A >= 0.0F) && (A <= 1.0F)) {
      u = A;
      guard1 = true;
    }
  } else if (D == 0.0F) {
    D = -B / (2.0F * A);
    if ((D >= 0.0F) && (D <= 1.0F)) {
      u = D;
      guard1 = true;
    }
  }

  if (guard1) {
    rty_P[0] = a_tmp * u + rtu_P1[0];
    rty_P[1] = b_a_tmp * u + rtu_P1[1];
  }

  *rty_u = u;
}

/*
 * Output and update for atomic system:
 *    '<S66>/OutRegionRegWP'
 *    '<S48>/OutRegionRegWP'
 *    '<S181>/OutRegionRegWP'
 */
void Formation_FMS_OutRegionRegWP_o(const real32_T rtu_P0[2], const real32_T
  rtu_P_Vehicle[2], real32_T rtu_R, real32_T rtu_L1, const real32_T rtu_n[2],
  real32_T rty_P[2])
{
  real32_T x;
  real32_T y_idx_0;
  x = rtu_P_Vehicle[0] - rtu_P0[0];
  y_idx_0 = x * x;
  x = rtu_P_Vehicle[1] - rtu_P0[1];
  if (x * x + y_idx_0 > rtu_R * rtu_R) {
    rty_P[0] = rtu_P0[0];
    rty_P[1] = rtu_P0[1];
  } else {
    rty_P[0] = rtu_n[0] * rtu_L1 + rtu_P_Vehicle[0];
    rty_P[1] = rtu_n[1] * rtu_L1 + rtu_P_Vehicle[1];
  }
}

/*
 * Output and update for atomic system:
 *    '<S66>/SearchL1RefWP'
 *    '<S48>/SearchL1RefWP'
 *    '<S181>/SearchL1RefWP'
 */
void Formation_FMS_SearchL1RefWP_i(const real32_T rtu_P_0[2], const real32_T
  rtu_P_Vehicle[2], real32_T rtu_R, real32_T rtu_L1, real32_T rty_P[2], real_T
  *rty_n)
{
  real32_T a;
  real32_T a_tmp;
  real32_T d;
  real32_T n0;
  real32_T n0_idx_0;
  rty_P[0] = 0.0F;
  rty_P[1] = 0.0F;
  if ((rtu_P_Vehicle[0] == rtu_P_0[0]) && (rtu_P_Vehicle[1] == rtu_P_0[1]) &&
      (rtu_R == rtu_L1)) {
    *rty_n = 0.0;
  } else {
    n0 = rtu_P_0[0] - rtu_P_Vehicle[0];
    a = n0 * n0;
    n0_idx_0 = n0;
    n0 = rtu_P_0[1] - rtu_P_Vehicle[1];
    d = (real32_T)sqrt(n0 * n0 + a);
    a_tmp = rtu_L1 * rtu_L1;
    a = ((d * d + a_tmp) - rtu_R * rtu_R) / (2.0F * d);
    n0_idx_0 /= d;
    n0 /= d;
    d = a * a;
    if (d > a_tmp) {
      *rty_n = 0.0;
    } else if (d == a_tmp) {
      *rty_n = 1.0;
      rty_P[0] = a * n0_idx_0 + rtu_P_Vehicle[0];
      rty_P[1] = a * n0 + rtu_P_Vehicle[1];
    } else {
      *rty_n = 2.0;
      a_tmp = (real32_T)sqrt(a_tmp - d);
      rty_P[0] = (0.0F * n0_idx_0 - n0) * a_tmp + (a * n0_idx_0 + rtu_P_Vehicle
        [0]);
      rty_P[1] = (0.0F * n0 + n0_idx_0) * a_tmp + (a * n0 + rtu_P_Vehicle[1]);
    }
  }
}

/*
 * System initialize for action system:
 *    '<S12>/Default'
 *    '<S10>/Default'
 */
void Formation_FMS_Default_Init(DW_Default_Formation_FMS_T *localDW)
{
  /* InitializeConditions for Delay: '<S22>/Delay' */
  localDW->icLoad = true;
}

/*
 * System reset for action system:
 *    '<S12>/Default'
 *    '<S10>/Default'
 */
void Formation_FMS_Default_Reset(DW_Default_Formation_FMS_T *localDW)
{
  /* InitializeConditions for Delay: '<S22>/Delay' */
  localDW->icLoad = true;
}

/*
 * Output and update for action system:
 *    '<S12>/Default'
 *    '<S10>/Default'
 */
void Formation_FMS_Default(const real32_T *rtu_FMS_In, const real32_T
  *rtu_FMS_In_h, FMS_Out_Bus *rty_FMS_Out, DW_Default_Formation_FMS_T *localDW)
{
  real32_T rtb_h_err_R_m;
  real32_T rtb_v_error_d;

  /* Sum: '<S21>/Sum' incorporates:
   *  Constant: '<S21>/Constant'
   */
  rtb_v_error_d = FMS_PARAM.FW_AIRSPD_TRIM - *rtu_FMS_In;

  /* Delay: '<S22>/Delay' */
  if (localDW->icLoad) {
    localDW->Delay_DSTATE = *rtu_FMS_In_h;
  }

  /* Sum: '<S22>/Sum' incorporates:
   *  Delay: '<S22>/Delay'
   */
  rtb_h_err_R_m = localDW->Delay_DSTATE - *rtu_FMS_In_h;
  memset(rty_FMS_Out, 0, sizeof(FMS_Out_Bus));

  /* BusAssignment: '<S16>/Bus Assignment' incorporates:
   *  Constant: '<S16>/Constant1'
   */
  rty_FMS_Out->state = VehicleState_None;
  rty_FMS_Out->ax_cmd = rtb_v_error_d;
  rty_FMS_Out->ay_cmd = 0.0F;

  /* Gain: '<S22>/Gain2' */
  rtb_v_error_d = FMS_PARAM.Z_P * rtb_h_err_R_m;

  /* Saturate: '<S22>/Saturation' */
  if (rtb_v_error_d > CONTROL_PARAM.FW_T_CLMB_MAX) {
    /* BusAssignment: '<S16>/Bus Assignment' */
    rty_FMS_Out->vh_cmd = CONTROL_PARAM.FW_T_CLMB_MAX;
  } else if (rtb_v_error_d < -CONTROL_PARAM.FW_T_SINK_MAX) {
    /* BusAssignment: '<S16>/Bus Assignment' */
    rty_FMS_Out->vh_cmd = -CONTROL_PARAM.FW_T_SINK_MAX;
  } else {
    /* BusAssignment: '<S16>/Bus Assignment' */
    rty_FMS_Out->vh_cmd = rtb_v_error_d;
  }

  /* End of Saturate: '<S22>/Saturation' */

  /* Update for Delay: '<S22>/Delay' */
  localDW->icLoad = false;
}

/* Function for Chart: '<Root>/FMS State Machine' */
static void Formation_exit_internal_Vehicle(DW_Formation_FMS_f_T *localDW)
{
  if (localDW->is_Vehicle == Formation_FMS_IN_Formation) {
    if (localDW->is_Formation == Formation_FMS_IN_FormAssemble) {
      localDW->is_FormAssemble = Formation_FM_IN_NO_ACTIVE_CHILD;
      localDW->is_Formation = Formation_FM_IN_NO_ACTIVE_CHILD;
    } else {
      localDW->is_FormMission = Formation_FM_IN_NO_ACTIVE_CHILD;
      localDW->is_Formation = Formation_FM_IN_NO_ACTIVE_CHILD;
    }

    localDW->is_Vehicle = Formation_FM_IN_NO_ACTIVE_CHILD;
  } else {
    localDW->is_Vehicle = Formation_FM_IN_NO_ACTIVE_CHILD;
  }
}

/* Function for MATLAB Function: '<S196>/MATLAB Function' */
static real_T Formation_FMS_atan3_9FIHpp9F(real_T x, real_T y, real_T x0, real_T
  b_y0, B_Formation_FMS_c_T *localB)
{
  real_T result;
  localB->deltax = x - x0;
  localB->deltay = y - b_y0;
  if ((localB->deltax > 0.0) && (localB->deltay > 0.0)) {
    result = atan(localB->deltay / localB->deltax);
  } else if ((localB->deltax > 0.0) && (localB->deltay < 0.0)) {
    result = atan(localB->deltay / localB->deltax) + 6.2831853071795862;
  } else if (localB->deltax < 0.0) {
    result = atan(localB->deltay / localB->deltax) + 3.1415926535897931;
  } else if (localB->deltax == 0.0) {
    if (localB->deltay > 0.0) {
      result = 1.5707963267948966;
    } else {
      result = 4.71238898038469;
    }
  } else if (localB->deltax > 0.0) {
    result = 0.0;
  } else {
    result = 3.1415926535897931;
  }

  return result;
}

/* Function for MATLAB Function: '<S196>/MATLAB Function' */
static void Formation_FMS_Dubins(const captured_var_Formation_FMS_T *PhiMaximum,
  const captured_var_Formation_FMS_T *rad2deg, real_T xs, real_T ys, real_T
  psi_s, real_T xf, real_T yf, real_T psi_f, real_T v, real_T xts[4], real_T
  yts[4], real_T xtf[4], real_T ytf[4], real_T cs[4], real_T cf[4], real_T lt[4],
  real_T l[4], real_T pos[4], real_T pof[4], real_T *r, real_T m_data[], int32_T
  *m_size, real_T n_data[], int32_T *n_size, B_Formation_FMS_c_T *localB)
{
  int32_T idx;
  int32_T ii;
  int32_T jj;
  int8_T j_data[4];
  boolean_T exitg1;
  boolean_T guard1;
  *r = v * v / 9.81 / tan(PhiMaximum->contents / rad2deg->contents);
  localB->l_best = psi_s / rad2deg->contents;
  localB->pos_tmp = sin(localB->l_best) * *r;
  pos[0] = localB->pos_tmp + xs;
  localB->l_best = cos(localB->l_best) * *r;
  pos[2] = ys - localB->l_best;
  pos[1] = xs - localB->pos_tmp;
  pos[3] = localB->l_best + ys;
  localB->l_best = psi_f / rad2deg->contents;
  localB->pos_tmp = sin(localB->l_best) * *r;
  pof[0] = localB->pos_tmp + xf;
  localB->l_best = cos(localB->l_best) * *r;
  pof[2] = yf - localB->l_best;
  pof[1] = xf - localB->pos_tmp;
  pof[3] = localB->l_best + yf;
  localB->pos_tmp = Formation_FMS_atan3_9FIHpp9F(pof[0], pof[2], pos[0], pos[2],
    localB);
  localB->xts_tmp = *r * sin(localB->pos_tmp);
  localB->l_best = pos[0] - localB->xts_tmp;
  localB->yts_tmp = *r * cos(localB->pos_tmp);
  localB->pos_tmp = localB->yts_tmp + pos[2];
  localB->xts_tmp = pof[0] - localB->xts_tmp;
  localB->yts_tmp += pof[2];
  localB->cs_tmp = xs - localB->l_best;
  localB->cs_tmp_c = ys - localB->pos_tmp;
  localB->cs_tmp = asin(sqrt(localB->cs_tmp * localB->cs_tmp + localB->cs_tmp_c *
    localB->cs_tmp_c) / 2.0 / *r) * (*r * 2.0);
  localB->cs_tmp_c = xf - localB->xts_tmp;
  localB->cf_tmp = yf - localB->yts_tmp;
  localB->cs_tmp_c = asin(sqrt(localB->cs_tmp_c * localB->cs_tmp_c +
    localB->cf_tmp * localB->cf_tmp) / 2.0 / *r) * (*r * 2.0);
  localB->cf_tmp = localB->xts_tmp - localB->l_best;
  localB->lt_tmp = localB->yts_tmp - localB->pos_tmp;
  localB->cf_tmp = sqrt(localB->cf_tmp * localB->cf_tmp + localB->lt_tmp *
                        localB->lt_tmp);
  l[0] = (localB->cs_tmp + localB->cs_tmp_c) + localB->cf_tmp;
  xts[0] = localB->l_best;
  yts[0] = localB->pos_tmp;
  xtf[0] = localB->xts_tmp;
  ytf[0] = localB->yts_tmp;
  cs[0] = localB->cs_tmp;
  cf[0] = localB->cs_tmp_c;
  lt[0] = localB->cf_tmp;
  localB->l_best = pof[3] - pos[2];
  localB->pos_tmp = pof[1] - pos[0];
  localB->pos_tmp = acos(*r * 2.0 / sqrt(localB->l_best * localB->l_best +
    localB->pos_tmp * localB->pos_tmp)) + Formation_FMS_atan3_9FIHpp9F(pof[1],
    pof[3], pos[0], pos[2], localB);
  localB->xts_tmp = *r * cos(localB->pos_tmp);
  localB->l_best = localB->xts_tmp + pos[0];
  localB->yts_tmp = *r * sin(localB->pos_tmp);
  localB->pos_tmp = localB->yts_tmp + pos[2];
  localB->xts_tmp = pof[1] - localB->xts_tmp;
  localB->yts_tmp = pof[3] - localB->yts_tmp;
  localB->cs_tmp = xs - localB->l_best;
  localB->cs_tmp_c = ys - localB->pos_tmp;
  localB->cs_tmp = asin(sqrt(localB->cs_tmp * localB->cs_tmp + localB->cs_tmp_c *
    localB->cs_tmp_c) / 2.0 / *r) * (*r * 2.0);
  localB->cs_tmp_c = xf - localB->xts_tmp;
  localB->cf_tmp = yf - localB->yts_tmp;
  localB->cs_tmp_c = asin(sqrt(localB->cs_tmp_c * localB->cs_tmp_c +
    localB->cf_tmp * localB->cf_tmp) / 2.0 / *r) * (*r * 2.0);
  localB->cf_tmp = localB->xts_tmp - localB->l_best;
  localB->lt_tmp = localB->yts_tmp - localB->pos_tmp;
  localB->cf_tmp = sqrt(localB->cf_tmp * localB->cf_tmp + localB->lt_tmp *
                        localB->lt_tmp);
  l[2] = (localB->cs_tmp + localB->cs_tmp_c) + localB->cf_tmp;
  xts[2] = localB->l_best;
  yts[2] = localB->pos_tmp;
  xtf[2] = localB->xts_tmp;
  ytf[2] = localB->yts_tmp;
  cs[2] = localB->cs_tmp;
  cf[2] = localB->cs_tmp_c;
  lt[2] = localB->cf_tmp;
  localB->l_best = pof[2] - pos[3];
  localB->pos_tmp = pof[0] - pos[1];
  localB->pos_tmp = Formation_FMS_atan3_9FIHpp9F(pof[0], pof[2], pos[1], pos[3],
    localB) - acos(*r * 2.0 / sqrt(localB->l_best * localB->l_best +
    localB->pos_tmp * localB->pos_tmp));
  localB->xts_tmp = *r * cos(localB->pos_tmp);
  localB->l_best = localB->xts_tmp + pos[1];
  localB->yts_tmp = *r * sin(localB->pos_tmp);
  localB->pos_tmp = localB->yts_tmp + pos[3];
  localB->xts_tmp = pof[0] - localB->xts_tmp;
  localB->yts_tmp = pof[2] - localB->yts_tmp;
  localB->cs_tmp = xs - localB->l_best;
  localB->cs_tmp_c = ys - localB->pos_tmp;
  localB->cs_tmp = asin(sqrt(localB->cs_tmp * localB->cs_tmp + localB->cs_tmp_c *
    localB->cs_tmp_c) / 2.0 / *r) * (*r * 2.0);
  localB->cs_tmp_c = xf - localB->xts_tmp;
  localB->cf_tmp = yf - localB->yts_tmp;
  localB->cs_tmp_c = asin(sqrt(localB->cs_tmp_c * localB->cs_tmp_c +
    localB->cf_tmp * localB->cf_tmp) / 2.0 / *r) * (*r * 2.0);
  localB->cf_tmp = localB->xts_tmp - localB->l_best;
  localB->lt_tmp = localB->yts_tmp - localB->pos_tmp;
  localB->cf_tmp = sqrt(localB->cf_tmp * localB->cf_tmp + localB->lt_tmp *
                        localB->lt_tmp);
  l[1] = (localB->cs_tmp + localB->cs_tmp_c) + localB->cf_tmp;
  xts[1] = localB->l_best;
  yts[1] = localB->pos_tmp;
  xtf[1] = localB->xts_tmp;
  ytf[1] = localB->yts_tmp;
  cs[1] = localB->cs_tmp;
  cf[1] = localB->cs_tmp_c;
  lt[1] = localB->cf_tmp;
  localB->pos_tmp = Formation_FMS_atan3_9FIHpp9F(pof[1], pof[3], pos[1], pos[3],
    localB);
  localB->xts_tmp = *r * sin(localB->pos_tmp);
  localB->l_best = localB->xts_tmp + pos[1];
  localB->yts_tmp = *r * cos(localB->pos_tmp);
  localB->pos_tmp = pos[3] - localB->yts_tmp;
  localB->xts_tmp += pof[1];
  localB->yts_tmp = pof[3] - localB->yts_tmp;
  localB->cs_tmp = xs - localB->l_best;
  localB->cs_tmp_c = ys - localB->pos_tmp;
  localB->cs_tmp = asin(sqrt(localB->cs_tmp * localB->cs_tmp + localB->cs_tmp_c *
    localB->cs_tmp_c) / 2.0 / *r) * (*r * 2.0);
  localB->cs_tmp_c = xf - localB->xts_tmp;
  localB->cf_tmp = yf - localB->yts_tmp;
  localB->cs_tmp_c = asin(sqrt(localB->cs_tmp_c * localB->cs_tmp_c +
    localB->cf_tmp * localB->cf_tmp) / 2.0 / *r) * (*r * 2.0);
  localB->cf_tmp = localB->xts_tmp - localB->l_best;
  localB->lt_tmp = localB->yts_tmp - localB->pos_tmp;
  localB->cf_tmp = sqrt(localB->cf_tmp * localB->cf_tmp + localB->lt_tmp *
                        localB->lt_tmp);
  l[3] = (localB->cs_tmp + localB->cs_tmp_c) + localB->cf_tmp;
  xts[3] = localB->l_best;
  yts[3] = localB->pos_tmp;
  xtf[3] = localB->xts_tmp;
  ytf[3] = localB->yts_tmp;
  cs[3] = localB->cs_tmp;
  cf[3] = localB->cs_tmp_c;
  lt[3] = localB->cf_tmp;
  localB->l_best = (rtInf);
  if ((l[0] >= 0.0) && (l[0] < (rtInf))) {
    localB->l_best = l[0];
  }

  if ((l[1] >= 0.0) && (l[1] < localB->l_best)) {
    localB->l_best = l[1];
  }

  if ((l[2] >= 0.0) && (l[2] < localB->l_best)) {
    localB->l_best = l[2];
  }

  if ((l[3] >= 0.0) && (l[3] < localB->l_best)) {
    localB->l_best = l[3];
  }

  idx = -1;
  ii = 1;
  jj = 1;
  exitg1 = false;
  while ((!exitg1) && (jj <= 2)) {
    guard1 = false;
    if (l[(((jj - 1) << 1) + ii) - 1] == localB->l_best) {
      idx++;
      localB->i_data[idx] = ii;
      j_data[idx] = (int8_T)jj;
      if (idx + 1 >= 4) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      ii++;
      if (ii > 2) {
        ii = 1;
        jj++;
      }
    }
  }

  if (idx + 1 < 1) {
    jj = -1;
  } else {
    jj = idx;
  }

  *m_size = jj + 1;
  for (ii = 0; ii <= jj; ii++) {
    m_data[ii] = localB->i_data[ii];
  }

  if (idx + 1 < 1) {
    idx = -1;
  }

  *n_size = idx + 1;
  for (ii = 0; ii <= idx; ii++) {
    n_data[ii] = j_data[ii];
  }

  if ((jj + 1 != 1) || (idx + 1 != 1)) {
    idx = (int32_T)m_data[0];
    *m_size = 1;
    m_data[0] = idx;
    idx = (int32_T)n_data[0];
    *n_size = 1;
    n_data[0] = idx;
  }
}

/* Function for Chart: '<Root>/FMS State Machine' */
static void Format_enter_internal_Formation(const INS_Out_Bus
  *BusConversion_InsertedFor_FMS_c, const Formation_Cross_Bus
  *BusConversion_InsertedFor_FMS_p, B_Formation_FMS_c_T *localB,
  DW_Formation_FMS_f_T *localDW)
{
  int32_T i;
  boolean_T out;
  static const real_T d[8] = { 1000.0, 100.0, 10.0, 1.0, 0.1, 0.01, 0.001,
    0.0001 };

  int32_T b_size;
  int32_T c_size;
  boolean_T exitg1;

  /* Delay: '<S5>/Delay' incorporates:
   *  MATLAB Function: '<S196>/MATLAB Function'
   */
  if (localDW->Delay_DSTATE_j == PilotMode_FormAssemble) {
    localDW->is_Formation = Formation_FMS_IN_FormAssemble;
    localB->Cmd_In.sp_waypoint[0] = BusConversion_InsertedFor_FMS_c->x_R;
    localB->Cmd_In.sp_waypoint[1] = BusConversion_InsertedFor_FMS_c->y_R;
    localB->Cmd_In.sp_waypoint[2] = BusConversion_InsertedFor_FMS_c->h_R;
    if (FORMATION_PARAM.UAV_ID == 1U) {
      /* Outputs for Function Call SubSystem: '<S3>/Vehicle.Formation.FormAssemble.dubinsPath' */
      /* MATLAB Function: '<S196>/MATLAB Function' */
      for (localB->k = 0; localB->k < 9; localB->k++) {
        localB->Goal[localB->k] = FORMATION_PARAM.FORM_POINT[localB->k];
      }

      localB->PhiMaximum.contents = 40.0;
      localB->rad2deg.contents = 57.295779513082323;
      localB->obj.xs = 0.0;
      localB->obj.ys = 0.0;
      localB->obj.psi_s = 0.0;
      localB->obj.xf = 0.0;
      localB->obj.yf = 0.0;
      localB->obj.psi_f = 0.0;
      localB->obj.v = 0.0;
      localB->obj.r = 0.0;
      localB->obj.pos[0] = 0.0;
      localB->obj.pof[0] = 0.0;
      localB->obj.xts[0] = 0.0;
      localB->obj.yts[0] = 0.0;
      localB->obj.xtf[0] = 0.0;
      localB->obj.ytf[0] = 0.0;
      localB->obj.cs[0] = 0.0;
      localB->obj.cf[0] = 0.0;
      localB->obj.lt[0] = 0.0;
      localB->obj.l[0] = 0.0;
      localB->obj.pos[1] = 0.0;
      localB->obj.pof[1] = 0.0;
      localB->obj.xts[1] = 0.0;
      localB->obj.yts[1] = 0.0;
      localB->obj.xtf[1] = 0.0;
      localB->obj.ytf[1] = 0.0;
      localB->obj.cs[1] = 0.0;
      localB->obj.cf[1] = 0.0;
      localB->obj.lt[1] = 0.0;
      localB->obj.l[1] = 0.0;
      localB->obj.pos[2] = 0.0;
      localB->obj.pof[2] = 0.0;
      localB->obj.xts[2] = 0.0;
      localB->obj.yts[2] = 0.0;
      localB->obj.xtf[2] = 0.0;
      localB->obj.ytf[2] = 0.0;
      localB->obj.cs[2] = 0.0;
      localB->obj.cf[2] = 0.0;
      localB->obj.lt[2] = 0.0;
      localB->obj.l[2] = 0.0;
      localB->obj.pos[3] = 0.0;
      localB->obj.pof[3] = 0.0;
      localB->obj.xts[3] = 0.0;
      localB->obj.yts[3] = 0.0;
      localB->obj.xtf[3] = 0.0;
      localB->obj.ytf[3] = 0.0;
      localB->obj.cs[3] = 0.0;
      localB->obj.cf[3] = 0.0;
      localB->obj.lt[3] = 0.0;
      localB->obj.l[3] = 0.0;
      localB->obj.index_dubins[0] = 0.0;
      localB->obj.index_dubins[1] = 0.0;
      localB->obj.l_ad = 0.0;
      localB->obj.precision_flag = 0.0;
      localB->obj.xm = 0.0;
      localB->obj.ym = 0.0;
      localB->object[0] = localB->obj;
      localB->object[1] = localB->obj;
      localB->object[2] = localB->obj;
      localB->object[0].xs = BusConversion_InsertedFor_FMS_p->x_R[0];
      localB->object[0].ys = BusConversion_InsertedFor_FMS_p->y_R[0];
      localB->object[0].psi_s = rt_atan2f_snf
        (BusConversion_InsertedFor_FMS_p->ve[0],
         BusConversion_InsertedFor_FMS_p->vn[0]);
      localB->object[0].xf = localB->Goal[0];
      localB->object[0].yf = localB->Goal[3];
      localB->object[0].psi_f = localB->Goal[6];
      localB->object[0].v = 25.0;
      localB->object[1].xs = BusConversion_InsertedFor_FMS_p->x_R[1];
      localB->object[1].ys = BusConversion_InsertedFor_FMS_p->y_R[1];
      localB->object[1].psi_s = rt_atan2f_snf
        (BusConversion_InsertedFor_FMS_p->ve[1],
         BusConversion_InsertedFor_FMS_p->vn[1]);
      localB->object[1].xf = localB->Goal[1];
      localB->object[1].yf = localB->Goal[4];
      localB->object[1].psi_f = localB->Goal[7];
      localB->object[1].v = 25.0;
      localB->object[2].xs = BusConversion_InsertedFor_FMS_p->x_R[2];
      localB->object[2].ys = BusConversion_InsertedFor_FMS_p->y_R[2];
      localB->object[2].psi_s = rt_atan2f_snf
        (BusConversion_InsertedFor_FMS_p->ve[2],
         BusConversion_InsertedFor_FMS_p->vn[2]);
      localB->object[2].xf = localB->Goal[2];
      localB->object[2].yf = localB->Goal[5];
      localB->object[2].psi_f = localB->Goal[8];
      localB->object[2].v = 25.0;
      localB->target = -1;
      localB->l_ref = -1.0;
      localB->obj = localB->object[0];

      /* MATLAB Function: '<S196>/MATLAB Function' */
      Formation_FMS_Dubins(&localB->PhiMaximum, &localB->rad2deg, (real_T)
                           BusConversion_InsertedFor_FMS_p->x_R[0], (real_T)
                           BusConversion_InsertedFor_FMS_p->y_R[0], (real_T)
                           rt_atan2f_snf(BusConversion_InsertedFor_FMS_p->ve[0],
        BusConversion_InsertedFor_FMS_p->vn[0]), localB->Goal[0], localB->Goal[3],
                           localB->Goal[6], 25.0, localB->obj.xts,
                           localB->obj.yts, localB->obj.xtf, localB->obj.ytf,
                           localB->obj.cs, localB->obj.cf, localB->obj.lt,
                           localB->obj.l, localB->obj.pos, localB->obj.pof,
                           &localB->obj.r, localB->b_data, &b_size,
                           localB->c_data, &c_size, localB);
      localB->obj.index_dubins[0] = localB->b_data[0];
      localB->obj.index_dubins[1] = localB->c_data[0];
      localB->goal = localB->obj.l[((((int32_T)localB->c_data[0] - 1) << 1) +
        (int32_T)localB->b_data[0]) - 1];
      if (localB->goal > -1.0) {
        localB->l_ref = localB->goal;
        localB->target = 1;
      }

      localB->object[0] = localB->obj;
      localB->obj = localB->object[1];

      /* MATLAB Function: '<S196>/MATLAB Function' */
      Formation_FMS_Dubins(&localB->PhiMaximum, &localB->rad2deg, (real_T)
                           BusConversion_InsertedFor_FMS_p->x_R[1], (real_T)
                           BusConversion_InsertedFor_FMS_p->y_R[1], (real_T)
                           rt_atan2f_snf(BusConversion_InsertedFor_FMS_p->ve[1],
        BusConversion_InsertedFor_FMS_p->vn[1]), localB->Goal[1], localB->Goal[4],
                           localB->Goal[7], 25.0, localB->obj.xts,
                           localB->obj.yts, localB->obj.xtf, localB->obj.ytf,
                           localB->obj.cs, localB->obj.cf, localB->obj.lt,
                           localB->obj.l, localB->obj.pos, localB->obj.pof,
                           &localB->obj.r, localB->b_data, &b_size,
                           localB->c_data, &c_size, localB);
      localB->obj.index_dubins[0] = localB->b_data[0];
      localB->obj.index_dubins[1] = localB->c_data[0];
      localB->goal = localB->obj.l[((((int32_T)localB->c_data[0] - 1) << 1) +
        (int32_T)localB->b_data[0]) - 1];
      if (localB->goal > localB->l_ref) {
        localB->l_ref = localB->goal;
        localB->target = 2;
      }

      localB->object[1] = localB->obj;
      localB->obj = localB->object[2];

      /* MATLAB Function: '<S196>/MATLAB Function' */
      Formation_FMS_Dubins(&localB->PhiMaximum, &localB->rad2deg, (real_T)
                           BusConversion_InsertedFor_FMS_p->x_R[2], (real_T)
                           BusConversion_InsertedFor_FMS_p->y_R[2], (real_T)
                           rt_atan2f_snf(BusConversion_InsertedFor_FMS_p->ve[2],
        BusConversion_InsertedFor_FMS_p->vn[2]), localB->Goal[2], localB->Goal[5],
                           localB->Goal[8], 25.0, localB->obj.xts,
                           localB->obj.yts, localB->obj.xtf, localB->obj.ytf,
                           localB->obj.cs, localB->obj.cf, localB->obj.lt,
                           localB->obj.l, localB->obj.pos, localB->obj.pof,
                           &localB->obj.r, localB->b_data, &b_size,
                           localB->c_data, &c_size, localB);
      localB->obj.index_dubins[0] = localB->b_data[0];
      localB->obj.index_dubins[1] = localB->c_data[0];
      localB->goal = localB->obj.l[((((int32_T)localB->c_data[0] - 1) << 1) +
        (int32_T)localB->b_data[0]) - 1];
      if (localB->goal > localB->l_ref) {
        localB->l_ref = localB->goal;
        localB->target = 3;
      }

      localB->object[2] = localB->obj;

      /* MATLAB Function: '<S196>/MATLAB Function' */
      for (localB->k = 0; localB->k < 3; localB->k++) {
        if (localB->k + 1 != localB->target) {
          localB->search_floor = 0.0;
          localB->search = 0.0;
          localB->stop_flag = 1;
          for (i = 0; i < 7; i++) {
            exitg1 = false;
            while ((!exitg1) && (localB->stop_flag != 0)) {
              localB->search += d[i];
              localB->goal = localB->object[localB->k].psi_f /
                57.295779513082323;
              localB->object[localB->k].xm = localB->object[localB->k].xf - cos
                (localB->goal) * localB->search;
              localB->object[localB->k].ym = localB->object[localB->k].yf - sin
                (localB->goal) * localB->search;
              Formation_FMS_Dubins(&localB->PhiMaximum, &localB->rad2deg,
                                   localB->object[localB->k].xs, localB->
                                   object[localB->k].ys, localB->object
                                   [localB->k].psi_s, localB->object[localB->k].
                                   xm, localB->object[localB->k].ym,
                                   localB->object[localB->k].psi_f,
                                   localB->object[localB->k].v, localB->
                                   object[localB->k].xts, localB->object
                                   [localB->k].yts, localB->object[localB->k].
                                   xtf, localB->object[localB->k].ytf,
                                   localB->object[localB->k].cs, localB->
                                   object[localB->k].cf, localB->object
                                   [localB->k].lt, localB->object[localB->k].l,
                                   localB->object[localB->k].pos, localB->
                                   object[localB->k].pof, &localB->goal,
                                   localB->b_data, &b_size, localB->c_data,
                                   &c_size, localB);
              localB->object[localB->k].index_dubins[0] = localB->b_data[0];
              localB->object[localB->k].index_dubins[1] = localB->c_data[0];
              localB->goal = (localB->object[localB->k].l[((((int32_T)
                localB->object[localB->k].index_dubins[1] - 1) << 1) + (int32_T)
                localB->object[localB->k].index_dubins[0]) - 1] + localB->search)
                - localB->l_ref;
              if (localB->goal > 0.0) {
                localB->object[localB->k].l_ad = localB->search_floor;
                localB->object[localB->k].precision_flag = d[i];
                localB->search = localB->search_floor;
                exitg1 = true;
              } else if (localB->goal < 0.0) {
                localB->search_floor = localB->search;
              } else {
                localB->object[localB->k].l_ad = localB->search;
                localB->object[localB->k].precision_flag = 0.0;
                localB->stop_flag = 0;
              }
            }
          }
        } else {
          localB->object[localB->k].xm = localB->object[localB->k].xf;
          localB->object[localB->k].ym = localB->object[localB->k].yf;
        }
      }

      if (localB->target == 1) {
        localB->result[0] = localB->object[0].xs;
        localB->result[15] = localB->object[0].ys;
        localB->result[30] = localB->object[0].psi_s;
        localB->result[45] = localB->object[0].r;
        localB->k = ((((int32_T)localB->object[0].index_dubins[1] - 1) << 1) +
                     (int32_T)localB->object[0].index_dubins[0]) - 1;
        localB->result[60] = localB->object[0].cs[localB->k];
        localB->result[3] = localB->object[0].xts[localB->k];
        localB->result[18] = localB->object[0].yts[localB->k];
        localB->l_ref = Formation_FMS_atan3_9FIHpp9F(localB->object[0]
          .xtf[localB->k], localB->object[0].ytf[localB->k], localB->object[0].
          xts[localB->k], localB->object[0].yts[localB->k], localB) *
          57.295779513082323;
        localB->result[33] = localB->l_ref;
        localB->result[48] = 0.0;
        localB->result[63] = localB->object[0].lt[localB->k];
        localB->result[6] = localB->object[0].xtf[localB->k];
        localB->result[21] = localB->object[0].ytf[localB->k];
        localB->result[36] = localB->l_ref;
        localB->result[51] = localB->object[0].r;
        localB->result[66] = localB->object[0].cf[localB->k];
        localB->result[9] = localB->object[0].xm;
        localB->result[24] = localB->object[0].ym;
        localB->result[39] = localB->object[0].psi_f;
        localB->result[54] = 0.0;
        localB->result[69] = 0.0;
        localB->result[12] = localB->object[0].xf;
        localB->result[27] = localB->object[0].yf;
        localB->result[42] = localB->object[0].psi_f;
        localB->result[57] = 0.0;
        localB->result[72] = localB->object[0].l[localB->k];
      } else {
        localB->result[0] = localB->object[0].xs;
        localB->result[15] = localB->object[0].ys;
        localB->result[30] = localB->object[0].psi_s;
        localB->result[45] = localB->object[0].r;
        localB->k = ((((int32_T)localB->object[0].index_dubins[1] - 1) << 1) +
                     (int32_T)localB->object[0].index_dubins[0]) - 1;
        localB->result[60] = localB->object[0].cs[localB->k];
        localB->l_ref = localB->object[0].xts[localB->k];
        localB->result[3] = localB->l_ref;
        localB->search_floor = localB->object[0].yts[localB->k];
        localB->result[18] = localB->search_floor;
        localB->search = localB->object[0].xtf[localB->k];
        localB->goal = localB->object[0].ytf[localB->k];
        localB->l_ref = Formation_FMS_atan3_9FIHpp9F(localB->search,
          localB->goal, localB->l_ref, localB->search_floor, localB) *
          57.295779513082323;
        localB->result[33] = localB->l_ref;
        localB->result[48] = 0.0;
        localB->result[63] = localB->object[0].lt[localB->k];
        localB->result[6] = localB->search;
        localB->result[21] = localB->goal;
        localB->result[36] = localB->l_ref;
        localB->result[51] = localB->object[0].r;
        localB->result[66] = localB->object[0].cf[localB->k];
        localB->result[9] = localB->object[0].xm;
        localB->result[24] = localB->object[0].ym;
        localB->result[39] = localB->object[0].psi_f;
        localB->result[54] = 0.0;
        localB->result[69] = localB->object[0].l_ad;
        localB->result[12] = localB->object[0].xf;
        localB->result[27] = localB->object[0].yf;
        localB->result[42] = localB->object[0].psi_f;
        localB->result[57] = 0.0;
        localB->result[72] = localB->object[0].l[localB->k] + localB->object[0].
          l_ad;
      }

      if (localB->target == 2) {
        localB->result[1] = localB->object[1].xs;
        localB->result[16] = localB->object[1].ys;
        localB->result[31] = localB->object[1].psi_s;
        localB->result[46] = localB->object[1].r;
        localB->k = ((((int32_T)localB->object[1].index_dubins[1] - 1) << 1) +
                     (int32_T)localB->object[1].index_dubins[0]) - 1;
        localB->result[61] = localB->object[1].cs[localB->k];
        localB->result[4] = localB->object[1].xts[localB->k];
        localB->result[19] = localB->object[1].yts[localB->k];
        localB->l_ref = Formation_FMS_atan3_9FIHpp9F(localB->object[1]
          .xtf[localB->k], localB->object[1].ytf[localB->k], localB->object[1].
          xts[localB->k], localB->object[1].yts[localB->k], localB) *
          57.295779513082323;
        localB->result[34] = localB->l_ref;
        localB->result[49] = 0.0;
        localB->result[64] = localB->object[1].lt[localB->k];
        localB->result[7] = localB->object[1].xtf[localB->k];
        localB->result[22] = localB->object[1].ytf[localB->k];
        localB->result[37] = localB->l_ref;
        localB->result[52] = localB->object[1].r;
        localB->result[67] = localB->object[1].cf[localB->k];
        localB->result[10] = localB->object[1].xm;
        localB->result[25] = localB->object[1].ym;
        localB->result[40] = localB->object[1].psi_f;
        localB->result[55] = 0.0;
        localB->result[70] = 0.0;
        localB->result[13] = localB->object[1].xf;
        localB->result[28] = localB->object[1].yf;
        localB->result[43] = localB->object[1].psi_f;
        localB->result[58] = 0.0;
        localB->result[73] = localB->object[1].l[localB->k];
      } else {
        localB->result[1] = localB->object[1].xs;
        localB->result[16] = localB->object[1].ys;
        localB->result[31] = localB->object[1].psi_s;
        localB->result[46] = localB->object[1].r;
        localB->k = ((((int32_T)localB->object[1].index_dubins[1] - 1) << 1) +
                     (int32_T)localB->object[1].index_dubins[0]) - 1;
        localB->result[61] = localB->object[1].cs[localB->k];
        localB->l_ref = localB->object[1].xts[localB->k];
        localB->result[4] = localB->l_ref;
        localB->search_floor = localB->object[1].yts[localB->k];
        localB->result[19] = localB->search_floor;
        localB->search = localB->object[1].xtf[localB->k];
        localB->goal = localB->object[1].ytf[localB->k];
        localB->l_ref = Formation_FMS_atan3_9FIHpp9F(localB->search,
          localB->goal, localB->l_ref, localB->search_floor, localB) *
          57.295779513082323;
        localB->result[34] = localB->l_ref;
        localB->result[49] = 0.0;
        localB->result[64] = localB->object[1].lt[localB->k];
        localB->result[7] = localB->search;
        localB->result[22] = localB->goal;
        localB->result[37] = localB->l_ref;
        localB->result[52] = localB->object[1].r;
        localB->result[67] = localB->object[1].cf[localB->k];
        localB->result[10] = localB->object[1].xm;
        localB->result[25] = localB->object[1].ym;
        localB->result[40] = localB->object[1].psi_f;
        localB->result[55] = 0.0;
        localB->result[70] = localB->object[1].l_ad;
        localB->result[13] = localB->object[1].xf;
        localB->result[28] = localB->object[1].yf;
        localB->result[43] = localB->object[1].psi_f;
        localB->result[58] = 0.0;
        localB->result[73] = localB->object[1].l[localB->k] + localB->object[1].
          l_ad;
      }

      if (localB->target == 3) {
        localB->result[2] = localB->object[2].xs;
        localB->result[17] = localB->object[2].ys;
        localB->result[32] = localB->object[2].psi_s;
        localB->result[47] = localB->object[2].r;
        localB->k = ((((int32_T)localB->object[2].index_dubins[1] - 1) << 1) +
                     (int32_T)localB->object[2].index_dubins[0]) - 1;
        localB->result[62] = localB->object[2].cs[localB->k];
        localB->result[5] = localB->object[2].xts[localB->k];
        localB->result[20] = localB->object[2].yts[localB->k];
        localB->l_ref = Formation_FMS_atan3_9FIHpp9F(localB->object[2]
          .xtf[localB->k], localB->object[2].ytf[localB->k], localB->object[2].
          xts[localB->k], localB->object[2].yts[localB->k], localB) *
          57.295779513082323;
        localB->result[35] = localB->l_ref;
        localB->result[50] = 0.0;
        localB->result[65] = localB->object[2].lt[localB->k];
        localB->result[8] = localB->object[2].xtf[localB->k];
        localB->result[23] = localB->object[2].ytf[localB->k];
        localB->result[38] = localB->l_ref;
        localB->result[53] = localB->object[2].r;
        localB->result[68] = localB->object[2].cf[localB->k];
        localB->result[11] = localB->object[2].xm;
        localB->result[26] = localB->object[2].ym;
        localB->result[41] = localB->object[2].psi_f;
        localB->result[56] = 0.0;
        localB->result[71] = 0.0;
        localB->result[14] = localB->object[2].xf;
        localB->result[29] = localB->object[2].yf;
        localB->result[44] = localB->object[2].psi_f;
        localB->result[59] = 0.0;
        localB->result[74] = localB->object[2].l[localB->k];
      } else {
        localB->result[2] = localB->object[2].xs;
        localB->result[17] = localB->object[2].ys;
        localB->result[32] = localB->object[2].psi_s;
        localB->result[47] = localB->object[2].r;
        localB->k = ((((int32_T)localB->object[2].index_dubins[1] - 1) << 1) +
                     (int32_T)localB->object[2].index_dubins[0]) - 1;
        localB->result[62] = localB->object[2].cs[localB->k];
        localB->l_ref = localB->object[2].xts[localB->k];
        localB->result[5] = localB->l_ref;
        localB->search_floor = localB->object[2].yts[localB->k];
        localB->result[20] = localB->search_floor;
        localB->search = localB->object[2].xtf[localB->k];
        localB->goal = localB->object[2].ytf[localB->k];
        localB->l_ref = Formation_FMS_atan3_9FIHpp9F(localB->search,
          localB->goal, localB->l_ref, localB->search_floor, localB) *
          57.295779513082323;
        localB->result[35] = localB->l_ref;
        localB->result[50] = 0.0;
        localB->result[65] = localB->object[2].lt[localB->k];
        localB->result[8] = localB->search;
        localB->result[23] = localB->goal;
        localB->result[38] = localB->l_ref;
        localB->result[53] = localB->object[2].r;
        localB->result[68] = localB->object[2].cf[localB->k];
        localB->result[11] = localB->object[2].xm;
        localB->result[26] = localB->object[2].ym;
        localB->result[41] = localB->object[2].psi_f;
        localB->result[56] = 0.0;
        localB->result[71] = localB->object[2].l_ad;
        localB->result[14] = localB->object[2].xf;
        localB->result[29] = localB->object[2].yf;
        localB->result[44] = localB->object[2].psi_f;
        localB->result[59] = 0.0;
        localB->result[74] = localB->object[2].l[localB->k] + localB->object[2].
          l_ad;
      }

      memset(&localB->Other_Mission_Data.x[0], 0, 24U * sizeof(real32_T));
      memset(&localB->Other_Mission_Data.y[0], 0, 24U * sizeof(real32_T));
      memset(&localB->Other_Mission_Data.z[0], 0, 24U * sizeof(real32_T));
      memset(&localB->Other_Mission_Data.heading[0], 0, 24U * sizeof(real32_T));
      memset(&localB->Other_Mission_Data.ext1[0], 0, 24U * sizeof(real32_T));
      memset(&localB->Other_Mission_Data.ext2[0], 0, 24U * sizeof(real32_T));
      for (localB->k = 0; localB->k < 5; localB->k++) {
        localB->Other_Mission_Data.x[3 * localB->k] = (real32_T)localB->result[3
          * localB->k];
        localB->Other_Mission_Data.y[3 * localB->k] = (real32_T)localB->result[3
          * localB->k + 15];
        i = 3 * localB->k + 1;
        localB->Other_Mission_Data.x[i] = (real32_T)localB->result[i];
        localB->Other_Mission_Data.y[i] = (real32_T)localB->result[3 * localB->k
          + 16];
        i = 3 * localB->k + 2;
        localB->Other_Mission_Data.x[i] = (real32_T)localB->result[i];
        localB->Other_Mission_Data.y[i] = (real32_T)localB->result[3 * localB->k
          + 17];
      }

      for (i = 0; i < 3; i++) {
        for (localB->k = 0; localB->k < 5; localB->k++) {
          localB->Other_Mission_Data.z[i + 3 * localB->k] =
            BusConversion_InsertedFor_FMS_p->h_R[i];
        }
      }

      for (localB->k = 0; localB->k < 5; localB->k++) {
        localB->Other_Mission_Data.heading[3 * localB->k] = (real32_T)
          localB->result[3 * localB->k + 30];
        localB->Other_Mission_Data.ext1[3 * localB->k] = (real32_T)
          localB->result[3 * localB->k + 45];
        localB->Other_Mission_Data.ext2[3 * localB->k] = (real32_T)
          localB->result[3 * localB->k + 60];
        i = 3 * localB->k + 1;
        localB->Other_Mission_Data.heading[i] = (real32_T)localB->result[3 *
          localB->k + 31];
        localB->Other_Mission_Data.ext1[i] = (real32_T)localB->result[3 *
          localB->k + 46];
        localB->Other_Mission_Data.ext2[i] = (real32_T)localB->result[3 *
          localB->k + 61];
        i = 3 * localB->k + 2;
        localB->Other_Mission_Data.heading[i] = (real32_T)localB->result[3 *
          localB->k + 32];
        localB->Other_Mission_Data.ext1[i] = (real32_T)localB->result[3 *
          localB->k + 47];
        localB->Other_Mission_Data.ext2[i] = (real32_T)localB->result[3 *
          localB->k + 62];
      }

      localB->Other_Mission_Data.timestamp = 9999U;
      localB->Other_Mission_Data.type[0] = 1U;
      localB->Other_Mission_Data.valid_items[0] = 5U;
      localB->Other_Mission_Data.type[1] = 1U;
      localB->Other_Mission_Data.valid_items[1] = 5U;
      localB->Other_Mission_Data.type[2] = 1U;
      localB->Other_Mission_Data.valid_items[2] = 5U;

      /* End of Outputs for SubSystem: '<S3>/Vehicle.Formation.FormAssemble.dubinsPath' */
    }

    localDW->is_FormAssemble = Formation_FMS_IN_WaitForUpdate;
    localDW->temporalCounter_i1 = 0U;
    localB->Cmd_In.cur_waypoint[0] = BusConversion_InsertedFor_FMS_c->x_R;
    localB->Cmd_In.cur_waypoint[1] = BusConversion_InsertedFor_FMS_c->y_R;
    localB->Cmd_In.cur_waypoint[2] = BusConversion_InsertedFor_FMS_c->h_R;
    localB->state = VehicleState_FormHold;
  } else {
    if (localDW->Delay_DSTATE_j == PilotMode_FormMission) {
      localB->b[0] = (uint32_T)((localB->Cmd_In.form_valid & 1U) != 0U);
      localB->b[1] = (uint32_T)((localB->Cmd_In.form_valid & 2U) != 0U);
      localB->b[2] = (uint32_T)((localB->Cmd_In.form_valid & 4U) != 0U);
      out = true;
      localB->k = 0;
      exitg1 = false;
      while ((!exitg1) && (localB->k < 3)) {
        if (localB->b[localB->k] <= 0U) {
          out = false;
          exitg1 = true;
        } else {
          localB->k++;
        }
      }
    } else {
      out = false;
    }

    if (out) {
      localDW->is_Formation = Formation_FMS_IN_FormMission;
      localB->Cmd_In.sp_waypoint[0] = BusConversion_InsertedFor_FMS_c->x_R;
      localB->Cmd_In.sp_waypoint[1] = BusConversion_InsertedFor_FMS_c->y_R;
      localB->Cmd_In.sp_waypoint[2] = BusConversion_InsertedFor_FMS_c->h_R;
      localB->state = VehicleState_FormMission;
      localDW->is_FormMission = Formation_FMS_IN_NextWP;
    } else if (localDW->Delay_DSTATE_j == PilotMode_FormDisband) {
      localDW->is_Formation = Formation_FMS_IN_FormDisband;
      localB->Cmd_In.cur_waypoint[0] = FORMATION_PARAM.DISBAND_POINT[(int32_T)
        FORMATION_PARAM.UAV_ID - 1];
      localB->Cmd_In.cur_waypoint[1] = FORMATION_PARAM.DISBAND_POINT[(int32_T)
        FORMATION_PARAM.UAV_ID + 2];
      localB->Cmd_In.cur_waypoint[2] = FORMATION_PARAM.DISBAND_POINT[(int32_T)
        FORMATION_PARAM.UAV_ID + 5];
      localB->state = VehicleState_FormDisband;
    } else {
      localDW->is_Formation = Formation_FMS_IN_InvalidMode;
    }
  }

  /* End of Delay: '<S5>/Delay' */
}

/* System initialize for referenced model: 'Formation_FMS' */
void Formation_FMS_Init(uint32_T *rty_FMS_Out_timestamp, VehicleState
  *rty_FMS_Out_state, real32_T *rty_FMS_Out_ax_cmd, real32_T *rty_FMS_Out_ay_cmd,
  real32_T *rty_FMS_Out_vh_cmd, uint32_T *rty_Other_Mission_Data_timestam,
  uint32_T rty_Other_Mission_Data_type[3], uint8_T
  rty_Other_Mission_Data_valid_it[3], real32_T rty_Other_Mission_Data_x[24],
  real32_T rty_Other_Mission_Data_y[24], real32_T rty_Other_Mission_Data_z[24],
  real32_T rty_Other_Mission_Data_heading[24], real32_T
  rty_Other_Mission_Data_ext1[24], real32_T rty_Other_Mission_Data_ext2[24],
  B_Formation_FMS_c_T *localB, DW_Formation_FMS_f_T *localDW)
{
  int32_T i;

  /* Start for SwitchCase: '<S10>/Switch Case' */
  localDW->SwitchCase_ActiveSubsystem = -1;

  /* SystemInitialize for Chart: '<Root>/FMS State Machine' */
  localB->Cmd_In.form_valid = 0U;
  localB->Other_Mission_Data.timestamp = 0U;
  localB->Cmd_In.sp_waypoint[0] = 0.0F;
  localB->Cmd_In.cur_waypoint[0] = 0.0F;
  localB->Other_Mission_Data.type[0] = 0U;
  localB->Other_Mission_Data.valid_items[0] = 0U;
  localB->Cmd_In.sp_waypoint[1] = 0.0F;
  localB->Cmd_In.cur_waypoint[1] = 0.0F;
  localB->Other_Mission_Data.type[1] = 0U;
  localB->Other_Mission_Data.valid_items[1] = 0U;
  localB->Cmd_In.sp_waypoint[2] = 0.0F;
  localB->Cmd_In.cur_waypoint[2] = 0.0F;
  localB->Other_Mission_Data.type[2] = 0U;
  localB->Other_Mission_Data.valid_items[2] = 0U;
  for (i = 0; i < 24; i++) {
    localB->Other_Mission_Data.x[i] = 0.0F;
    localB->Other_Mission_Data.y[i] = 0.0F;
    localB->Other_Mission_Data.z[i] = 0.0F;
    localB->Other_Mission_Data.heading[i] = 0.0F;
    localB->Other_Mission_Data.ext1[i] = 0.0F;
    localB->Other_Mission_Data.ext2[i] = 0.0F;

    /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
    rty_Other_Mission_Data_x[i] = localB->Other_Mission_Data.x[i];

    /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
    rty_Other_Mission_Data_y[i] = localB->Other_Mission_Data.y[i];

    /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
    rty_Other_Mission_Data_z[i] = localB->Other_Mission_Data.z[i];

    /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
    rty_Other_Mission_Data_heading[i] = localB->Other_Mission_Data.heading[i];

    /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
    rty_Other_Mission_Data_ext1[i] = localB->Other_Mission_Data.ext1[i];

    /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
    rty_Other_Mission_Data_ext2[i] = localB->Other_Mission_Data.ext2[i];
  }

  /* End of SystemInitialize for Chart: '<Root>/FMS State Machine' */

  /* SystemInitialize for IfAction SubSystem: '<S10>/Form_Subsystem' */
  /* Start for SwitchCase: '<S12>/Switch Case' */
  localDW->SwitchCase_ActiveSubsystem_c = -1;

  /* SystemInitialize for IfAction SubSystem: '<S12>/FormAssemble' */
  /* SystemInitialize for Resettable SubSystem: '<S17>/Mission_SubSystem' */
  /* InitializeConditions for Delay: '<S26>/Delay' */
  localDW->icLoad_h = true;

  /* End of SystemInitialize for SubSystem: '<S17>/Mission_SubSystem' */
  /* End of SystemInitialize for SubSystem: '<S12>/FormAssemble' */

  /* SystemInitialize for IfAction SubSystem: '<S12>/FormHold' */
  /* InitializeConditions for Delay: '<S68>/start_vel' */
  localDW->icLoad_l3 = true;

  /* InitializeConditions for Delay: '<S64>/Delay' */
  localDW->icLoad_p = true;

  /* End of SystemInitialize for SubSystem: '<S12>/FormHold' */

  /* SystemInitialize for IfAction SubSystem: '<S12>/FormMission' */
  /* SystemInitialize for Resettable SubSystem: '<S20>/FormMission_SubSystem' */
  /* InitializeConditions for Delay: '<S156>/Delay' */
  localDW->icLoad_k = true;

  /* End of SystemInitialize for SubSystem: '<S20>/FormMission_SubSystem' */
  /* End of SystemInitialize for SubSystem: '<S12>/FormMission' */

  /* SystemInitialize for IfAction SubSystem: '<S12>/FormDisband' */
  /* InitializeConditions for Delay: '<S50>/start_vel' */
  localDW->icLoad_a = true;

  /* InitializeConditions for Delay: '<S46>/Delay' */
  localDW->icLoad_l = true;

  /* End of SystemInitialize for SubSystem: '<S12>/FormDisband' */

  /* SystemInitialize for IfAction SubSystem: '<S12>/Default' */
  Formation_FMS_Default_Init(&localDW->Default_d);

  /* End of SystemInitialize for SubSystem: '<S12>/Default' */
  /* End of SystemInitialize for SubSystem: '<S10>/Form_Subsystem' */

  /* SystemInitialize for IfAction SubSystem: '<S10>/Hold' */
  /* InitializeConditions for Delay: '<S183>/start_vel' */
  localDW->icLoad = true;

  /* InitializeConditions for Delay: '<S179>/Delay' */
  localDW->icLoad_j = true;

  /* End of SystemInitialize for SubSystem: '<S10>/Hold' */

  /* SystemInitialize for IfAction SubSystem: '<S10>/Default' */
  Formation_FMS_Default_Init(&localDW->Default);

  /* End of SystemInitialize for SubSystem: '<S10>/Default' */

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_type[0] = localB->Other_Mission_Data.type[0];

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_valid_it[0] = localB->Other_Mission_Data.valid_items[0];

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_type[1] = localB->Other_Mission_Data.type[1];

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_valid_it[1] = localB->Other_Mission_Data.valid_items[1];

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_type[2] = localB->Other_Mission_Data.type[2];

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_valid_it[2] = localB->Other_Mission_Data.valid_items[2];

  /* SystemInitialize for SignalConversion generated from: '<Root>/Other_Mission_Data' */
  *rty_Other_Mission_Data_timestam = localB->Other_Mission_Data.timestamp;

  /* SystemInitialize for Merge: '<S10>/Merge' */
  memset(&localB->Merge, 0, sizeof(FMS_Out_Bus));

  /* SystemInitialize for SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_timestamp = localB->Merge.timestamp;

  /* SystemInitialize for SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_state = localB->Merge.state;

  /* SystemInitialize for SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_ax_cmd = localB->Merge.ax_cmd;

  /* SystemInitialize for SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_ay_cmd = localB->Merge.ay_cmd;

  /* SystemInitialize for SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_vh_cmd = localB->Merge.vh_cmd;
}

/* Disable for referenced model: 'Formation_FMS' */
void Formation_FMS_Disable(DW_Formation_FMS_f_T *localDW)
{
  /* Disable for SwitchCase: '<S10>/Switch Case' */
  if (localDW->SwitchCase_ActiveSubsystem == 0) {
    /* Disable for SwitchCase: '<S12>/Switch Case' */
    localDW->SwitchCase_ActiveSubsystem_c = -1;
  }

  localDW->SwitchCase_ActiveSubsystem = -1;

  /* End of Disable for SwitchCase: '<S10>/Switch Case' */
}

/* Output and update for referenced model: 'Formation_FMS' */
void Formation_FMS(const uint32_T *rtu_Pilot_Cmd_timestamp, const uint32_T
                   *rtu_Pilot_Cmd_mode, const uint32_T
                   *rtu_Mission_Data_timestamp, const uint32_T
                   *rtu_Mission_Data_type, const uint8_T
                   *rtu_Mission_Data_valid_items, const real32_T
                   rtu_Mission_Data_x[8], const real32_T rtu_Mission_Data_y[8],
                   const real32_T rtu_Mission_Data_z[8], const uint32_T
                   *rtu_INS_Out_timestamp, const real32_T *rtu_INS_Out_phi,
                   const real32_T *rtu_INS_Out_theta, const real32_T
                   *rtu_INS_Out_psi, const real32_T *rtu_INS_Out_p, const
                   real32_T *rtu_INS_Out_q, const real32_T *rtu_INS_Out_r, const
                   real32_T rtu_INS_Out_quat[4], const real32_T *rtu_INS_Out_x_R,
                   const real32_T *rtu_INS_Out_y_R, const real32_T
                   *rtu_INS_Out_h_R, const real32_T *rtu_INS_Out_airspeed, const
                   real32_T *rtu_INS_Out_ax, const real32_T *rtu_INS_Out_ay,
                   const real32_T *rtu_INS_Out_az, const real32_T
                   *rtu_INS_Out_vn, const real32_T *rtu_INS_Out_ve, const
                   real32_T *rtu_INS_Out_vd, const uint32_T
                   rtu_Formation_Cross_timestamp[3], const real32_T
                   rtu_Formation_Cross_x_R[3], const real32_T
                   rtu_Formation_Cross_y_R[3], const real32_T
                   rtu_Formation_Cross_h_R[3], const real32_T
                   rtu_Formation_Cross_vn[3], const real32_T
                   rtu_Formation_Cross_ve[3], const real32_T
                   rtu_Formation_Cross_vd[3], const real32_T
                   rtu_Formation_Cross_ay[3], uint32_T *rty_FMS_Out_timestamp,
                   VehicleState *rty_FMS_Out_state, real32_T *rty_FMS_Out_ax_cmd,
                   real32_T *rty_FMS_Out_ay_cmd, real32_T *rty_FMS_Out_vh_cmd,
                   uint32_T *rty_Other_Mission_Data_timestam, uint32_T
                   rty_Other_Mission_Data_type[3], uint8_T
                   rty_Other_Mission_Data_valid_it[3], real32_T
                   rty_Other_Mission_Data_x[24], real32_T
                   rty_Other_Mission_Data_y[24], real32_T
                   rty_Other_Mission_Data_z[24], real32_T
                   rty_Other_Mission_Data_heading[24], real32_T
                   rty_Other_Mission_Data_ext1[24], real32_T
                   rty_Other_Mission_Data_ext2[24], B_Formation_FMS_c_T *localB,
                   DW_Formation_FMS_f_T *localDW, ZCE_Formation_FMS_T *localZCE)
{
  int8_T rtPrevAction;
  boolean_T x[3];
  boolean_T b_x[2];
  boolean_T exitg1;
  boolean_T rtb_FixPtRelationalOperator_c;

  /* BusCreator generated from: '<Root>/FMS State Machine' */
  localB->BusConversion_InsertedFor_FMS_c.timestamp = *rtu_INS_Out_timestamp;
  localB->BusConversion_InsertedFor_FMS_c.phi = *rtu_INS_Out_phi;
  localB->BusConversion_InsertedFor_FMS_c.theta = *rtu_INS_Out_theta;
  localB->BusConversion_InsertedFor_FMS_c.psi = *rtu_INS_Out_psi;
  localB->BusConversion_InsertedFor_FMS_c.p = *rtu_INS_Out_p;
  localB->BusConversion_InsertedFor_FMS_c.q = *rtu_INS_Out_q;
  localB->BusConversion_InsertedFor_FMS_c.r = *rtu_INS_Out_r;
  localB->BusConversion_InsertedFor_FMS_c.quat[0] = rtu_INS_Out_quat[0];
  localB->BusConversion_InsertedFor_FMS_c.quat[1] = rtu_INS_Out_quat[1];
  localB->BusConversion_InsertedFor_FMS_c.quat[2] = rtu_INS_Out_quat[2];
  localB->BusConversion_InsertedFor_FMS_c.quat[3] = rtu_INS_Out_quat[3];
  localB->BusConversion_InsertedFor_FMS_c.x_R = *rtu_INS_Out_x_R;
  localB->BusConversion_InsertedFor_FMS_c.y_R = *rtu_INS_Out_y_R;
  localB->BusConversion_InsertedFor_FMS_c.h_R = *rtu_INS_Out_h_R;
  localB->BusConversion_InsertedFor_FMS_c.airspeed = *rtu_INS_Out_airspeed;
  localB->BusConversion_InsertedFor_FMS_c.ax = *rtu_INS_Out_ax;
  localB->BusConversion_InsertedFor_FMS_c.ay = *rtu_INS_Out_ay;
  localB->BusConversion_InsertedFor_FMS_c.az = *rtu_INS_Out_az;
  localB->BusConversion_InsertedFor_FMS_c.vn = *rtu_INS_Out_vn;
  localB->BusConversion_InsertedFor_FMS_c.ve = *rtu_INS_Out_ve;
  localB->BusConversion_InsertedFor_FMS_c.vd = *rtu_INS_Out_vd;

  /* BusCreator generated from: '<Root>/FMS State Machine' */
  localB->BusConversion_InsertedFor_FMS_p.timestamp[0] =
    rtu_Formation_Cross_timestamp[0];
  localB->BusConversion_InsertedFor_FMS_p.x_R[0] = rtu_Formation_Cross_x_R[0];
  localB->BusConversion_InsertedFor_FMS_p.y_R[0] = rtu_Formation_Cross_y_R[0];
  localB->BusConversion_InsertedFor_FMS_p.h_R[0] = rtu_Formation_Cross_h_R[0];
  localB->BusConversion_InsertedFor_FMS_p.vn[0] = rtu_Formation_Cross_vn[0];
  localB->BusConversion_InsertedFor_FMS_p.ve[0] = rtu_Formation_Cross_ve[0];
  localB->BusConversion_InsertedFor_FMS_p.vd[0] = rtu_Formation_Cross_vd[0];
  localB->BusConversion_InsertedFor_FMS_p.ay[0] = rtu_Formation_Cross_ay[0];
  localB->BusConversion_InsertedFor_FMS_p.timestamp[1] =
    rtu_Formation_Cross_timestamp[1];
  localB->BusConversion_InsertedFor_FMS_p.x_R[1] = rtu_Formation_Cross_x_R[1];
  localB->BusConversion_InsertedFor_FMS_p.y_R[1] = rtu_Formation_Cross_y_R[1];
  localB->BusConversion_InsertedFor_FMS_p.h_R[1] = rtu_Formation_Cross_h_R[1];
  localB->BusConversion_InsertedFor_FMS_p.vn[1] = rtu_Formation_Cross_vn[1];
  localB->BusConversion_InsertedFor_FMS_p.ve[1] = rtu_Formation_Cross_ve[1];
  localB->BusConversion_InsertedFor_FMS_p.vd[1] = rtu_Formation_Cross_vd[1];
  localB->BusConversion_InsertedFor_FMS_p.ay[1] = rtu_Formation_Cross_ay[1];
  localB->BusConversion_InsertedFor_FMS_p.timestamp[2] =
    rtu_Formation_Cross_timestamp[2];
  localB->BusConversion_InsertedFor_FMS_p.x_R[2] = rtu_Formation_Cross_x_R[2];
  localB->BusConversion_InsertedFor_FMS_p.y_R[2] = rtu_Formation_Cross_y_R[2];
  localB->BusConversion_InsertedFor_FMS_p.h_R[2] = rtu_Formation_Cross_h_R[2];
  localB->BusConversion_InsertedFor_FMS_p.vn[2] = rtu_Formation_Cross_vn[2];
  localB->BusConversion_InsertedFor_FMS_p.ve[2] = rtu_Formation_Cross_ve[2];
  localB->BusConversion_InsertedFor_FMS_p.vd[2] = rtu_Formation_Cross_vd[2];
  localB->BusConversion_InsertedFor_FMS_p.ay[2] = rtu_Formation_Cross_ay[2];

  /* BusCreator generated from: '<Root>/FMS State Machine' */
  localB->BusConversion_InsertedFor_FM_pb = *rtu_Mission_Data_timestamp;
  localB->BusConversion_InsertedFor_FM_cv = *rtu_Mission_Data_type;
  for (localB->i = 0; localB->i < 8; localB->i++) {
    localB->BusConversion_InsertedFor_FMSSt[localB->i] =
      rtu_Mission_Data_x[localB->i];
    localB->BusConversion_InsertedFor_FM_cl[localB->i] =
      rtu_Mission_Data_y[localB->i];
    localB->BusConversion_InsertedFor_FMS_k[localB->i] =
      rtu_Mission_Data_z[localB->i];
  }

  /* End of BusCreator generated from: '<Root>/FMS State Machine' */

  /* Outputs for Atomic SubSystem: '<Root>/CommandProcess' */
  /* RelationalOperator: '<S7>/FixPt Relational Operator' incorporates:
   *  UnitDelay: '<S7>/Delay Input1'
   *
   * Block description for '<S7>/Delay Input1':
   *
   *  Store in Global RAM
   */
  rtb_FixPtRelationalOperator_c = (*rtu_Pilot_Cmd_timestamp !=
    localDW->DelayInput1_DSTATE);

  /* DiscreteIntegrator: '<S4>/Discrete-Time Integrator1' */
  if (rtb_FixPtRelationalOperator_c) {
    localDW->DiscreteTimeIntegrator1_DSTATE = 0.0F;
  }

  /* RelationalOperator: '<S8>/Compare' incorporates:
   *  Constant: '<S8>/Constant'
   */
  rtb_FixPtRelationalOperator_c = (*rtu_Pilot_Cmd_mode != 0U);

  /* Switch: '<S5>/Switch' incorporates:
   *  Constant: '<S6>/Constant'
   *  DiscreteIntegrator: '<S4>/Discrete-Time Integrator1'
   *  Logic: '<S5>/Logical Operator1'
   *  RelationalOperator: '<S6>/Compare'
   */
  if (rtb_FixPtRelationalOperator_c && (localDW->DiscreteTimeIntegrator1_DSTATE <
       0.5F)) {
    /* Delay: '<S5>/Delay' incorporates:
     *  DataTypeConversion: '<S5>/Data Type Conversion2'
     */
    localDW->Delay_DSTATE_j = (PilotMode)*rtu_Pilot_Cmd_mode;
  }

  /* End of Switch: '<S5>/Switch' */

  /* Update for UnitDelay: '<S7>/Delay Input1'
   *
   * Block description for '<S7>/Delay Input1':
   *
   *  Store in Global RAM
   */
  localDW->DelayInput1_DSTATE = *rtu_Pilot_Cmd_timestamp;

  /* Update for DiscreteIntegrator: '<S4>/Discrete-Time Integrator1' */
  localDW->DiscreteTimeIntegrator1_DSTATE += 0.04F;

  /* End of Outputs for SubSystem: '<Root>/CommandProcess' */

  /* Chart: '<Root>/FMS State Machine' incorporates:
   *  Constant: '<Root>/ACCEPT_R'
   *  Delay: '<S5>/Delay'
   *  MATLAB Function: '<S197>/MATLAB Function'
   */
  if (localDW->temporalCounter_i1 < 255U) {
    localDW->temporalCounter_i1++;
  }

  localDW->Mission_Data_timestamp_prev = localDW->Mission_Data_timestamp_start;
  localDW->Mission_Data_timestamp_start =
    localB->BusConversion_InsertedFor_FM_pb;
  localDW->mode_prev = localDW->mode_start;
  localDW->mode_start = localDW->Delay_DSTATE_j;
  if (localDW->is_active_c3_Formation_FMS == 0U) {
    localDW->Mission_Data_timestamp_prev =
      localB->BusConversion_InsertedFor_FM_pb;
    localDW->mode_prev = localDW->Delay_DSTATE_j;
    localDW->is_active_c3_Formation_FMS = 1U;
    x[0] = (localDW->Delay_DSTATE_j == PilotMode_FormAssemble);
    x[1] = (localDW->Delay_DSTATE_j == PilotMode_FormMission);
    x[2] = (localDW->Delay_DSTATE_j == PilotMode_FormDisband);
    rtb_FixPtRelationalOperator_c = false;
    localB->i = 0;
    exitg1 = false;
    while ((!exitg1) && (localB->i < 3)) {
      if (x[localB->i]) {
        rtb_FixPtRelationalOperator_c = true;
        exitg1 = true;
      } else {
        localB->i++;
      }
    }

    if (rtb_FixPtRelationalOperator_c) {
      localDW->is_Vehicle = Formation_FMS_IN_Formation;
      Format_enter_internal_Formation(&localB->BusConversion_InsertedFor_FMS_c,
        &localB->BusConversion_InsertedFor_FMS_p, localB, localDW);
    } else {
      localDW->is_Vehicle = Formation_FMS_IN_Standby;
      localB->Cmd_In.cur_waypoint[0] =
        localB->BusConversion_InsertedFor_FMS_c.x_R;
      localB->Cmd_In.cur_waypoint[1] =
        localB->BusConversion_InsertedFor_FMS_c.y_R;
      localB->Cmd_In.cur_waypoint[2] =
        localB->BusConversion_InsertedFor_FMS_c.h_R;
      localB->state = VehicleState_Hold;
    }
  } else if ((localDW->mode_prev != localDW->mode_start) &&
             (localDW->Delay_DSTATE_j != PilotMode_None)) {
    x[0] = (localDW->Delay_DSTATE_j == PilotMode_FormAssemble);
    x[1] = (localDW->Delay_DSTATE_j == PilotMode_FormMission);
    x[2] = (localDW->Delay_DSTATE_j == PilotMode_FormDisband);
    rtb_FixPtRelationalOperator_c = false;
    localB->i = 0;
    exitg1 = false;
    while ((!exitg1) && (localB->i < 3)) {
      if (x[localB->i]) {
        rtb_FixPtRelationalOperator_c = true;
        exitg1 = true;
      } else {
        localB->i++;
      }
    }

    if (rtb_FixPtRelationalOperator_c) {
      Formation_exit_internal_Vehicle(localDW);
      localDW->is_Vehicle = Formation_FMS_IN_Formation;
      Format_enter_internal_Formation(&localB->BusConversion_InsertedFor_FMS_c,
        &localB->BusConversion_InsertedFor_FMS_p, localB, localDW);
    } else {
      Formation_exit_internal_Vehicle(localDW);
      localDW->is_Vehicle = Formation_FMS_IN_Standby;
      localB->Cmd_In.cur_waypoint[0] =
        localB->BusConversion_InsertedFor_FMS_c.x_R;
      localB->Cmd_In.cur_waypoint[1] =
        localB->BusConversion_InsertedFor_FMS_c.y_R;
      localB->Cmd_In.cur_waypoint[2] =
        localB->BusConversion_InsertedFor_FMS_c.h_R;
      localB->state = VehicleState_Hold;
    }
  } else if (localDW->is_Vehicle == Formation_FMS_IN_Formation) {
    /* Outputs for Function Call SubSystem: '<S3>/Vehicle.Formation.check_form_valid' */
    /* MATLAB Function: '<S197>/MATLAB Function' */
    localB->BusConversion_InsertedFor_FM_pb = 0U;
    for (localB->i = 0; localB->i < 3; localB->i++) {
      localB->unit_center_to_pose[0] =
        localB->BusConversion_InsertedFor_FMS_p.x_R[localB->i] -
        localB->BusConversion_InsertedFor_FMS_c.x_R;
      localB->unit_center_to_pose[1] =
        localB->BusConversion_InsertedFor_FMS_p.y_R[localB->i] -
        localB->BusConversion_InsertedFor_FMS_c.y_R;
      localB->scale = 1.29246971E-26F;
      localB->absxk = (real32_T)fabs(localB->unit_center_to_pose[0]);
      if (localB->absxk > 1.29246971E-26F) {
        localB->rtb_vd_idx_2 = 1.0F;
        localB->scale = localB->absxk;
      } else {
        localB->t = localB->absxk / 1.29246971E-26F;
        localB->rtb_vd_idx_2 = localB->t * localB->t;
      }

      localB->absxk = (real32_T)fabs(localB->unit_center_to_pose[1]);
      if (localB->absxk > localB->scale) {
        localB->t = localB->scale / localB->absxk;
        localB->rtb_vd_idx_2 = localB->rtb_vd_idx_2 * localB->t * localB->t +
          1.0F;
        localB->scale = localB->absxk;
      } else {
        localB->t = localB->absxk / localB->scale;
        localB->rtb_vd_idx_2 += localB->t * localB->t;
      }

      if (localB->scale * (real32_T)sqrt(localB->rtb_vd_idx_2) <=
          FORMATION_PARAM.FORM_RADIUS) {
        localB->BusConversion_InsertedFor_FM_pb |= 1U << localB->i;
      } else {
        localB->BusConversion_InsertedFor_FM_pb &= ~(1U << localB->i);
      }
    }

    localB->Cmd_In.form_valid = localB->BusConversion_InsertedFor_FM_pb;

    /* End of Outputs for SubSystem: '<S3>/Vehicle.Formation.check_form_valid' */
    switch (localDW->is_Formation) {
     case Formation_FMS_IN_FormAssemble:
      switch (localDW->is_FormAssemble) {
       case Formation_FMS_IN_Hold:
        localB->state = VehicleState_FormHold;
        if ((FORMATION_PARAM.UAV_ID != 1U) && ((localB->Cmd_In.form_valid & 1U)
             != 0U)) {
          rtb_FixPtRelationalOperator_c = true;
        } else if (FORMATION_PARAM.UAV_ID == 1U) {
          localB->e[0] = (uint32_T)((localB->Cmd_In.form_valid & 1U) != 0U);
          localB->e[1] = (uint32_T)((localB->Cmd_In.form_valid & 2U) != 0U);
          localB->e[2] = (uint32_T)((localB->Cmd_In.form_valid & 4U) != 0U);
          rtb_FixPtRelationalOperator_c = true;
          localB->i = 0;
          exitg1 = false;
          while ((!exitg1) && (localB->i < 3)) {
            if (localB->e[localB->i] <= 0U) {
              rtb_FixPtRelationalOperator_c = false;
              exitg1 = true;
            } else {
              localB->i++;
            }
          }
        } else {
          rtb_FixPtRelationalOperator_c = false;
        }

        if (rtb_FixPtRelationalOperator_c) {
          localDW->is_FormAssemble = Formation_FM_IN_NO_ACTIVE_CHILD;
          localDW->is_Formation = Formation_FMS_IN_FormMission;
          localB->Cmd_In.sp_waypoint[0] =
            localB->BusConversion_InsertedFor_FMS_c.x_R;
          localB->Cmd_In.sp_waypoint[1] =
            localB->BusConversion_InsertedFor_FMS_c.y_R;
          localB->Cmd_In.sp_waypoint[2] =
            localB->BusConversion_InsertedFor_FMS_c.h_R;
          localB->state = VehicleState_FormMission;
          localDW->is_FormMission = Formation_FMS_IN_NextWP;
        }
        break;

       case Formation_FMS_IN_NextWP:
        localB->state = VehicleState_FormAssemble;
        if (localB->wp_index <= *rtu_Mission_Data_valid_items) {
          localDW->is_FormAssemble = Formation_FMS_IN_Waypoint;
          localB->Cmd_In.cur_waypoint[0] = localB->Cmd_In.sp_waypoint[0];
          localB->Cmd_In.cur_waypoint[1] = localB->Cmd_In.sp_waypoint[1];
          localB->Cmd_In.cur_waypoint[2] = localB->Cmd_In.sp_waypoint[2];
          localB->Cmd_In.sp_waypoint[0] =
            localB->BusConversion_InsertedFor_FMSSt[localB->wp_index - 1];
          localB->Cmd_In.sp_waypoint[1] =
            localB->BusConversion_InsertedFor_FM_cl[localB->wp_index - 1];
          localB->Cmd_In.sp_waypoint[2] =
            localB->BusConversion_InsertedFor_FMS_k[localB->wp_index - 1];
          b_x[0] = (localB->Cmd_In.cur_waypoint[0] == localB->
                    Cmd_In.sp_waypoint[0]);
          b_x[1] = (localB->Cmd_In.cur_waypoint[1] == localB->
                    Cmd_In.sp_waypoint[1]);
          rtb_FixPtRelationalOperator_c = true;
          localB->i = 0;
          exitg1 = false;
          while ((!exitg1) && (localB->i < 2)) {
            if (!b_x[localB->i]) {
              rtb_FixPtRelationalOperator_c = false;
              exitg1 = true;
            } else {
              localB->i++;
            }
          }

          if (rtb_FixPtRelationalOperator_c) {
            localB->Cmd_In.sp_waypoint[0] += 0.01F;
            localB->Cmd_In.sp_waypoint[1] += 0.01F;
          }
        } else {
          localDW->is_FormAssemble = Formation_FMS_IN_Hold;
          localB->Cmd_In.cur_waypoint[0] = localB->Cmd_In.sp_waypoint[0];
          localB->Cmd_In.cur_waypoint[1] = localB->Cmd_In.sp_waypoint[1];
          localB->Cmd_In.cur_waypoint[2] = localB->Cmd_In.sp_waypoint[2];
          localB->state = VehicleState_FormHold;
        }
        break;

       case Formation_FMS_IN_WaitForUpdate:
        localB->state = VehicleState_FormHold;
        if ((localDW->Mission_Data_timestamp_prev !=
             localDW->Mission_Data_timestamp_start) &&
            (localB->BusConversion_InsertedFor_FM_cv == 1U)) {
          localB->wp_index = 1U;
          localDW->is_FormAssemble = Formation_FMS_IN_NextWP;
          localB->state = VehicleState_FormAssemble;
        } else if (localDW->temporalCounter_i1 >= 250U) {
          localDW->is_FormAssemble = Formation_FM_IN_NO_ACTIVE_CHILD;
          localDW->is_Formation = Formation_FM_IN_NO_ACTIVE_CHILD;
          localDW->is_Vehicle = Formation_FMS_IN_Standby;
          localB->Cmd_In.cur_waypoint[0] =
            localB->BusConversion_InsertedFor_FMS_c.x_R;
          localB->Cmd_In.cur_waypoint[1] =
            localB->BusConversion_InsertedFor_FMS_c.y_R;
          localB->Cmd_In.cur_waypoint[2] =
            localB->BusConversion_InsertedFor_FMS_c.h_R;
          localB->state = VehicleState_Hold;
        }
        break;

       default:
        /* case IN_Waypoint: */
        localB->unit_center_to_pose[0] =
          localB->BusConversion_InsertedFor_FMS_c.x_R -
          localB->Cmd_In.sp_waypoint[0];
        localB->unit_center_to_pose[1] =
          localB->BusConversion_InsertedFor_FMS_c.y_R -
          localB->Cmd_In.sp_waypoint[1];
        if (norm_FqL2z8vo(localB->unit_center_to_pose) <= FMS_PARAM.ACCEPT_R) {
          localB->BusConversion_InsertedFor_FM_cv = localB->wp_index + 1U;
          if (localB->wp_index + 1U > 65535U) {
            localB->BusConversion_InsertedFor_FM_cv = 65535U;
          }

          localB->wp_index = (uint16_T)localB->BusConversion_InsertedFor_FM_cv;
          localDW->is_FormAssemble = Formation_FMS_IN_NextWP;
          localB->state = VehicleState_FormAssemble;
        }
        break;
      }
      break;

     case Formation_FMS_IN_FormDisband:
      localB->state = VehicleState_FormDisband;
      break;

     case Formation_FMS_IN_FormMission:
      localB->state = VehicleState_FormMission;
      if ((localDW->Mission_Data_timestamp_prev !=
           localDW->Mission_Data_timestamp_start) &&
          (localB->BusConversion_InsertedFor_FM_cv == 3U)) {
        localB->wp_index = 1U;
        localDW->is_FormMission = Formation_FMS_IN_NextWP;
      } else {
        switch (localDW->is_FormMission) {
         case Formation_FMS_IN_Follower:
          break;

         case Formation_FMS_IN_NextWP:
          if (FORMATION_PARAM.UAV_ID != 1U) {
            localDW->is_FormMission = Formation_FMS_IN_Follower;
          } else if (localB->wp_index <= *rtu_Mission_Data_valid_items) {
            localDW->is_FormMission = Formation_FMS_IN_Waypoint_g;
            localB->Cmd_In.cur_waypoint[0] = localB->Cmd_In.sp_waypoint[0];
            localB->Cmd_In.cur_waypoint[1] = localB->Cmd_In.sp_waypoint[1];
            localB->Cmd_In.cur_waypoint[2] = localB->Cmd_In.sp_waypoint[2];
            localB->Cmd_In.sp_waypoint[0] =
              localB->BusConversion_InsertedFor_FMSSt[localB->wp_index - 1];
            localB->Cmd_In.sp_waypoint[1] =
              localB->BusConversion_InsertedFor_FM_cl[localB->wp_index - 1];
            localB->Cmd_In.sp_waypoint[2] =
              localB->BusConversion_InsertedFor_FMS_k[localB->wp_index - 1];
          }
          break;

         default:
          /* case IN_Waypoint: */
          localB->unit_center_to_pose[0] =
            localB->BusConversion_InsertedFor_FMS_c.x_R -
            localB->Cmd_In.sp_waypoint[0];
          localB->unit_center_to_pose[1] =
            localB->BusConversion_InsertedFor_FMS_c.y_R -
            localB->Cmd_In.sp_waypoint[1];
          if (norm_FqL2z8vo(localB->unit_center_to_pose) <= FMS_PARAM.ACCEPT_R)
          {
            localB->BusConversion_InsertedFor_FM_cv = localB->wp_index + 1U;
            if (localB->wp_index + 1U > 65535U) {
              localB->BusConversion_InsertedFor_FM_cv = 65535U;
            }

            localB->wp_index = (uint16_T)localB->BusConversion_InsertedFor_FM_cv;
            localDW->is_FormMission = Formation_FMS_IN_NextWP;
          }
          break;
        }
      }
      break;

     default:
      /* case IN_InvalidMode: */
      localDW->is_Formation = Formation_FM_IN_NO_ACTIVE_CHILD;
      localDW->is_Vehicle = Formation_FMS_IN_Standby;
      localB->Cmd_In.cur_waypoint[0] =
        localB->BusConversion_InsertedFor_FMS_c.x_R;
      localB->Cmd_In.cur_waypoint[1] =
        localB->BusConversion_InsertedFor_FMS_c.y_R;
      localB->Cmd_In.cur_waypoint[2] =
        localB->BusConversion_InsertedFor_FMS_c.h_R;
      localB->state = VehicleState_Hold;
      break;
    }
  } else {
    /* case IN_Standby: */
    localB->state = VehicleState_Hold;
  }

  /* End of Chart: '<Root>/FMS State Machine' */

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  memcpy(&rty_Other_Mission_Data_x[0], &localB->Other_Mission_Data.x[0], 24U *
         sizeof(real32_T));

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  memcpy(&rty_Other_Mission_Data_y[0], &localB->Other_Mission_Data.y[0], 24U *
         sizeof(real32_T));

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  memcpy(&rty_Other_Mission_Data_z[0], &localB->Other_Mission_Data.z[0], 24U *
         sizeof(real32_T));

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  memcpy(&rty_Other_Mission_Data_heading[0], &localB->
         Other_Mission_Data.heading[0], 24U * sizeof(real32_T));

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  memcpy(&rty_Other_Mission_Data_ext1[0], &localB->Other_Mission_Data.ext1[0],
         24U * sizeof(real32_T));

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  memcpy(&rty_Other_Mission_Data_ext2[0], &localB->Other_Mission_Data.ext2[0],
         24U * sizeof(real32_T));

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_vd_idx_0 = rtu_Formation_Cross_vd[0];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_ve_idx_0 = rtu_Formation_Cross_ve[0];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->vn[0] = rtu_Formation_Cross_vn[0];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->h_R[0] = rtu_Formation_Cross_h_R[0];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_y_R_idx_0 = rtu_Formation_Cross_y_R[0];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_x_R_idx_0 = rtu_Formation_Cross_x_R[0];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_vd_idx_1 = rtu_Formation_Cross_vd[1];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->u_n = rtu_Formation_Cross_ve[1];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->vn[1] = rtu_Formation_Cross_vn[1];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->h_R[1] = rtu_Formation_Cross_h_R[1];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_y_R_idx_1 = rtu_Formation_Cross_y_R[1];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_x_R_idx_1 = rtu_Formation_Cross_x_R[1];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->rtb_vd_idx_2 = rtu_Formation_Cross_vd[2];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->DotProduct = rtu_Formation_Cross_ve[2];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->vn[2] = rtu_Formation_Cross_vn[2];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->h_R[2] = rtu_Formation_Cross_h_R[2];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->scale = rtu_Formation_Cross_y_R[2];

  /* SignalConversion generated from: '<S9>/Signal Copy2' */
  localB->absxk = rtu_Formation_Cross_x_R[2];

  /* SwitchCase: '<S10>/Switch Case' incorporates:
   *  Merge: '<S10>/Merge'
   *  Product: '<S193>/Divide'
   *  Product: '<S194>/Divide'
   */
  rtPrevAction = localDW->SwitchCase_ActiveSubsystem;
  switch (localB->state) {
   case VehicleState_FormAssemble:
   case VehicleState_FormHold:
   case VehicleState_FormMission:
   case VehicleState_FormDisband:
    localDW->SwitchCase_ActiveSubsystem = 0;
    break;

   case VehicleState_Hold:
    localDW->SwitchCase_ActiveSubsystem = 1;
    break;

   default:
    localDW->SwitchCase_ActiveSubsystem = 2;
    break;
  }

  if ((rtPrevAction != localDW->SwitchCase_ActiveSubsystem) && (rtPrevAction ==
       0)) {
    /* SwitchCase: '<S12>/Switch Case' */
    localDW->SwitchCase_ActiveSubsystem_c = -1;
  }

  switch (localDW->SwitchCase_ActiveSubsystem) {
   case 0:
    /* Outputs for IfAction SubSystem: '<S10>/Form_Subsystem' incorporates:
     *  ActionPort: '<S12>/Action Port'
     */
    /* SwitchCase: '<S12>/Switch Case' incorporates:
     *  Gain: '<S92>/POS_KP'
     *  MATLAB Function: '<S89>/Calc_Position_Velocity_Setpoint'
     *  Merge: '<S10>/Merge'
     *  Product: '<S43>/Divide'
     *  Product: '<S60>/Divide'
     *  Product: '<S61>/Divide'
     *  Product: '<S78>/Divide'
     *  Product: '<S79>/Divide'
     *  Sum: '<S31>/Subtract'
     */
    rtPrevAction = localDW->SwitchCase_ActiveSubsystem_c;
    switch (localB->state) {
     case VehicleState_FormAssemble:
      localDW->SwitchCase_ActiveSubsystem_c = 0;
      break;

     case VehicleState_FormHold:
      localDW->SwitchCase_ActiveSubsystem_c = 1;
      break;

     case VehicleState_FormMission:
      localDW->SwitchCase_ActiveSubsystem_c = 2;
      break;

     case VehicleState_FormDisband:
      localDW->SwitchCase_ActiveSubsystem_c = 3;
      break;

     default:
      localDW->SwitchCase_ActiveSubsystem_c = 4;
      break;
    }

    switch (localDW->SwitchCase_ActiveSubsystem_c) {
     case 0:
      /* Outputs for IfAction SubSystem: '<S12>/FormAssemble' incorporates:
       *  ActionPort: '<S17>/Action Port'
       */
      /* RelationalOperator: '<S23>/FixPt Relational Operator' incorporates:
       *  UnitDelay: '<S23>/Delay Input1'
       *
       * Block description for '<S23>/Delay Input1':
       *
       *  Store in Global RAM
       */
      rtb_FixPtRelationalOperator_c = (localB->wp_index !=
        localDW->DelayInput1_DSTATE_h);

      /* Outputs for Resettable SubSystem: '<S17>/Mission_SubSystem' incorporates:
       *  ResetPort: '<S24>/Reset'
       */
      /* InitializeConditions for Delay: '<S26>/Delay' */
      localDW->icLoad_h = ((rtb_FixPtRelationalOperator_c &&
                            (localZCE->Mission_SubSystem_Reset_ZCE != POS_ZCSIG))
                           || localDW->icLoad_h);
      localZCE->Mission_SubSystem_Reset_ZCE = rtb_FixPtRelationalOperator_c;

      /* SignalConversion generated from: '<S39>/Square' */
      localB->TmpSignalConversionAtSqua_c[0] = *rtu_INS_Out_vn;
      localB->TmpSignalConversionAtSqua_c[1] = *rtu_INS_Out_ve;

      /* Sum: '<S40>/Sum of Elements' incorporates:
       *  Math: '<S40>/Math Function'
       *  Sum: '<S39>/Sum of Elements'
       */
      localB->rtb_y_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0] +
        localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1];

      /* Math: '<S40>/Math Function1' incorporates:
       *  Sum: '<S40>/Sum of Elements'
       *
       * About '<S40>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_y_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_0 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_y_R_idx_0));
      } else {
        localB->rtb_x_R_idx_0 = (real32_T)sqrt(localB->rtb_y_R_idx_0);
      }

      /* End of Math: '<S40>/Math Function1' */

      /* Switch: '<S40>/Switch' incorporates:
       *  Constant: '<S40>/Constant'
       *  Product: '<S40>/Product'
       */
      if (localB->rtb_x_R_idx_0 > 0.0F) {
        localB->rtb_vd_idx_0 = *rtu_INS_Out_vn;
        localB->rtb_vd_idx_1 = *rtu_INS_Out_ve;
        localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_0;
      } else {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->rtb_vd_idx_2 = 1.0F;
      }

      /* End of Switch: '<S40>/Switch' */

      /* Reshape: '<S29>/Reshape2' */
      localB->Reshape2_bq[0] = *rtu_INS_Out_x_R;
      localB->Reshape2_bq[1] = *rtu_INS_Out_y_R;

      /* MATLAB Function: '<S30>/NearbyRefWP' incorporates:
       *  Constant: '<S28>/L1'
       */
      Formation_FMS_NearbyRefWP(&localB->Cmd_In.sp_waypoint[0],
        localB->Reshape2_bq, FMS_PARAM.L1, localB->TmpSignalConversionAtSqua_c,
        &localB->rtb_ve_idx_0);

      /* MATLAB Function: '<S30>/SearchL1RefWP' incorporates:
       *  Constant: '<S28>/L1'
       */
      Formation_FMS_SearchL1RefWP(&localB->Cmd_In.cur_waypoint[0],
        &localB->Cmd_In.sp_waypoint[0], localB->Reshape2_bq, FMS_PARAM.L1,
        localB->P_d, &localB->u_n);

      /* MATLAB Function: '<S30>/OutRegionRegWP' */
      Formation_FMS_OutRegionRegWP(&localB->Cmd_In.cur_waypoint[0],
        &localB->Cmd_In.sp_waypoint[0], localB->Reshape2_bq, localB->P_i);

      /* Switch: '<S30>/Switch1' incorporates:
       *  Constant: '<S32>/Constant'
       *  Constant: '<S33>/Constant'
       *  Product: '<S43>/Divide'
       *  RelationalOperator: '<S32>/Compare'
       *  RelationalOperator: '<S33>/Compare'
       *  Switch: '<S30>/Switch'
       */
      if (localB->rtb_ve_idx_0 > 0.0F) {
        localB->P_d[0] = localB->TmpSignalConversionAtSqua_c[0];
        localB->P_d[1] = localB->TmpSignalConversionAtSqua_c[1];
      } else if (!(localB->u_n >= 0.0F)) {
        /* Product: '<S43>/Divide' incorporates:
         *  Switch: '<S30>/Switch'
         */
        localB->P_d[0] = localB->P_i[0];
        localB->P_d[1] = localB->P_i[1];
      }

      /* End of Switch: '<S30>/Switch1' */

      /* Sum: '<S31>/Subtract' incorporates:
       *  Product: '<S43>/Divide'
       *  Reshape: '<S29>/Reshape2'
       */
      localB->u_n = localB->P_d[0] - localB->Reshape2_bq[0];

      /* Math: '<S41>/Math Function' */
      localB->TmpSignalConversionAtSqua_c[0] = localB->u_n * localB->u_n;
      localB->P_d[0] = localB->u_n;

      /* Sum: '<S31>/Subtract' incorporates:
       *  Product: '<S43>/Divide'
       *  Reshape: '<S29>/Reshape2'
       */
      localB->u_n = localB->P_d[1] - localB->Reshape2_bq[1];

      /* Sum: '<S41>/Sum of Elements' incorporates:
       *  Math: '<S41>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->u_n * localB->u_n +
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S41>/Math Function1' incorporates:
       *  Sum: '<S41>/Sum of Elements'
       *
       * About '<S41>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_0 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_0 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S41>/Math Function1' */

      /* Switch: '<S41>/Switch' incorporates:
       *  Constant: '<S41>/Constant'
       *  Product: '<S41>/Product'
       *  Sum: '<S31>/Subtract'
       *  Switch: '<S44>/Switch'
       */
      if (localB->rtb_x_R_idx_0 > 0.0F) {
        localB->rtb_ve_idx_0 = localB->P_d[0];
        localB->DotProduct = localB->rtb_x_R_idx_0;
      } else {
        localB->rtb_ve_idx_0 = localB->P_d[0] * 0.0F;
        localB->u_n *= 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S41>/Switch' */

      /* Product: '<S40>/Divide' */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_vd_idx_0 /
        localB->rtb_vd_idx_2;
      localB->TmpSignalConversionAtSqua_c[1] = localB->rtb_vd_idx_1 /
        localB->rtb_vd_idx_2;

      /* Sum: '<S43>/Sum of Elements' incorporates:
       *  Math: '<S43>/Math Function'
       *  SignalConversion generated from: '<S43>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1] +
        localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S43>/Math Function1' incorporates:
       *  Sum: '<S43>/Sum of Elements'
       *
       * About '<S43>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_0 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_0 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S43>/Math Function1' */

      /* Switch: '<S43>/Switch' incorporates:
       *  Constant: '<S43>/Constant'
       *  Product: '<S43>/Product'
       *  SignalConversion generated from: '<S43>/Math Function'
       */
      if (localB->rtb_x_R_idx_0 > 0.0F) {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0];
        localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_0;
      } else {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->rtb_vd_idx_2 = 1.0F;
      }

      /* End of Switch: '<S43>/Switch' */

      /* Product: '<S41>/Divide' incorporates:
       *  Product: '<S44>/Divide'
       */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_ve_idx_0 /
        localB->DotProduct;
      localB->TmpSignalConversionAtSqua_c[1] = localB->u_n / localB->DotProduct;

      /* Sum: '<S44>/Sum of Elements' incorporates:
       *  Math: '<S44>/Math Function'
       *  SignalConversion generated from: '<S44>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1] +
        localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S44>/Math Function1' incorporates:
       *  Sum: '<S44>/Sum of Elements'
       *
       * About '<S44>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_0 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_0 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S44>/Math Function1' */

      /* Switch: '<S44>/Switch' incorporates:
       *  Constant: '<S44>/Constant'
       *  Product: '<S44>/Product'
       *  SignalConversion generated from: '<S44>/Math Function'
       */
      if (localB->rtb_x_R_idx_0 > 0.0F) {
        localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
        localB->u_n = localB->TmpSignalConversionAtSqua_c[0];
        localB->DotProduct = localB->rtb_x_R_idx_0;
      } else {
        localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S44>/Switch' */

      /* Product: '<S44>/Divide' */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_ve_idx_0 /
        localB->DotProduct;

      /* Product: '<S43>/Divide' */
      localB->P_d[0] = localB->rtb_vd_idx_0 / localB->rtb_vd_idx_2;

      /* Product: '<S44>/Divide' */
      localB->TmpSignalConversionAtSqua_c[1] = localB->u_n / localB->DotProduct;

      /* Product: '<S43>/Divide' */
      localB->P_d[1] = localB->rtb_vd_idx_1 / localB->rtb_vd_idx_2;

      /* Sqrt: '<S39>/Sqrt' */
      localB->rtb_x_R_idx_0 = (real32_T)sqrt(localB->rtb_y_R_idx_0);

      /* Math: '<S37>/Square' */
      localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_0 * localB->rtb_x_R_idx_0;

      /* Sum: '<S42>/Subtract' incorporates:
       *  Product: '<S42>/Multiply'
       *  Product: '<S42>/Multiply1'
       */
      localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * localB->P_d[1] -
        localB->P_d[0] * localB->TmpSignalConversionAtSqua_c[1];

      /* Signum: '<S38>/Sign1' */
      if (rtIsNaNF(localB->u_n)) {
        localB->rtb_x_R_idx_0 = (rtNaNF);
      } else if (localB->u_n < 0.0F) {
        localB->rtb_x_R_idx_0 = -1.0F;
      } else {
        localB->rtb_x_R_idx_0 = (real32_T)(localB->u_n > 0.0F);
      }

      /* End of Signum: '<S38>/Sign1' */

      /* Switch: '<S38>/Switch2' incorporates:
       *  Constant: '<S38>/Constant4'
       */
      if (localB->rtb_x_R_idx_0 != 0.0F) {
        localB->u_n = localB->rtb_x_R_idx_0;
      } else {
        localB->u_n = 1.0F;
      }

      /* End of Switch: '<S38>/Switch2' */

      /* Sum: '<S25>/Sum' incorporates:
       *  Constant: '<S25>/Constant'
       */
      localB->rtb_x_R_idx_0 = FMS_PARAM.FW_AIRSPD_TRIM - *rtu_INS_Out_airspeed;

      /* Delay: '<S26>/Delay' incorporates:
       *  Constant: '<S26>/Constant'
       */
      if (localDW->icLoad_h) {
        localDW->Delay_DSTATE_l = FMS_PARAM.FW_HEIGHT_TRIM;
      }

      /* Sum: '<S26>/Sum' incorporates:
       *  Delay: '<S26>/Delay'
       */
      localB->rtb_ve_idx_0 = localDW->Delay_DSTATE_l - *rtu_INS_Out_h_R;

      /* End of Outputs for SubSystem: '<S17>/Mission_SubSystem' */
      /* End of Outputs for SubSystem: '<S12>/FormAssemble' */
      memset(&localB->Merge, 0, sizeof(FMS_Out_Bus));

      /* Outputs for IfAction SubSystem: '<S12>/FormAssemble' incorporates:
       *  ActionPort: '<S17>/Action Port'
       */
      /* Outputs for Resettable SubSystem: '<S17>/Mission_SubSystem' incorporates:
       *  ResetPort: '<S24>/Reset'
       */
      /* BusAssignment: '<S24>/Bus Assignment' incorporates:
       *  Constant: '<S24>/Constant1'
       *  Merge: '<S10>/Merge'
       */
      localB->Merge.state = VehicleState_FormAssemble;
      localB->Merge.ax_cmd = localB->rtb_x_R_idx_0;

      /* DotProduct: '<S38>/Dot Product' */
      localB->rtb_x_R_idx_0 = localB->P_d[0] *
        localB->TmpSignalConversionAtSqua_c[0] + localB->P_d[1] *
        localB->TmpSignalConversionAtSqua_c[1];

      /* Trigonometry: '<S38>/Acos' incorporates:
       *  DotProduct: '<S38>/Dot Product'
       */
      if (localB->rtb_x_R_idx_0 > 1.0F) {
        localB->rtb_x_R_idx_0 = 1.0F;
      } else if (localB->rtb_x_R_idx_0 < -1.0F) {
        localB->rtb_x_R_idx_0 = -1.0F;
      }

      /* Product: '<S38>/Multiply' incorporates:
       *  Trigonometry: '<S38>/Acos'
       */
      localB->rtb_x_R_idx_0 = (real32_T)acos(localB->rtb_x_R_idx_0) *
        localB->u_n;

      /* Saturate: '<S37>/Saturation' */
      if (localB->rtb_x_R_idx_0 > 1.57079637F) {
        localB->rtb_x_R_idx_0 = 1.57079637F;
      } else if (localB->rtb_x_R_idx_0 < -1.57079637F) {
        localB->rtb_x_R_idx_0 = -1.57079637F;
      }

      /* Product: '<S37>/Divide' incorporates:
       *  Constant: '<S28>/L1'
       *  Gain: '<S37>/Gain'
       *  Product: '<S37>/Multiply1'
       *  Saturate: '<S37>/Saturation'
       *  Trigonometry: '<S37>/Sin'
       */
      localB->rtb_x_R_idx_0 = 2.0F * localB->rtb_vd_idx_2 * (real32_T)sin
        (localB->rtb_x_R_idx_0) / FMS_PARAM.L1;

      /* Saturate: '<S28>/Saturation1' */
      if (localB->rtb_x_R_idx_0 > FMS_PARAM.ACC_Y_LIM) {
        /* BusAssignment: '<S24>/Bus Assignment' */
        localB->Merge.ay_cmd = FMS_PARAM.ACC_Y_LIM;
      } else if (localB->rtb_x_R_idx_0 < -FMS_PARAM.ACC_Y_LIM) {
        /* BusAssignment: '<S24>/Bus Assignment' */
        localB->Merge.ay_cmd = -FMS_PARAM.ACC_Y_LIM;
      } else {
        /* BusAssignment: '<S24>/Bus Assignment' */
        localB->Merge.ay_cmd = localB->rtb_x_R_idx_0;
      }

      /* End of Saturate: '<S28>/Saturation1' */

      /* Gain: '<S26>/Gain2' */
      localB->rtb_x_R_idx_0 = FMS_PARAM.Z_P * localB->rtb_ve_idx_0;

      /* Saturate: '<S26>/Saturation' */
      if (localB->rtb_x_R_idx_0 > CONTROL_PARAM.FW_T_CLMB_MAX) {
        /* BusAssignment: '<S24>/Bus Assignment' */
        localB->Merge.vh_cmd = CONTROL_PARAM.FW_T_CLMB_MAX;
      } else if (localB->rtb_x_R_idx_0 < -CONTROL_PARAM.FW_T_SINK_MAX) {
        /* BusAssignment: '<S24>/Bus Assignment' */
        localB->Merge.vh_cmd = -CONTROL_PARAM.FW_T_SINK_MAX;
      } else {
        /* BusAssignment: '<S24>/Bus Assignment' */
        localB->Merge.vh_cmd = localB->rtb_x_R_idx_0;
      }

      /* End of Saturate: '<S26>/Saturation' */

      /* Update for Delay: '<S26>/Delay' */
      localDW->icLoad_h = false;

      /* End of Outputs for SubSystem: '<S17>/Mission_SubSystem' */

      /* Update for UnitDelay: '<S23>/Delay Input1'
       *
       * Block description for '<S23>/Delay Input1':
       *
       *  Store in Global RAM
       */
      localDW->DelayInput1_DSTATE_h = localB->wp_index;

      /* End of Outputs for SubSystem: '<S12>/FormAssemble' */
      break;

     case 1:
      if (localDW->SwitchCase_ActiveSubsystem_c != rtPrevAction) {
        /* InitializeConditions for IfAction SubSystem: '<S12>/FormHold' incorporates:
         *  ActionPort: '<S19>/Action Port'
         */
        /* InitializeConditions for SwitchCase: '<S12>/Switch Case' incorporates:
         *  Delay: '<S64>/Delay'
         *  Delay: '<S68>/start_vel'
         */
        localDW->icLoad_l3 = true;
        localDW->icLoad_p = true;

        /* End of InitializeConditions for SubSystem: '<S12>/FormHold' */
      }

      /* Outputs for IfAction SubSystem: '<S12>/FormHold' incorporates:
       *  ActionPort: '<S19>/Action Port'
       */
      /* SignalConversion generated from: '<S74>/Square' */
      localB->TmpSignalConversionAtSqua_c[0] = *rtu_INS_Out_vn;
      localB->TmpSignalConversionAtSqua_c[1] = *rtu_INS_Out_ve;

      /* Sum: '<S75>/Sum of Elements' incorporates:
       *  Math: '<S75>/Math Function'
       *  Sum: '<S74>/Sum of Elements'
       */
      localB->rtb_y_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0] +
        localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1];

      /* Math: '<S75>/Math Function1' incorporates:
       *  Sum: '<S75>/Sum of Elements'
       *
       * About '<S75>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_y_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_y_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_y_R_idx_0);
      }

      /* End of Math: '<S75>/Math Function1' */

      /* Switch: '<S75>/Switch' incorporates:
       *  Constant: '<S75>/Constant'
       *  Product: '<S75>/Product'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_vd_idx_0 = *rtu_INS_Out_vn;
        localB->rtb_vd_idx_1 = *rtu_INS_Out_ve;
        localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->rtb_vd_idx_2 = 1.0F;
      }

      /* End of Switch: '<S75>/Switch' */

      /* Delay: '<S68>/start_vel' */
      if (localDW->icLoad_l3) {
        localDW->start_vel_DSTATE_n[0] = localB->TmpSignalConversionAtSqua_c[0];
        localDW->start_vel_DSTATE_n[1] = localB->TmpSignalConversionAtSqua_c[1];
      }

      /* Math: '<S80>/Math Function' incorporates:
       *  Delay: '<S68>/start_vel'
       */
      localB->TmpSignalConversionAtSqua_c[0] = localDW->start_vel_DSTATE_n[0] *
        localDW->start_vel_DSTATE_n[0];
      localB->TmpSignalConversionAtSqua_c[1] = localDW->start_vel_DSTATE_n[1] *
        localDW->start_vel_DSTATE_n[1];

      /* Sum: '<S80>/Sum of Elements' */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] +
        localB->TmpSignalConversionAtSqua_c[1];

      /* Math: '<S80>/Math Function1' incorporates:
       *  Sum: '<S80>/Sum of Elements'
       *
       * About '<S80>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S80>/Math Function1' */

      /* Switch: '<S80>/Switch' incorporates:
       *  Constant: '<S80>/Constant'
       *  Delay: '<S68>/start_vel'
       *  Product: '<S80>/Product'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_ve_idx_0 = localDW->start_vel_DSTATE_n[0];
        localB->u_n = localDW->start_vel_DSTATE_n[1];
        localB->DotProduct = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_ve_idx_0 = localDW->start_vel_DSTATE_n[0] * 0.0F;
        localB->u_n = localDW->start_vel_DSTATE_n[1] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S80>/Switch' */

      /* Product: '<S75>/Divide' */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_vd_idx_0 /
        localB->rtb_vd_idx_2;
      localB->TmpSignalConversionAtSqua_c[1] = localB->rtb_vd_idx_1 /
        localB->rtb_vd_idx_2;

      /* SignalConversion generated from: '<S78>/Math Function' */
      localB->P_d[0] = localB->TmpSignalConversionAtSqua_c[1];
      localB->P_d[1] = localB->TmpSignalConversionAtSqua_c[0];

      /* Sum: '<S78>/Sum of Elements' incorporates:
       *  Math: '<S78>/Math Function'
       *  SignalConversion generated from: '<S78>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1] +
        localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S78>/Math Function1' incorporates:
       *  Sum: '<S78>/Sum of Elements'
       *
       * About '<S78>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S78>/Math Function1' */

      /* Switch: '<S78>/Switch' incorporates:
       *  Constant: '<S78>/Constant'
       *  Product: '<S78>/Product'
       *  SignalConversion generated from: '<S78>/Math Function'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0];
        localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->rtb_vd_idx_2 = 1.0F;
      }

      /* End of Switch: '<S78>/Switch' */

      /* Gain: '<S66>/Gain' incorporates:
       *  Constant: '<S65>/L1'
       */
      localB->rtb_x_R_idx_1 = 0.5F * FMS_PARAM.L1;

      /* MinMax: '<S66>/Max' incorporates:
       *  Constant: '<S65>/R'
       */
      if ((FMS_PARAM.LOITER_R >= localB->rtb_x_R_idx_1) || rtIsNaNF
          (localB->rtb_x_R_idx_1)) {
        localB->rtb_x_R_idx_1 = FMS_PARAM.LOITER_R;
      }

      /* End of MinMax: '<S66>/Max' */

      /* Reshape: '<S66>/Reshape2' */
      localB->Reshape2_bq[0] = *rtu_INS_Out_x_R;
      localB->Reshape2_bq[1] = *rtu_INS_Out_y_R;

      /* MATLAB Function: '<S66>/SearchL1RefWP' incorporates:
       *  Constant: '<S65>/L1'
       */
      Formation_FMS_SearchL1RefWP_i(&localB->Cmd_In.cur_waypoint[0],
        localB->Reshape2_bq, localB->rtb_x_R_idx_1, FMS_PARAM.L1,
        localB->TmpSignalConversionAtSqua_c, &localB->n_j);

      /* RelationalOperator: '<S69>/Compare' incorporates:
       *  Constant: '<S69>/Constant'
       */
      rtb_FixPtRelationalOperator_c = (localB->n_j > 0.0);

      /* Product: '<S80>/Divide' */
      localB->P_i[0] = localB->rtb_ve_idx_0 / localB->DotProduct;
      localB->P_i[1] = localB->u_n / localB->DotProduct;

      /* MATLAB Function: '<S66>/OutRegionRegWP' incorporates:
       *  Constant: '<S65>/L1'
       */
      Formation_FMS_OutRegionRegWP_o(&localB->Cmd_In.cur_waypoint[0],
        localB->Reshape2_bq, localB->rtb_x_R_idx_1, FMS_PARAM.L1, localB->P_i,
        localB->P_d);

      /* Switch: '<S66>/Switch' */
      if (rtb_FixPtRelationalOperator_c) {
        localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0];
      } else {
        localB->rtb_x_R_idx_0 = localB->P_d[0];
      }

      /* Sum: '<S67>/Subtract' incorporates:
       *  Switch: '<S66>/Switch'
       */
      localB->Reshape2_bq[0] = localB->rtb_x_R_idx_0 - *rtu_INS_Out_x_R;

      /* Switch: '<S66>/Switch' */
      if (rtb_FixPtRelationalOperator_c) {
        localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
      } else {
        localB->rtb_x_R_idx_0 = localB->P_d[1];
      }

      /* Sum: '<S67>/Subtract' incorporates:
       *  Switch: '<S66>/Switch'
       */
      localB->Reshape2_bq[1] = localB->rtb_x_R_idx_0 - *rtu_INS_Out_y_R;

      /* Sum: '<S76>/Sum of Elements' incorporates:
       *  Math: '<S76>/Math Function'
       *  Sum: '<S67>/Subtract'
       */
      localB->rtb_x_R_idx_0 = localB->Reshape2_bq[0] * localB->Reshape2_bq[0] +
        localB->Reshape2_bq[1] * localB->Reshape2_bq[1];

      /* Math: '<S76>/Math Function1' incorporates:
       *  Sum: '<S76>/Sum of Elements'
       *
       * About '<S76>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S76>/Math Function1' */

      /* Switch: '<S76>/Switch' incorporates:
       *  Constant: '<S76>/Constant'
       *  Product: '<S76>/Product'
       *  Sum: '<S67>/Subtract'
       *  Switch: '<S79>/Switch'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_ve_idx_0 = localB->Reshape2_bq[0];
        localB->u_n = localB->Reshape2_bq[1];
        localB->DotProduct = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_ve_idx_0 = localB->Reshape2_bq[0] * 0.0F;
        localB->u_n = localB->Reshape2_bq[1] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S76>/Switch' */

      /* Product: '<S76>/Divide' incorporates:
       *  Product: '<S79>/Divide'
       */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_ve_idx_0 /
        localB->DotProduct;
      localB->TmpSignalConversionAtSqua_c[1] = localB->u_n / localB->DotProduct;

      /* Sum: '<S79>/Sum of Elements' incorporates:
       *  Math: '<S79>/Math Function'
       *  SignalConversion generated from: '<S79>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1] +
        localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S79>/Math Function1' incorporates:
       *  Sum: '<S79>/Sum of Elements'
       *
       * About '<S79>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S79>/Math Function1' */

      /* Switch: '<S79>/Switch' incorporates:
       *  Constant: '<S79>/Constant'
       *  Product: '<S79>/Product'
       *  SignalConversion generated from: '<S79>/Math Function'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
        localB->u_n = localB->TmpSignalConversionAtSqua_c[0];
        localB->DotProduct = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S79>/Switch' */

      /* Product: '<S79>/Divide' */
      localB->rtb_y_R_idx_1 = localB->rtb_ve_idx_0 / localB->DotProduct;

      /* Product: '<S78>/Divide' */
      localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_0 / localB->rtb_vd_idx_2;

      /* DotProduct: '<S73>/Dot Product' */
      localB->rtb_ve_idx_0 = localB->rtb_x_R_idx_0 * localB->rtb_y_R_idx_1;
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_y_R_idx_1;
      localB->P_d[0] = localB->rtb_x_R_idx_0;

      /* Product: '<S79>/Divide' incorporates:
       *  Product: '<S78>/Divide'
       */
      localB->rtb_y_R_idx_1 = localB->u_n / localB->DotProduct;

      /* Product: '<S78>/Divide' */
      localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_1 / localB->rtb_vd_idx_2;

      /* DotProduct: '<S73>/Dot Product' */
      localB->rtb_ve_idx_0 += localB->rtb_x_R_idx_0 * localB->rtb_y_R_idx_1;

      /* Sqrt: '<S74>/Sqrt' */
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_y_R_idx_0);

      /* Math: '<S72>/Square' */
      localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1 * localB->rtb_x_R_idx_1;

      /* Sum: '<S77>/Subtract' incorporates:
       *  Product: '<S77>/Multiply'
       *  Product: '<S77>/Multiply1'
       */
      localB->u_n = localB->TmpSignalConversionAtSqua_c[0] *
        localB->rtb_x_R_idx_0 - localB->P_d[0] * localB->rtb_y_R_idx_1;

      /* Signum: '<S73>/Sign1' */
      if (rtIsNaNF(localB->u_n)) {
        localB->rtb_x_R_idx_1 = (rtNaNF);
      } else if (localB->u_n < 0.0F) {
        localB->rtb_x_R_idx_1 = -1.0F;
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)(localB->u_n > 0.0F);
      }

      /* End of Signum: '<S73>/Sign1' */

      /* Trigonometry: '<S73>/Acos' incorporates:
       *  DotProduct: '<S73>/Dot Product'
       */
      if (localB->rtb_ve_idx_0 > 1.0F) {
        localB->rtb_ve_idx_0 = 1.0F;
      } else if (localB->rtb_ve_idx_0 < -1.0F) {
        localB->rtb_ve_idx_0 = -1.0F;
      }

      /* Switch: '<S73>/Switch2' incorporates:
       *  Constant: '<S73>/Constant4'
       */
      if (!(localB->rtb_x_R_idx_1 != 0.0F)) {
        localB->rtb_x_R_idx_1 = 1.0F;
      }

      /* Product: '<S73>/Multiply' incorporates:
       *  Switch: '<S73>/Switch2'
       *  Trigonometry: '<S73>/Acos'
       */
      localB->rtb_x_R_idx_0 = (real32_T)acos(localB->rtb_ve_idx_0) *
        localB->rtb_x_R_idx_1;

      /* Sum: '<S63>/Sum' incorporates:
       *  Constant: '<S63>/Constant'
       */
      localB->rtb_x_R_idx_1 = FMS_PARAM.FW_AIRSPD_TRIM - *rtu_INS_Out_airspeed;

      /* Delay: '<S64>/Delay' incorporates:
       *  Constant: '<S64>/Constant'
       */
      if (localDW->icLoad_p) {
        localDW->Delay_DSTATE_a = FMS_PARAM.FW_HEIGHT_TRIM;
      }

      /* Sum: '<S64>/Sum' incorporates:
       *  Delay: '<S64>/Delay'
       */
      localB->rtb_ve_idx_0 = localDW->Delay_DSTATE_a - *rtu_INS_Out_h_R;

      /* End of Outputs for SubSystem: '<S12>/FormHold' */
      memset(&localB->Merge, 0, sizeof(FMS_Out_Bus));

      /* Outputs for IfAction SubSystem: '<S12>/FormHold' incorporates:
       *  ActionPort: '<S19>/Action Port'
       */
      /* BusAssignment: '<S19>/Bus Assignment' incorporates:
       *  Constant: '<S19>/Constant1'
       *  Merge: '<S10>/Merge'
       */
      localB->Merge.state = VehicleState_FormHold;
      localB->Merge.ax_cmd = localB->rtb_x_R_idx_1;

      /* Saturate: '<S72>/Saturation' */
      if (localB->rtb_x_R_idx_0 > 1.57079637F) {
        localB->rtb_x_R_idx_0 = 1.57079637F;
      } else if (localB->rtb_x_R_idx_0 < -1.57079637F) {
        localB->rtb_x_R_idx_0 = -1.57079637F;
      }

      /* BusAssignment: '<S19>/Bus Assignment' incorporates:
       *  Constant: '<S65>/L1'
       *  Gain: '<S72>/Gain'
       *  Product: '<S72>/Divide'
       *  Product: '<S72>/Multiply1'
       *  Saturate: '<S72>/Saturation'
       *  Trigonometry: '<S72>/Sin'
       */
      localB->Merge.ay_cmd = 2.0F * localB->rtb_vd_idx_2 * (real32_T)sin
        (localB->rtb_x_R_idx_0) / FMS_PARAM.L1;

      /* Gain: '<S64>/Gain2' */
      localB->rtb_x_R_idx_0 = FMS_PARAM.Z_P * localB->rtb_ve_idx_0;

      /* Saturate: '<S64>/Saturation' */
      if (localB->rtb_x_R_idx_0 > CONTROL_PARAM.FW_T_CLMB_MAX) {
        /* BusAssignment: '<S19>/Bus Assignment' */
        localB->Merge.vh_cmd = CONTROL_PARAM.FW_T_CLMB_MAX;
      } else if (localB->rtb_x_R_idx_0 < -CONTROL_PARAM.FW_T_SINK_MAX) {
        /* BusAssignment: '<S19>/Bus Assignment' */
        localB->Merge.vh_cmd = -CONTROL_PARAM.FW_T_SINK_MAX;
      } else {
        /* BusAssignment: '<S19>/Bus Assignment' */
        localB->Merge.vh_cmd = localB->rtb_x_R_idx_0;
      }

      /* End of Saturate: '<S64>/Saturation' */

      /* Update for Delay: '<S68>/start_vel' */
      localDW->icLoad_l3 = false;

      /* Update for Delay: '<S64>/Delay' */
      localDW->icLoad_p = false;

      /* End of Outputs for SubSystem: '<S12>/FormHold' */
      break;

     case 2:
      /* Outputs for IfAction SubSystem: '<S12>/FormMission' incorporates:
       *  ActionPort: '<S20>/Action Port'
       */
      /* RelationalOperator: '<S81>/FixPt Relational Operator' incorporates:
       *  UnitDelay: '<S81>/Delay Input1'
       *
       * Block description for '<S81>/Delay Input1':
       *
       *  Store in Global RAM
       */
      rtb_FixPtRelationalOperator_c = (localB->wp_index !=
        localDW->DelayInput1_DSTATE_d);

      /* Outputs for Resettable SubSystem: '<S20>/FormMission_SubSystem' incorporates:
       *  ResetPort: '<S82>/Reset'
       */
      if (rtb_FixPtRelationalOperator_c &&
          (localZCE->FormMission_SubSystem_Reset_ZCE != POS_ZCSIG)) {
        /* InitializeConditions for Delay: '<S156>/Delay' */
        localDW->icLoad_k = true;

        /* InitializeConditions for DiscreteIntegrator: '<S137>/Integrator' */
        localDW->Integrator_DSTATE = 0.0F;
      }

      localZCE->FormMission_SubSystem_Reset_ZCE = rtb_FixPtRelationalOperator_c;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[0] = localB->rtb_x_R_idx_0;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[3] = localB->rtb_y_R_idx_0;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[6] = localB->h_R[0];

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[0] = localB->vn[0];

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[3] = localB->rtb_ve_idx_0;

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[6] = localB->rtb_vd_idx_0;
      localB->rtb_ve_idx_0 = 0.0F;
      localB->rtb_vd_idx_0 = 0.0F;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' incorporates:
       *  MATLAB Function: '<S89>/Calc_Position_Velocity_Setpoint'
       */
      localB->xyz_O_nx3[1] = localB->rtb_x_R_idx_1;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[4] = localB->rtb_y_R_idx_1;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[7] = localB->h_R[1];

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[1] = localB->vn[1];

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[4] = localB->u_n;

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[7] = localB->rtb_vd_idx_1;
      localB->u_n = 0.0F;
      localB->rtb_vd_idx_1 = 0.0F;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' incorporates:
       *  MATLAB Function: '<S89>/Calc_Position_Velocity_Setpoint'
       */
      localB->xyz_O_nx3[2] = localB->absxk;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[5] = localB->scale;

      /* SignalConversion generated from: '<S89>/Vector Concatenate' */
      localB->xyz_O_nx3[8] = localB->h_R[2];

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[2] = localB->vn[2];

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[5] = localB->DotProduct;

      /* SignalConversion generated from: '<S89>/Vector Concatenate1' */
      localB->vNED_O_nx3[8] = localB->rtb_vd_idx_2;
      localB->rtb_vd_idx_2 = 0.0F;

      /* MATLAB Function: '<S89>/Calc_Position_Velocity_Setpoint' incorporates:
       *  Concatenate: '<S89>/Vector Concatenate'
       *  Concatenate: '<S89>/Vector Concatenate1'
       */
      localB->rtb_x_R_idx_0 = 0.0F;
      localB->DotProduct = 0.0F;
      for (localB->i = 0; localB->i < 3; localB->i++) {
        if ((uint32_T)(localB->i + 1) != FORMATION_PARAM.UAV_ID) {
          localB->idx = (3 * localB->i + (int32_T)FORMATION_PARAM.UAV_ID) - 1;
          localB->rtb_y_R_idx_0 = FORMATION_PARAM.ADJ_MARTIX[localB->idx];
          if (!(localB->rtb_y_R_idx_0 == 0.0F)) {
            localB->scale = 1.29246971E-26F;
            localB->absxk = (real32_T)fabs(localB->vNED_O_nx3[localB->i]);
            if (localB->absxk > 1.29246971E-26F) {
              localB->rtb_x_R_idx_1 = 1.0F;
              localB->scale = localB->absxk;
            } else {
              localB->t = localB->absxk / 1.29246971E-26F;
              localB->rtb_x_R_idx_1 = localB->t * localB->t;
            }

            localB->rtb_y_R_idx_1 = localB->vNED_O_nx3[localB->i + 3];
            localB->absxk = (real32_T)fabs(localB->rtb_y_R_idx_1);
            if (localB->absxk > localB->scale) {
              localB->t = localB->scale / localB->absxk;
              localB->rtb_x_R_idx_1 = localB->rtb_x_R_idx_1 * localB->t *
                localB->t + 1.0F;
              localB->scale = localB->absxk;
            } else {
              localB->t = localB->absxk / localB->scale;
              localB->rtb_x_R_idx_1 += localB->t * localB->t;
            }

            localB->rtb_x_R_idx_1 = localB->scale * (real32_T)sqrt
              (localB->rtb_x_R_idx_1);
            if ((real32_T)fabs(rtu_Formation_Cross_ay[localB->i]) < 0.1) {
              localB->scale = rt_atan2f_snf(localB->rtb_y_R_idx_1,
                localB->vNED_O_nx3[localB->i]);
              localB->rtb_y_R_idx_1 = (real32_T)sin(localB->scale);
              localB->scale = (real32_T)cos(localB->scale);
              localB->b_tmp[0] = localB->scale;
              localB->b_tmp[3] = -localB->rtb_y_R_idx_1;
              localB->b_tmp[6] = 0.0F;
              localB->b_tmp[1] = localB->rtb_y_R_idx_1;
              localB->b_tmp[4] = localB->scale;
              localB->b_tmp[7] = 0.0F;
              localB->b_tmp[2] = 0.0F;
              localB->b_tmp[5] = 0.0F;
              localB->b_tmp[8] = 1.0F;
              localB->rtb_y_R_idx_1 = FORMATION_PARAM.REL_X_MATRIX[localB->idx];
              localB->scale = FORMATION_PARAM.REL_Y_MATRIX[localB->idx];
              localB->absxk = FORMATION_PARAM.REL_Z_MATRIX[localB->idx];
              for (localB->idx = 0; localB->idx < 3; localB->idx++) {
                localB->rtb_vn_tmp = 3 * localB->idx + localB->i;
                localB->vn[localB->idx] = ((localB->b_tmp[localB->idx + 3] *
                  localB->scale + localB->b_tmp[localB->idx] *
                  localB->rtb_y_R_idx_1) + localB->b_tmp[localB->idx + 6] *
                  localB->absxk) + localB->xyz_O_nx3[localB->rtb_vn_tmp];
                localB->h_R[localB->idx] = localB->vNED_O_nx3[localB->rtb_vn_tmp];
              }

              localB->rtb_y_R_idx_1 = rtu_Formation_Cross_ay[localB->i] /
                localB->rtb_x_R_idx_1;
            } else {
              localB->absxk = localB->rtb_x_R_idx_1 * localB->rtb_x_R_idx_1 /
                rtu_Formation_Cross_ay[localB->i];
              localB->scale = localB->absxk -
                FORMATION_PARAM.REL_Y_MATRIX[localB->idx];
              localB->dpsi = FORMATION_PARAM.REL_X_MATRIX[localB->idx] /
                localB->absxk;
              localB->unit_center_to_pose[0] = localB->rtb_y_R_idx_1 /
                localB->rtb_x_R_idx_1;
              localB->unit_center_to_pose[1] = -localB->vNED_O_nx3[localB->i] /
                localB->rtb_x_R_idx_1;
              localB->t = (real32_T)sin(localB->dpsi);
              localB->dpsi = (real32_T)cos(localB->dpsi);
              localB->rtb_y_R_idx_1 = (localB->dpsi *
                localB->unit_center_to_pose[0] + -localB->t *
                localB->unit_center_to_pose[1]) * localB->scale;
              localB->rtb_TmpSignalConversionAtSqua_b = localB->rtb_y_R_idx_1 -
                localB->unit_center_to_pose[0] * localB->absxk;
              localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_y_R_idx_1;
              localB->rtb_y_R_idx_1 = (localB->t * localB->unit_center_to_pose[0]
                + localB->dpsi * localB->unit_center_to_pose[1]) * localB->scale;
              localB->vn[0] = localB->xyz_O_nx3[localB->i] +
                localB->rtb_TmpSignalConversionAtSqua_b;
              localB->vn[1] = (localB->rtb_y_R_idx_1 -
                               localB->unit_center_to_pose[1] * localB->absxk) +
                localB->xyz_O_nx3[localB->i + 3];
              localB->vn[2] = localB->xyz_O_nx3[localB->i + 6] +
                FORMATION_PARAM.REL_Z_MATRIX[localB->idx];
              localB->absxk = localB->scale / localB->absxk *
                localB->rtb_x_R_idx_1;
              localB->h_R[0] = localB->absxk * -localB->rtb_y_R_idx_1 /
                localB->scale;
              localB->h_R[1] = localB->absxk *
                localB->TmpSignalConversionAtSqua_c[0] / localB->scale;
              localB->rtb_y_R_idx_1 = rtu_Formation_Cross_ay[localB->i] /
                localB->rtb_x_R_idx_1;
            }

            localB->rtb_vd_idx_0 += localB->rtb_y_R_idx_0 * localB->vn[0];
            localB->rtb_ve_idx_0 += localB->rtb_y_R_idx_0 * localB->h_R[0];
            localB->rtb_vd_idx_1 += localB->rtb_y_R_idx_0 * localB->vn[1];
            localB->u_n += localB->rtb_y_R_idx_0 * localB->h_R[1];
            localB->rtb_vd_idx_2 += localB->rtb_y_R_idx_0 * localB->vn[2];
            localB->DotProduct += localB->rtb_y_R_idx_0 * localB->rtb_y_R_idx_1;
            localB->rtb_x_R_idx_0++;
          }
        }
      }

      if (localB->rtb_x_R_idx_0 != 0.0F) {
        localB->rtb_vd_idx_0 = localB->rtb_vd_idx_0 / localB->rtb_x_R_idx_0 -
          localB->xyz_O_nx3[(int32_T)FORMATION_PARAM.UAV_ID - 1];
        localB->rtb_ve_idx_0 /= localB->rtb_x_R_idx_0;
        localB->rtb_vd_idx_1 = localB->rtb_vd_idx_1 / localB->rtb_x_R_idx_0 -
          localB->xyz_O_nx3[(int32_T)FORMATION_PARAM.UAV_ID + 2];
        localB->u_n /= localB->rtb_x_R_idx_0;
        localB->rtb_vd_idx_2 = localB->rtb_vd_idx_2 / localB->rtb_x_R_idx_0 -
          localB->xyz_O_nx3[(int32_T)FORMATION_PARAM.UAV_ID + 5];
        localB->DotProduct /= localB->rtb_x_R_idx_0;
      } else {
        localB->rtb_ve_idx_0 = localB->vNED_O_nx3[(int32_T)
          FORMATION_PARAM.UAV_ID - 1];
        localB->u_n = localB->vNED_O_nx3[(int32_T)FORMATION_PARAM.UAV_ID + 2];
      }

      /* SignalConversion generated from: '<S94>/Square' */
      localB->TmpSignalConversionAtSqua_c[0] = *rtu_INS_Out_vn;
      localB->TmpSignalConversionAtSqua_c[1] = *rtu_INS_Out_ve;

      /* Reshape: '<S162>/Reshape2' */
      localB->Reshape2_bq[0] = *rtu_INS_Out_x_R;
      localB->Reshape2_bq[1] = *rtu_INS_Out_y_R;

      /* MATLAB Function: '<S163>/NearbyRefWP' incorporates:
       *  Constant: '<S161>/L1'
       */
      Formation_FMS_NearbyRefWP(&localB->Cmd_In.sp_waypoint[0],
        localB->Reshape2_bq, FMS_PARAM.L1, localB->unit_center_to_pose,
        &localB->rtb_x_R_idx_1);

      /* MATLAB Function: '<S163>/SearchL1RefWP' incorporates:
       *  Constant: '<S161>/L1'
       */
      Formation_FMS_SearchL1RefWP(&localB->Cmd_In.cur_waypoint[0],
        &localB->Cmd_In.sp_waypoint[0], localB->Reshape2_bq, FMS_PARAM.L1,
        localB->P_i, &localB->dpsi);

      /* MATLAB Function: '<S163>/OutRegionRegWP' */
      Formation_FMS_OutRegionRegWP(&localB->Cmd_In.cur_waypoint[0],
        &localB->Cmd_In.sp_waypoint[0], localB->Reshape2_bq, localB->P_b);

      /* MATLAB Function: '<S157>/min_radius' */
      if (!rtIsNaNF(FORMATION_PARAM.REL_Y_MATRIX[0])) {
        localB->idx = 1;
      } else {
        localB->idx = 0;
        localB->i = 2;
        exitg1 = false;
        while ((!exitg1) && (localB->i < 4)) {
          if (!rtIsNaNF(FORMATION_PARAM.REL_Y_MATRIX[(localB->i - 1) * 3])) {
            localB->idx = localB->i;
            exitg1 = true;
          } else {
            localB->i++;
          }
        }
      }

      if (localB->idx == 0) {
        localB->rtb_y_R_idx_1 = FORMATION_PARAM.REL_Y_MATRIX[0];
      } else {
        localB->rtb_y_R_idx_1 = FORMATION_PARAM.REL_Y_MATRIX[(localB->idx - 1) *
          3];
        for (localB->i = localB->idx + 1; localB->i < 4; localB->i++) {
          localB->rtb_y_R_idx_0 = FORMATION_PARAM.REL_Y_MATRIX[(localB->i - 1) *
            3];
          if (localB->rtb_y_R_idx_1 < localB->rtb_y_R_idx_0) {
            localB->rtb_y_R_idx_1 = localB->rtb_y_R_idx_0;
          }
        }
      }

      /* Delay: '<S156>/Delay' incorporates:
       *  Constant: '<S156>/Constant'
       */
      if (localDW->icLoad_k) {
        localDW->Delay_DSTATE_o = FMS_PARAM.FW_HEIGHT_TRIM;
      }

      /* MATLAB Function: '<S92>/Pout_Saturation' */
      localB->scale = 1.29246971E-26F;

      /* Gain: '<S92>/POS_KP' */
      localB->rtb_x_R_idx_0 = FORMATION_PARAM.FORM_POS_KP * localB->rtb_vd_idx_0;

      /* MATLAB Function: '<S92>/Pout_Saturation' */
      localB->absxk = (real32_T)fabs(localB->rtb_x_R_idx_0);
      if (localB->absxk > 1.29246971E-26F) {
        localB->rtb_y_R_idx_0 = 1.0F;
        localB->scale = localB->absxk;
      } else {
        localB->t = localB->absxk / 1.29246971E-26F;
        localB->rtb_y_R_idx_0 = localB->t * localB->t;
      }

      localB->P_d[0] = localB->rtb_x_R_idx_0;

      /* Gain: '<S92>/POS_KP' */
      localB->rtb_x_R_idx_0 = FORMATION_PARAM.FORM_POS_KP * localB->rtb_vd_idx_1;

      /* MATLAB Function: '<S92>/Pout_Saturation' */
      localB->absxk = (real32_T)fabs(localB->rtb_x_R_idx_0);
      if (localB->absxk > localB->scale) {
        localB->t = localB->scale / localB->absxk;
        localB->rtb_y_R_idx_0 = localB->rtb_y_R_idx_0 * localB->t * localB->t +
          1.0F;
        localB->scale = localB->absxk;
      } else {
        localB->t = localB->absxk / localB->scale;
        localB->rtb_y_R_idx_0 += localB->t * localB->t;
      }

      localB->P_d[1] = localB->rtb_x_R_idx_0;

      /* MATLAB Function: '<S92>/Pout_Saturation' incorporates:
       *  Gain: '<S92>/POS_KP'
       */
      localB->rtb_y_R_idx_0 = localB->scale * (real32_T)sqrt
        (localB->rtb_y_R_idx_0);
      if (localB->rtb_y_R_idx_0 > FORMATION_PARAM.FORM_POUT_LIM) {
        localB->P_d[0] = localB->P_d[0] / localB->rtb_y_R_idx_0 *
          FORMATION_PARAM.FORM_POUT_LIM;
        localB->P_d[1] = localB->rtb_x_R_idx_0 / localB->rtb_y_R_idx_0 *
          FORMATION_PARAM.FORM_POUT_LIM;
      }

      /* Switch: '<S82>/Switch' incorporates:
       *  Constant: '<S82>/Constant2'
       *  Constant: '<S83>/Constant'
       *  Constant: '<S95>/Constant'
       *  DiscreteIntegrator: '<S137>/Integrator'
       *  Gain: '<S142>/Proportional Gain'
       *  Gain: '<S155>/Gain'
       *  MATLAB Function: '<S89>/Calc_Position_Velocity_Setpoint'
       *  Math: '<S96>/Square'
       *  Product: '<S176>/Divide'
       *  Product: '<S91>/Product'
       *  Product: '<S91>/Product1'
       *  RelationalOperator: '<S83>/Compare'
       *  RelationalOperator: '<S95>/Compare'
       *  Saturate: '<S91>/Saturation'
       *  Sqrt: '<S96>/Sqrt'
       *  Sum: '<S146>/Sum'
       *  Sum: '<S91>/Sum1'
       *  Sum: '<S92>/Add1'
       *  Sum: '<S96>/Sum of Elements'
       */
      if (FORMATION_PARAM.UAV_ID == 1U) {
        /* Switch: '<S163>/Switch1' incorporates:
         *  Constant: '<S166>/Constant'
         *  RelationalOperator: '<S166>/Compare'
         */
        if (!(localB->rtb_x_R_idx_1 > 0.0F)) {
          /* Switch: '<S163>/Switch' incorporates:
           *  Constant: '<S165>/Constant'
           *  RelationalOperator: '<S165>/Compare'
           *  Switch: '<S163>/Switch1'
           */
          if (localB->dpsi >= 0.0F) {
            localB->unit_center_to_pose[0] = localB->P_i[0];
            localB->unit_center_to_pose[1] = localB->P_i[1];
          } else {
            localB->unit_center_to_pose[0] = localB->P_b[0];
            localB->unit_center_to_pose[1] = localB->P_b[1];
          }

          /* End of Switch: '<S163>/Switch' */
        }

        /* End of Switch: '<S163>/Switch1' */

        /* Sum: '<S164>/Subtract' incorporates:
         *  Reshape: '<S162>/Reshape2'
         *  Switch: '<S163>/Switch1'
         */
        localB->rtb_vd_idx_1 = localB->unit_center_to_pose[0] -
          localB->Reshape2_bq[0];
        localB->rtb_vd_idx_0 = localB->unit_center_to_pose[1] -
          localB->Reshape2_bq[1];

        /* Sum: '<S173>/Sum of Elements' incorporates:
         *  Math: '<S173>/Math Function'
         *  Sum: '<S172>/Sum of Elements'
         */
        localB->rtb_y_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] *
          localB->TmpSignalConversionAtSqua_c[0] +
          localB->TmpSignalConversionAtSqua_c[1] *
          localB->TmpSignalConversionAtSqua_c[1];

        /* Math: '<S173>/Math Function1' incorporates:
         *  Sum: '<S173>/Sum of Elements'
         *
         * About '<S173>/Math Function1':
         *  Operator: sqrt
         */
        if (localB->rtb_y_R_idx_0 < 0.0F) {
          localB->DotProduct = -(real32_T)sqrt((real32_T)fabs
            (localB->rtb_y_R_idx_0));
        } else {
          localB->DotProduct = (real32_T)sqrt(localB->rtb_y_R_idx_0);
        }

        /* End of Math: '<S173>/Math Function1' */

        /* Switch: '<S173>/Switch' incorporates:
         *  Constant: '<S173>/Constant'
         *  Product: '<S173>/Product'
         */
        if (localB->DotProduct > 0.0F) {
          localB->rtb_ve_idx_0 = *rtu_INS_Out_vn;
          localB->u_n = *rtu_INS_Out_ve;
        } else {
          localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
          localB->u_n = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
          localB->DotProduct = 1.0F;
        }

        /* End of Switch: '<S173>/Switch' */

        /* Product: '<S173>/Divide' */
        localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_ve_idx_0 /
          localB->DotProduct;
        localB->TmpSignalConversionAtSqua_c[1] = localB->u_n /
          localB->DotProduct;

        /* Sum: '<S174>/Sum of Elements' incorporates:
         *  Math: '<S174>/Math Function'
         */
        localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_1 * localB->rtb_vd_idx_1 +
          localB->rtb_vd_idx_0 * localB->rtb_vd_idx_0;

        /* Math: '<S174>/Math Function1' incorporates:
         *  Sum: '<S174>/Sum of Elements'
         *
         * About '<S174>/Math Function1':
         *  Operator: sqrt
         */
        if (localB->rtb_x_R_idx_0 < 0.0F) {
          localB->DotProduct = -(real32_T)sqrt((real32_T)fabs
            (localB->rtb_x_R_idx_0));
        } else {
          localB->DotProduct = (real32_T)sqrt(localB->rtb_x_R_idx_0);
        }

        /* End of Math: '<S174>/Math Function1' */

        /* Switch: '<S174>/Switch' incorporates:
         *  Constant: '<S174>/Constant'
         *  Product: '<S174>/Product'
         *  Switch: '<S176>/Switch'
         */
        if (localB->DotProduct > 0.0F) {
          localB->rtb_ve_idx_0 = localB->rtb_vd_idx_1;
          localB->u_n = localB->rtb_vd_idx_0;
        } else {
          localB->rtb_ve_idx_0 = localB->rtb_vd_idx_1 * 0.0F;
          localB->u_n = localB->rtb_vd_idx_0 * 0.0F;
          localB->DotProduct = 1.0F;
        }

        /* End of Switch: '<S174>/Switch' */

        /* Sum: '<S155>/Sum' incorporates:
         *  Constant: '<S155>/Constant'
         */
        localB->rtb_vd_idx_1 = FMS_PARAM.FW_AIRSPD_TRIM - *rtu_INS_Out_airspeed;

        /* MATLAB Function: '<S157>/min_radius' */
        localB->rtb_x_R_idx_0 = CONTROL_PARAM.FW_AIRSPD_MAX -
          CONTROL_PARAM.FW_AIRSPD_TRIM;
        localB->rtb_vd_idx_0 = CONTROL_PARAM.FW_AIRSPD_TRIM -
          CONTROL_PARAM.FW_AIRSPD_MIN;
        if ((localB->rtb_x_R_idx_0 <= localB->rtb_vd_idx_0) || rtIsNaNF
            (localB->rtb_vd_idx_0)) {
          localB->rtb_vd_idx_0 = localB->rtb_x_R_idx_0;
        }

        /* MinMax: '<S157>/Min' incorporates:
         *  Gain: '<S157>/Satefy'
         *  MATLAB Function: '<S157>/min_radius'
         */
        localB->n_j = localB->rtb_y_R_idx_1 * CONTROL_PARAM.FW_AIRSPD_TRIM /
          localB->rtb_vd_idx_0 * FMS_PARAM.FW_RADIUS_RATIO;
        if (!(localB->n_j >= 0.1)) {
          localB->n_j = 0.1;
        }

        /* Product: '<S157>/Divide' incorporates:
         *  Constant: '<S157>/v^2'
         *  MinMax: '<S157>/Min'
         */
        localB->n_j = CONTROL_PARAM.FW_AIRSPD_TRIM *
          CONTROL_PARAM.FW_AIRSPD_TRIM / localB->n_j;

        /* Product: '<S174>/Divide' incorporates:
         *  Product: '<S177>/Divide'
         */
        localB->unit_center_to_pose[0] = localB->rtb_ve_idx_0 /
          localB->DotProduct;
        localB->unit_center_to_pose[1] = localB->u_n / localB->DotProduct;

        /* Sqrt: '<S172>/Sqrt' */
        localB->DotProduct = (real32_T)sqrt(localB->rtb_y_R_idx_0);

        /* Math: '<S170>/Square' */
        localB->rtb_vd_idx_0 = localB->DotProduct * localB->DotProduct;

        /* Sum: '<S177>/Sum of Elements' incorporates:
         *  Math: '<S177>/Math Function'
         *  SignalConversion generated from: '<S177>/Math Function'
         */
        localB->rtb_y_R_idx_0 = localB->unit_center_to_pose[1] *
          localB->unit_center_to_pose[1] + localB->unit_center_to_pose[0] *
          localB->unit_center_to_pose[0];

        /* Math: '<S177>/Math Function1' incorporates:
         *  Sum: '<S177>/Sum of Elements'
         *
         * About '<S177>/Math Function1':
         *  Operator: sqrt
         */
        if (localB->rtb_y_R_idx_0 < 0.0F) {
          localB->DotProduct = -(real32_T)sqrt((real32_T)fabs
            (localB->rtb_y_R_idx_0));
        } else {
          localB->DotProduct = (real32_T)sqrt(localB->rtb_y_R_idx_0);
        }

        /* End of Math: '<S177>/Math Function1' */

        /* Switch: '<S177>/Switch' incorporates:
         *  Constant: '<S177>/Constant'
         *  Product: '<S177>/Product'
         *  SignalConversion generated from: '<S177>/Math Function'
         */
        if (localB->DotProduct > 0.0F) {
          localB->rtb_ve_idx_0 = localB->unit_center_to_pose[1];
          localB->u_n = localB->unit_center_to_pose[0];
        } else {
          localB->rtb_ve_idx_0 = localB->unit_center_to_pose[1] * 0.0F;
          localB->u_n = localB->unit_center_to_pose[0] * 0.0F;
          localB->DotProduct = 1.0F;
        }

        /* End of Switch: '<S177>/Switch' */

        /* Product: '<S177>/Divide' */
        localB->unit_center_to_pose[0] = localB->rtb_ve_idx_0 /
          localB->DotProduct;
        localB->unit_center_to_pose[1] = localB->u_n / localB->DotProduct;

        /* Sum: '<S176>/Sum of Elements' incorporates:
         *  Math: '<S176>/Math Function'
         *  SignalConversion generated from: '<S176>/Math Function'
         */
        localB->rtb_y_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
          localB->TmpSignalConversionAtSqua_c[1] +
          localB->TmpSignalConversionAtSqua_c[0] *
          localB->TmpSignalConversionAtSqua_c[0];

        /* Math: '<S176>/Math Function1' incorporates:
         *  Sum: '<S176>/Sum of Elements'
         *
         * About '<S176>/Math Function1':
         *  Operator: sqrt
         */
        if (localB->rtb_y_R_idx_0 < 0.0F) {
          localB->DotProduct = -(real32_T)sqrt((real32_T)fabs
            (localB->rtb_y_R_idx_0));
        } else {
          localB->DotProduct = (real32_T)sqrt(localB->rtb_y_R_idx_0);
        }

        /* End of Math: '<S176>/Math Function1' */

        /* Switch: '<S176>/Switch' incorporates:
         *  Constant: '<S176>/Constant'
         *  Product: '<S176>/Product'
         *  SignalConversion generated from: '<S176>/Math Function'
         */
        if (localB->DotProduct > 0.0F) {
          localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
          localB->u_n = localB->TmpSignalConversionAtSqua_c[0];
        } else {
          localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
          localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
          localB->DotProduct = 1.0F;
        }

        /* End of Switch: '<S176>/Switch' */

        /* Product: '<S176>/Divide' */
        localB->rtb_y_R_idx_1 = localB->rtb_ve_idx_0 / localB->DotProduct;

        /* DotProduct: '<S171>/Dot Product' */
        localB->rtb_ve_idx_0 = localB->rtb_y_R_idx_1 *
          localB->unit_center_to_pose[0];
        localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_y_R_idx_1;

        /* Product: '<S176>/Divide' */
        localB->rtb_y_R_idx_1 = localB->u_n / localB->DotProduct;

        /* DotProduct: '<S171>/Dot Product' */
        localB->rtb_ve_idx_0 += localB->rtb_y_R_idx_1 *
          localB->unit_center_to_pose[1];

        /* Sum: '<S175>/Subtract' incorporates:
         *  Product: '<S175>/Multiply'
         *  Product: '<S175>/Multiply1'
         */
        localB->u_n = localB->unit_center_to_pose[0] * localB->rtb_y_R_idx_1 -
          localB->TmpSignalConversionAtSqua_c[0] * localB->unit_center_to_pose[1];

        /* Signum: '<S171>/Sign1' */
        if (rtIsNaNF(localB->u_n)) {
          localB->u_n = (rtNaNF);
        } else if (localB->u_n < 0.0F) {
          localB->u_n = -1.0F;
        } else {
          localB->u_n = (real32_T)(localB->u_n > 0.0F);
        }

        /* End of Signum: '<S171>/Sign1' */

        /* Trigonometry: '<S171>/Acos' incorporates:
         *  DotProduct: '<S171>/Dot Product'
         */
        if (localB->rtb_ve_idx_0 > 1.0F) {
          localB->rtb_ve_idx_0 = 1.0F;
        } else if (localB->rtb_ve_idx_0 < -1.0F) {
          localB->rtb_ve_idx_0 = -1.0F;
        }

        /* Switch: '<S171>/Switch2' incorporates:
         *  Constant: '<S171>/Constant4'
         */
        if (!(localB->u_n != 0.0F)) {
          localB->u_n = 1.0F;
        }

        /* Product: '<S171>/Multiply' incorporates:
         *  Switch: '<S171>/Switch2'
         *  Trigonometry: '<S171>/Acos'
         */
        localB->rtb_x_R_idx_0 = (real32_T)acos(localB->rtb_ve_idx_0) *
          localB->u_n;

        /* Saturate: '<S170>/Saturation' */
        if (localB->rtb_x_R_idx_0 > 1.57079637F) {
          localB->rtb_x_R_idx_0 = 1.57079637F;
        } else if (localB->rtb_x_R_idx_0 < -1.57079637F) {
          localB->rtb_x_R_idx_0 = -1.57079637F;
        }

        /* Product: '<S170>/Divide' incorporates:
         *  Constant: '<S161>/L1'
         *  Gain: '<S170>/Gain'
         *  Product: '<S170>/Multiply1'
         *  Saturate: '<S170>/Saturation'
         *  Trigonometry: '<S170>/Sin'
         */
        localB->u_n = 2.0F * localB->rtb_vd_idx_0 * (real32_T)sin
          (localB->rtb_x_R_idx_0) / FMS_PARAM.L1;

        /* Saturate: '<S161>/Saturation1' */
        if (localB->u_n > FMS_PARAM.ACC_Y_LIM) {
          localB->u_n = FMS_PARAM.ACC_Y_LIM;
        } else if (localB->u_n < -FMS_PARAM.ACC_Y_LIM) {
          localB->u_n = -FMS_PARAM.ACC_Y_LIM;
        }

        /* End of Saturate: '<S161>/Saturation1' */

        /* Switch: '<S159>/Switch2' incorporates:
         *  Gain: '<S157>/Gain'
         *  RelationalOperator: '<S159>/LowerRelop1'
         *  RelationalOperator: '<S159>/UpperRelop'
         *  Switch: '<S159>/Switch'
         */
        if (localB->u_n > localB->n_j) {
          localB->u_n = (real32_T)localB->n_j;
        } else if (localB->u_n < -localB->n_j) {
          /* Switch: '<S159>/Switch' incorporates:
           *  Gain: '<S157>/Gain'
           */
          localB->u_n = (real32_T)-localB->n_j;
        }

        /* End of Switch: '<S159>/Switch2' */

        /* Sum: '<S156>/Sum' incorporates:
         *  Delay: '<S156>/Delay'
         */
        localB->rtb_vd_idx_0 = localDW->Delay_DSTATE_o - *rtu_INS_Out_h_R;
        localB->rtb_ve_idx_0 = FMS_PARAM.AIRSPD_P * localB->rtb_vd_idx_1;

        /* Gain: '<S156>/Gain2' incorporates:
         *  Gain: '<S155>/Gain'
         */
        localB->DotProduct = FMS_PARAM.Z_P * localB->rtb_vd_idx_0;

        /* Saturate: '<S156>/Saturation' */
        if (localB->DotProduct > CONTROL_PARAM.FW_T_CLMB_MAX) {
          localB->DotProduct = CONTROL_PARAM.FW_T_CLMB_MAX;
        } else if (localB->DotProduct < -CONTROL_PARAM.FW_T_SINK_MAX) {
          localB->DotProduct = -CONTROL_PARAM.FW_T_SINK_MAX;
        }

        /* End of Saturate: '<S156>/Saturation' */
      } else {
        /* Sum: '<S92>/Add1' */
        localB->rtb_x_R_idx_0 = localB->P_d[0] + localB->rtb_ve_idx_0;

        /* Math: '<S93>/Square' incorporates:
         *  Math: '<S97>/Square'
         */
        localB->unit_center_to_pose[0] = localB->rtb_x_R_idx_0 *
          localB->rtb_x_R_idx_0;
        localB->P_d[0] = localB->rtb_x_R_idx_0;

        /* Sum: '<S92>/Add1' incorporates:
         *  Math: '<S96>/Square'
         */
        localB->rtb_x_R_idx_0 = localB->P_d[1] + localB->u_n;

        /* Trigonometry: '<S98>/Atan2' */
        localB->rtb_ve_idx_0 = rt_atan2f_snf(*rtu_INS_Out_ve, *rtu_INS_Out_vn);

        /* Sum: '<S91>/Sum' incorporates:
         *  Trigonometry: '<S99>/Atan2'
         */
        localB->u_n = rt_atan2f_snf(localB->rtb_x_R_idx_0, localB->P_d[0]) -
          localB->rtb_ve_idx_0;

        /* Abs: '<S100>/Abs' */
        localB->rtb_ve_idx_0 = (real32_T)fabs(localB->u_n);

        /* Switch: '<S100>/Switch' incorporates:
         *  Constant: '<S100>/Constant'
         *  Constant: '<S101>/Constant'
         *  Product: '<S100>/Multiply'
         *  RelationalOperator: '<S101>/Compare'
         *  Signum: '<S100>/Sign'
         *  Sum: '<S100>/Subtract'
         */
        if (localB->rtb_ve_idx_0 > 3.14159274F) {
          /* Signum: '<S100>/Sign' */
          if (rtIsNaNF(localB->u_n)) {
            localB->rtb_y_R_idx_0 = (rtNaNF);
          } else if (localB->u_n < 0.0F) {
            localB->rtb_y_R_idx_0 = -1.0F;
          } else {
            localB->rtb_y_R_idx_0 = (real32_T)(localB->u_n > 0.0F);
          }

          localB->u_n = (localB->rtb_ve_idx_0 - 6.28318548F) *
            localB->rtb_y_R_idx_0;
        }

        /* End of Switch: '<S100>/Switch' */

        /* Sqrt: '<S93>/Sqrt' incorporates:
         *  Math: '<S93>/Square'
         *  Sqrt: '<S97>/Sqrt'
         *  Sum: '<S93>/Sum of Elements'
         */
        localB->rtb_y_R_idx_0 = (real32_T)sqrt(localB->rtb_x_R_idx_0 *
          localB->rtb_x_R_idx_0 + localB->unit_center_to_pose[0]);

        /* Gain: '<S90>/Vel_Kp' incorporates:
         *  Math: '<S94>/Square'
         *  Sqrt: '<S93>/Sqrt'
         *  Sqrt: '<S94>/Sqrt'
         *  Sum: '<S90>/Sum1'
         *  Sum: '<S94>/Sum of Elements'
         */
        localB->rtb_ve_idx_0 = (localB->rtb_y_R_idx_0 - (real32_T)sqrt
          (localB->TmpSignalConversionAtSqua_c[0] *
           localB->TmpSignalConversionAtSqua_c[0] +
           localB->TmpSignalConversionAtSqua_c[1] *
           localB->TmpSignalConversionAtSqua_c[1])) *
          FORMATION_PARAM.FORM_VEL_KP;

        /* Saturate: '<S90>/Saturation' */
        if (localB->rtb_ve_idx_0 > FORMATION_PARAM.FORM_VOUT_LIM) {
          localB->rtb_ve_idx_0 = FORMATION_PARAM.FORM_VOUT_LIM;
        } else if (localB->rtb_ve_idx_0 < -FORMATION_PARAM.FORM_VOUT_LIM) {
          localB->rtb_ve_idx_0 = -FORMATION_PARAM.FORM_VOUT_LIM;
        }

        /* End of Saturate: '<S90>/Saturation' */

        /* Gain: '<S91>/Gain' */
        localB->rtb_x_R_idx_0 = FORMATION_PARAM.FORM_HEAD_KP * localB->u_n;

        /* Saturate: '<S91>/Saturation' */
        if (localB->rtb_x_R_idx_0 > FORMATION_PARAM.FORM_HOUT_LIM) {
          localB->rtb_x_R_idx_0 = FORMATION_PARAM.FORM_HOUT_LIM;
        } else if (localB->rtb_x_R_idx_0 < -FORMATION_PARAM.FORM_HOUT_LIM) {
          localB->rtb_x_R_idx_0 = -FORMATION_PARAM.FORM_HOUT_LIM;
        }

        localB->u_n = ((real32_T)((real32_T)sqrt(localB->rtb_vd_idx_0 *
          localB->rtb_vd_idx_0 + localB->rtb_vd_idx_1 * localB->rtb_vd_idx_1) <=
          FORMATION_PARAM.FORM_HRFF_THRE) * localB->DotProduct +
                       localB->rtb_x_R_idx_0) * localB->rtb_y_R_idx_0;
        localB->DotProduct = FORMATION_PARAM.FORM_HGT_KP * localB->rtb_vd_idx_2
          + localDW->Integrator_DSTATE;
      }

      /* End of Switch: '<S82>/Switch' */
      /* End of Outputs for SubSystem: '<S20>/FormMission_SubSystem' */
      /* End of Outputs for SubSystem: '<S12>/FormMission' */
      memset(&localB->Merge, 0, sizeof(FMS_Out_Bus));

      /* Outputs for IfAction SubSystem: '<S12>/FormMission' incorporates:
       *  ActionPort: '<S20>/Action Port'
       */
      /* Outputs for Resettable SubSystem: '<S20>/FormMission_SubSystem' incorporates:
       *  ResetPort: '<S82>/Reset'
       */
      /* BusAssignment: '<S82>/Bus Assignment' incorporates:
       *  Constant: '<S82>/Constant1'
       *  Merge: '<S10>/Merge'
       */
      localB->Merge.state = VehicleState_FormMission;
      localB->Merge.ax_cmd = localB->rtb_ve_idx_0;
      localB->Merge.ay_cmd = localB->u_n;
      localB->Merge.vh_cmd = localB->DotProduct;

      /* Update for Delay: '<S156>/Delay' */
      localDW->icLoad_k = false;

      /* Update for DiscreteIntegrator: '<S137>/Integrator' incorporates:
       *  Gain: '<S134>/Integral Gain'
       */
      localDW->Integrator_DSTATE += 0.0F * localB->rtb_vd_idx_2 * 0.04F;

      /* End of Outputs for SubSystem: '<S20>/FormMission_SubSystem' */

      /* Update for UnitDelay: '<S81>/Delay Input1'
       *
       * Block description for '<S81>/Delay Input1':
       *
       *  Store in Global RAM
       */
      localDW->DelayInput1_DSTATE_d = localB->wp_index;

      /* End of Outputs for SubSystem: '<S12>/FormMission' */
      break;

     case 3:
      if (localDW->SwitchCase_ActiveSubsystem_c != rtPrevAction) {
        /* InitializeConditions for IfAction SubSystem: '<S12>/FormDisband' incorporates:
         *  ActionPort: '<S18>/Action Port'
         */
        /* InitializeConditions for SwitchCase: '<S12>/Switch Case' incorporates:
         *  Delay: '<S46>/Delay'
         *  Delay: '<S50>/start_vel'
         */
        localDW->icLoad_a = true;
        localDW->icLoad_l = true;

        /* End of InitializeConditions for SubSystem: '<S12>/FormDisband' */
      }

      /* Outputs for IfAction SubSystem: '<S12>/FormDisband' incorporates:
       *  ActionPort: '<S18>/Action Port'
       */
      /* SignalConversion generated from: '<S56>/Square' */
      localB->TmpSignalConversionAtSqua_c[0] = *rtu_INS_Out_vn;
      localB->TmpSignalConversionAtSqua_c[1] = *rtu_INS_Out_ve;

      /* Sum: '<S57>/Sum of Elements' incorporates:
       *  Math: '<S57>/Math Function'
       *  Sum: '<S56>/Sum of Elements'
       */
      localB->rtb_y_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0] +
        localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1];

      /* Math: '<S57>/Math Function1' incorporates:
       *  Sum: '<S57>/Sum of Elements'
       *
       * About '<S57>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_y_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_y_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_y_R_idx_0);
      }

      /* End of Math: '<S57>/Math Function1' */

      /* Switch: '<S57>/Switch' incorporates:
       *  Constant: '<S57>/Constant'
       *  Product: '<S57>/Product'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_vd_idx_0 = *rtu_INS_Out_vn;
        localB->rtb_vd_idx_1 = *rtu_INS_Out_ve;
        localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->rtb_vd_idx_2 = 1.0F;
      }

      /* End of Switch: '<S57>/Switch' */

      /* Delay: '<S50>/start_vel' */
      if (localDW->icLoad_a) {
        localDW->start_vel_DSTATE_m[0] = localB->TmpSignalConversionAtSqua_c[0];
        localDW->start_vel_DSTATE_m[1] = localB->TmpSignalConversionAtSqua_c[1];
      }

      /* Math: '<S62>/Math Function' incorporates:
       *  Delay: '<S50>/start_vel'
       */
      localB->TmpSignalConversionAtSqua_c[0] = localDW->start_vel_DSTATE_m[0] *
        localDW->start_vel_DSTATE_m[0];
      localB->TmpSignalConversionAtSqua_c[1] = localDW->start_vel_DSTATE_m[1] *
        localDW->start_vel_DSTATE_m[1];

      /* Sum: '<S62>/Sum of Elements' */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] +
        localB->TmpSignalConversionAtSqua_c[1];

      /* Math: '<S62>/Math Function1' incorporates:
       *  Sum: '<S62>/Sum of Elements'
       *
       * About '<S62>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S62>/Math Function1' */

      /* Switch: '<S62>/Switch' incorporates:
       *  Constant: '<S62>/Constant'
       *  Delay: '<S50>/start_vel'
       *  Product: '<S62>/Product'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_ve_idx_0 = localDW->start_vel_DSTATE_m[0];
        localB->u_n = localDW->start_vel_DSTATE_m[1];
        localB->DotProduct = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_ve_idx_0 = localDW->start_vel_DSTATE_m[0] * 0.0F;
        localB->u_n = localDW->start_vel_DSTATE_m[1] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S62>/Switch' */

      /* Product: '<S57>/Divide' */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_vd_idx_0 /
        localB->rtb_vd_idx_2;
      localB->TmpSignalConversionAtSqua_c[1] = localB->rtb_vd_idx_1 /
        localB->rtb_vd_idx_2;

      /* SignalConversion generated from: '<S60>/Math Function' */
      localB->P_d[0] = localB->TmpSignalConversionAtSqua_c[1];
      localB->P_d[1] = localB->TmpSignalConversionAtSqua_c[0];

      /* Sum: '<S60>/Sum of Elements' incorporates:
       *  Math: '<S60>/Math Function'
       *  SignalConversion generated from: '<S60>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1] +
        localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S60>/Math Function1' incorporates:
       *  Sum: '<S60>/Sum of Elements'
       *
       * About '<S60>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S60>/Math Function1' */

      /* Switch: '<S60>/Switch' incorporates:
       *  Constant: '<S60>/Constant'
       *  Product: '<S60>/Product'
       *  SignalConversion generated from: '<S60>/Math Function'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0];
        localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->rtb_vd_idx_2 = 1.0F;
      }

      /* End of Switch: '<S60>/Switch' */

      /* Gain: '<S48>/Gain' incorporates:
       *  Constant: '<S47>/L1'
       */
      localB->rtb_x_R_idx_1 = 0.5F * FMS_PARAM.L1;

      /* MinMax: '<S48>/Max' incorporates:
       *  Constant: '<S47>/R'
       */
      if ((FMS_PARAM.LOITER_R >= localB->rtb_x_R_idx_1) || rtIsNaNF
          (localB->rtb_x_R_idx_1)) {
        localB->rtb_x_R_idx_1 = FMS_PARAM.LOITER_R;
      }

      /* End of MinMax: '<S48>/Max' */

      /* Reshape: '<S48>/Reshape2' */
      localB->Reshape2_bq[0] = *rtu_INS_Out_x_R;
      localB->Reshape2_bq[1] = *rtu_INS_Out_y_R;

      /* MATLAB Function: '<S48>/SearchL1RefWP' incorporates:
       *  Constant: '<S47>/L1'
       */
      Formation_FMS_SearchL1RefWP_i(&localB->Cmd_In.cur_waypoint[0],
        localB->Reshape2_bq, localB->rtb_x_R_idx_1, FMS_PARAM.L1,
        localB->TmpSignalConversionAtSqua_c, &localB->n_j);

      /* RelationalOperator: '<S51>/Compare' incorporates:
       *  Constant: '<S51>/Constant'
       */
      rtb_FixPtRelationalOperator_c = (localB->n_j > 0.0);

      /* Product: '<S62>/Divide' */
      localB->P_i[0] = localB->rtb_ve_idx_0 / localB->DotProduct;
      localB->P_i[1] = localB->u_n / localB->DotProduct;

      /* MATLAB Function: '<S48>/OutRegionRegWP' incorporates:
       *  Constant: '<S47>/L1'
       */
      Formation_FMS_OutRegionRegWP_o(&localB->Cmd_In.cur_waypoint[0],
        localB->Reshape2_bq, localB->rtb_x_R_idx_1, FMS_PARAM.L1, localB->P_i,
        localB->P_d);

      /* Switch: '<S48>/Switch' */
      if (rtb_FixPtRelationalOperator_c) {
        localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0];
      } else {
        localB->rtb_x_R_idx_0 = localB->P_d[0];
      }

      /* Sum: '<S49>/Subtract' incorporates:
       *  Switch: '<S48>/Switch'
       */
      localB->Reshape2_bq[0] = localB->rtb_x_R_idx_0 - *rtu_INS_Out_x_R;

      /* Switch: '<S48>/Switch' */
      if (rtb_FixPtRelationalOperator_c) {
        localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
      } else {
        localB->rtb_x_R_idx_0 = localB->P_d[1];
      }

      /* Sum: '<S49>/Subtract' incorporates:
       *  Switch: '<S48>/Switch'
       */
      localB->Reshape2_bq[1] = localB->rtb_x_R_idx_0 - *rtu_INS_Out_y_R;

      /* Sum: '<S58>/Sum of Elements' incorporates:
       *  Math: '<S58>/Math Function'
       *  Sum: '<S49>/Subtract'
       */
      localB->rtb_x_R_idx_0 = localB->Reshape2_bq[0] * localB->Reshape2_bq[0] +
        localB->Reshape2_bq[1] * localB->Reshape2_bq[1];

      /* Math: '<S58>/Math Function1' incorporates:
       *  Sum: '<S58>/Sum of Elements'
       *
       * About '<S58>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S58>/Math Function1' */

      /* Switch: '<S58>/Switch' incorporates:
       *  Constant: '<S58>/Constant'
       *  Product: '<S58>/Product'
       *  Sum: '<S49>/Subtract'
       *  Switch: '<S61>/Switch'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_ve_idx_0 = localB->Reshape2_bq[0];
        localB->u_n = localB->Reshape2_bq[1];
        localB->DotProduct = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_ve_idx_0 = localB->Reshape2_bq[0] * 0.0F;
        localB->u_n = localB->Reshape2_bq[1] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S58>/Switch' */

      /* Product: '<S58>/Divide' incorporates:
       *  Product: '<S61>/Divide'
       */
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_ve_idx_0 /
        localB->DotProduct;
      localB->TmpSignalConversionAtSqua_c[1] = localB->u_n / localB->DotProduct;

      /* Sum: '<S61>/Sum of Elements' incorporates:
       *  Math: '<S61>/Math Function'
       *  SignalConversion generated from: '<S61>/Math Function'
       */
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
        localB->TmpSignalConversionAtSqua_c[1] +
        localB->TmpSignalConversionAtSqua_c[0] *
        localB->TmpSignalConversionAtSqua_c[0];

      /* Math: '<S61>/Math Function1' incorporates:
       *  Sum: '<S61>/Sum of Elements'
       *
       * About '<S61>/Math Function1':
       *  Operator: sqrt
       */
      if (localB->rtb_x_R_idx_0 < 0.0F) {
        localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
          (localB->rtb_x_R_idx_0));
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
      }

      /* End of Math: '<S61>/Math Function1' */

      /* Switch: '<S61>/Switch' incorporates:
       *  Constant: '<S61>/Constant'
       *  Product: '<S61>/Product'
       *  SignalConversion generated from: '<S61>/Math Function'
       */
      if (localB->rtb_x_R_idx_1 > 0.0F) {
        localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
        localB->u_n = localB->TmpSignalConversionAtSqua_c[0];
        localB->DotProduct = localB->rtb_x_R_idx_1;
      } else {
        localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
        localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
        localB->DotProduct = 1.0F;
      }

      /* End of Switch: '<S61>/Switch' */

      /* Product: '<S61>/Divide' */
      localB->rtb_y_R_idx_1 = localB->rtb_ve_idx_0 / localB->DotProduct;

      /* Product: '<S60>/Divide' */
      localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_0 / localB->rtb_vd_idx_2;

      /* DotProduct: '<S55>/Dot Product' */
      localB->rtb_ve_idx_0 = localB->rtb_x_R_idx_0 * localB->rtb_y_R_idx_1;
      localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_y_R_idx_1;
      localB->P_d[0] = localB->rtb_x_R_idx_0;

      /* Product: '<S61>/Divide' incorporates:
       *  Product: '<S60>/Divide'
       */
      localB->rtb_y_R_idx_1 = localB->u_n / localB->DotProduct;

      /* Product: '<S60>/Divide' */
      localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_1 / localB->rtb_vd_idx_2;

      /* DotProduct: '<S55>/Dot Product' */
      localB->rtb_ve_idx_0 += localB->rtb_x_R_idx_0 * localB->rtb_y_R_idx_1;

      /* Sqrt: '<S56>/Sqrt' */
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_y_R_idx_0);

      /* Math: '<S54>/Square' */
      localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1 * localB->rtb_x_R_idx_1;

      /* Sum: '<S59>/Subtract' incorporates:
       *  Product: '<S59>/Multiply'
       *  Product: '<S59>/Multiply1'
       */
      localB->u_n = localB->TmpSignalConversionAtSqua_c[0] *
        localB->rtb_x_R_idx_0 - localB->P_d[0] * localB->rtb_y_R_idx_1;

      /* Signum: '<S55>/Sign1' */
      if (rtIsNaNF(localB->u_n)) {
        localB->rtb_x_R_idx_1 = (rtNaNF);
      } else if (localB->u_n < 0.0F) {
        localB->rtb_x_R_idx_1 = -1.0F;
      } else {
        localB->rtb_x_R_idx_1 = (real32_T)(localB->u_n > 0.0F);
      }

      /* End of Signum: '<S55>/Sign1' */

      /* Trigonometry: '<S55>/Acos' incorporates:
       *  DotProduct: '<S55>/Dot Product'
       */
      if (localB->rtb_ve_idx_0 > 1.0F) {
        localB->rtb_ve_idx_0 = 1.0F;
      } else if (localB->rtb_ve_idx_0 < -1.0F) {
        localB->rtb_ve_idx_0 = -1.0F;
      }

      /* Switch: '<S55>/Switch2' incorporates:
       *  Constant: '<S55>/Constant4'
       */
      if (!(localB->rtb_x_R_idx_1 != 0.0F)) {
        localB->rtb_x_R_idx_1 = 1.0F;
      }

      /* Product: '<S55>/Multiply' incorporates:
       *  Switch: '<S55>/Switch2'
       *  Trigonometry: '<S55>/Acos'
       */
      localB->rtb_x_R_idx_0 = (real32_T)acos(localB->rtb_ve_idx_0) *
        localB->rtb_x_R_idx_1;

      /* Sum: '<S45>/Sum' incorporates:
       *  Constant: '<S45>/Constant'
       */
      localB->rtb_x_R_idx_1 = FMS_PARAM.FW_AIRSPD_TRIM - *rtu_INS_Out_airspeed;

      /* Delay: '<S46>/Delay' incorporates:
       *  Constant: '<S46>/Constant'
       */
      if (localDW->icLoad_l) {
        localDW->Delay_DSTATE_i = FMS_PARAM.FW_HEIGHT_TRIM;
      }

      /* Sum: '<S46>/Sum' incorporates:
       *  Delay: '<S46>/Delay'
       */
      localB->rtb_ve_idx_0 = localDW->Delay_DSTATE_i - *rtu_INS_Out_h_R;

      /* End of Outputs for SubSystem: '<S12>/FormDisband' */
      memset(&localB->Merge, 0, sizeof(FMS_Out_Bus));

      /* Outputs for IfAction SubSystem: '<S12>/FormDisband' incorporates:
       *  ActionPort: '<S18>/Action Port'
       */
      /* BusAssignment: '<S18>/Bus Assignment' incorporates:
       *  Constant: '<S18>/Constant1'
       *  Merge: '<S10>/Merge'
       */
      localB->Merge.state = VehicleState_FormDisband;
      localB->Merge.ax_cmd = localB->rtb_x_R_idx_1;

      /* Saturate: '<S54>/Saturation' */
      if (localB->rtb_x_R_idx_0 > 1.57079637F) {
        localB->rtb_x_R_idx_0 = 1.57079637F;
      } else if (localB->rtb_x_R_idx_0 < -1.57079637F) {
        localB->rtb_x_R_idx_0 = -1.57079637F;
      }

      /* BusAssignment: '<S18>/Bus Assignment' incorporates:
       *  Constant: '<S47>/L1'
       *  Gain: '<S54>/Gain'
       *  Product: '<S54>/Divide'
       *  Product: '<S54>/Multiply1'
       *  Saturate: '<S54>/Saturation'
       *  Trigonometry: '<S54>/Sin'
       */
      localB->Merge.ay_cmd = 2.0F * localB->rtb_vd_idx_2 * (real32_T)sin
        (localB->rtb_x_R_idx_0) / FMS_PARAM.L1;

      /* Gain: '<S46>/Gain2' */
      localB->rtb_x_R_idx_0 = FMS_PARAM.Z_P * localB->rtb_ve_idx_0;

      /* Saturate: '<S46>/Saturation' */
      if (localB->rtb_x_R_idx_0 > CONTROL_PARAM.FW_T_CLMB_MAX) {
        /* BusAssignment: '<S18>/Bus Assignment' */
        localB->Merge.vh_cmd = CONTROL_PARAM.FW_T_CLMB_MAX;
      } else if (localB->rtb_x_R_idx_0 < -CONTROL_PARAM.FW_T_SINK_MAX) {
        /* BusAssignment: '<S18>/Bus Assignment' */
        localB->Merge.vh_cmd = -CONTROL_PARAM.FW_T_SINK_MAX;
      } else {
        /* BusAssignment: '<S18>/Bus Assignment' */
        localB->Merge.vh_cmd = localB->rtb_x_R_idx_0;
      }

      /* End of Saturate: '<S46>/Saturation' */

      /* Update for Delay: '<S50>/start_vel' */
      localDW->icLoad_a = false;

      /* Update for Delay: '<S46>/Delay' */
      localDW->icLoad_l = false;

      /* End of Outputs for SubSystem: '<S12>/FormDisband' */
      break;

     default:
      if (localDW->SwitchCase_ActiveSubsystem_c != rtPrevAction) {
        /* SystemReset for IfAction SubSystem: '<S12>/Default' incorporates:
         *  ActionPort: '<S16>/Action Port'
         */
        /* SystemReset for SwitchCase: '<S12>/Switch Case' */
        Formation_FMS_Default_Reset(&localDW->Default_d);

        /* End of SystemReset for SubSystem: '<S12>/Default' */
      }

      /* Outputs for IfAction SubSystem: '<S12>/Default' incorporates:
       *  ActionPort: '<S16>/Action Port'
       */
      Formation_FMS_Default(rtu_INS_Out_airspeed, rtu_INS_Out_h_R,
                            &localB->Merge, &localDW->Default_d);

      /* End of Outputs for SubSystem: '<S12>/Default' */
      break;
    }

    /* End of Outputs for SubSystem: '<S10>/Form_Subsystem' */
    break;

   case 1:
    if (localDW->SwitchCase_ActiveSubsystem != rtPrevAction) {
      /* InitializeConditions for IfAction SubSystem: '<S10>/Hold' incorporates:
       *  ActionPort: '<S13>/Action Port'
       */
      /* InitializeConditions for SwitchCase: '<S10>/Switch Case' incorporates:
       *  Delay: '<S179>/Delay'
       *  Delay: '<S183>/start_vel'
       */
      localDW->icLoad = true;
      localDW->icLoad_j = true;

      /* End of InitializeConditions for SubSystem: '<S10>/Hold' */
    }

    /* Outputs for IfAction SubSystem: '<S10>/Hold' incorporates:
     *  ActionPort: '<S13>/Action Port'
     */
    /* SignalConversion generated from: '<S189>/Square' */
    localB->TmpSignalConversionAtSqua_c[0] = *rtu_INS_Out_vn;
    localB->TmpSignalConversionAtSqua_c[1] = *rtu_INS_Out_ve;

    /* Sum: '<S190>/Sum of Elements' incorporates:
     *  Math: '<S190>/Math Function'
     *  Sum: '<S189>/Sum of Elements'
     */
    localB->rtb_y_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] *
      localB->TmpSignalConversionAtSqua_c[0] +
      localB->TmpSignalConversionAtSqua_c[1] *
      localB->TmpSignalConversionAtSqua_c[1];

    /* Math: '<S190>/Math Function1' incorporates:
     *  Sum: '<S190>/Sum of Elements'
     *
     * About '<S190>/Math Function1':
     *  Operator: sqrt
     */
    if (localB->rtb_y_R_idx_0 < 0.0F) {
      localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
        (localB->rtb_y_R_idx_0));
    } else {
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_y_R_idx_0);
    }

    /* End of Math: '<S190>/Math Function1' */

    /* Switch: '<S190>/Switch' incorporates:
     *  Constant: '<S190>/Constant'
     *  Product: '<S190>/Product'
     */
    if (localB->rtb_x_R_idx_1 > 0.0F) {
      localB->rtb_vd_idx_0 = *rtu_INS_Out_vn;
      localB->rtb_vd_idx_1 = *rtu_INS_Out_ve;
      localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1;
    } else {
      localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
      localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
      localB->rtb_vd_idx_2 = 1.0F;
    }

    /* End of Switch: '<S190>/Switch' */

    /* Delay: '<S183>/start_vel' */
    if (localDW->icLoad) {
      localDW->start_vel_DSTATE[0] = localB->TmpSignalConversionAtSqua_c[0];
      localDW->start_vel_DSTATE[1] = localB->TmpSignalConversionAtSqua_c[1];
    }

    /* Math: '<S195>/Math Function' incorporates:
     *  Delay: '<S183>/start_vel'
     */
    localB->TmpSignalConversionAtSqua_c[0] = localDW->start_vel_DSTATE[0] *
      localDW->start_vel_DSTATE[0];
    localB->TmpSignalConversionAtSqua_c[1] = localDW->start_vel_DSTATE[1] *
      localDW->start_vel_DSTATE[1];

    /* Sum: '<S195>/Sum of Elements' */
    localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0] +
      localB->TmpSignalConversionAtSqua_c[1];

    /* Math: '<S195>/Math Function1' incorporates:
     *  Sum: '<S195>/Sum of Elements'
     *
     * About '<S195>/Math Function1':
     *  Operator: sqrt
     */
    if (localB->rtb_x_R_idx_0 < 0.0F) {
      localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
        (localB->rtb_x_R_idx_0));
    } else {
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
    }

    /* End of Math: '<S195>/Math Function1' */

    /* Switch: '<S195>/Switch' incorporates:
     *  Constant: '<S195>/Constant'
     *  Delay: '<S183>/start_vel'
     *  Product: '<S195>/Product'
     */
    if (localB->rtb_x_R_idx_1 > 0.0F) {
      localB->rtb_ve_idx_0 = localDW->start_vel_DSTATE[0];
      localB->u_n = localDW->start_vel_DSTATE[1];
      localB->DotProduct = localB->rtb_x_R_idx_1;
    } else {
      localB->rtb_ve_idx_0 = localDW->start_vel_DSTATE[0] * 0.0F;
      localB->u_n = localDW->start_vel_DSTATE[1] * 0.0F;
      localB->DotProduct = 1.0F;
    }

    /* End of Switch: '<S195>/Switch' */

    /* Product: '<S190>/Divide' */
    localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_vd_idx_0 /
      localB->rtb_vd_idx_2;
    localB->TmpSignalConversionAtSqua_c[1] = localB->rtb_vd_idx_1 /
      localB->rtb_vd_idx_2;

    /* SignalConversion generated from: '<S193>/Math Function' */
    localB->P_d[0] = localB->TmpSignalConversionAtSqua_c[1];
    localB->P_d[1] = localB->TmpSignalConversionAtSqua_c[0];

    /* Sum: '<S193>/Sum of Elements' incorporates:
     *  Math: '<S193>/Math Function'
     *  SignalConversion generated from: '<S193>/Math Function'
     */
    localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
      localB->TmpSignalConversionAtSqua_c[1] +
      localB->TmpSignalConversionAtSqua_c[0] *
      localB->TmpSignalConversionAtSqua_c[0];

    /* Math: '<S193>/Math Function1' incorporates:
     *  Sum: '<S193>/Sum of Elements'
     *
     * About '<S193>/Math Function1':
     *  Operator: sqrt
     */
    if (localB->rtb_x_R_idx_0 < 0.0F) {
      localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
        (localB->rtb_x_R_idx_0));
    } else {
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
    }

    /* End of Math: '<S193>/Math Function1' */

    /* Switch: '<S193>/Switch' incorporates:
     *  Constant: '<S193>/Constant'
     *  Product: '<S193>/Product'
     *  SignalConversion generated from: '<S193>/Math Function'
     */
    if (localB->rtb_x_R_idx_1 > 0.0F) {
      localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
      localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0];
      localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1;
    } else {
      localB->rtb_vd_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
      localB->rtb_vd_idx_1 = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
      localB->rtb_vd_idx_2 = 1.0F;
    }

    /* End of Switch: '<S193>/Switch' */

    /* Gain: '<S181>/Gain' incorporates:
     *  Constant: '<S180>/L1'
     */
    localB->rtb_x_R_idx_1 = 0.5F * FMS_PARAM.L1;

    /* MinMax: '<S181>/Max' incorporates:
     *  Constant: '<S180>/R'
     */
    if ((FMS_PARAM.LOITER_R >= localB->rtb_x_R_idx_1) || rtIsNaNF
        (localB->rtb_x_R_idx_1)) {
      localB->rtb_x_R_idx_1 = FMS_PARAM.LOITER_R;
    }

    /* End of MinMax: '<S181>/Max' */

    /* Reshape: '<S181>/Reshape2' */
    localB->Reshape2_bq[0] = *rtu_INS_Out_x_R;
    localB->Reshape2_bq[1] = *rtu_INS_Out_y_R;

    /* MATLAB Function: '<S181>/SearchL1RefWP' incorporates:
     *  Constant: '<S180>/L1'
     */
    Formation_FMS_SearchL1RefWP_i(&localB->Cmd_In.cur_waypoint[0],
      localB->Reshape2_bq, localB->rtb_x_R_idx_1, FMS_PARAM.L1,
      localB->TmpSignalConversionAtSqua_c, &localB->n_j);

    /* RelationalOperator: '<S184>/Compare' incorporates:
     *  Constant: '<S184>/Constant'
     */
    rtb_FixPtRelationalOperator_c = (localB->n_j > 0.0);

    /* Product: '<S195>/Divide' */
    localB->P_i[0] = localB->rtb_ve_idx_0 / localB->DotProduct;
    localB->P_i[1] = localB->u_n / localB->DotProduct;

    /* MATLAB Function: '<S181>/OutRegionRegWP' incorporates:
     *  Constant: '<S180>/L1'
     */
    Formation_FMS_OutRegionRegWP_o(&localB->Cmd_In.cur_waypoint[0],
      localB->Reshape2_bq, localB->rtb_x_R_idx_1, FMS_PARAM.L1, localB->P_i,
      localB->P_d);

    /* Switch: '<S181>/Switch' */
    if (rtb_FixPtRelationalOperator_c) {
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[0];
    } else {
      localB->rtb_x_R_idx_0 = localB->P_d[0];
    }

    /* Sum: '<S182>/Subtract' incorporates:
     *  Switch: '<S181>/Switch'
     */
    localB->Reshape2_bq[0] = localB->rtb_x_R_idx_0 - *rtu_INS_Out_x_R;

    /* Switch: '<S181>/Switch' */
    if (rtb_FixPtRelationalOperator_c) {
      localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
    } else {
      localB->rtb_x_R_idx_0 = localB->P_d[1];
    }

    /* Sum: '<S182>/Subtract' incorporates:
     *  Switch: '<S181>/Switch'
     */
    localB->Reshape2_bq[1] = localB->rtb_x_R_idx_0 - *rtu_INS_Out_y_R;

    /* Sum: '<S191>/Sum of Elements' incorporates:
     *  Math: '<S191>/Math Function'
     *  Sum: '<S182>/Subtract'
     */
    localB->rtb_x_R_idx_0 = localB->Reshape2_bq[0] * localB->Reshape2_bq[0] +
      localB->Reshape2_bq[1] * localB->Reshape2_bq[1];

    /* Math: '<S191>/Math Function1' incorporates:
     *  Sum: '<S191>/Sum of Elements'
     *
     * About '<S191>/Math Function1':
     *  Operator: sqrt
     */
    if (localB->rtb_x_R_idx_0 < 0.0F) {
      localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
        (localB->rtb_x_R_idx_0));
    } else {
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
    }

    /* End of Math: '<S191>/Math Function1' */

    /* Switch: '<S191>/Switch' incorporates:
     *  Constant: '<S191>/Constant'
     *  Product: '<S191>/Product'
     *  Sum: '<S182>/Subtract'
     *  Switch: '<S194>/Switch'
     */
    if (localB->rtb_x_R_idx_1 > 0.0F) {
      localB->rtb_ve_idx_0 = localB->Reshape2_bq[0];
      localB->u_n = localB->Reshape2_bq[1];
      localB->DotProduct = localB->rtb_x_R_idx_1;
    } else {
      localB->rtb_ve_idx_0 = localB->Reshape2_bq[0] * 0.0F;
      localB->u_n = localB->Reshape2_bq[1] * 0.0F;
      localB->DotProduct = 1.0F;
    }

    /* End of Switch: '<S191>/Switch' */

    /* Product: '<S191>/Divide' incorporates:
     *  Product: '<S194>/Divide'
     */
    localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_ve_idx_0 /
      localB->DotProduct;
    localB->TmpSignalConversionAtSqua_c[1] = localB->u_n / localB->DotProduct;

    /* Sum: '<S194>/Sum of Elements' incorporates:
     *  Math: '<S194>/Math Function'
     *  SignalConversion generated from: '<S194>/Math Function'
     */
    localB->rtb_x_R_idx_0 = localB->TmpSignalConversionAtSqua_c[1] *
      localB->TmpSignalConversionAtSqua_c[1] +
      localB->TmpSignalConversionAtSqua_c[0] *
      localB->TmpSignalConversionAtSqua_c[0];

    /* Math: '<S194>/Math Function1' incorporates:
     *  Sum: '<S194>/Sum of Elements'
     *
     * About '<S194>/Math Function1':
     *  Operator: sqrt
     */
    if (localB->rtb_x_R_idx_0 < 0.0F) {
      localB->rtb_x_R_idx_1 = -(real32_T)sqrt((real32_T)fabs
        (localB->rtb_x_R_idx_0));
    } else {
      localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_x_R_idx_0);
    }

    /* End of Math: '<S194>/Math Function1' */

    /* Switch: '<S194>/Switch' incorporates:
     *  Constant: '<S194>/Constant'
     *  Product: '<S194>/Product'
     *  SignalConversion generated from: '<S194>/Math Function'
     */
    if (localB->rtb_x_R_idx_1 > 0.0F) {
      localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1];
      localB->u_n = localB->TmpSignalConversionAtSqua_c[0];
      localB->DotProduct = localB->rtb_x_R_idx_1;
    } else {
      localB->rtb_ve_idx_0 = localB->TmpSignalConversionAtSqua_c[1] * 0.0F;
      localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * 0.0F;
      localB->DotProduct = 1.0F;
    }

    /* End of Switch: '<S194>/Switch' */

    /* Product: '<S194>/Divide' */
    localB->rtb_y_R_idx_1 = localB->rtb_ve_idx_0 / localB->DotProduct;

    /* Product: '<S193>/Divide' */
    localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_0 / localB->rtb_vd_idx_2;

    /* DotProduct: '<S188>/Dot Product' */
    localB->rtb_ve_idx_0 = localB->rtb_x_R_idx_0 * localB->rtb_y_R_idx_1;
    localB->TmpSignalConversionAtSqua_c[0] = localB->rtb_y_R_idx_1;
    localB->P_d[0] = localB->rtb_x_R_idx_0;

    /* Product: '<S194>/Divide' incorporates:
     *  Product: '<S193>/Divide'
     */
    localB->rtb_y_R_idx_1 = localB->u_n / localB->DotProduct;

    /* Product: '<S193>/Divide' */
    localB->rtb_x_R_idx_0 = localB->rtb_vd_idx_1 / localB->rtb_vd_idx_2;

    /* DotProduct: '<S188>/Dot Product' */
    localB->rtb_ve_idx_0 += localB->rtb_x_R_idx_0 * localB->rtb_y_R_idx_1;

    /* Sqrt: '<S189>/Sqrt' */
    localB->rtb_x_R_idx_1 = (real32_T)sqrt(localB->rtb_y_R_idx_0);

    /* Math: '<S187>/Square' */
    localB->rtb_vd_idx_2 = localB->rtb_x_R_idx_1 * localB->rtb_x_R_idx_1;

    /* Sum: '<S192>/Subtract' incorporates:
     *  Product: '<S192>/Multiply'
     *  Product: '<S192>/Multiply1'
     */
    localB->u_n = localB->TmpSignalConversionAtSqua_c[0] * localB->rtb_x_R_idx_0
      - localB->P_d[0] * localB->rtb_y_R_idx_1;

    /* Signum: '<S188>/Sign1' */
    if (rtIsNaNF(localB->u_n)) {
      localB->rtb_x_R_idx_1 = (rtNaNF);
    } else if (localB->u_n < 0.0F) {
      localB->rtb_x_R_idx_1 = -1.0F;
    } else {
      localB->rtb_x_R_idx_1 = (real32_T)(localB->u_n > 0.0F);
    }

    /* End of Signum: '<S188>/Sign1' */

    /* Trigonometry: '<S188>/Acos' incorporates:
     *  DotProduct: '<S188>/Dot Product'
     */
    if (localB->rtb_ve_idx_0 > 1.0F) {
      localB->rtb_ve_idx_0 = 1.0F;
    } else if (localB->rtb_ve_idx_0 < -1.0F) {
      localB->rtb_ve_idx_0 = -1.0F;
    }

    /* Switch: '<S188>/Switch2' incorporates:
     *  Constant: '<S188>/Constant4'
     */
    if (!(localB->rtb_x_R_idx_1 != 0.0F)) {
      localB->rtb_x_R_idx_1 = 1.0F;
    }

    /* Product: '<S188>/Multiply' incorporates:
     *  Switch: '<S188>/Switch2'
     *  Trigonometry: '<S188>/Acos'
     */
    localB->rtb_x_R_idx_0 = (real32_T)acos(localB->rtb_ve_idx_0) *
      localB->rtb_x_R_idx_1;

    /* Sum: '<S178>/Sum' incorporates:
     *  Constant: '<S178>/Constant'
     */
    localB->rtb_x_R_idx_1 = FMS_PARAM.FW_AIRSPD_TRIM - *rtu_INS_Out_airspeed;

    /* Delay: '<S179>/Delay' incorporates:
     *  Constant: '<S179>/Constant'
     */
    if (localDW->icLoad_j) {
      localDW->Delay_DSTATE = FMS_PARAM.FW_HEIGHT_TRIM;
    }

    /* Sum: '<S179>/Sum' incorporates:
     *  Delay: '<S179>/Delay'
     */
    localB->rtb_ve_idx_0 = localDW->Delay_DSTATE - *rtu_INS_Out_h_R;

    /* End of Outputs for SubSystem: '<S10>/Hold' */
    memset(&localB->Merge, 0, sizeof(FMS_Out_Bus));

    /* Outputs for IfAction SubSystem: '<S10>/Hold' incorporates:
     *  ActionPort: '<S13>/Action Port'
     */
    /* BusAssignment: '<S13>/Bus Assignment' incorporates:
     *  Constant: '<S13>/Constant1'
     *  Merge: '<S10>/Merge'
     */
    localB->Merge.state = VehicleState_Hold;
    localB->Merge.ax_cmd = localB->rtb_x_R_idx_1;

    /* Saturate: '<S187>/Saturation' */
    if (localB->rtb_x_R_idx_0 > 1.57079637F) {
      localB->rtb_x_R_idx_0 = 1.57079637F;
    } else if (localB->rtb_x_R_idx_0 < -1.57079637F) {
      localB->rtb_x_R_idx_0 = -1.57079637F;
    }

    /* BusAssignment: '<S13>/Bus Assignment' incorporates:
     *  Constant: '<S180>/L1'
     *  Gain: '<S187>/Gain'
     *  Product: '<S187>/Divide'
     *  Product: '<S187>/Multiply1'
     *  Saturate: '<S187>/Saturation'
     *  Trigonometry: '<S187>/Sin'
     */
    localB->Merge.ay_cmd = 2.0F * localB->rtb_vd_idx_2 * (real32_T)sin
      (localB->rtb_x_R_idx_0) / FMS_PARAM.L1;

    /* Gain: '<S179>/Gain2' */
    localB->rtb_x_R_idx_0 = FMS_PARAM.Z_P * localB->rtb_ve_idx_0;

    /* Saturate: '<S179>/Saturation' */
    if (localB->rtb_x_R_idx_0 > CONTROL_PARAM.FW_T_CLMB_MAX) {
      /* BusAssignment: '<S13>/Bus Assignment' */
      localB->Merge.vh_cmd = CONTROL_PARAM.FW_T_CLMB_MAX;
    } else if (localB->rtb_x_R_idx_0 < -CONTROL_PARAM.FW_T_SINK_MAX) {
      /* BusAssignment: '<S13>/Bus Assignment' */
      localB->Merge.vh_cmd = -CONTROL_PARAM.FW_T_SINK_MAX;
    } else {
      /* BusAssignment: '<S13>/Bus Assignment' */
      localB->Merge.vh_cmd = localB->rtb_x_R_idx_0;
    }

    /* End of Saturate: '<S179>/Saturation' */

    /* Update for Delay: '<S183>/start_vel' */
    localDW->icLoad = false;

    /* Update for Delay: '<S179>/Delay' */
    localDW->icLoad_j = false;

    /* End of Outputs for SubSystem: '<S10>/Hold' */
    break;

   default:
    if (localDW->SwitchCase_ActiveSubsystem != rtPrevAction) {
      /* SystemReset for IfAction SubSystem: '<S10>/Default' incorporates:
       *  ActionPort: '<S11>/Action Port'
       */
      /* SystemReset for SwitchCase: '<S10>/Switch Case' */
      Formation_FMS_Default_Reset(&localDW->Default);

      /* End of SystemReset for SubSystem: '<S10>/Default' */
    }

    /* Outputs for IfAction SubSystem: '<S10>/Default' incorporates:
     *  ActionPort: '<S11>/Action Port'
     */
    Formation_FMS_Default(rtu_INS_Out_airspeed, rtu_INS_Out_h_R, &localB->Merge,
                          &localDW->Default);

    /* End of Outputs for SubSystem: '<S10>/Default' */
    break;
  }

  /* End of SwitchCase: '<S10>/Switch Case' */

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_type[0] = localB->Other_Mission_Data.type[0];

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_valid_it[0] = localB->Other_Mission_Data.valid_items[0];

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_type[1] = localB->Other_Mission_Data.type[1];

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_valid_it[1] = localB->Other_Mission_Data.valid_items[1];

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_type[2] = localB->Other_Mission_Data.type[2];

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  rty_Other_Mission_Data_valid_it[2] = localB->Other_Mission_Data.valid_items[2];

  /* SignalConversion generated from: '<Root>/Other_Mission_Data' */
  *rty_Other_Mission_Data_timestam = localB->Other_Mission_Data.timestamp;

  /* SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_timestamp = localB->Merge.timestamp;

  /* SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_state = localB->Merge.state;

  /* SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_ax_cmd = localB->Merge.ax_cmd;

  /* SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_ay_cmd = localB->Merge.ay_cmd;

  /* SignalConversion generated from: '<Root>/FMS_Out' */
  *rty_FMS_Out_vh_cmd = localB->Merge.vh_cmd;
}

/* Model initialize function */
void Formation_FMS_initialize(const char_T **rt_errorStatus,
  RT_MODEL_Formation_FMS_T *const Formation_FMS_M, ZCE_Formation_FMS_T *localZCE)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize error status */
  rtmSetErrorStatusPointer(Formation_FMS_M, rt_errorStatus);
  localZCE->Mission_SubSystem_Reset_ZCE = POS_ZCSIG;
  localZCE->FormMission_SubSystem_Reset_ZCE = POS_ZCSIG;
}

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
