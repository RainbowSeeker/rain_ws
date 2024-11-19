/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: Formation_FMS.h
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

#ifndef RTW_HEADER_Formation_FMS_h_
#define RTW_HEADER_Formation_FMS_h_
#ifndef Formation_FMS_COMMON_INCLUDES_
#define Formation_FMS_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* Formation_FMS_COMMON_INCLUDES_ */

#include "rtGetInf.h"
#include "rt_nonfinite.h"
#include "zero_crossing_types.h"

/* Forward declaration for rtModel */
typedef struct tag_RTM_Formation_FMS_T RT_MODEL_Formation_FMS_T;

#ifndef DEFINED_TYPEDEF_FOR_Pilot_Cmd_Bus_
#define DEFINED_TYPEDEF_FOR_Pilot_Cmd_Bus_

typedef struct {
  uint32_T timestamp;
  uint32_T mode;
} Pilot_Cmd_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_Mission_Data_Bus_
#define DEFINED_TYPEDEF_FOR_Mission_Data_Bus_

typedef struct {
  uint32_T timestamp;

  /* FormAssemble(1),
     FormDisband(2),
     FormMission(3),
     Mission(4), */
  uint32_T type;
  uint8_T valid_items;
  real32_T x[8];
  real32_T y[8];
  real32_T z[8];
  real32_T heading[8];
  real32_T ext1[8];
  real32_T ext2[8];
} Mission_Data_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_INS_Out_Bus_
#define DEFINED_TYPEDEF_FOR_INS_Out_Bus_

typedef struct {
  uint32_T timestamp;
  real32_T phi;
  real32_T theta;
  real32_T psi;
  real32_T p;
  real32_T q;
  real32_T r;

  /* Quaternion */
  real32_T quat[4];
  real32_T x_R;
  real32_T y_R;
  real32_T h_R;
  real32_T airspeed;
  real32_T ax;
  real32_T ay;
  real32_T az;
  real32_T vn;
  real32_T ve;
  real32_T vd;
} INS_Out_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_Formation_Cross_Bus_
#define DEFINED_TYPEDEF_FOR_Formation_Cross_Bus_

/* Supports up to 3 drone. */
typedef struct {
  uint32_T timestamp[3];
  real32_T x_R[3];
  real32_T y_R[3];
  real32_T h_R[3];
  real32_T vn[3];
  real32_T ve[3];
  real32_T vd[3];
  real32_T ay[3];
} Formation_Cross_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_PilotMode_
#define DEFINED_TYPEDEF_FOR_PilotMode_

/* enumeration of pilot mode */
typedef enum {
  PilotMode_None = 0,                  /* Default value */
  PilotMode_Hold,
  PilotMode_FormAssemble,
  PilotMode_FormHold,
  PilotMode_FormMission,
  PilotMode_FormDisband
} PilotMode;

#endif

#ifndef DEFINED_TYPEDEF_FOR_VehicleState_
#define DEFINED_TYPEDEF_FOR_VehicleState_

/* enumeration to track active leaf state of FMS/FMS State Machine/Vehicle */
typedef enum {
  VehicleState_None = 0,               /* Default value */
  VehicleState_Hold,
  VehicleState_FormAssemble,
  VehicleState_FormHold,
  VehicleState_FormMission,
  VehicleState_FormDisband
} VehicleState;

#endif

#ifndef DEFINED_TYPEDEF_FOR_FMS_Out_Bus_
#define DEFINED_TYPEDEF_FOR_FMS_Out_Bus_

typedef struct {
  uint32_T timestamp;
  VehicleState state;
  real32_T ax_cmd;
  real32_T ay_cmd;
  real32_T vh_cmd;
} FMS_Out_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_Commander_In_Bus_
#define DEFINED_TYPEDEF_FOR_Commander_In_Bus_

typedef struct {
  uint32_T form_valid;
  real32_T sp_waypoint[3];
  real32_T cur_waypoint[3];
} Commander_In_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_Other_Mission_Data_Bus_
#define DEFINED_TYPEDEF_FOR_Other_Mission_Data_Bus_

typedef struct {
  uint32_T timestamp;

  /* FormAssemble(1),
     FormDisband(2),
     FormMission(3),
     Mission(4), */
  uint32_T type[3];
  uint8_T valid_items[3];
  real32_T x[24];
  real32_T y[24];
  real32_T z[24];
  real32_T heading[24];
  real32_T ext1[24];
  real32_T ext2[24];
} Other_Mission_Data_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_0NI0yCUy2KXHIXsCBP30IB_
#define DEFINED_TYPEDEF_FOR_struct_0NI0yCUy2KXHIXsCBP30IB_

typedef struct {
  uint32_T UAV_ID;
  uint32_T NUM_UAV;
  real32_T ADJ_MARTIX[9];
  real32_T REL_X_MATRIX[9];
  real32_T REL_Y_MATRIX[9];
  real32_T REL_Z_MATRIX[9];
  real32_T FORM_POINT[9];
  real32_T DISBAND_POINT[9];
  real32_T FORM_RADIUS;
  real32_T FORM_HGT_KP;
  real32_T FORM_POS_KP;
  real32_T FORM_POUT_LIM;
  real32_T FORM_VEL_KP;
  real32_T FORM_VOUT_LIM;
  real32_T FORM_HEAD_KP;
  real32_T FORM_HOUT_LIM;
  real32_T FORM_HRFF_THRE;
} struct_0NI0yCUy2KXHIXsCBP30IB;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_FebwIpiU9Ih55vl7WG22GB_
#define DEFINED_TYPEDEF_FOR_struct_FebwIpiU9Ih55vl7WG22GB_

typedef struct {
  real32_T ROLL_P;
  real32_T PITCH_P;
  real32_T ROLL_RATE_P;
  real32_T PITCH_RATE_P;
  real32_T YAW_RATE_P;
  real32_T ROLL_RATE_I;
  real32_T PITCH_RATE_I;
  real32_T YAW_RATE_I;
  real32_T RATE_I_MIN;
  real32_T RATE_I_MAX;
  real32_T TRIM_ROLL;
  real32_T TRIM_PITCH;
  real32_T TRIM_YAW;
  real32_T FW_PSP_OFF;
  real32_T FW_AIRSPD_MIN;
  real32_T FW_AIRSPD_MAX;
  real32_T FW_AIRSPD_TRIM;
  real32_T FW_AIRSPD_STALL;
  int32_T FW_ARSP_MODE;
  int32_T FW_ARSP_SCALE_EN;
  real32_T FW_T_TAS_TC;
  real32_T FW_T_I_GAIN_PIT;
  real32_T FW_T_I_GAIN_THR;
  real32_T FW_T_THR_DAMP;
  real32_T FW_T_SPDWEIGHT;
  real32_T FW_T_CLMB_MAX;
  real32_T FW_T_SINK_MIN;
  real32_T FW_T_SINK_MAX;
  real32_T FW_T_CLMB_R_SP;
  real32_T FW_T_SINK_R_SP;
  real32_T FW_P_LIM_MAX;
  real32_T FW_P_LIM_MIN;
  real32_T FW_R_LIM;
  real32_T FW_T_VERT_ACC;
  real32_T FW_T_ALT_TC;
  real32_T FW_T_HRATE_FF;
  real32_T FW_T_RLL2THR;
  real32_T FW_RR_FF;
  real32_T FW_PR_FF;
  real32_T FW_YR_FF;
  real32_T FW_R_RMAX;
  real32_T FW_P_RMAX;
  real32_T FW_Y_RMAX;
  real32_T FW_THR_MAX;
  real32_T FW_THR_MIN;
  real32_T FW_THR_TRIM;
  real32_T FW_T_SEB_R_FF;
  real32_T FW_T_PTCH_DAMP;
} struct_FebwIpiU9Ih55vl7WG22GB;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_krKDLIJ9OZPHOAucYn9ayF_
#define DEFINED_TYPEDEF_FOR_struct_krKDLIJ9OZPHOAucYn9ayF_

typedef struct {
  real32_T FW_AIRSPD_TRIM;
  real32_T FW_HEIGHT_TRIM;
  real32_T FW_RADIUS_RATIO;
  real32_T AIRSPD_P;
  real32_T Z_P;
  real32_T L1;
  real32_T ACCEPT_R;
  real32_T LOITER_R;
  real32_T ACC_Y_LIM;
} struct_krKDLIJ9OZPHOAucYn9ayF;

#endif

/* Custom Type definition for MATLAB Function: '<S196>/MATLAB Function' */
#ifndef struct_tag_lliwpDezBRCIPtkwTpZvmC
#define struct_tag_lliwpDezBRCIPtkwTpZvmC

struct tag_lliwpDezBRCIPtkwTpZvmC
{
  real_T contents;
};

#endif                                 /* struct_tag_lliwpDezBRCIPtkwTpZvmC */

#ifndef typedef_captured_var_Formation_FMS_T
#define typedef_captured_var_Formation_FMS_T

typedef struct tag_lliwpDezBRCIPtkwTpZvmC captured_var_Formation_FMS_T;

#endif                                /* typedef_captured_var_Formation_FMS_T */

#ifndef struct_tag_sGXMO9PmoeW2JqDMJvYKnvG
#define struct_tag_sGXMO9PmoeW2JqDMJvYKnvG

struct tag_sGXMO9PmoeW2JqDMJvYKnvG
{
  real_T xs;
  real_T ys;
  real_T psi_s;
  real_T xf;
  real_T yf;
  real_T psi_f;
  real_T v;
  real_T r;
  real_T pos[4];
  real_T pof[4];
  real_T xts[4];
  real_T yts[4];
  real_T xtf[4];
  real_T ytf[4];
  real_T cs[4];
  real_T cf[4];
  real_T lt[4];
  real_T l[4];
  real_T index_dubins[2];
  real_T l_ad;
  real_T precision_flag;
  real_T xm;
  real_T ym;
};

#endif                                 /* struct_tag_sGXMO9PmoeW2JqDMJvYKnvG */

#ifndef typedef_sGXMO9PmoeW2JqDMJvYKnvG_Forma_T
#define typedef_sGXMO9PmoeW2JqDMJvYKnvG_Forma_T

typedef struct tag_sGXMO9PmoeW2JqDMJvYKnvG sGXMO9PmoeW2JqDMJvYKnvG_Forma_T;

#endif                             /* typedef_sGXMO9PmoeW2JqDMJvYKnvG_Forma_T */

/* Block states (default storage) for system '<S12>/Default' */
typedef struct {
  real32_T Delay_DSTATE;               /* '<S22>/Delay' */
  boolean_T icLoad;                    /* '<S22>/Delay' */
} DW_Default_Formation_FMS_T;

/* Block signals for model 'Formation_FMS' */
typedef struct {
  sGXMO9PmoeW2JqDMJvYKnvG_Forma_T object[3];
  real_T result[75];
  Other_Mission_Data_Bus Other_Mission_Data;/* '<Root>/FMS State Machine' */
  sGXMO9PmoeW2JqDMJvYKnvG_Forma_T obj;
  Formation_Cross_Bus BusConversion_InsertedFor_FMS_p;
  INS_Out_Bus BusConversion_InsertedFor_FMS_c;
  real_T Goal[9];
  real32_T xyz_O_nx3[9];               /* '<S89>/Vector Concatenate' */
  real32_T vNED_O_nx3[9];              /* '<S89>/Vector Concatenate1' */
  real32_T b_tmp[9];
  real32_T BusConversion_InsertedFor_FMSSt[8];
  real32_T BusConversion_InsertedFor_FM_cl[8];
  real32_T BusConversion_InsertedFor_FMS_k[8];
  Commander_In_Bus Cmd_In;             /* '<Root>/FMS State Machine' */
  real_T b_data[4];
  real_T c_data[4];
  FMS_Out_Bus Merge;                   /* '<S10>/Merge' */
  int32_T i_data[4];
  uint32_T e[3];
  uint32_T b[3];
  real32_T vn[3];                      /* '<S9>/Signal Copy2' */
  real32_T h_R[3];                     /* '<S9>/Signal Copy2' */
  real_T SFunction_o9;                 /* '<Root>/FMS State Machine' */
  real_T SFunction_o10;                /* '<Root>/FMS State Machine' */
  real_T SFunction_o11;                /* '<Root>/FMS State Machine' */
  real_T SFunction_o12;                /* '<Root>/FMS State Machine' */
  real_T n_j;                          /* '<S66>/SearchL1RefWP' */
  real_T l_ref;
  real_T search_floor;
  real_T search;
  real_T goal;
  real_T l_best;
  real_T pos_tmp;
  real_T xts_tmp;
  real_T yts_tmp;
  real_T cs_tmp;
  real_T cs_tmp_c;
  real_T cf_tmp;
  real_T lt_tmp;
  real_T deltax;
  real_T deltay;
  real32_T unit_center_to_pose[2];
  real32_T TmpSignalConversionAtSqua_c[2];
  real32_T Reshape2_bq[2];             /* '<S29>/Reshape2' */
  real32_T P_d[2];                     /* '<S30>/SearchL1RefWP' */
  real32_T P_i[2];                     /* '<S30>/OutRegionRegWP' */
  real32_T P_b[2];                     /* '<S163>/OutRegionRegWP' */
  captured_var_Formation_FMS_T PhiMaximum;
  captured_var_Formation_FMS_T rad2deg;
  real32_T dpsi;
  real32_T scale;
  real32_T absxk;
  real32_T t;
  real32_T u_n;                        /* '<S30>/SearchL1RefWP' */
  real32_T DotProduct;                 /* '<S171>/Dot Product' */
  real32_T rtb_vd_idx_0;
  real32_T rtb_vd_idx_1;
  real32_T rtb_vd_idx_2;
  real32_T rtb_ve_idx_0;
  real32_T rtb_TmpSignalConversionAtSqua_b;
  real32_T rtb_x_R_idx_0;
  real32_T rtb_y_R_idx_0;
  real32_T rtb_x_R_idx_1;
  real32_T rtb_y_R_idx_1;
  int32_T idx;
  int32_T i;
  int32_T rtb_vn_tmp;
  int32_T k;
  int32_T target;
  int32_T stop_flag;
  uint32_T BusConversion_InsertedFor_FM_pb;
  uint32_T BusConversion_InsertedFor_FM_cv;
  VehicleState state;                  /* '<Root>/FMS State Machine' */
  uint16_T wp_index;                   /* '<Root>/FMS State Machine' */
} B_Formation_FMS_c_T;

/* Block states (default storage) for model 'Formation_FMS' */
typedef struct {
  real32_T start_vel_DSTATE[2];        /* '<S183>/start_vel' */
  real32_T Delay_DSTATE;               /* '<S179>/Delay' */
  real32_T start_vel_DSTATE_m[2];      /* '<S50>/start_vel' */
  real32_T Delay_DSTATE_i;             /* '<S46>/Delay' */
  real32_T Delay_DSTATE_o;             /* '<S156>/Delay' */
  real32_T Integrator_DSTATE;          /* '<S137>/Integrator' */
  real32_T start_vel_DSTATE_n[2];      /* '<S68>/start_vel' */
  real32_T Delay_DSTATE_a;             /* '<S64>/Delay' */
  real32_T Delay_DSTATE_l;             /* '<S26>/Delay' */
  real32_T DiscreteTimeIntegrator1_DSTATE;/* '<S4>/Discrete-Time Integrator1' */
  uint32_T DelayInput1_DSTATE;         /* '<S7>/Delay Input1' */
  PilotMode Delay_DSTATE_j;            /* '<S5>/Delay' */
  uint32_T Mission_Data_timestamp_prev;/* '<Root>/FMS State Machine' */
  uint32_T Mission_Data_timestamp_start;/* '<Root>/FMS State Machine' */
  PilotMode mode_prev;                 /* '<Root>/FMS State Machine' */
  PilotMode mode_start;                /* '<Root>/FMS State Machine' */
  uint16_T DelayInput1_DSTATE_d;       /* '<S81>/Delay Input1' */
  uint16_T DelayInput1_DSTATE_h;       /* '<S23>/Delay Input1' */
  int8_T SwitchCase_ActiveSubsystem;   /* '<S10>/Switch Case' */
  int8_T SwitchCase_ActiveSubsystem_c; /* '<S12>/Switch Case' */
  uint8_T is_Vehicle;                  /* '<Root>/FMS State Machine' */
  uint8_T is_Formation;                /* '<Root>/FMS State Machine' */
  uint8_T is_FormAssemble;             /* '<Root>/FMS State Machine' */
  uint8_T is_FormMission;              /* '<Root>/FMS State Machine' */
  uint8_T is_active_c3_Formation_FMS;  /* '<Root>/FMS State Machine' */
  uint8_T temporalCounter_i1;          /* '<Root>/FMS State Machine' */
  boolean_T icLoad;                    /* '<S183>/start_vel' */
  boolean_T icLoad_j;                  /* '<S179>/Delay' */
  boolean_T icLoad_a;                  /* '<S50>/start_vel' */
  boolean_T icLoad_l;                  /* '<S46>/Delay' */
  boolean_T icLoad_k;                  /* '<S156>/Delay' */
  boolean_T icLoad_l3;                 /* '<S68>/start_vel' */
  boolean_T icLoad_p;                  /* '<S64>/Delay' */
  boolean_T icLoad_h;                  /* '<S26>/Delay' */
  DW_Default_Formation_FMS_T Default;  /* '<S10>/Default' */
  DW_Default_Formation_FMS_T Default_d;/* '<S12>/Default' */
} DW_Formation_FMS_f_T;

/* Zero-crossing (trigger) state for model 'Formation_FMS' */
typedef struct {
  real_T FormMission_SubSystem_Reset_ZC;/* '<S20>/FormMission_SubSystem' */
  real_T Mission_SubSystem_Reset_ZC;   /* '<S17>/Mission_SubSystem' */
} ZCV_Formation_FMS_g_T;

/* Zero-crossing (trigger) state for model 'Formation_FMS' */
typedef struct {
  ZCSigState FormMission_SubSystem_Reset_ZCE;/* '<S20>/FormMission_SubSystem' */
  ZCSigState Mission_SubSystem_Reset_ZCE;/* '<S17>/Mission_SubSystem' */
} ZCE_Formation_FMS_T;

/* Real-time Model Data Structure */
struct tag_RTM_Formation_FMS_T {
  const char_T **errorStatus;
};

typedef struct {
  B_Formation_FMS_c_T rtb;
  DW_Formation_FMS_f_T rtdw;
  RT_MODEL_Formation_FMS_T rtm;
  ZCE_Formation_FMS_T rtzce;
} MdlrefDW_Formation_FMS_T;

/*
 * Exported Global Parameters
 *
 * Note: Exported global parameters are tunable parameters with an exported
 * global storage class designation.  Code generation will declare the memory for
 * these parameters and exports their symbols.
 *
 */
extern struct_0NI0yCUy2KXHIXsCBP30IB FORMATION_PARAM;/* Variable: FORMATION_PARAM
                                                      * Referenced by:
                                                      *   '<Root>/FMS State Machine'
                                                      *   '<S196>/MATLAB Function'
                                                      *   '<S197>/MATLAB Function'
                                                      *   '<S82>/Constant2'
                                                      *   '<S157>/min_radius'
                                                      *   '<S89>/Calc_Position_Velocity_Setpoint'
                                                      *   '<S90>/Vel_Kp'
                                                      *   '<S90>/Saturation'
                                                      *   '<S91>/Gain'
                                                      *   '<S91>/Saturation'
                                                      *   '<S92>/Pout_Saturation'
                                                      *   '<S92>/POS_KP'
                                                      *   '<S95>/Constant'
                                                      *   '<S142>/Proportional Gain'
                                                      */
extern struct_FebwIpiU9Ih55vl7WG22GB CONTROL_PARAM;/* Variable: CONTROL_PARAM
                                                    * Referenced by:
                                                    *   '<S15>/Saturation'
                                                    *   '<S179>/Saturation'
                                                    *   '<S22>/Saturation'
                                                    *   '<S46>/Saturation'
                                                    *   '<S64>/Saturation'
                                                    *   '<S26>/Saturation'
                                                    *   '<S156>/Saturation'
                                                    *   '<S157>/min_radius'
                                                    *   '<S157>/v^2'
                                                    */
extern struct_krKDLIJ9OZPHOAucYn9ayF FMS_PARAM;/* Variable: FMS_PARAM
                                                * Referenced by:
                                                *   '<Root>/ACCEPT_R'
                                                *   '<S14>/Constant'
                                                *   '<S15>/Gain2'
                                                *   '<S178>/Constant'
                                                *   '<S179>/Constant'
                                                *   '<S179>/Gain2'
                                                *   '<S180>/L1'
                                                *   '<S180>/R'
                                                *   '<S21>/Constant'
                                                *   '<S22>/Gain2'
                                                *   '<S45>/Constant'
                                                *   '<S46>/Constant'
                                                *   '<S46>/Gain2'
                                                *   '<S47>/L1'
                                                *   '<S47>/R'
                                                *   '<S63>/Constant'
                                                *   '<S64>/Constant'
                                                *   '<S64>/Gain2'
                                                *   '<S65>/L1'
                                                *   '<S65>/R'
                                                *   '<S25>/Constant'
                                                *   '<S26>/Constant'
                                                *   '<S26>/Gain2'
                                                *   '<S28>/L1'
                                                *   '<S28>/Saturation1'
                                                *   '<S155>/Constant'
                                                *   '<S155>/Gain'
                                                *   '<S156>/Constant'
                                                *   '<S156>/Gain2'
                                                *   '<S157>/Satefy'
                                                *   '<S161>/L1'
                                                *   '<S161>/Saturation1'
                                                */

/* Model reference registration function */
extern void Formation_FMS_initialize(const char_T **rt_errorStatus,
  RT_MODEL_Formation_FMS_T *const Formation_FMS_M, ZCE_Formation_FMS_T *localZCE);
extern void Formation_FMS_NearbyRefWP(const real32_T rtu_P2[2], const real32_T
  rtu_P3[2], real32_T rtu_L1, real32_T rty_P[2], real32_T *rty_d);
extern void Formation_FMS_OutRegionRegWP(const real32_T rtu_P1[2], const
  real32_T rtu_P2[2], const real32_T rtu_P3[2], real32_T rty_P[2]);
extern void Formation_FMS_SearchL1RefWP(const real32_T rtu_P1[2], const real32_T
  rtu_P2[2], const real32_T rtu_P3[2], real32_T rtu_L1, real32_T rty_P[2],
  real32_T *rty_u);
extern void Formation_FMS_OutRegionRegWP_o(const real32_T rtu_P0[2], const
  real32_T rtu_P_Vehicle[2], real32_T rtu_R, real32_T rtu_L1, const real32_T
  rtu_n[2], real32_T rty_P[2]);
extern void Formation_FMS_SearchL1RefWP_i(const real32_T rtu_P_0[2], const
  real32_T rtu_P_Vehicle[2], real32_T rtu_R, real32_T rtu_L1, real32_T rty_P[2],
  real_T *rty_n);
extern void Formation_FMS_Default_Init(DW_Default_Formation_FMS_T *localDW);
extern void Formation_FMS_Default_Reset(DW_Default_Formation_FMS_T *localDW);
extern void Formation_FMS_Default(const real32_T *rtu_FMS_In, const real32_T
  *rtu_FMS_In_h, FMS_Out_Bus *rty_FMS_Out, DW_Default_Formation_FMS_T *localDW);
extern void Formation_FMS_Init(uint32_T *rty_FMS_Out_timestamp, VehicleState
  *rty_FMS_Out_state, real32_T *rty_FMS_Out_ax_cmd, real32_T *rty_FMS_Out_ay_cmd,
  real32_T *rty_FMS_Out_vh_cmd, uint32_T *rty_Other_Mission_Data_timestam,
  uint32_T rty_Other_Mission_Data_type[3], uint8_T
  rty_Other_Mission_Data_valid_it[3], real32_T rty_Other_Mission_Data_x[24],
  real32_T rty_Other_Mission_Data_y[24], real32_T rty_Other_Mission_Data_z[24],
  real32_T rty_Other_Mission_Data_heading[24], real32_T
  rty_Other_Mission_Data_ext1[24], real32_T rty_Other_Mission_Data_ext2[24],
  B_Formation_FMS_c_T *localB, DW_Formation_FMS_f_T *localDW);
extern void Formation_FMS_Disable(DW_Formation_FMS_f_T *localDW);
extern void Formation_FMS(const uint32_T *rtu_Pilot_Cmd_timestamp, const
  uint32_T *rtu_Pilot_Cmd_mode, const uint32_T *rtu_Mission_Data_timestamp,
  const uint32_T *rtu_Mission_Data_type, const uint8_T
  *rtu_Mission_Data_valid_items, const real32_T rtu_Mission_Data_x[8], const
  real32_T rtu_Mission_Data_y[8], const real32_T rtu_Mission_Data_z[8], const
  uint32_T *rtu_INS_Out_timestamp, const real32_T *rtu_INS_Out_phi, const
  real32_T *rtu_INS_Out_theta, const real32_T *rtu_INS_Out_psi, const real32_T
  *rtu_INS_Out_p, const real32_T *rtu_INS_Out_q, const real32_T *rtu_INS_Out_r,
  const real32_T rtu_INS_Out_quat[4], const real32_T *rtu_INS_Out_x_R, const
  real32_T *rtu_INS_Out_y_R, const real32_T *rtu_INS_Out_h_R, const real32_T
  *rtu_INS_Out_airspeed, const real32_T *rtu_INS_Out_ax, const real32_T
  *rtu_INS_Out_ay, const real32_T *rtu_INS_Out_az, const real32_T
  *rtu_INS_Out_vn, const real32_T *rtu_INS_Out_ve, const real32_T
  *rtu_INS_Out_vd, const uint32_T rtu_Formation_Cross_timestamp[3], const
  real32_T rtu_Formation_Cross_x_R[3], const real32_T rtu_Formation_Cross_y_R[3],
  const real32_T rtu_Formation_Cross_h_R[3], const real32_T
  rtu_Formation_Cross_vn[3], const real32_T rtu_Formation_Cross_ve[3], const
  real32_T rtu_Formation_Cross_vd[3], const real32_T rtu_Formation_Cross_ay[3],
  uint32_T *rty_FMS_Out_timestamp, VehicleState *rty_FMS_Out_state, real32_T
  *rty_FMS_Out_ax_cmd, real32_T *rty_FMS_Out_ay_cmd, real32_T
  *rty_FMS_Out_vh_cmd, uint32_T *rty_Other_Mission_Data_timestam, uint32_T
  rty_Other_Mission_Data_type[3], uint8_T rty_Other_Mission_Data_valid_it[3],
  real32_T rty_Other_Mission_Data_x[24], real32_T rty_Other_Mission_Data_y[24],
  real32_T rty_Other_Mission_Data_z[24], real32_T
  rty_Other_Mission_Data_heading[24], real32_T rty_Other_Mission_Data_ext1[24],
  real32_T rty_Other_Mission_Data_ext2[24], B_Formation_FMS_c_T *localB,
  DW_Formation_FMS_f_T *localDW, ZCE_Formation_FMS_T *localZCE);

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S43>/Data Type Duplicate' : Unused code path elimination
 * Block '<S44>/Data Type Duplicate' : Unused code path elimination
 * Block '<S40>/Data Type Duplicate' : Unused code path elimination
 * Block '<S41>/Data Type Duplicate' : Unused code path elimination
 * Block '<S60>/Data Type Duplicate' : Unused code path elimination
 * Block '<S61>/Data Type Duplicate' : Unused code path elimination
 * Block '<S57>/Data Type Duplicate' : Unused code path elimination
 * Block '<S58>/Data Type Duplicate' : Unused code path elimination
 * Block '<S62>/Data Type Duplicate' : Unused code path elimination
 * Block '<S78>/Data Type Duplicate' : Unused code path elimination
 * Block '<S79>/Data Type Duplicate' : Unused code path elimination
 * Block '<S75>/Data Type Duplicate' : Unused code path elimination
 * Block '<S76>/Data Type Duplicate' : Unused code path elimination
 * Block '<S80>/Data Type Duplicate' : Unused code path elimination
 * Block '<S103>/Sqrt' : Unused code path elimination
 * Block '<S103>/Square' : Unused code path elimination
 * Block '<S103>/Sum of Elements' : Unused code path elimination
 * Block '<S104>/Sqrt' : Unused code path elimination
 * Block '<S104>/Square' : Unused code path elimination
 * Block '<S104>/Sum of Elements' : Unused code path elimination
 * Block '<S159>/Data Type Duplicate' : Unused code path elimination
 * Block '<S159>/Data Type Propagation' : Unused code path elimination
 * Block '<S176>/Data Type Duplicate' : Unused code path elimination
 * Block '<S177>/Data Type Duplicate' : Unused code path elimination
 * Block '<S173>/Data Type Duplicate' : Unused code path elimination
 * Block '<S174>/Data Type Duplicate' : Unused code path elimination
 * Block '<S193>/Data Type Duplicate' : Unused code path elimination
 * Block '<S194>/Data Type Duplicate' : Unused code path elimination
 * Block '<S190>/Data Type Duplicate' : Unused code path elimination
 * Block '<S191>/Data Type Duplicate' : Unused code path elimination
 * Block '<S195>/Data Type Duplicate' : Unused code path elimination
 * Block '<S9>/Signal Copy3' : Eliminate redundant signal conversion block
 * Block '<S9>/Signal Copy5' : Eliminate redundant signal conversion block
 * Block '<S9>/Signal Copy6' : Eliminate redundant signal conversion block
 * Block '<S14>/Gain' : Eliminated nontunable gain of 1
 * Block '<S21>/Gain' : Eliminated nontunable gain of 1
 * Block '<S25>/Gain' : Eliminated nontunable gain of 1
 * Block '<S30>/Reshape' : Reshape block reduction
 * Block '<S29>/Reshape' : Reshape block reduction
 * Block '<S29>/Reshape1' : Reshape block reduction
 * Block '<S45>/Gain' : Eliminated nontunable gain of 1
 * Block '<S48>/Reshape' : Reshape block reduction
 * Block '<S48>/Reshape1' : Reshape block reduction
 * Block '<S48>/Reshape3' : Reshape block reduction
 * Block '<S63>/Gain' : Eliminated nontunable gain of 1
 * Block '<S66>/Reshape' : Reshape block reduction
 * Block '<S66>/Reshape1' : Reshape block reduction
 * Block '<S66>/Reshape3' : Reshape block reduction
 * Block '<S91>/Gain1' : Eliminated nontunable gain of 1
 * Block '<S163>/Reshape' : Reshape block reduction
 * Block '<S162>/Reshape' : Reshape block reduction
 * Block '<S162>/Reshape1' : Reshape block reduction
 * Block '<S178>/Gain' : Eliminated nontunable gain of 1
 * Block '<S181>/Reshape' : Reshape block reduction
 * Block '<S181>/Reshape1' : Reshape block reduction
 * Block '<S181>/Reshape3' : Reshape block reduction
 */

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Formation_FMS'
 * '<S1>'   : 'Formation_FMS/CommandProcess'
 * '<S2>'   : 'Formation_FMS/FMS Commander'
 * '<S3>'   : 'Formation_FMS/FMS State Machine'
 * '<S4>'   : 'Formation_FMS/CommandProcess/Check Valid'
 * '<S5>'   : 'Formation_FMS/CommandProcess/Mode PreProcess'
 * '<S6>'   : 'Formation_FMS/CommandProcess/Check Valid/Compare To Constant1'
 * '<S7>'   : 'Formation_FMS/CommandProcess/Check Valid/Detect Change1'
 * '<S8>'   : 'Formation_FMS/CommandProcess/Mode PreProcess/Compare To Zero1'
 * '<S9>'   : 'Formation_FMS/FMS Commander/Bus Creator'
 * '<S10>'  : 'Formation_FMS/FMS Commander/Commander'
 * '<S11>'  : 'Formation_FMS/FMS Commander/Commander/Default'
 * '<S12>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem'
 * '<S13>'  : 'Formation_FMS/FMS Commander/Commander/Hold'
 * '<S14>'  : 'Formation_FMS/FMS Commander/Commander/Default/Airspeed Hold'
 * '<S15>'  : 'Formation_FMS/FMS Commander/Commander/Default/Altitude_Hold'
 * '<S16>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/Default'
 * '<S17>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble'
 * '<S18>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband'
 * '<S19>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold'
 * '<S20>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission'
 * '<S21>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/Default/Airspeed Hold'
 * '<S22>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/Default/Altitude_Hold'
 * '<S23>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Detect Change'
 * '<S24>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem'
 * '<S25>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Airspeed Hold'
 * '<S26>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Altitude_Hold'
 * '<S27>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command'
 * '<S28>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control'
 * '<S29>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/WayPoints'
 * '<S30>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1 Reference WP'
 * '<S31>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core'
 * '<S32>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1 Reference WP/Compare To Constant'
 * '<S33>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1 Reference WP/Compare To Constant1'
 * '<S34>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1 Reference WP/NearbyRefWP'
 * '<S35>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1 Reference WP/OutRegionRegWP'
 * '<S36>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1 Reference WP/SearchL1RefWP'
 * '<S37>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration'
 * '<S38>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle'
 * '<S39>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Vector Modulus'
 * '<S40>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize'
 * '<S41>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize1'
 * '<S42>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle/2D Cross Product'
 * '<S43>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize'
 * '<S44>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormAssemble/Mission_SubSystem/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize1'
 * '<S45>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Airspeed Hold'
 * '<S46>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Altitude_Hold'
 * '<S47>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control'
 * '<S48>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1 Reference WP'
 * '<S49>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core'
 * '<S50>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/Reference_Point'
 * '<S51>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1 Reference WP/Compare To Constant'
 * '<S52>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1 Reference WP/OutRegionRegWP'
 * '<S53>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1 Reference WP/SearchL1RefWP'
 * '<S54>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration'
 * '<S55>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle'
 * '<S56>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Modulus'
 * '<S57>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize'
 * '<S58>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize1'
 * '<S59>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/2D Cross Product'
 * '<S60>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize'
 * '<S61>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize1'
 * '<S62>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormDisband/Loiter_Control/Reference_Point/Vector Normalize'
 * '<S63>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Airspeed Hold'
 * '<S64>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Altitude_Hold'
 * '<S65>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control'
 * '<S66>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1 Reference WP'
 * '<S67>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core'
 * '<S68>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/Reference_Point'
 * '<S69>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1 Reference WP/Compare To Constant'
 * '<S70>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1 Reference WP/OutRegionRegWP'
 * '<S71>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1 Reference WP/SearchL1RefWP'
 * '<S72>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration'
 * '<S73>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle'
 * '<S74>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Modulus'
 * '<S75>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize'
 * '<S76>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize1'
 * '<S77>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/2D Cross Product'
 * '<S78>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize'
 * '<S79>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize1'
 * '<S80>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormHold/Loiter_Control/Reference_Point/Vector Normalize'
 * '<S81>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/Detect Change'
 * '<S82>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem'
 * '<S83>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Compare To Constant'
 * '<S84>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller'
 * '<S85>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller'
 * '<S86>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller'
 * '<S87>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control'
 * '<S88>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller'
 * '<S89>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Vector_Mapping'
 * '<S90>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Forward_Control'
 * '<S91>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control'
 * '<S92>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Vector_Combination'
 * '<S93>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Forward_Control/Vector Modulus'
 * '<S94>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Forward_Control/Vector Modulus1'
 * '<S95>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/Compare To Constant'
 * '<S96>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/Vector Modulus'
 * '<S97>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/Vector Modulus1'
 * '<S98>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/Vector_Angle'
 * '<S99>'  : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/Vector_Angle1'
 * '<S100>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/psi_err_saturation'
 * '<S101>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Lateral_Control/psi_err_saturation/Compare To Constant'
 * '<S102>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Vector_Combination/Pout_Saturation'
 * '<S103>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Vector_Combination/Vector Modulus'
 * '<S104>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Forward_Lateral_control/Vector_Combination/Vector Modulus1'
 * '<S105>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller'
 * '<S106>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Anti-windup'
 * '<S107>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/D Gain'
 * '<S108>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Filter'
 * '<S109>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Filter ICs'
 * '<S110>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/I Gain'
 * '<S111>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Ideal P Gain'
 * '<S112>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Ideal P Gain Fdbk'
 * '<S113>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Integrator'
 * '<S114>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Integrator ICs'
 * '<S115>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/N Copy'
 * '<S116>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/N Gain'
 * '<S117>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/P Copy'
 * '<S118>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Parallel P Gain'
 * '<S119>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Reset Signal'
 * '<S120>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Saturation'
 * '<S121>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Saturation Fdbk'
 * '<S122>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Sum'
 * '<S123>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Sum Fdbk'
 * '<S124>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tracking Mode'
 * '<S125>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tracking Mode Sum'
 * '<S126>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tsamp - Integral'
 * '<S127>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tsamp - Ngain'
 * '<S128>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/postSat Signal'
 * '<S129>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/preSat Signal'
 * '<S130>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Anti-windup/Passthrough'
 * '<S131>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/D Gain/Disabled'
 * '<S132>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Filter/Disabled'
 * '<S133>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Filter ICs/Disabled'
 * '<S134>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/I Gain/Internal Parameters'
 * '<S135>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Ideal P Gain/Passthrough'
 * '<S136>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Ideal P Gain Fdbk/Disabled'
 * '<S137>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Integrator/Discrete'
 * '<S138>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Integrator ICs/Internal IC'
 * '<S139>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/N Copy/Disabled wSignal Specification'
 * '<S140>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/N Gain/Disabled'
 * '<S141>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/P Copy/Disabled'
 * '<S142>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Parallel P Gain/Internal Parameters'
 * '<S143>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Reset Signal/Disabled'
 * '<S144>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Saturation/Passthrough'
 * '<S145>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Saturation Fdbk/Disabled'
 * '<S146>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Sum/Sum_PI'
 * '<S147>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Sum Fdbk/Disabled'
 * '<S148>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tracking Mode/Disabled'
 * '<S149>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tracking Mode Sum/Passthrough'
 * '<S150>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tsamp - Integral/Passthrough'
 * '<S151>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/Tsamp - Ngain/Passthrough'
 * '<S152>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/postSat Signal/Forward_Path'
 * '<S153>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Height_Controller/Discrete PID Controller/preSat Signal/Forward_Path'
 * '<S154>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Formation_Controller/Consensus Controller/Vector_Mapping/Calc_Position_Velocity_Setpoint'
 * '<S155>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Airspeed Hold'
 * '<S156>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Altitude_Hold'
 * '<S157>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/FormRadiusLimit'
 * '<S158>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command'
 * '<S159>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/FormRadiusLimit/Saturation Dynamic'
 * '<S160>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/FormRadiusLimit/min_radius'
 * '<S161>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control'
 * '<S162>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/WayPoints'
 * '<S163>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1 Reference WP'
 * '<S164>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core'
 * '<S165>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1 Reference WP/Compare To Constant'
 * '<S166>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1 Reference WP/Compare To Constant1'
 * '<S167>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1 Reference WP/NearbyRefWP'
 * '<S168>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1 Reference WP/OutRegionRegWP'
 * '<S169>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1 Reference WP/SearchL1RefWP'
 * '<S170>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration'
 * '<S171>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle'
 * '<S172>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Vector Modulus'
 * '<S173>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize'
 * '<S174>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize1'
 * '<S175>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle/2D Cross Product'
 * '<S176>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize'
 * '<S177>' : 'Formation_FMS/FMS Commander/Commander/Form_Subsystem/FormMission/FormMission_SubSystem/Leader_Controller/Position Command/L1 Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize1'
 * '<S178>' : 'Formation_FMS/FMS Commander/Commander/Hold/Airspeed Hold'
 * '<S179>' : 'Formation_FMS/FMS Commander/Commander/Hold/Altitude_Hold'
 * '<S180>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control'
 * '<S181>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1 Reference WP'
 * '<S182>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core'
 * '<S183>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/Reference_Point'
 * '<S184>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1 Reference WP/Compare To Constant'
 * '<S185>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1 Reference WP/OutRegionRegWP'
 * '<S186>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1 Reference WP/SearchL1RefWP'
 * '<S187>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration'
 * '<S188>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle'
 * '<S189>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Modulus'
 * '<S190>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize'
 * '<S191>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Vector Normalize1'
 * '<S192>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/2D Cross Product'
 * '<S193>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize'
 * '<S194>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/L1_Core/L1_Lateral_Acceleration/Included Angle/Vector Normalize1'
 * '<S195>' : 'Formation_FMS/FMS Commander/Commander/Hold/Loiter_Control/Reference_Point/Vector Normalize'
 * '<S196>' : 'Formation_FMS/FMS State Machine/Vehicle.Formation.FormAssemble.dubinsPath'
 * '<S197>' : 'Formation_FMS/FMS State Machine/Vehicle.Formation.check_form_valid'
 * '<S198>' : 'Formation_FMS/FMS State Machine/Vehicle.Formation.FormAssemble.dubinsPath/MATLAB Function'
 * '<S199>' : 'Formation_FMS/FMS State Machine/Vehicle.Formation.check_form_valid/MATLAB Function'
 */
#endif                                 /* RTW_HEADER_Formation_FMS_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
