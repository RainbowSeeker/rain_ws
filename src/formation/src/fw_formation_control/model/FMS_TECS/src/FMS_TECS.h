/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: FMS_TECS.h
 *
 * Code generated for Simulink model 'FMS_TECS'.
 *
 * Model version                  : 1.43
 * Simulink Coder version         : 9.8 (R2022b) 13-May-2022
 * C/C++ source code generated on : Fri Aug 16 11:40:32 2024
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Intel->x86-64 (Linux 64)
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_FMS_TECS_h_
#define RTW_HEADER_FMS_TECS_h_
#ifndef FMS_TECS_COMMON_INCLUDES_
#define FMS_TECS_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif                                 /* FMS_TECS_COMMON_INCLUDES_ */

#include "Formation_FMS.h"
#include "zero_crossing_types.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetErrorStatusPointer
#define rtmGetErrorStatusPointer(rtm)  ((const char_T **)(&((rtm)->errorStatus)))
#endif

/* Forward declaration for rtModel */
typedef struct tag_RTM_FMS_TECS_T RT_MODEL_FMS_TECS_T;

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

/* Block signals (default storage) */
typedef struct {
  real32_T ske_rate;                   /* '<S11>/Product3' */
  real32_T out_d;                      /* '<S17>/Integrator_Clamping' */
  real32_T Product3;                   /* '<S14>/Product3' */
  real32_T Add_l;                      /* '<S13>/Add' */
  real32_T Product2;                   /* '<S14>/Product2' */
  real32_T DiscreteTimeIntegrator_h;   /* '<S23>/Discrete-Time Integrator' */
  real32_T Product;                    /* '<S14>/Product' */
  real32_T out;                        /* '<S7>/Integrator_Clamping' */
  real32_T Product_h;                  /* '<S21>/Product' */
  real32_T rtb_new_state_idx_0;
  real32_T innovation_idx_0;
  real32_T ay_cmd;                     /* '<Root>/Model' */
  real32_T ax_cmd;                     /* '<Root>/Model' */
  real32_T vh_cmd;                     /* '<Root>/Model' */
} B_FMS_TECS_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  real32_T DiscreteTimeIntegrator_DSTATE;/* '<S17>/Discrete-Time Integrator' */
  real32_T Delay_DSTATE;               /* '<S12>/Delay' */
  real32_T DiscreteTimeIntegrator5_DSTATE;/* '<S24>/Discrete-Time Integrator5' */
  real32_T DiscreteTimeIntegrator_DSTATE_l;/* '<S23>/Discrete-Time Integrator' */
  real32_T Delay_DSTATE_n;             /* '<S13>/Delay' */
  real32_T kalman_gain[4];             /* '<S6>/MATLAB Function' */
  real32_T airspeed_state[2];          /* '<S6>/MATLAB Function' */
  int8_T DiscreteTimeIntegrator_PrevRese;/* '<S17>/Discrete-Time Integrator' */
  int8_T DiscreteTimeIntegrator_PrevRe_m;/* '<S23>/Discrete-Time Integrator' */
  uint8_T DiscreteTimeIntegrator5_IC_LOAD;/* '<S24>/Discrete-Time Integrator5' */
  MdlrefDW_Formation_FMS_T Model_InstanceData;/* '<Root>/Model' */
} DW_FMS_TECS_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  Pilot_Cmd_Bus Pilot_Cmd;             /* '<Root>/Pilot_Cmd' */
  Mission_Data_Bus Mission_Data;       /* '<Root>/Mission_Data' */
  INS_Out_Bus INS_Out;                 /* '<Root>/INS_Out' */
  Formation_Cross_Bus Formation_Cross; /* '<Root>/Formation_Cross' */
} ExtU_FMS_TECS_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  FMS_Out_Bus FMS_Out;                 /* '<Root>/FMS_Out' */
  real32_T att_cmd[2];                 /* '<Root>/att_cmd' */
  real32_T throttle_cmd;               /* '<Root>/throttle_cmd' */
  Other_Mission_Data_Bus Other_Mission_Data;/* '<Root>/Other_Mission_Data' */
} ExtY_FMS_TECS_T;

/* Real-time Model Data Structure */
struct tag_RTM_FMS_TECS_T {
  const char_T *errorStatus;
};

/* Block signals (default storage) */
extern B_FMS_TECS_T FMS_TECS_B;

/* Block states (default storage) */
extern DW_FMS_TECS_T FMS_TECS_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_FMS_TECS_T FMS_TECS_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_FMS_TECS_T FMS_TECS_Y;

/*
 * Exported Global Parameters
 *
 * Note: Exported global parameters are tunable parameters with an exported
 * global storage class designation.  Code generation will declare the memory for
 * these parameters and exports their symbols.
 *
 */
extern struct_0NI0yCUy2KXHIXsCBP30IB FORMATION_PARAM;/* Variable: FORMATION_PARAM
                                                      * Referenced by: '<Root>/Model'
                                                      */
extern struct_FebwIpiU9Ih55vl7WG22GB CONTROL_PARAM;/* Variable: CONTROL_PARAM
                                                    * Referenced by:
                                                    *   '<Root>/Model'
                                                    *   '<S3>/Saturation'
                                                    *   '<S4>/Saturation'
                                                    *   '<S7>/Integrator_Clamping'
                                                    *   '<S7>/Constant'
                                                    *   '<S7>/Constant1'
                                                    *   '<S13>/Saturation'
                                                    *   '<S14>/Constant'
                                                    *   '<S16>/pitch_damping_gain'
                                                    *   '<S16>/seb_rate_ff'
                                                    *   '<S17>/Integrator_Clamping'
                                                    *   '<S17>/I_Gain'
                                                    *   '<S18>/Pitch_Saturation'
                                                    *   '<S18>/Saturation'
                                                    *   '<S21>/throttle_above_trim'
                                                    *   '<S21>/throttle_below_trim'
                                                    *   '<S21>/throttle_trim_adjusted'
                                                    *   '<S21>/Saturation'
                                                    *   '<S22>/STE_rate_to_throttle'
                                                    *   '<S22>/throttle_damping_gain'
                                                    *   '<S23>/Integrator_Clamping'
                                                    *   '<S23>/STE_rate_to_throttle'
                                                    *   '<S23>/I_Gain'
                                                    *   '<S25>/load_factor_correction'
                                                    */
extern struct_krKDLIJ9OZPHOAucYn9ayF FMS_PARAM;/* Variable: FMS_PARAM
                                                * Referenced by: '<Root>/Model'
                                                */

/* Model entry point functions */
extern void FMS_TECS_initialize(void);
extern void FMS_TECS_step(void);
extern void FMS_TECS_terminate(void);

/* Real-time Model object */
extern RT_MODEL_FMS_TECS_T *const FMS_TECS_M;

/*-
 * These blocks were eliminated from the model due to optimizations:
 *
 * Block '<S10>/Data Type Duplicate' : Unused code path elimination
 * Block '<S10>/Data Type Propagation' : Unused code path elimination
 * Block '<S3>/Constant' : Unused code path elimination
 * Block '<S24>/Data Type Conversion' : Eliminate redundant data type conversion
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
 * '<Root>' : 'FMS_TECS'
 * '<S1>'   : 'FMS_TECS/Lateral_Longitudinal_Control'
 * '<S2>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Airspeed_Filter'
 * '<S3>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Altitude_Control'
 * '<S4>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Lateral_Control'
 * '<S5>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control'
 * '<S6>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Airspeed_Filter/Airspeed_Filter'
 * '<S7>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Airspeed_Filter/Calc_TAS_Rate_Setpoint'
 * '<S8>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Airspeed_Filter/Airspeed_Filter/MATLAB Function'
 * '<S9>'   : 'FMS_TECS/Lateral_Longitudinal_Control/Airspeed_Filter/Calc_TAS_Rate_Setpoint/Integrator_Clamping'
 * '<S10>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Airspeed_Filter/Calc_TAS_Rate_Setpoint/Saturation Dynamic'
 * '<S11>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Calc_Specific_Energy_Rates'
 * '<S12>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control'
 * '<S13>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control'
 * '<S14>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/Calc_Seb_Rate'
 * '<S15>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/Calc_climb_angle_to_SEB_rate'
 * '<S16>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/FF & Damp'
 * '<S17>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/I'
 * '<S18>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/Pitch_Saturation'
 * '<S19>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/I/Integrator_Clamping'
 * '<S20>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Pitch_Control/Pitch_Saturation/Pitch_Saturation'
 * '<S21>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control/Calc_STE_Rate'
 * '<S22>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control/Damp'
 * '<S23>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control/I'
 * '<S24>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control/Calc_STE_Rate/First Order LPF'
 * '<S25>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control/Calc_STE_Rate/induced_drag_scales'
 * '<S26>'  : 'FMS_TECS/Lateral_Longitudinal_Control/Longitudinal_Control/Throttle_Control/I/Integrator_Clamping'
 */
#endif                                 /* RTW_HEADER_FMS_TECS_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
