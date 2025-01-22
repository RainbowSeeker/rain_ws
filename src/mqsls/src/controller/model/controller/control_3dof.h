/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * File: control_3dof.h
 *
 * Code generated for Simulink model 'control_3dof'.
 *
 * Model version                  : 1.755
 * Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
 * C/C++ source code generated on : Tue Jan 14 10:43:17 2025
 *
 * Target selection: ert.tlc
 * Embedded hardware selection: Custom Processor->Custom Processor
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef control_3dof_h_
#define control_3dof_h_
#ifndef control_3dof_COMMON_INCLUDES_
#define control_3dof_COMMON_INCLUDES_
#include <stdbool.h>
#include <stdint.h>
#include "math.h"
#endif                                 /* control_3dof_COMMON_INCLUDES_ */

/* Macros for accessing real-time model data structure */
#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

/* Forward declaration for rtModel */
typedef struct tag_RTM_control_3dof_T RT_MODEL_control_3dof_T;

#ifndef DEFINED_TYPEDEF_FOR_Payload_Out_Bus_
#define DEFINED_TYPEDEF_FOR_Payload_Out_Bus_

typedef struct {
  uint64_t timestamp;
  double pL[3];
  double vL[3];
  double p_1[3];
  double p_2[3];
  double p_3[3];
  double v_1[3];
  double v_2[3];
  double v_3[3];
} Payload_Out_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_TrajCmd_In_Bus_
#define DEFINED_TYPEDEF_FOR_TrajCmd_In_Bus_

typedef struct {
  uint64_t timestamp;
  double pos_sp[3];
  double vel_sp[3];
  double acc_ff[3];
} TrajCmd_In_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_DirCmd_In_Bus_
#define DEFINED_TYPEDEF_FOR_DirCmd_In_Bus_

/* 3x Cable Direction Cmd. */
typedef struct {
  uint64_t timestamp;
  double q_sp1[3];
  double q_sp2[3];
  double q_sp3[3];
} DirCmd_In_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_Control_In_Bus_
#define DEFINED_TYPEDEF_FOR_Control_In_Bus_

typedef struct {
  uint64_t timestamp;

  /* Used by 6dof */
  double vforce[3];

  /* Used by 3dof */
  double q_sp[3];

  /* unit vector for direction of cable i */
  double q[3];
  double w[3];
  double ai_sp[3];
  double pi[3];
  double vi[3];
} Control_In_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_Control_State_Bus_
#define DEFINED_TYPEDEF_FOR_Control_State_Bus_

typedef struct {
  uint64_t timestamp;

  /* Payload Disturbance Force */
  double dL[3];
  double margin;
} Control_State_Bus;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_SXZiMrrc3I4047EGQQvBKH_
#define DEFINED_TYPEDEF_FOR_struct_SXZiMrrc3I4047EGQQvBKH_

typedef struct {
  double MASS_LOAD;
  double MASS_UAV;
  double CABLE_LEN;
  double TENSION_MIN;
  double TENSION_MAX;
  double ESO_PL[3];
  double ESO_VL[3];
  double ESO_PI[3];
  double ESO_VI[3];
  double KP;
  double KV;
  double KQ;
  double KW;
  double KQI;
} struct_SXZiMrrc3I4047EGQQvBKH;

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_s1X7MNNw1KKuxbD9HCjUEH_
#define DEFINED_TYPEDEF_FOR_struct_s1X7MNNw1KKuxbD9HCjUEH_

typedef struct {
  int32_t period;
} struct_s1X7MNNw1KKuxbD9HCjUEH;

#endif

#ifndef struct_tag_7XesmMQLqBKi3FDUE8hRQE
#define struct_tag_7XesmMQLqBKi3FDUE8hRQE

struct tag_7XesmMQLqBKi3FDUE8hRQE
{
  int32_t isInitialized;
  double K_p[3];
  double K_v[3];
  double Z_1[3];
  double Z_2[3];
  double Z_3[3];
  bool is_init;
};

#endif                                 /* struct_tag_7XesmMQLqBKi3FDUE8hRQE */

#ifndef typedef_pos_2nd_eso_control_3dof_T
#define typedef_pos_2nd_eso_control_3dof_T

typedef struct tag_7XesmMQLqBKi3FDUE8hRQE pos_2nd_eso_control_3dof_T;

#endif                                 /* typedef_pos_2nd_eso_control_3dof_T */

/* Block signals for system '<S1>/MATLAB System' */
typedef struct {
  double MATLABSystem_o3[3];           /* '<S1>/MATLAB System' */
} B_MATLABSystem_control_3dof_T;

/* Block states (default storage) for system '<S1>/MATLAB System' */
typedef struct {
  pos_2nd_eso_control_3dof_T obj;      /* '<S1>/MATLAB System' */
  bool objisempty;                     /* '<S1>/MATLAB System' */
} DW_MATLABSystem_control_3dof_T;

/* Block signals for system '<S5>/Parallel Control' */
typedef struct {
  double rtu_q[9];
  double rtp_mi[9];
} B_ParallelControl_control_3do_T;

/* Block signals for system '<S6>/Vertical Control' */
typedef struct {
  double Sqi[9];
  double Sqi_m[9];
  double e_q[3];
  double e_w[3];
  double w_sp_tmp[3];
  double e_w_tmp[3];
  double b_a[3];
  double w_sp_tmp_c[3];
  double a;
  double rtu_q;
  double e_w_tmp_k;
  double e_q_idx_1;
  double e_q_idx_2;
  double Sqi_c;
  double Sqi_b;
} B_VerticalControl_control_3do_T;

/* Block signals (default storage) */
typedef struct {
  double A_wrench_tmp[18];
  double state_f[12];                  /* '<S6>/Vertical Control' */
  double Q_new[9];                     /* '<S20>/Vector Concatenate' */
  double rtb_Q_new_m[9];
  double dv[9];
  double V[9];
  double U[9];
  double b_A[9];
  double Vf[9];
  double b_wrench[6];
  double B[6];
  double B_tmp[6];
  double q_2[3];
  double q_3[3];
  double F_trim[3];
  double xv[3];
  double e_1[3];
  double e_2[3];
  double rtb_control_in1_w[3];
  double rtb_control_in2_vforce[3];
  double rtb_control_in2_w[3];
  double rtb_control_in3_w[3];
  double f_vertical[3];                /* '<S6>/Vertical Control' */
  double dot_err_g[3];                 /* '<S6>/Vertical Control' */
  double f_parallel[3];                /* '<S5>/Parallel Control' */
  double dot_err_c[3];                 /* '<S10>/Vertical Control' */
  double dot_err[3];                   /* '<S14>/Vertical Control' */
  double s[3];
  double b_s[3];
  double e[3];
  double work[3];
  double alpha;
  double upper;
  double lower;
  double A_wrench_tmp_c;
  double e_1_k;
  double e_2_c;
  double e_2_b;
  double e_2_p;
  double UnitDelay_DSTATE;
  double UnitDelay_DSTATE_c;
  double UnitDelay_DSTATE_f;
  double F_trim_g;
  double F_trim_g1;
  double absx;
  double cscale;
  double anrm;
  double nrm;
  double rt;
  double r;
  double smm1;
  double emm1;
  double sqds;
  double shift;
  double cfromc;
  double ctoc;
  double cfrom1;
  double cto1;
  double mul;
  double cfromc_m;
  double ctoc_n;
  double cfrom1_p;
  double cto1_l;
  double mul_j;
  double roe;
  double absa;
  double absb;
  double scale;
  double ads;
  double bds;
  double scale_d;
  double absxk;
  double t;
  double scale_g;
  double absxk_l;
  double t_d;
  double scale_dy;
  double absxk_lx;
  double t_o;
  double temp;
  double temp_tmp;
  double absxk_b;
  double temp_n;
  int32_t c_k;
  int32_t i;
  int32_t r_b;
  int32_t vcol;
  int32_t j;
  int32_t ar;
  int32_t b;
  int32_t ib;
  int32_t b_ic;
  int32_t qp1;
  int32_t qq;
  int32_t qjj;
  int32_t m;
  int32_t qs;
  int32_t e_k;
  int32_t offset;
  int32_t j_l;
  int32_t b_i;
  int32_t i_h;
  int32_t offset_b;
  int32_t j_d;
  int32_t b_i_e;
  int32_t i1;
  int32_t kend;
  int32_t k;
  int32_t kend_b;
  int32_t k_j;
  int32_t k_f;
  int32_t k_a;
  int32_t i2;
  int32_t k_ju;
  int32_t i3;
  int32_t k_jz;
  int32_t i4;
  int32_t k_o;
  bool FixPtRelationalOperator;        /* '<S22>/FixPt Relational Operator' */
  bool p;
  bool doscale;
  bool apply_transform;
  bool notdone;
  bool notdone_n;
  B_VerticalControl_control_3do_T sf_VerticalControl_h;/* '<S14>/Vertical Control' */
  B_ParallelControl_control_3do_T sf_ParallelControl_p;/* '<S13>/Parallel Control' */
  B_MATLABSystem_control_3dof_T MATLABSystem_a;/* '<S1>/MATLAB System' */
  B_VerticalControl_control_3do_T sf_VerticalControl_m;/* '<S10>/Vertical Control' */
  B_ParallelControl_control_3do_T sf_ParallelControl_k;/* '<S9>/Parallel Control' */
  B_MATLABSystem_control_3dof_T MATLABSystem_g;/* '<S1>/MATLAB System' */
  B_VerticalControl_control_3do_T sf_VerticalControl;/* '<S6>/Vertical Control' */
  B_ParallelControl_control_3do_T sf_ParallelControl;/* '<S5>/Parallel Control' */
  B_MATLABSystem_control_3dof_T MATLABSystem;/* '<S1>/MATLAB System' */
} B_control_3dof_T;

/* Block states (default storage) for system '<Root>' */
typedef struct {
  pos_2nd_eso_control_3dof_T obj;      /* '<S4>/Position 2nd ESO' */
  double Delay_DSTATE[9];              /* '<S20>/Delay' */
  double UnitDelay_DSTATE[3];          /* '<S4>/Unit Delay' */
  double DiscreteTimeIntegrator_DSTATE[3];/* '<S14>/Discrete-Time Integrator' */
  double UnitDelay_DSTATE_l[3];        /* '<S3>/Unit Delay' */
  double DiscreteTimeIntegrator_DSTATE_l[3];/* '<S10>/Discrete-Time Integrator' */
  double UnitDelay_DSTATE_a[3];        /* '<S2>/Unit Delay' */
  double DiscreteTimeIntegrator_DSTATE_p[3];/* '<S6>/Discrete-Time Integrator' */
  double UnitDelay_DSTATE_f[3];        /* '<S1>/Unit Delay' */
  uint64_t DelayInput1_DSTATE;         /* '<S22>/Delay Input1' */
  uint64_t Accumulator_DSTATE;         /* '<S21>/Accumulator' */
  bool icLoad;                         /* '<S20>/Delay' */
  DW_MATLABSystem_control_3dof_T MATLABSystem_a;/* '<S1>/MATLAB System' */
  DW_MATLABSystem_control_3dof_T MATLABSystem_g;/* '<S1>/MATLAB System' */
  DW_MATLABSystem_control_3dof_T MATLABSystem;/* '<S1>/MATLAB System' */
} DW_control_3dof_T;

/* External inputs (root inport signals with default storage) */
typedef struct {
  Payload_Out_Bus Payload_Out;         /* '<Root>/Payload_Out' */
  TrajCmd_In_Bus Traj_sp;              /* '<Root>/Traj_sp' */
  DirCmd_In_Bus Dir_sp;                /* '<Root>/Dir_sp' */
} ExtU_control_3dof_T;

/* External outputs (root outports fed by signals with default storage) */
typedef struct {
  double force_sp1[3];                 /* '<Root>/force_sp1' */
  double force_sp2[3];                 /* '<Root>/force_sp2' */
  double force_sp3[3];                 /* '<Root>/force_sp3' */
  Control_State_Bus state;             /* '<Root>/state' */
} ExtY_control_3dof_T;

/* Real-time Model Data Structure */
struct tag_RTM_control_3dof_T {
  const char * volatile errorStatus;
};

/* Block signals (default storage) */
extern B_control_3dof_T control_3dof_B;

/* Block states (default storage) */
extern DW_control_3dof_T control_3dof_DW;

/* External inputs (root inport signals with default storage) */
extern ExtU_control_3dof_T control_3dof_U;

/* External outputs (root outports fed by signals with default storage) */
extern ExtY_control_3dof_T control_3dof_Y;

/*
 * Exported Global Parameters
 *
 * Note: Exported global parameters are tunable parameters with an exported
 * global storage class designation.  Code generation will declare the memory for
 * these parameters and exports their symbols.
 *
 */
extern struct_SXZiMrrc3I4047EGQQvBKH CONTROL_PARAM;/* Variable: CONTROL_PARAM
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
                                                    *   '<S18>/KP'
                                                    *   '<S18>/KV'
                                                    */
extern struct_s1X7MNNw1KKuxbD9HCjUEH CONTROL_EXPORT;/* Variable: CONTROL_EXPORT
                                                     * Referenced by: '<S21>/Period [ms]'
                                                     */

/* Model entry point functions */
extern void control_3dof_initialize(void);
extern void control_3dof_step(void);
extern void control_3dof_terminate(void);

/* Real-time Model object */
extern RT_MODEL_control_3dof_T *const control_3dof_M;

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
 * '<Root>' : 'control_3dof'
 * '<S1>'   : 'control_3dof/Cable Controller'
 * '<S2>'   : 'control_3dof/Cable Controller1'
 * '<S3>'   : 'control_3dof/Cable Controller2'
 * '<S4>'   : 'control_3dof/Payload Controller'
 * '<S5>'   : 'control_3dof/Cable Controller/Parallel Control'
 * '<S6>'   : 'control_3dof/Cable Controller/Vertical Control'
 * '<S7>'   : 'control_3dof/Cable Controller/Parallel Control/Parallel Control'
 * '<S8>'   : 'control_3dof/Cable Controller/Vertical Control/Vertical Control'
 * '<S9>'   : 'control_3dof/Cable Controller1/Parallel Control'
 * '<S10>'  : 'control_3dof/Cable Controller1/Vertical Control'
 * '<S11>'  : 'control_3dof/Cable Controller1/Parallel Control/Parallel Control'
 * '<S12>'  : 'control_3dof/Cable Controller1/Vertical Control/Vertical Control'
 * '<S13>'  : 'control_3dof/Cable Controller2/Parallel Control'
 * '<S14>'  : 'control_3dof/Cable Controller2/Vertical Control'
 * '<S15>'  : 'control_3dof/Cable Controller2/Parallel Control/Parallel Control'
 * '<S16>'  : 'control_3dof/Cable Controller2/Vertical Control/Vertical Control'
 * '<S17>'  : 'control_3dof/Payload Controller/Force Saturation && Disturbution'
 * '<S18>'  : 'control_3dof/Payload Controller/Load Control Force'
 * '<S19>'  : 'control_3dof/Payload Controller/State Bus Creater'
 * '<S20>'  : 'control_3dof/Payload Controller/Update Cable Direction'
 * '<S21>'  : 'control_3dof/Payload Controller/State Bus Creater/timestamp'
 * '<S22>'  : 'control_3dof/Payload Controller/Update Cable Direction/Detect Change'
 */
#endif                                 /* control_3dof_h_ */

/*
 * File trailer for generated code.
 *
 * [EOF]
 */
