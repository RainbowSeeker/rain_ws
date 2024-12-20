//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: control_3dof.h
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
#ifndef control_3dof_h_
#define control_3dof_h_
#include <cmath>
#include "rtwtypes.h"
#ifndef DEFINED_TYPEDEF_FOR_Payload_Out_Bus_
#define DEFINED_TYPEDEF_FOR_Payload_Out_Bus_

struct Payload_Out_Bus
{
  uint64_T timestamp;
  real_T pL[3];
  real_T vL[3];
  real_T p_1[3];
  real_T p_2[3];
  real_T p_3[3];
  real_T v_1[3];
  real_T v_2[3];
  real_T v_3[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_TrajCmd_In_Bus_
#define DEFINED_TYPEDEF_FOR_TrajCmd_In_Bus_

struct TrajCmd_In_Bus
{
  uint64_T timestamp;
  real_T pos_sp[3];
  real_T vel_sp[3];
  real_T acc_ff[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_DirCmd_In_Bus_
#define DEFINED_TYPEDEF_FOR_DirCmd_In_Bus_

// 3x Cable Direction Cmd.
struct DirCmd_In_Bus
{
  uint64_T timestamp;
  real_T q_sp1[3];
  real_T q_sp2[3];
  real_T q_sp3[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_Control_In_Bus_
#define DEFINED_TYPEDEF_FOR_Control_In_Bus_

struct Control_In_Bus
{
  uint64_T timestamp;

  // Used by 6dof
  real_T vforce[3];

  // Used by 3dof
  real_T q_sp[3];

  // unit vector for direction of cable i
  real_T q[3];
  real_T w[3];
  real_T ai_sp[3];
  real_T pi[3];
  real_T vi[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_Control_State_Bus_
#define DEFINED_TYPEDEF_FOR_Control_State_Bus_

struct Control_State_Bus
{
  uint64_T timestamp;

  // Payload Disturbance Force
  real_T dL[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_VO89MlSMmM5pi10TnqAIZ_
#define DEFINED_TYPEDEF_FOR_struct_VO89MlSMmM5pi10TnqAIZ_

struct struct_VO89MlSMmM5pi10TnqAIZ
{
  real_T MASS_LOAD;
  real_T MASS_UAV;
  real_T CABLE_LEN;
  real_T TENSION_MIN;
  real_T TENSION_MAX;
  real_T ESO_PL[3];
  real_T ESO_VL[3];
  real_T KP;
  real_T KV;
  real_T KQ;
  real_T KW;
  real_T KQI;
};

#endif

#ifndef struct_pos_2nd_eso_control_3dof_T
#define struct_pos_2nd_eso_control_3dof_T

struct pos_2nd_eso_control_3dof_T
{
  int32_T isInitialized;
  real_T K_p[3];
  real_T K_v[3];
  real_T Z_1[3];
  real_T Z_2[3];
  real_T Z_3[3];
  boolean_T is_init;
};

#endif                                 // struct_pos_2nd_eso_control_3dof_T

//
//  Exported Global Parameters
//
//  Note: Exported global parameters are tunable parameters with an exported
//  global storage class designation.  Code generation will declare the memory for
//  these parameters and exports their symbols.
//

extern struct_VO89MlSMmM5pi10TnqAIZ CONTROL_PARAM;// Variable: CONTROL_PARAM
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
  static real_T rtGetNaN(void);
  static real32_T rtGetNaNF(void);
}                                      // extern "C"

extern "C"
{
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
}                                      // extern "C"

// Class declaration for model control_3dof
class control_3dof final
{
  // public data and function members
 public:
  // Block signals for system '<S6>/Parallel Control'
  struct B_ParallelControl_control_3do_T {
    real_T rtu_q[9];
    real_T rtp_mi[9];
    real_T a;
    real_T scale;
    real_T absxk;
    real_T t;
    real_T rtu_ai_sp;
    real_T rtu_ai_sp_m;
    real_T rtu_ai_sp_c;
    int32_T i;
    int32_T rtu_q_tmp;
  };

  // Block signals for system '<S7>/Vertical Control'
  struct B_VerticalControl_control_3do_T {
    real_T Sqi[9];
    real_T Sqi_m[9];
    real_T e_q[3];
    real_T e_w[3];
    real_T w_sp_tmp[3];
    real_T e_w_tmp[3];
    real_T b_a[3];
    real_T w_sp_tmp_c[3];
    real_T dv[2];
    real_T a;
    real_T rtu_q;
    real_T e_q_idx_1;
    real_T e_q_idx_2;
    real_T Sqi_k;
    real_T Sqi_c;
    real_T Sqi_b;
    real_T w_sp_tmp_p;
    real_T w_sp_tmp_cv;
    real_T w_sp_tmp_f;
    real_T rtu_ai_sp;
    real_T rtu_ai_sp_g;
    real_T rtu_ai_sp_g1;
    real_T e_w_tmp_tmp;
    int32_T i;
    int32_T i_m;
  };

  // Block signals (default storage)
  struct B_control_3dof_T {
    real_T state_f[12];                // '<S7>/Vertical Control'
    real_T Q_new[9];                   // '<S20>/Vector Concatenate'
    real_T rtb_vu_tmp[9];
    real_T Q_new_m[9];
    real_T dv[9];
    real_T V[9];
    real_T U[9];
    real_T b_A[9];
    real_T Vf[9];
    real_T q_2[3];
    real_T q_3[3];
    real_T F_trim[3];
    real_T F_way[3];
    real_T e_1[3];
    real_T e_2[3];
    real_T acc_limit[3];               // '<S19>/acc_limit'
    real_T rtb_control_in2_w[3];
    real_T rtb_control_in3_w[3];
    real_T vu[3];                    // '<S4>/Force Saturation && Disturbution'
    real_T f_vertical[3];              // '<S7>/Vertical Control'
    real_T dot_err_g[3];               // '<S7>/Vertical Control'
    real_T dot_err_c[3];               // '<S11>/Vertical Control'
    real_T dot_err[3];                 // '<S15>/Vertical Control'
    real_T vu_c[3];
    real_T s[3];
    real_T b_s[3];
    real_T e[3];
    real_T work[3];
    real_T dv1[2];
    real_T dv2[2];
    real_T dv3[2];
    real_T q_3_k;
    real_T F_trim_c;
    real_T e_2_b;
    real_T u0;
    real_T u0_p;
    real_T UnitDelay_DSTATE;
    real_T UnitDelay_DSTATE_c;
    real_T UnitDelay_DSTATE_f;
    real_T rtb_vu_tmp_g;
    real_T UnitDelay_DSTATE_tmp;
    real_T UnitDelay_DSTATE_tmp_g;
    real_T UnitDelay_DSTATE_tmp_m;
    real_T result;
    real_T F_d_idx_0;
    real_T F_d_idx_1;
    real_T F_d_idx_2;
    real_T absx;
    real_T cscale;
    real_T anrm;
    real_T nrm;
    real_T rt;
    real_T r;
    real_T smm1;
    real_T c;
    real_T shift;
    real_T cfromc;
    real_T ctoc;
    real_T cfrom1;
    real_T cto1;
    real_T mul;
    real_T cfromc_n;
    real_T ctoc_p;
    real_T cfrom1_l;
    real_T cto1_j;
    real_T mul_d;
    real_T roe;
    real_T absa;
    real_T absb;
    real_T scale;
    real_T ads;
    real_T bds;
    real_T scale_g;
    real_T absxk;
    real_T t;
    real_T scale_l;
    real_T absxk_d;
    real_T t_d;
    real_T temp;
    real_T temp_tmp;
    int32_T i;
    int32_T i_l;
    int32_T b_i;
    int32_T r_o;
    int32_T vcol;
    int32_T j;
    int32_T exponent;
    int32_T ar;
    int32_T b;
    int32_T ib;
    int32_T b_ic;
    int32_T vectorUB;
    int32_T qp1;
    int32_T qq;
    int32_T qjj;
    int32_T m;
    int32_T kase;
    int32_T d;
    int32_T e_k;
    int32_T qq_tmp;
    int32_T offset;
    int32_T j_b;
    int32_T b_i_n;
    int32_T scalarLB;
    int32_T vectorUB_b;
    int32_T i_ln;
    int32_T offset_h;
    int32_T j_bn;
    int32_T b_i_d;
    int32_T scalarLB_e;
    int32_T vectorUB_bj;
    int32_T i1;
    int32_T k;
    int32_T scalarLB_j;
    int32_T vectorUB_f;
    int32_T i2;
    int32_T k_a;
    int32_T scalarLB_ju;
    int32_T vectorUB_j;
    int32_T i3;
    int32_T kend;
    int32_T k_o;
    int32_T kend_n;
    int32_T k_i;
    boolean_T x[3];
    boolean_T found[3];
    boolean_T FixPtRelationalOperator; // '<S21>/FixPt Relational Operator'
    boolean_T is_valid;
    boolean_T p;
    boolean_T doscale;
    boolean_T apply_transform;
    boolean_T notdone;
    boolean_T notdone_o;
    B_VerticalControl_control_3do_T sf_VerticalControl_h;// '<S15>/Vertical Control' 
    B_ParallelControl_control_3do_T sf_ParallelControl_p;// '<S14>/Parallel Control' 
    B_VerticalControl_control_3do_T sf_VerticalControl_m;// '<S11>/Vertical Control' 
    B_ParallelControl_control_3do_T sf_ParallelControl_k;// '<S10>/Parallel Control' 
    B_VerticalControl_control_3do_T sf_VerticalControl;// '<S7>/Vertical Control' 
    B_ParallelControl_control_3do_T sf_ParallelControl;// '<S6>/Parallel Control' 
  };

  // Block states (default storage) for system '<Root>'
  struct DW_control_3dof_T {
    pos_2nd_eso_control_3dof_T obj;    // '<S4>/Position 2nd ESO'
    real_T Delay_DSTATE[9];            // '<S20>/Delay'
    real_T UnitDelay_DSTATE[3];        // '<S4>/Unit Delay'
    real_T DiscreteTimeIntegrator_DSTATE[3];// '<S15>/Discrete-Time Integrator'
    real_T DiscreteTimeIntegrator_DSTATE_l[3];// '<S11>/Discrete-Time Integrator' 
    real_T DiscreteTimeIntegrator_DSTATE_p[3];// '<S7>/Discrete-Time Integrator' 
    uint64_T DelayInput1_DSTATE;       // '<S21>/Delay Input1'
    uint64_T DiscreteTimeIntegrator_DSTAT_ps;// '<S22>/Discrete-Time Integrator' 
    boolean_T icLoad;                  // '<S20>/Delay'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_control_3dof_T {
    Payload_Out_Bus Payload_Out;       // '<Root>/Payload_Out'
    TrajCmd_In_Bus Traj_sp;            // '<Root>/Traj_sp'
    DirCmd_In_Bus Dir_sp;              // '<Root>/Dir_sp'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_control_3dof_T {
    real_T force_sp1[3];               // '<Root>/force_sp1'
    real_T force_sp2[3];               // '<Root>/force_sp2'
    real_T force_sp3[3];               // '<Root>/force_sp3'
    Control_State_Bus state;           // '<Root>/state'
  };

  // Real-time Model Data Structure
  struct RT_MODEL_control_3dof_T {
    const char_T * volatile errorStatus;
    const char_T* getErrorStatus() const;
    void setErrorStatus(const char_T* const volatile aErrorStatus);
  };

  // Copy Constructor
  control_3dof(control_3dof const&) = delete;

  // Assignment Operator
  control_3dof& operator= (control_3dof const&) & = delete;

  // Move Constructor
  control_3dof(control_3dof &&) = delete;

  // Move Assignment Operator
  control_3dof& operator= (control_3dof &&) = delete;

  // Real-Time Model get method
  control_3dof::RT_MODEL_control_3dof_T * getRTM();

  // Root inports set method
  void setExternalInputs(const ExtU_control_3dof_T *pExtU_control_3dof_T)
  {
    control_3dof_U = *pExtU_control_3dof_T;
  }

  // Root outports get method
  const ExtY_control_3dof_T &getExternalOutputs() const
  {
    return control_3dof_Y;
  }

  // model initialize function
  void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  control_3dof();

  // Destructor
  ~control_3dof();

  // private data and function members
 private:
  // External inputs
  ExtU_control_3dof_T control_3dof_U;

  // External outputs
  ExtY_control_3dof_T control_3dof_Y;

  // Block signals
  B_control_3dof_T control_3dof_B;

  // Block states
  DW_control_3dof_T control_3dof_DW;

  // private member function(s) for subsystem '<S6>/Parallel Control'
  static void control_3dof_ParallelControl(const real_T rtu_q[3], const real_T
    rtu_w[3], const real_T rtu_vforce[3], const real_T rtu_ai_sp[3], real_T
    rty_f_parallel[3], real_T rtp_li, real_T rtp_mi,
    B_ParallelControl_control_3do_T *localB);

  // private member function(s) for subsystem '<S7>/Vertical Control'
  static void control_3dof_VerticalControl(const real_T rtu_q_sp[3], const
    real_T rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
    rty_f_vertical[3], real_T rty_dot_err[3], real_T rty_state[12], real_T
    rtp_Kq, real_T rtp_Kw, real_T rtp_li, real_T rtp_mi,
    B_VerticalControl_control_3do_T *localB);

  // private member function(s) for subsystem '<Root>'
  real_T control_3dof_norm(const real_T x[3]);
  real_T control_3dof_xzlangeM(const real_T x[9]);
  void control_3dof_xzlascl(real_T cfrom, real_T cto, int32_T m, int32_T n,
    real_T A[9], int32_T iA0, int32_T lda);
  real_T control_3dof_xnrm2(int32_T n, const real_T x[9], int32_T ix0);
  real_T control_3dof_xdotc(int32_T n, const real_T x[9], int32_T ix0, const
    real_T y[9], int32_T iy0);
  void control_3dof_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[9], int32_T
    iy0);
  real_T control_3dof_xnrm2_n(int32_T n, const real_T x[3], int32_T ix0);
  void control_3dof_xaxpy_f(int32_T n, real_T a, const real_T x[9], int32_T ix0,
    real_T y[3], int32_T iy0);
  void control_3dof_xaxpy_fh(int32_T n, real_T a, const real_T x[3], int32_T ix0,
    real_T y[9], int32_T iy0);
  void control_3dof_xzlascl_l(real_T cfrom, real_T cto, int32_T m, int32_T n,
    real_T A[3], int32_T iA0, int32_T lda);
  void control_3dof_xswap(real_T x[9], int32_T ix0, int32_T iy0);
  void control_3dof_xrotg(real_T *a, real_T *b, real_T *c, real_T *s);
  void control_3dof_xrot(real_T x[9], int32_T ix0, int32_T iy0, real_T c, real_T
    s);
  void control_3dof_svd(const real_T A[9], real_T U[9], real_T s[3], real_T V[9]);
  void control_3dof_pinv(const real_T A[9], real_T X[9]);
  void control_3dof_TCISolver(const real_T F_sp[3], const real_T F_fixed[3],
    const real_T Q[9], real_T b_min_tension, real_T b_max_tension, real_T
    result[3], real_T F_actual[3]);

  // Real-Time Model
  RT_MODEL_control_3dof_T control_3dof_M;
};

//-
//  The generated code includes comments that allow you to trace directly
//  back to the appropriate location in the model.  The basic format
//  is <system>/block_name, where system is the system number (uniquely
//  assigned by Simulink) and block_name is the name of the block.
//
//  Use the MATLAB hilite_system command to trace the generated code back
//  to the model.  For example,
//
//  hilite_system('<S3>')    - opens system 3
//  hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
//
//  Here is the system hierarchy for this model
//
//  '<Root>' : 'control_3dof'
//  '<S1>'   : 'control_3dof/Cable Controller'
//  '<S2>'   : 'control_3dof/Cable Controller1'
//  '<S3>'   : 'control_3dof/Cable Controller2'
//  '<S4>'   : 'control_3dof/Payload Controller'
//  '<S5>'   : 'control_3dof/State Bus Creater'
//  '<S6>'   : 'control_3dof/Cable Controller/Parallel Control'
//  '<S7>'   : 'control_3dof/Cable Controller/Vertical Control'
//  '<S8>'   : 'control_3dof/Cable Controller/Parallel Control/Parallel Control'
//  '<S9>'   : 'control_3dof/Cable Controller/Vertical Control/Vertical Control'
//  '<S10>'  : 'control_3dof/Cable Controller1/Parallel Control'
//  '<S11>'  : 'control_3dof/Cable Controller1/Vertical Control'
//  '<S12>'  : 'control_3dof/Cable Controller1/Parallel Control/Parallel Control'
//  '<S13>'  : 'control_3dof/Cable Controller1/Vertical Control/Vertical Control'
//  '<S14>'  : 'control_3dof/Cable Controller2/Parallel Control'
//  '<S15>'  : 'control_3dof/Cable Controller2/Vertical Control'
//  '<S16>'  : 'control_3dof/Cable Controller2/Parallel Control/Parallel Control'
//  '<S17>'  : 'control_3dof/Cable Controller2/Vertical Control/Vertical Control'
//  '<S18>'  : 'control_3dof/Payload Controller/Force Saturation && Disturbution'
//  '<S19>'  : 'control_3dof/Payload Controller/Load Control Force'
//  '<S20>'  : 'control_3dof/Payload Controller/Update Cable Direction'
//  '<S21>'  : 'control_3dof/Payload Controller/Update Cable Direction/Detect Change'
//  '<S22>'  : 'control_3dof/State Bus Creater/timestamp'

#endif                                 // control_3dof_h_

//
// File trailer for generated code.
//
// [EOF]
//
