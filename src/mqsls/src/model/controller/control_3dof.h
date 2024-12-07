//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: control_3dof.h
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
#ifndef control_3dof_h_
#define control_3dof_h_
#include <cmath>
#include "rtwtypes.h"
#ifndef DEFINED_TYPEDEF_FOR_Payload_Out_Bus_
#define DEFINED_TYPEDEF_FOR_Payload_Out_Bus_

struct Payload_Out_Bus
{
  uint32_T timestamp;
  real_T pL[3];
  real_T vL[3];
  real_T q_1[3];
  real_T q_2[3];
  real_T q_3[3];
  real_T w_1[3];
  real_T w_2[3];
  real_T w_3[3];
  real_T euler[3];
  real_T omega[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_Control_In_Bus_
#define DEFINED_TYPEDEF_FOR_Control_In_Bus_

struct Control_In_Bus
{
  uint32_T timestamp;

  // Used by 6dof
  real_T vforce[3];

  // Used by 3dof
  real_T q_sp[3];

  // unit vector for direction of cable i
  real_T q[3];
  real_T w[3];
  real_T ai_sp[3];
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_KOZZ2f41NjtyX9PsCFQIaF_
#define DEFINED_TYPEDEF_FOR_struct_KOZZ2f41NjtyX9PsCFQIaF_

struct struct_KOZZ2f41NjtyX9PsCFQIaF
{
  real_T MASS_LOAD;
  real_T MASS_UAV;
  real_T CABLE_LEN;
  real_T KP;
  real_T KV;
  real_T KQ;
  real_T KW;
  real_T RHO[9];
  real_T KR;
  real_T KOMEGA;
};

#endif

//
//  Exported Global Parameters
//
//  Note: Exported global parameters are tunable parameters with an exported
//  global storage class designation.  Code generation will declare the memory for
//  these parameters and exports their symbols.
//

extern struct_KOZZ2f41NjtyX9PsCFQIaF CONTROL_PARAM;// Variable: CONTROL_PARAM
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
    real_T rtp_mi[9];
    real_T e_q[3];
    real_T w_sp[3];
    real_T e_w[3];
    real_T rtp_Kw[3];
    real_T a;
    real_T e_w_c;
    real_T rtp_Kw_k;
    real_T rtu_ai_sp;
    real_T rtu_ai_sp_c;
    real_T rtu_ai_sp_b;
    real_T e_w_tmp;
    int32_T i;
    int32_T i_p;
    int32_T e_w_tmp_tmp;
  };

  // Block signals (default storage)
  struct B_control_3dof_T {
    real_T state_e[12];                // '<S7>/Vertical Control'
    real_T a[9];
    real_T A[9];
    real_T V[9];
    real_T U[9];
    real_T y[9];                       // '<S4>/Q sp'
    real_T b_A[9];
    real_T Vf[9];
    real_T vu[3];
    real_T aL_sp[3];                   // '<S19>/g'
    real_T rtb_control_in1_vforce[3];
    real_T rtb_control_in2_vforce[3];
    real_T y_n[3];                     // '<S7>/Vertical Control'
    real_T b_s[3];
    real_T e[3];
    real_T work[3];
    real_T dv[2];
    real_T absx;
    real_T rtb_aL_sp_m;
    real_T rtb_aL_sp_c;
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
    real_T cfromc_k;
    real_T ctoc_c;
    real_T cfrom1_b;
    real_T cto1_p;
    real_T mul_c;
    real_T roe;
    real_T absa;
    real_T absb;
    real_T scale;
    real_T ads;
    real_T bds;
    real_T scale_f;
    real_T absxk;
    real_T t;
    real_T scale_g;
    real_T absxk_g;
    real_T t_m;
    real_T temp;
    real_T temp_tmp;
    int32_T r_n;
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
    int32_T j_p;
    int32_T b_i;
    int32_T scalarLB;
    int32_T vectorUB_l;
    int32_T i;
    int32_T offset_j;
    int32_T j_d;
    int32_T b_i_g;
    int32_T scalarLB_l;
    int32_T vectorUB_d;
    int32_T i1;
    int32_T k;
    int32_T scalarLB_d;
    int32_T vectorUB_lx;
    int32_T i2;
    int32_T k_o;
    int32_T scalarLB_b;
    int32_T vectorUB_n;
    int32_T i3;
    int32_T kend;
    int32_T k_b;
    int32_T kend_l;
    int32_T k_h;
    boolean_T p;
    boolean_T doscale;
    boolean_T apply_transform;
    boolean_T notdone;
    boolean_T notdone_b;
    B_VerticalControl_control_3do_T sf_VerticalControl_i;// '<S15>/Vertical Control' 
    B_ParallelControl_control_3do_T sf_ParallelControl_k;// '<S14>/Parallel Control' 
    B_VerticalControl_control_3do_T sf_VerticalControl_o;// '<S11>/Vertical Control' 
    B_ParallelControl_control_3do_T sf_ParallelControl_g;// '<S10>/Parallel Control' 
    B_VerticalControl_control_3do_T sf_VerticalControl;// '<S7>/Vertical Control' 
    B_ParallelControl_control_3do_T sf_ParallelControl;// '<S6>/Parallel Control' 
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_control_3dof_T {
    Payload_Out_Bus Payload_Out;       // '<Root>/Payload_Out'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_control_3dof_T {
    real_T force_sp1[3];               // '<Root>/force_sp1'
    real_T force_sp2[3];               // '<Root>/force_sp2'
    real_T force_sp3[3];               // '<Root>/force_sp3'
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
  static void initialize();

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

  // private member function(s) for subsystem '<S6>/Parallel Control'
  static void control_3dof_ParallelControl(const real_T rtu_q[3], const real_T
    rtu_w[3], const real_T rtu_vforce[3], const real_T rtu_ai_sp[3], real_T
    rty_y[3], real_T rtp_li, real_T rtp_mi, B_ParallelControl_control_3do_T
    *localB);

  // private member function(s) for subsystem '<S7>/Vertical Control'
  static void control_3dof_VerticalControl(const real_T rtu_q_sp[3], const
    real_T rtu_q[3], const real_T rtu_w[3], const real_T rtu_ai_sp[3], real_T
    rty_y[3], real_T rty_state[12], real_T rtp_Kq, real_T rtp_Kw, real_T rtp_li,
    real_T rtp_mi, B_VerticalControl_control_3do_T *localB);

  // private member function(s) for subsystem '<Root>'
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

  // Real-Time Model
  RT_MODEL_control_3dof_T control_3dof_M;
};

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S19>/Saturation' : Eliminated Saturate block


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
//  '<S5>'   : 'control_3dof/Trajectory Generator'
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
//  '<S18>'  : 'control_3dof/Payload Controller/Force Disturbution'
//  '<S19>'  : 'control_3dof/Payload Controller/Load Control Force'
//  '<S20>'  : 'control_3dof/Payload Controller/Q sp'

#endif                                 // control_3dof_h_

//
// File trailer for generated code.
//
// [EOF]
//
