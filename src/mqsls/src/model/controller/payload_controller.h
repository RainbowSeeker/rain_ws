//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: payload_controller.h
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
#ifndef payload_controller_h_
#define payload_controller_h_
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
};

#endif

#ifndef DEFINED_TYPEDEF_FOR_struct_yHOZp7GlS2mevFwObUORyE_
#define DEFINED_TYPEDEF_FOR_struct_yHOZp7GlS2mevFwObUORyE_

struct struct_yHOZp7GlS2mevFwObUORyE
{
  real_T MASS_LOAD;
  real_T MASS_UAV;
  real_T CABLE_LEN;
  real_T KP;
  real_T KV;
  real_T KQ;
  real_T KW;
};

#endif

//
//  Exported Global Parameters
//
//  Note: Exported global parameters are tunable parameters with an exported
//  global storage class designation.  Code generation will declare the memory for
//  these parameters and exports their symbols.
//

extern struct_yHOZp7GlS2mevFwObUORyE CONTROL_PARAM;// Variable: CONTROL_PARAM
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

// Class declaration for model payload_controller
class payload_controller final
{
  // public data and function members
 public:
  // Block signals for system '<S11>/Parallel Control'
  struct B_ParallelControl_payload_con_T {
    real_T rtp_mi[9];
    real_T a;
    real_T scale;
    real_T absxk;
    real_T t;
    real_T rtu_q;
    int32_T i;
  };

  // Block signals for system '<S12>/Vertical Control'
  struct B_VerticalControl_payload_con_T {
    real_T Sqi[9];
    real_T e_q[3];
    real_T e_w[3];
    real_T rtp_Kw[3];
  };

  // Block signals (default storage)
  struct B_payload_controller_T {
    real_T V[9];
    real_T U[9];
    real_T Q[9];                       // '<S23>/Vector Concatenate'
    real_T y[9];                       // '<S23>/Moore–Penrose inverse'
    real_T b_A[9];
    real_T Vf[9];
    real_T state_e[6];                 // '<S12>/Vertical Control'
    real_T aL_sp[3];                   // '<S24>/Add'
    real_T y_f[3];                     // '<S16>/Vertical Control'
    real_T MatrixMultiply[3];          // '<S23>/Matrix Multiply'
    real_T b_s[3];
    real_T e[3];
    real_T work[3];
    real_T dv[2];
    real_T absx;
    real_T force_sp3;
    real_T force_sp3_m;
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
    real_T cfromc_c;
    real_T ctoc_k;
    real_T cfrom1_c;
    real_T cto1_b;
    real_T mul_p;
    real_T roe;
    real_T absa;
    real_T absb;
    real_T scale;
    real_T ads;
    real_T bds;
    real_T scale_c;
    real_T absxk;
    real_T t;
    real_T scale_f;
    real_T absxk_g;
    real_T t_g;
    real_T temp;
    real_T temp_tmp;
    int32_T vcol;
    int32_T j;
    int32_T exponent;
    int32_T ar;
    int32_T b;
    int32_T ib;
    int32_T b_ic;
    int32_T c_k;
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
    int32_T j_m;
    int32_T b_i;
    int32_T scalarLB;
    int32_T vectorUB_n;
    int32_T i;
    int32_T offset_p;
    int32_T j_l;
    int32_T b_i_j;
    int32_T scalarLB_d;
    int32_T vectorUB_g;
    int32_T i1;
    int32_T k;
    int32_T scalarLB_l;
    int32_T vectorUB_d;
    int32_T i2;
    int32_T k_d;
    int32_T scalarLB_lx;
    int32_T vectorUB_o;
    int32_T i3;
    int32_T kend;
    int32_T k_b;
    int32_T kend_n;
    int32_T k_bs;
    boolean_T p;
    boolean_T doscale;
    boolean_T apply_transform;
    boolean_T notdone;
    boolean_T notdone_l;
    B_VerticalControl_payload_con_T sf_VerticalControl_i;// '<S20>/Vertical Control' 
    B_ParallelControl_payload_con_T sf_ParallelControl_k;// '<S19>/Parallel Control' 
    B_VerticalControl_payload_con_T sf_VerticalControl_o;// '<S16>/Vertical Control' 
    B_ParallelControl_payload_con_T sf_ParallelControl_g;// '<S15>/Parallel Control' 
    B_VerticalControl_payload_con_T sf_VerticalControl;// '<S12>/Vertical Control' 
    B_ParallelControl_payload_con_T sf_ParallelControl;// '<S11>/Parallel Control' 
  };

  // Invariant block signals (default storage)
  struct ConstB_payload_controller_T {
    real_T g[3];                       // '<S24>/-g'
  };

  // External inputs (root inport signals with default storage)
  struct ExtU_payload_controller_T {
    Payload_Out_Bus Payload_Out1;      // '<Root>/Payload_Out1'
  };

  // External outputs (root outports fed by signals with default storage)
  struct ExtY_payload_controller_T {
    real_T force_sp1[3];               // '<Root>/force_sp1'
    real_T force_sp2[3];               // '<Root>/force_sp2'
    real_T force_sp3[3];               // '<Root>/force_sp3'
  };

  // Real-time Model Data Structure
  struct RT_MODEL_payload_controller_T {
    const char_T * volatile errorStatus;
    const char_T* getErrorStatus() const;
    void setErrorStatus(const char_T* const volatile aErrorStatus);
  };

  // Copy Constructor
  payload_controller(payload_controller const&) = delete;

  // Assignment Operator
  payload_controller& operator= (payload_controller const&) & = delete;

  // Move Constructor
  payload_controller(payload_controller &&) = delete;

  // Move Assignment Operator
  payload_controller& operator= (payload_controller &&) = delete;

  // Real-Time Model get method
  payload_controller::RT_MODEL_payload_controller_T * getRTM();

  // Root inports set method
  void setExternalInputs(const ExtU_payload_controller_T
    *pExtU_payload_controller_T)
  {
    payload_controller_U = *pExtU_payload_controller_T;
  }

  // Root outports get method
  const ExtY_payload_controller_T &getExternalOutputs() const
  {
    return payload_controller_Y;
  }

  // model initialize function
  static void initialize();

  // model step function
  void step();

  // model terminate function
  static void terminate();

  // Constructor
  payload_controller();

  // Destructor
  ~payload_controller();

  // private data and function members
 private:
  // External inputs
  ExtU_payload_controller_T payload_controller_U;

  // External outputs
  ExtY_payload_controller_T payload_controller_Y;

  // Block signals
  B_payload_controller_T payload_controller_B;

  // private member function(s) for subsystem '<S1>/q1'
  static void payload_controller_q1(real_T rtu_theta, real_T rtu_psi, real_T
    rty_q[3]);

  // private member function(s) for subsystem '<S11>/Parallel Control'
  static void payload_control_ParallelControl(const real_T rtu_q[3], const
    real_T rtu_w[3], real_T rtu_vforce, const real_T rtu_aL[3], real_T rty_y[3],
    real_T rtp_li, real_T rtp_mi, B_ParallelControl_payload_con_T *localB);

  // private member function(s) for subsystem '<S12>/Vertical Control'
  static void payload_control_VerticalControl(const real_T rtu_q_sp[3], const
    real_T rtu_q[3], const real_T rtu_w[3], const real_T rtu_aL_sp[3], real_T
    rty_y[3], real_T rty_state[6], real_T rtp_Kq, real_T rtp_Kw, real_T rtp_li,
    real_T rtp_mi, B_VerticalControl_payload_con_T *localB);

  // private member function(s) for subsystem '<Root>'
  real_T payload_controller_xzlangeM(const real_T x[9]);
  void payload_controller_xzlascl(real_T cfrom, real_T cto, int32_T m, int32_T n,
    real_T A[9], int32_T iA0, int32_T lda);
  real_T payload_controller_xnrm2(int32_T n, const real_T x[9], int32_T ix0);
  real_T payload_controller_xdotc(int32_T n, const real_T x[9], int32_T ix0,
    const real_T y[9], int32_T iy0);
  void payload_controller_xaxpy(int32_T n, real_T a, int32_T ix0, real_T y[9],
    int32_T iy0);
  real_T payload_controller_xnrm2_d(int32_T n, const real_T x[3], int32_T ix0);
  void payload_controller_xaxpy_g(int32_T n, real_T a, const real_T x[9],
    int32_T ix0, real_T y[3], int32_T iy0);
  void payload_controller_xaxpy_g3(int32_T n, real_T a, const real_T x[3],
    int32_T ix0, real_T y[9], int32_T iy0);
  void payload_controller_xzlascl_a(real_T cfrom, real_T cto, int32_T m, int32_T
    n, real_T A[3], int32_T iA0, int32_T lda);
  void payload_controller_xswap(real_T x[9], int32_T ix0, int32_T iy0);
  void payload_controller_xrotg(real_T *a, real_T *b, real_T *c, real_T *s);
  void payload_controller_xrot(real_T x[9], int32_T ix0, int32_T iy0, real_T c,
    real_T s);
  void payload_controller_svd(const real_T A[9], real_T U[9], real_T s[3],
    real_T V[9]);

  // Real-Time Model
  RT_MODEL_payload_controller_T payload_controller_M;
};

extern const payload_controller::ConstB_payload_controller_T
  payload_controller_ConstB;           // constant block i/o

//-
//  These blocks were eliminated from the model due to optimizations:
//
//  Block '<S5>/Constant' : Unused code path elimination
//  Block '<S5>/Data Type Conversion' : Eliminate redundant data type conversion
//  Block '<S5>/Data Type Conversion1' : Eliminate redundant data type conversion
//  Block '<S5>/Data Type Conversion2' : Eliminate redundant data type conversion
//  Block '<S5>/Data Type Conversion3' : Eliminate redundant data type conversion
//  Block '<S5>/Data Type Conversion4' : Eliminate redundant data type conversion
//  Block '<S5>/Data Type Conversion5' : Eliminate redundant data type conversion


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
//  '<Root>' : 'payload_controller'
//  '<S1>'   : 'payload_controller/Cable Configuration Planning'
//  '<S2>'   : 'payload_controller/Cable Controller'
//  '<S3>'   : 'payload_controller/Cable Controller1'
//  '<S4>'   : 'payload_controller/Cable Controller2'
//  '<S5>'   : 'payload_controller/Control Bus Creater'
//  '<S6>'   : 'payload_controller/Payload Controller'
//  '<S7>'   : 'payload_controller/Trajectory Generator'
//  '<S8>'   : 'payload_controller/Cable Configuration Planning/q1'
//  '<S9>'   : 'payload_controller/Cable Configuration Planning/q2'
//  '<S10>'  : 'payload_controller/Cable Configuration Planning/q3'
//  '<S11>'  : 'payload_controller/Cable Controller/Parallel Control'
//  '<S12>'  : 'payload_controller/Cable Controller/Vertical Control'
//  '<S13>'  : 'payload_controller/Cable Controller/Parallel Control/Parallel Control'
//  '<S14>'  : 'payload_controller/Cable Controller/Vertical Control/Vertical Control'
//  '<S15>'  : 'payload_controller/Cable Controller1/Parallel Control'
//  '<S16>'  : 'payload_controller/Cable Controller1/Vertical Control'
//  '<S17>'  : 'payload_controller/Cable Controller1/Parallel Control/Parallel Control'
//  '<S18>'  : 'payload_controller/Cable Controller1/Vertical Control/Vertical Control'
//  '<S19>'  : 'payload_controller/Cable Controller2/Parallel Control'
//  '<S20>'  : 'payload_controller/Cable Controller2/Vertical Control'
//  '<S21>'  : 'payload_controller/Cable Controller2/Parallel Control/Parallel Control'
//  '<S22>'  : 'payload_controller/Cable Controller2/Vertical Control/Vertical Control'
//  '<S23>'  : 'payload_controller/Payload Controller/Force Distribution'
//  '<S24>'  : 'payload_controller/Payload Controller/Payload Position Controller'
//  '<S25>'  : 'payload_controller/Payload Controller/Force Distribution/Moore–Penrose inverse'

#endif                                 // payload_controller_h_

//
// File trailer for generated code.
//
// [EOF]
//
