/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: optimize_forces.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2024-12-19 20:42:18
 */

/* Include Files */
#include "optimize_forces.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  bool gradOK;
  bool fevalOK;
  bool done;
  bool stepAccepted;
  bool failedLineSearch;
  int stepType;
} struct_T;

#endif                                 /* typedef_struct_T */

#ifndef typedef_b_struct_T
#define typedef_b_struct_T

typedef struct {
  double center[3];
  double A_alpha[18];
  double b_alpha[6];
  double deltaT[3];
  double T_min[3];
} b_struct_T;

#endif                                 /* typedef_b_struct_T */

#ifndef typedef_anonymous_function
#define typedef_anonymous_function

typedef struct {
  b_struct_T workspace;
} anonymous_function;

#endif                                 /* typedef_anonymous_function */

#ifndef c_typedef_coder_internal_sticky
#define c_typedef_coder_internal_sticky

typedef struct {
  anonymous_function value;
} coder_internal_stickyStruct;

#endif                                 /* c_typedef_coder_internal_sticky */

#ifndef c_typedef_b_coder_internal_stic
#define c_typedef_b_coder_internal_stic

typedef struct {
  coder_internal_stickyStruct next;
} b_coder_internal_stickyStruct;

#endif                                 /* c_typedef_b_coder_internal_stic */

#ifndef c_typedef_c_coder_internal_stic
#define c_typedef_c_coder_internal_stic

typedef struct {
  b_coder_internal_stickyStruct next;
} c_coder_internal_stickyStruct;

#endif                                 /* c_typedef_c_coder_internal_stic */

#ifndef c_typedef_d_coder_internal_stic
#define c_typedef_d_coder_internal_stic

typedef struct {
  c_coder_internal_stickyStruct next;
} d_coder_internal_stickyStruct;

#endif                                 /* c_typedef_d_coder_internal_stic */

#ifndef c_typedef_e_coder_internal_stic
#define c_typedef_e_coder_internal_stic

typedef struct {
  d_coder_internal_stickyStruct next;
} e_coder_internal_stickyStruct;

#endif                                 /* c_typedef_e_coder_internal_stic */

#ifndef c_typedef_f_coder_internal_stic
#define c_typedef_f_coder_internal_stic

typedef struct {
  e_coder_internal_stickyStruct next;
} f_coder_internal_stickyStruct;

#endif                                 /* c_typedef_f_coder_internal_stic */

#ifndef c_typedef_g_coder_internal_stic
#define c_typedef_g_coder_internal_stic

typedef struct {
  f_coder_internal_stickyStruct next;
} g_coder_internal_stickyStruct;

#endif                                 /* c_typedef_g_coder_internal_stic */

#ifndef c_typedef_h_coder_internal_stic
#define c_typedef_h_coder_internal_stic

typedef struct {
  g_coder_internal_stickyStruct next;
} h_coder_internal_stickyStruct;

#endif                                 /* c_typedef_h_coder_internal_stic */

#ifndef c_typedef_i_coder_internal_stic
#define c_typedef_i_coder_internal_stic

typedef struct {
  h_coder_internal_stickyStruct next;
} i_coder_internal_stickyStruct;

#endif                                 /* c_typedef_i_coder_internal_stic */

#ifndef struct_emxArray_real_T_6
#define struct_emxArray_real_T_6

struct emxArray_real_T_6
{
  double data[6];
  int size[1];
};

#endif                                 /* struct_emxArray_real_T_6 */

#ifndef typedef_emxArray_real_T_6
#define typedef_emxArray_real_T_6

typedef struct emxArray_real_T_6 emxArray_real_T_6;

#endif                                 /* typedef_emxArray_real_T_6 */

#ifndef struct_emxArray_real_T_22
#define struct_emxArray_real_T_22

struct emxArray_real_T_22
{
  double data[22];
  int size[1];
};

#endif                                 /* struct_emxArray_real_T_22 */

#ifndef typedef_emxArray_real_T_22
#define typedef_emxArray_real_T_22

typedef struct emxArray_real_T_22 emxArray_real_T_22;

#endif                                 /* typedef_emxArray_real_T_22 */

#ifndef struct_emxArray_real_T_40
#define struct_emxArray_real_T_40

struct emxArray_real_T_40
{
  double data[40];
  int size[1];
};

#endif                                 /* struct_emxArray_real_T_40 */

#ifndef typedef_emxArray_real_T_40
#define typedef_emxArray_real_T_40

typedef struct emxArray_real_T_40 emxArray_real_T_40;

#endif                                 /* typedef_emxArray_real_T_40 */

#ifndef struct_emxArray_real_T_40x40
#define struct_emxArray_real_T_40x40

struct emxArray_real_T_40x40
{
  double data[1600];
  int size[2];
};

#endif                                 /* struct_emxArray_real_T_40x40 */

#ifndef typedef_emxArray_real_T_40x40
#define typedef_emxArray_real_T_40x40

typedef struct emxArray_real_T_40x40 emxArray_real_T_40x40;

#endif                                 /* typedef_emxArray_real_T_40x40 */

#ifndef struct_emxArray_real_T_40x22
#define struct_emxArray_real_T_40x22

struct emxArray_real_T_40x22
{
  double data[880];
  int size[2];
};

#endif                                 /* struct_emxArray_real_T_40x22 */

#ifndef typedef_emxArray_real_T_40x22
#define typedef_emxArray_real_T_40x22

typedef struct emxArray_real_T_40x22 emxArray_real_T_40x22;

#endif                                 /* typedef_emxArray_real_T_40x22 */

#ifndef struct_emxArray_real_T_132
#define struct_emxArray_real_T_132

struct emxArray_real_T_132
{
  double data[132];
  int size[1];
};

#endif                                 /* struct_emxArray_real_T_132 */

#ifndef typedef_emxArray_real_T_132
#define typedef_emxArray_real_T_132

typedef struct emxArray_real_T_132 emxArray_real_T_132;

#endif                                 /* typedef_emxArray_real_T_132 */

#ifndef struct_emxArray_int32_T_40
#define struct_emxArray_int32_T_40

struct emxArray_int32_T_40
{
  int data[40];
};

#endif                                 /* struct_emxArray_int32_T_40 */

#ifndef typedef_emxArray_int32_T_40
#define typedef_emxArray_int32_T_40

typedef struct emxArray_int32_T_40 emxArray_int32_T_40;

#endif                                 /* typedef_emxArray_int32_T_40 */

#ifndef struct_emxArray_real_T_22x3
#define struct_emxArray_real_T_22x3

struct emxArray_real_T_22x3
{
  double data[66];
};

#endif                                 /* struct_emxArray_real_T_22x3 */

#ifndef typedef_emxArray_real_T_22x3
#define typedef_emxArray_real_T_22x3

typedef struct emxArray_real_T_22x3 emxArray_real_T_22x3;

#endif                                 /* typedef_emxArray_real_T_22x3 */

#ifndef typedef_c_struct_T
#define typedef_c_struct_T

typedef struct {
  int mIneq;
  int iNonIneq0;
  double sqpFval;
  double sqpFval_old;
  double xstarsqp[9];
  double xstarsqp_old[9];
  emxArray_real_T_6 cIneq;
  emxArray_real_T_6 cIneq_old;
  double cEq[3];
  double cEq_old[3];
  emxArray_real_T_22 grad;
  emxArray_real_T_22 grad_old;
  int FunctionEvaluations;
  int sqpIterations;
  int sqpExitFlag;
  emxArray_real_T_40 lambdasqp;
  emxArray_real_T_40 lambdaStopTest;
  emxArray_real_T_40 lambdaStopTestPrev;
  double steplength;
  emxArray_real_T_22 delta_x;
  emxArray_real_T_22 socDirection;
  emxArray_int32_T_40 workingset_old;
  emxArray_real_T_22x3 JacCeqTrans_old;
  emxArray_real_T_22 gradLag;
  emxArray_real_T_22 delta_gradLag;
  emxArray_real_T_40 xstar;
  double fstar;
  emxArray_real_T_40 lambda;
  int state;
  double maxConstr;
  int iterations;
  emxArray_real_T_40 searchDir;
} c_struct_T;

#endif                                 /* typedef_c_struct_T */

#ifndef typedef_d_struct_T
#define typedef_d_struct_T

typedef struct {
  int ldq;
  emxArray_real_T_40x40 QR;
  emxArray_real_T_40x40 Q;
  emxArray_int32_T_40 jpvt;
  int mrows;
  int ncols;
  emxArray_real_T_40 tau;
  int minRowCol;
  bool usedPivoting;
} d_struct_T;

#endif                                 /* typedef_d_struct_T */

#ifndef typedef_e_struct_T
#define typedef_e_struct_T

typedef struct {
  emxArray_real_T_40x40 FMat;
  int ldm;
  int ndims;
  int info;
  bool ConvexCheck;
  double workspace_;
} e_struct_T;

#endif                                 /* typedef_e_struct_T */

#ifndef struct_emxArray_real_T_21
#define struct_emxArray_real_T_21

struct emxArray_real_T_21
{
  double data[21];
};

#endif                                 /* struct_emxArray_real_T_21 */

#ifndef typedef_emxArray_real_T_21
#define typedef_emxArray_real_T_21

typedef struct emxArray_real_T_21 emxArray_real_T_21;

#endif                                 /* typedef_emxArray_real_T_21 */

#ifndef typedef_f_struct_T
#define typedef_f_struct_T

typedef struct {
  emxArray_real_T_22 grad;
  emxArray_real_T_21 Hx;
  bool hasLinear;
  int nvar;
  int maxVar;
  double beta;
  double rho;
  int objtype;
  int prev_objtype;
  int prev_nvar;
  bool prev_hasLinear;
  double gammaScalar;
} f_struct_T;

#endif                                 /* typedef_f_struct_T */

#ifndef typedef_g_struct_T
#define typedef_g_struct_T

typedef struct {
  emxArray_real_T_40x22 workspace_float;
  emxArray_int32_T_40 workspace_int;
  emxArray_int32_T_40 workspace_sort;
} g_struct_T;

#endif                                 /* typedef_g_struct_T */

#ifndef struct_emxArray_real_T_66
#define struct_emxArray_real_T_66

struct emxArray_real_T_66
{
  double data[66];
};

#endif                                 /* struct_emxArray_real_T_66 */

#ifndef typedef_emxArray_real_T_66
#define typedef_emxArray_real_T_66

typedef struct emxArray_real_T_66 emxArray_real_T_66;

#endif                                 /* typedef_emxArray_real_T_66 */

#ifndef struct_emxArray_int32_T_22
#define struct_emxArray_int32_T_22

struct emxArray_int32_T_22
{
  int data[22];
};

#endif                                 /* struct_emxArray_int32_T_22 */

#ifndef typedef_emxArray_int32_T_22
#define typedef_emxArray_int32_T_22

typedef struct emxArray_int32_T_22 emxArray_int32_T_22;

#endif                                 /* typedef_emxArray_int32_T_22 */

#ifndef struct_emxArray_real_T_880
#define struct_emxArray_real_T_880

struct emxArray_real_T_880
{
  double data[880];
};

#endif                                 /* struct_emxArray_real_T_880 */

#ifndef typedef_emxArray_real_T_880
#define typedef_emxArray_real_T_880

typedef struct emxArray_real_T_880 emxArray_real_T_880;

#endif                                 /* typedef_emxArray_real_T_880 */

#ifndef struct_emxArray_boolean_T_40
#define struct_emxArray_boolean_T_40

struct emxArray_boolean_T_40
{
  bool data[40];
};

#endif                                 /* struct_emxArray_boolean_T_40 */

#ifndef typedef_emxArray_boolean_T_40
#define typedef_emxArray_boolean_T_40

typedef struct emxArray_boolean_T_40 emxArray_boolean_T_40;

#endif                                 /* typedef_emxArray_boolean_T_40 */

#ifndef typedef_h_struct_T
#define typedef_h_struct_T

typedef struct {
  int mConstr;
  int mConstrOrig;
  int mConstrMax;
  int nVar;
  int nVarMax;
  int ldA;
  emxArray_real_T_132 Aineq;
  emxArray_real_T_6 bineq;
  emxArray_real_T_66 Aeq;
  double beq[3];
  emxArray_real_T_22 lb;
  emxArray_real_T_22 ub;
  emxArray_int32_T_22 indexLB;
  emxArray_int32_T_22 indexUB;
  emxArray_int32_T_22 indexFixed;
  int mEqRemoved;
  int indexEqRemoved[3];
  emxArray_real_T_880 ATwset;
  emxArray_real_T_40 bwset;
  int nActiveConstr;
  emxArray_real_T_40 maxConstrWorkspace;
  int sizes[5];
  int sizesNormal[5];
  int sizesPhaseOne[5];
  int sizesRegularized[5];
  int sizesRegPhaseOne[5];
  int isActiveIdx[6];
  int isActiveIdxNormal[6];
  int isActiveIdxPhaseOne[6];
  int isActiveIdxRegularized[6];
  int isActiveIdxRegPhaseOne[6];
  emxArray_boolean_T_40 isActiveConstr;
  emxArray_int32_T_40 Wid;
  emxArray_int32_T_40 Wlocalidx;
  int nWConstr[5];
  int probType;
} h_struct_T;

#endif                                 /* typedef_h_struct_T */

#ifndef typedef_i_struct_T
#define typedef_i_struct_T

typedef struct {
  anonymous_function objfun;
  double f_1;
  double cEq_1[3];
  int numEvals;
  bool hasBounds;
} i_struct_T;

#endif                                 /* typedef_i_struct_T */

#ifndef typedef_j_struct_T
#define typedef_j_struct_T

typedef struct {
  double penaltyParam;
  double threshold;
  int nPenaltyDecreases;
  double linearizedConstrViol;
  double initFval;
  double initConstrViolationEq;
  double initConstrViolationIneq;
  double phi;
  double phiPrimePlus;
  double phiFullStep;
  double feasRelativeFactor;
  double nlpPrimalFeasError;
  double nlpDualFeasError;
  double nlpComplError;
  double firstOrderOpt;
} j_struct_T;

#endif                                 /* typedef_j_struct_T */

#ifndef typedef_k_struct_T
#define typedef_k_struct_T

typedef struct {
  char SolverName[7];
  int MaxIterations;
} k_struct_T;

#endif                                 /* typedef_k_struct_T */

/* Variable Definitions */
static const signed char iv[9] = { 1, 1, 0, 1, 1, 0, 1, 1, 0 };

/* Function Declarations */
static bool BFGSUpdate(int nvar, double Bk[81], const double sk_data[], double
  yk_data[], double workspace_data[]);
static void PresolveWorkingSet(c_struct_T *solution, g_struct_T *memspace,
  h_struct_T *workingset, d_struct_T *qrmanager);
static int RemoveDependentEq_(g_struct_T *memspace, h_struct_T *workingset,
  d_struct_T *qrmanager);
static void addAeqConstr(h_struct_T *obj, int idx_local);
static void addBoundToActiveSetMatrix_(h_struct_T *obj, int TYPE, int idx_local);
static bool all(const bool x[3]);
static void b_driver(const double H[81], const double f_data[], c_struct_T
                     *solution, g_struct_T *memspace, h_struct_T *workingset,
                     d_struct_T *qrmanager, e_struct_T *cholmanager, f_struct_T *
                     objective, const k_struct_T *options, const k_struct_T
                     *runTimeOptions);
static void b_factoryConstruct(int mIneqMax, int nVarMax, int mConstrMax,
  h_struct_T *obj);
static double b_maxConstraintViolation(h_struct_T *obj, const double x_data[]);
static void b_test_exit(struct_T *Flags, g_struct_T *memspace, j_struct_T
  *MeritFunction, int fscales_lineq_constraint_size, h_struct_T *WorkingSet,
  c_struct_T *TrialState, d_struct_T *QRManager);
static void b_updateWorkingSetForNewQP(const double xk[9], h_struct_T
  *WorkingSet, int mIneq, const double cIneq_data[], const double cEq[3]);
static void b_vecnorm(const double x[18], double y[6]);
static void b_xaxpy(double a, const double x[9], int ix0, double y[3]);
static void b_xgemm(int m, int n, int k, const double A_data[], int ia0, int lda,
                    const double B_data[], int ldb, double C_data[], int ldc);
static void b_xgemv(int m, int n, const double A_data[], int lda, const double
                    x_data[], int ix0, double y_data[]);
static double b_xnrm2(int n, const double x_data[], int ix0);
static void b_xzlascl(double cfrom, double cto, double A[3]);
static double c_maxConstraintViolation_AMats_(h_struct_T *obj, const double
  x_data[]);
static void c_xaxpy(double a, const double x[3], double y[9], int iy0);
static double c_xnrm2(int n, const double x_data[]);
static double computeComplError(int fscales_lineq_constraint_size, const double
  xCurrent[9], int mIneq, const double cIneq_data[], const int finiteLB_data[],
  int mLB, const int finiteUB_data[], int mUB, const double lambda_data[], int
  iL0);
static int computeConstraints_(const double x[9], const double
  Cineq_workspace_data[], int ineq0, double Ceq_workspace[3]);
static bool computeDualFeasError(int nVar, const double gradLag_data[], double
  *val);
static bool computeFiniteDifferences(i_struct_T *obj, double fCurrent, const
  double cEqCurrent[3], double xk[9], double gradf_data[], double
  JacCeqTrans_data[], int ldJE);
static double computeFval_ReuseHx(const f_struct_T *obj, double workspace_data[],
  const double f_data[], const double x_data[]);
static void computeGradLag(double workspace_data[], int ldA, int nVar, const
  double grad_data[], int mIneq, const double AineqTrans_data[], const double
  AeqTrans_data[], const int finiteFixed_data[], int mFixed, const int
  finiteLB_data[], int mLB, const int finiteUB_data[], int mUB, const double
  lambda_data[]);
static void computeGrad_StoreHx(f_struct_T *obj, const double H[81], const
  double f_data[], const double x_data[]);
static void computeLinearResiduals(const double x[9], int nVar, double
  workspaceIneq_data[], const int *workspaceIneq_size, int mLinIneq, const
  double AineqT_data[], const double bineq_data[], int ldAi);
static double computeMeritFcn(double obj_penaltyParam, double fval, const double
  Cineq_workspace_data[], int mIneq, const double Ceq_workspace[3], bool
  evalWellDefined);
static double computePrimalFeasError(const double x[9], int mLinIneq, const
  double cIneq_data[], const double cEq[3], const int finiteLB_data[], int mLB,
  const int finiteUB_data[], int mUB);
static void computeQ_(d_struct_T *obj, int nrows);
static void compute_deltax(const double H[81], c_struct_T *solution, g_struct_T *
  memspace, const d_struct_T *qrmanager, e_struct_T *cholmanager, const
  f_struct_T *objective, bool alwaysPositiveDef);
static void compute_lambda(double workspace_data[], c_struct_T *solution, const
  f_struct_T *objective, const d_struct_T *qrmanager);
static void countsort(int x_data[], int xLen, int workspace_data[], int xMin,
                      int xMax);
static double d_maxConstraintViolation_AMats_(h_struct_T *obj, const double
  x_data[]);
static double d_xnrm2(const double x[3]);
static void deleteColMoveEnd(d_struct_T *obj, int idx);
static void diag(const double v[3], double d[9]);
static void driver(const double bineq_data[], c_struct_T *TrialState, j_struct_T
                   *MeritFunction, const i_coder_internal_stickyStruct
                   *FcnEvaluator, i_struct_T *FiniteDifferences, g_struct_T
                   *memspace, h_struct_T *WorkingSet, d_struct_T *QRManager,
                   e_struct_T *CholManager, f_struct_T *QPObjective, int
                   fscales_lineq_constraint_size, double Hessian[81]);
static double evalObjAndConstr(const b_struct_T *c_obj_next_next_next_next_next_,
  const double x[9], const double Cineq_workspace_data[], int ineq0, double
  Ceq_workspace[3], int *status);
static void factorQR(d_struct_T *obj, const double A_data[], int mrows, int
                     ncols, int ldA);
static void factorQRE(d_struct_T *obj, const double A_data[], int mrows, int
                      ncols, int ldA);
static void factoryConstruct(int nVarMax, int mConstrMax, int mIneq, c_struct_T *
  obj);
static bool feasibleX0ForWorkingSet(double workspace_data[], const int
  workspace_size[2], double xCurrent_data[], h_struct_T *workingset, d_struct_T *
  qrmanager);
static double feasibleratiotest(const double solution_xstar_data[], const double
  solution_searchDir_data[], double workspace_data[], const int workspace_size[2],
  int workingset_nVar, int workingset_ldA, const double workingset_Aineq_data[],
  const double workingset_bineq_data[], const double workingset_lb_data[], const
  double workingset_ub_data[], const int workingset_indexLB_data[], const int
  workingset_indexUB_data[], const int workingset_sizes[5], const int
  workingset_isActiveIdx[6], const bool workingset_isActiveConstr_data[], const
  int workingset_nWConstr[5], bool isPhaseOne, bool *newBlocking, int
  *constrType, int *constrIdx);
static bool finDiffEvalAndChkErr(const double obj_objfun_workspace_center[3],
  const double obj_objfun_workspace_A_alpha[18], const double
  obj_objfun_workspace_b_alpha[6], const double obj_objfun_workspace_deltaT[3],
  const double obj_objfun_workspace_T_min[3], double cEqPlus[3], int dim, double
  delta, double xk[9], double *fplus);
static double fmincon(const double fun_workspace_center[3], const double
                      fun_workspace_A_alpha[18], const double
                      fun_workspace_b_alpha[6], const double
                      fun_workspace_deltaT[3], const double fun_workspace_T_min
                      [3], const double x0[9], const double Aineq_data[], const
                      double bineq_data[], const int bineq_size[2], double x[9],
                      double *exitflag, double *output_iterations, double
                      *output_funcCount, char output_algorithm[3], double
                      *output_constrviolation, double *output_stepsize, double
                      *output_lssteplength, double *output_firstorderopt);
static void fullColLDL2_(e_struct_T *obj, int NColsRemain);
static void initActiveSet(h_struct_T *obj);
static void linearForm_(bool obj_hasLinear, int obj_nvar, double workspace_data[],
  const double H[81], const double f_data[], const double x_data[]);
static double maxConstraintViolation(h_struct_T *obj, const double x_data[], int
  ix0);
static double minimum(const double x[6]);
static void modifyOverheadPhaseOne_(h_struct_T *obj);
static double optimize_forces_anonFcn1(const double center[3], const double
  A_alpha[18], const double b_alpha[6], const double deltaT[3], const double
  T_min[3], const double W[9]);
static void phaseone(const double H[81], const double f_data[], c_struct_T
                     *solution, g_struct_T *memspace, h_struct_T *workingset,
                     d_struct_T *qrmanager, e_struct_T *cholmanager, f_struct_T *
                     objective, const char options_SolverName[7], const
                     k_struct_T *runTimeOptions);
static void pinv(const double A[9], double X[9]);
static void qrf(double A_data[], const int A_size[2], int m, int n, int nfxd,
                double tau_data[]);
static void relaxed(const double Hessian[81], const double grad_data[],
                    c_struct_T *TrialState, j_struct_T *MeritFunction,
                    g_struct_T *memspace, h_struct_T *WorkingSet, d_struct_T
                    *QRManager, e_struct_T *CholManager, f_struct_T *QPObjective,
                    k_struct_T *qpoptions);
static void removeConstr(h_struct_T *obj, int idx_global);
static void setProblemType(h_struct_T *obj, int PROBLEM_TYPE);
static bool soc(const double Hessian[81], const double grad_data[], c_struct_T
                *TrialState, g_struct_T *memspace, h_struct_T *WorkingSet,
                d_struct_T *QRManager, e_struct_T *CholManager, f_struct_T
                *QPObjective, const k_struct_T *qpoptions);
static void solve(const e_struct_T *obj, double rhs_data[]);
static void sortLambdaQP(double lambda_data[], int WorkingSet_nActiveConstr,
  const int WorkingSet_sizes[5], const int WorkingSet_isActiveIdx[6], const int
  WorkingSet_Wid_data[], const int WorkingSet_Wlocalidx_data[], double
  workspace_data[]);
static void squareQ_appendCol(d_struct_T *obj, const double vec_data[], int iv0);
static void step(struct_T *stepFlags, double Hessian[81], c_struct_T *TrialState,
                 j_struct_T *MeritFunction, g_struct_T *memspace, h_struct_T
                 *WorkingSet, d_struct_T *QRManager, e_struct_T *CholManager,
                 f_struct_T *QPObjective, k_struct_T *qpoptions);
static void svd(const double A[9], double U[9], double s[3], double V[9]);
static bool test_exit(j_struct_T *MeritFunction, const h_struct_T *WorkingSet,
                      c_struct_T *TrialState, bool *Flags_fevalOK, bool
                      *Flags_done, bool *Flags_stepAccepted, bool
                      *Flags_failedLineSearch, int *Flags_stepType);
static void updateWorkingSetForNewQP(const double xk[9], h_struct_T *WorkingSet,
  int mIneq, const double cIneq_data[], const double cEq[3], int mLB, int mUB,
  int mFixed);
static void vecnorm(const double x[9], double y[3]);
static void xaxpy(int n, double a, int ix0, double y[9], int iy0);
static double xdotc(int n, const double x[9], int ix0, const double y[9], int
                    iy0);
static void xgemm(int m, int n, int k, const double A[81], int lda, const double
                  B_data[], int ib0, int ldb, double C_data[], int ldc);
static void xgemv(int m, int n, const double A_data[], int lda, const double
                  x_data[], double y_data[]);
static int xgeqp3(double A_data[], const int A_size[2], int m, int n, int
                  jpvt_data[], double tau_data[]);
static void xgerc(int m, int n, double alpha1, int ix0, const double y_data[],
                  double A_data[], int ia0, int lda);
static double xnrm2(int n, const double x[9], int ix0);
static int xpotrf(int n, double A_data[], int lda);
static void xrot(double x[9], int ix0, int iy0, double c, double s);
static double xrotg(double *a, double *b, double *s);
static void xswap(double x[9], int ix0, int iy0);
static double xzlangeM(const double x[9]);
static void xzlarf(int m, int n, int iv0, double tau, double C_data[], int ic0,
                   int ldc, double work_data[]);
static double xzlarfg(int n, double *alpha1, double x_data[], int ix0);
static void xzlascl(double cfrom, double cto, double A[9]);

/* Function Definitions */
/*
 * Arguments    : int nvar
 *                double Bk[81]
 *                const double sk_data[]
 *                double yk_data[]
 *                double workspace_data[]
 * Return Type  : bool
 */
static bool BFGSUpdate(int nvar, double Bk[81], const double sk_data[], double
  yk_data[], double workspace_data[])
{
  double curvatureS;
  double dotSY;
  double theta;
  int i;
  int i1;
  int ia;
  int iac;
  int ijA;
  int ix;
  int j;
  bool success;
  dotSY = 0.0;
  if (nvar >= 1) {
    for (ix = 0; ix < nvar; ix++) {
      dotSY += sk_data[ix] * yk_data[ix];
    }
  }

  i = (unsigned char)nvar;
  memset(&workspace_data[0], 0, (unsigned int)i * sizeof(double));
  ix = 0;
  i1 = 9 * (nvar - 1) + 1;
  for (iac = 1; iac <= i1; iac += 9) {
    j = iac + nvar;
    for (ia = iac; ia < j; ia++) {
      ijA = ia - iac;
      workspace_data[ijA] += Bk[ia - 1] * sk_data[ix];
    }

    ix++;
  }

  curvatureS = 0.0;
  if (nvar >= 1) {
    for (ix = 0; ix < nvar; ix++) {
      curvatureS += sk_data[ix] * workspace_data[ix];
    }
  }

  if (dotSY < 0.2 * curvatureS) {
    theta = 0.8 * curvatureS / (curvatureS - dotSY);
    for (ix = 0; ix < i; ix++) {
      yk_data[ix] *= theta;
    }

    if (!(1.0 - theta == 0.0)) {
      for (ix = 0; ix < nvar; ix++) {
        yk_data[ix] += (1.0 - theta) * workspace_data[ix];
      }
    }

    dotSY = 0.0;
    if (nvar >= 1) {
      for (ix = 0; ix < nvar; ix++) {
        dotSY += sk_data[ix] * yk_data[ix];
      }
    }
  }

  if ((curvatureS > 2.2204460492503131E-16) && (dotSY > 2.2204460492503131E-16))
  {
    success = true;
  } else {
    success = false;
  }

  if (success) {
    curvatureS = -1.0 / curvatureS;
    if (!(curvatureS == 0.0)) {
      ix = 1;
      for (j = 0; j < i; j++) {
        theta = workspace_data[j];
        if (theta != 0.0) {
          theta *= curvatureS;
          i1 = nvar + ix;
          for (ijA = ix; ijA < i1; ijA++) {
            Bk[ijA - 1] += workspace_data[ijA - ix] * theta;
          }
        }

        ix += 9;
      }
    }

    curvatureS = 1.0 / dotSY;
    if (!(curvatureS == 0.0)) {
      ix = 1;
      for (j = 0; j < i; j++) {
        theta = yk_data[j];
        if (theta != 0.0) {
          theta *= curvatureS;
          i1 = nvar + ix;
          for (ijA = ix; ijA < i1; ijA++) {
            Bk[ijA - 1] += yk_data[ijA - ix] * theta;
          }
        }

        ix += 9;
      }
    }
  }

  return success;
}

/*
 * Arguments    : c_struct_T *solution
 *                g_struct_T *memspace
 *                h_struct_T *workingset
 *                d_struct_T *qrmanager
 * Return Type  : void
 */
static void PresolveWorkingSet(c_struct_T *solution, g_struct_T *memspace,
  h_struct_T *workingset, d_struct_T *qrmanager)
{
  int i;
  int idxEndIneq;
  int idx_col;
  int idx_global;
  int k;
  solution->state = 82;
  i = RemoveDependentEq_(memspace, workingset, qrmanager);
  if ((i != -1) && (workingset->nActiveConstr <= qrmanager->ldq)) {
    double maxDiag;
    double tol;
    int idxDiag;
    int idxStartIneq;
    int nFixedConstr_tmp;
    int nVar_tmp;
    bool guard1;
    bool okWorkingSet;
    idxEndIneq = workingset->nActiveConstr;
    nFixedConstr_tmp = workingset->nWConstr[0] + workingset->nWConstr[1];
    nVar_tmp = workingset->nVar;
    if ((workingset->nWConstr[2] + workingset->nWConstr[3]) +
        workingset->nWConstr[4] > 0) {
      idxDiag = workingset->nVar;
      idxStartIneq = workingset->nActiveConstr;
      if (idxDiag >= idxStartIneq) {
        idxStartIneq = idxDiag;
      }

      tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)
                 idxStartIneq);
      for (idx_global = 0; idx_global < nFixedConstr_tmp; idx_global++) {
        qrmanager->jpvt.data[idx_global] = 1;
      }

      i = nFixedConstr_tmp + 1;
      if (i <= idxEndIneq) {
        memset(&qrmanager->jpvt.data[i + -1], 0, (unsigned int)((idxEndIneq - i)
                + 1) * sizeof(int));
      }

      i = workingset->nActiveConstr;
      for (idx_col = 0; idx_col < i; idx_col++) {
        idxEndIneq = qrmanager->ldq * idx_col;
        idxDiag = workingset->ldA * idx_col;
        idx_global = (unsigned char)nVar_tmp;
        for (k = 0; k < idx_global; k++) {
          qrmanager->QR.data[idxEndIneq + k] = workingset->ATwset.data[idxDiag +
            k];
        }
      }

      if (workingset->nVar * workingset->nActiveConstr == 0) {
        qrmanager->mrows = workingset->nVar;
        qrmanager->ncols = workingset->nActiveConstr;
        qrmanager->minRowCol = 0;
      } else {
        qrmanager->usedPivoting = true;
        qrmanager->mrows = workingset->nVar;
        qrmanager->ncols = workingset->nActiveConstr;
        idxDiag = workingset->nVar;
        idxStartIneq = workingset->nActiveConstr;
        if (idxDiag <= idxStartIneq) {
          idxStartIneq = idxDiag;
        }

        qrmanager->minRowCol = idxStartIneq;
        qrmanager->tau.size[0] = xgeqp3(qrmanager->QR.data, qrmanager->QR.size,
          workingset->nVar, workingset->nActiveConstr, qrmanager->jpvt.data,
          qrmanager->tau.data);
      }

      idxStartIneq = 0;
      for (idx_global = workingset->nActiveConstr - 1; idx_global + 1 > nVar_tmp;
           idx_global--) {
        idxStartIneq++;
        memspace->workspace_int.data[idxStartIneq - 1] = qrmanager->
          jpvt.data[idx_global];
      }

      maxDiag = fabs(qrmanager->QR.data[0]);
      for (idxEndIneq = 0; idxEndIneq < idx_global; idxEndIneq++) {
        maxDiag = fmax(maxDiag, fabs(qrmanager->QR.data[(qrmanager->ldq *
          (idxEndIneq + 1) + idxEndIneq) + 1]));
      }

      if (idx_global + 1 <= workingset->nVar) {
        idxDiag = idx_global + qrmanager->ldq * idx_global;
        while ((idx_global + 1 > nFixedConstr_tmp) && (fabs(qrmanager->
                 QR.data[idxDiag]) < tol * maxDiag)) {
          idxStartIneq++;
          memspace->workspace_int.data[idxStartIneq - 1] = qrmanager->
            jpvt.data[idx_global];
          idx_global--;
          idxDiag = (idxDiag - qrmanager->ldq) - 1;
        }
      }

      countsort(memspace->workspace_int.data, idxStartIneq,
                memspace->workspace_sort.data, nFixedConstr_tmp + 1,
                workingset->nActiveConstr);
      for (idx_global = idxStartIneq; idx_global >= 1; idx_global--) {
        removeConstr(workingset, memspace->workspace_int.data[idx_global - 1]);
      }
    }

    okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_float.data,
      memspace->workspace_float.size, solution->xstar.data, workingset,
      qrmanager);
    guard1 = false;
    if (!okWorkingSet) {
      idxEndIneq = workingset->nActiveConstr;
      i = workingset->nWConstr[0] + workingset->nWConstr[1];
      if ((workingset->nWConstr[2] + workingset->nWConstr[3]) +
          workingset->nWConstr[4] > 0) {
        idxDiag = workingset->nVar;
        idxStartIneq = workingset->nActiveConstr;
        if (idxDiag >= idxStartIneq) {
          idxStartIneq = idxDiag;
        }

        tol = 10.0 * fmin(1.4901161193847656E-8, 2.2204460492503131E-15 *
                          (double)idxStartIneq);
        for (idx_global = 0; idx_global < i; idx_global++) {
          qrmanager->jpvt.data[idx_global] = 1;
        }

        idx_global = i + 1;
        if (idx_global <= idxEndIneq) {
          memset(&qrmanager->jpvt.data[idx_global + -1], 0, (unsigned int)
                 ((idxEndIneq - idx_global) + 1) * sizeof(int));
        }

        idx_global = workingset->nActiveConstr;
        for (idx_col = 0; idx_col < idx_global; idx_col++) {
          idxEndIneq = qrmanager->ldq * idx_col;
          idxDiag = workingset->ldA * idx_col;
          idxStartIneq = (unsigned char)nVar_tmp;
          for (k = 0; k < idxStartIneq; k++) {
            qrmanager->QR.data[idxEndIneq + k] = workingset->ATwset.data[idxDiag
              + k];
          }
        }

        if (workingset->nVar * workingset->nActiveConstr == 0) {
          qrmanager->mrows = workingset->nVar;
          qrmanager->ncols = workingset->nActiveConstr;
          qrmanager->minRowCol = 0;
        } else {
          qrmanager->usedPivoting = true;
          qrmanager->mrows = workingset->nVar;
          qrmanager->ncols = workingset->nActiveConstr;
          idxDiag = workingset->nVar;
          idxStartIneq = workingset->nActiveConstr;
          if (idxDiag <= idxStartIneq) {
            idxStartIneq = idxDiag;
          }

          qrmanager->minRowCol = idxStartIneq;
          qrmanager->tau.size[0] = xgeqp3(qrmanager->QR.data, qrmanager->QR.size,
            workingset->nVar, workingset->nActiveConstr, qrmanager->jpvt.data,
            qrmanager->tau.data);
        }

        idxStartIneq = 0;
        for (idx_global = workingset->nActiveConstr - 1; idx_global + 1 >
             nVar_tmp; idx_global--) {
          idxStartIneq++;
          memspace->workspace_int.data[idxStartIneq - 1] = qrmanager->
            jpvt.data[idx_global];
        }

        maxDiag = fabs(qrmanager->QR.data[0]);
        for (idxEndIneq = 0; idxEndIneq < idx_global; idxEndIneq++) {
          maxDiag = fmax(maxDiag, fabs(qrmanager->QR.data[(qrmanager->ldq *
            (idxEndIneq + 1) + idxEndIneq) + 1]));
        }

        if (idx_global + 1 <= workingset->nVar) {
          idxDiag = idx_global + qrmanager->ldq * idx_global;
          while ((idx_global + 1 > i) && (fabs(qrmanager->QR.data[idxDiag]) <
                  tol * maxDiag)) {
            idxStartIneq++;
            memspace->workspace_int.data[idxStartIneq - 1] =
              qrmanager->jpvt.data[idx_global];
            idx_global--;
            idxDiag = (idxDiag - qrmanager->ldq) - 1;
          }
        }

        countsort(memspace->workspace_int.data, idxStartIneq,
                  memspace->workspace_sort.data, i + 1,
                  workingset->nActiveConstr);
        for (idx_global = idxStartIneq; idx_global >= 1; idx_global--) {
          removeConstr(workingset, memspace->workspace_int.data[idx_global - 1]);
        }
      }

      okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_float.data,
        memspace->workspace_float.size, solution->xstar.data, workingset,
        qrmanager);
      if (!okWorkingSet) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1 && (workingset->nWConstr[0] + workingset->nWConstr[1] ==
                   workingset->nVar)) {
      tol = b_maxConstraintViolation(workingset, solution->xstar.data);
      if (tol > 1.0E-6) {
        solution->state = -2;
      }
    }
  } else {
    int idxDiag;
    int idxStartIneq;
    solution->state = -3;
    idxDiag = workingset->nWConstr[0] + workingset->nWConstr[1];
    idxStartIneq = idxDiag + 1;
    idxEndIneq = workingset->nActiveConstr;
    for (idx_global = idxStartIneq; idx_global <= idxEndIneq; idx_global++) {
      workingset->isActiveConstr.data[(workingset->isActiveIdx
        [workingset->Wid.data[idx_global - 1] - 1] + workingset->
        Wlocalidx.data[idx_global - 1]) - 2] = false;
    }

    workingset->nWConstr[2] = 0;
    workingset->nWConstr[3] = 0;
    workingset->nWConstr[4] = 0;
    workingset->nActiveConstr = idxDiag;
  }
}

/*
 * Arguments    : g_struct_T *memspace
 *                h_struct_T *workingset
 *                d_struct_T *qrmanager
 * Return Type  : int
 */
static int RemoveDependentEq_(g_struct_T *memspace, h_struct_T *workingset,
  d_struct_T *qrmanager)
{
  int idx;
  int idxDiag;
  int idx_col;
  int k;
  int mTotalWorkingEq_tmp;
  int mWorkingFixed;
  int nDepInd;
  mWorkingFixed = workingset->nWConstr[0];
  mTotalWorkingEq_tmp = workingset->nWConstr[0] + workingset->nWConstr[1];
  nDepInd = 0;
  if (mTotalWorkingEq_tmp > 0) {
    double tol;
    int i;
    int i1;
    int u0;
    i = (unsigned char)workingset->nVar;
    for (idxDiag = 0; idxDiag < mTotalWorkingEq_tmp; idxDiag++) {
      for (idx_col = 0; idx_col < i; idx_col++) {
        qrmanager->QR.data[idxDiag + qrmanager->ldq * idx_col] =
          workingset->ATwset.data[idx_col + workingset->ldA * idxDiag];
      }
    }

    nDepInd = mTotalWorkingEq_tmp - workingset->nVar;
    if (nDepInd <= 0) {
      nDepInd = 0;
    }

    memset(&qrmanager->jpvt.data[0], 0, (unsigned int)i * sizeof(int));
    i1 = mTotalWorkingEq_tmp * workingset->nVar;
    if (i1 == 0) {
      qrmanager->mrows = mTotalWorkingEq_tmp;
      qrmanager->ncols = workingset->nVar;
      qrmanager->minRowCol = 0;
    } else {
      qrmanager->usedPivoting = true;
      qrmanager->mrows = mTotalWorkingEq_tmp;
      qrmanager->ncols = workingset->nVar;
      idxDiag = workingset->nVar;
      if (mTotalWorkingEq_tmp <= idxDiag) {
        idxDiag = mTotalWorkingEq_tmp;
      }

      qrmanager->minRowCol = idxDiag;
      qrmanager->tau.size[0] = xgeqp3(qrmanager->QR.data, qrmanager->QR.size,
        mTotalWorkingEq_tmp, workingset->nVar, qrmanager->jpvt.data,
        qrmanager->tau.data);
    }

    idxDiag = workingset->nVar;
    if (mTotalWorkingEq_tmp >= idxDiag) {
      idxDiag = mTotalWorkingEq_tmp;
    }

    tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)idxDiag);
    u0 = workingset->nVar;
    if (u0 > mTotalWorkingEq_tmp) {
      u0 = mTotalWorkingEq_tmp;
    }

    idxDiag = u0 + qrmanager->ldq * (u0 - 1);
    while ((idxDiag > 0) && (fabs(qrmanager->QR.data[idxDiag - 1]) < tol * fabs
            (qrmanager->QR.data[0]))) {
      idxDiag = (idxDiag - qrmanager->ldq) - 1;
      nDepInd++;
    }

    if (nDepInd > 0) {
      bool exitg1;
      computeQ_(qrmanager, qrmanager->mrows);
      idx = 0;
      exitg1 = false;
      while ((!exitg1) && (idx <= nDepInd - 1)) {
        double qtb;
        idxDiag = qrmanager->ldq * ((mTotalWorkingEq_tmp - idx) - 1);
        qtb = 0.0;
        for (k = 0; k < mTotalWorkingEq_tmp; k++) {
          qtb += qrmanager->Q.data[idxDiag + k] * workingset->bwset.data[k];
        }

        if (fabs(qtb) >= tol) {
          nDepInd = -1;
          exitg1 = true;
        } else {
          idx++;
        }
      }
    }

    if (nDepInd > 0) {
      int ix0;
      for (idx_col = 0; idx_col < mTotalWorkingEq_tmp; idx_col++) {
        idxDiag = qrmanager->ldq * idx_col;
        ix0 = workingset->ldA * idx_col;
        for (k = 0; k < i; k++) {
          qrmanager->QR.data[idxDiag + k] = workingset->ATwset.data[ix0 + k];
        }
      }

      for (idx = 0; idx < mWorkingFixed; idx++) {
        qrmanager->jpvt.data[idx] = 1;
      }

      idx_col = workingset->nWConstr[0] + 1;
      if (idx_col <= mTotalWorkingEq_tmp) {
        memset(&qrmanager->jpvt.data[idx_col + -1], 0, (unsigned int)
               ((mTotalWorkingEq_tmp - idx_col) + 1) * sizeof(int));
      }

      if (i1 == 0) {
        qrmanager->mrows = workingset->nVar;
        qrmanager->ncols = mTotalWorkingEq_tmp;
        qrmanager->minRowCol = 0;
      } else {
        qrmanager->usedPivoting = true;
        qrmanager->mrows = workingset->nVar;
        qrmanager->ncols = mTotalWorkingEq_tmp;
        qrmanager->minRowCol = u0;
        qrmanager->tau.size[0] = xgeqp3(qrmanager->QR.data, qrmanager->QR.size,
          workingset->nVar, mTotalWorkingEq_tmp, qrmanager->jpvt.data,
          qrmanager->tau.data);
      }

      for (idx = 0; idx < nDepInd; idx++) {
        memspace->workspace_int.data[idx] = qrmanager->jpvt.data
          [(mTotalWorkingEq_tmp - nDepInd) + idx];
      }

      countsort(memspace->workspace_int.data, nDepInd,
                memspace->workspace_sort.data, 1, mTotalWorkingEq_tmp);
      for (idx = nDepInd; idx >= 1; idx--) {
        i1 = workingset->nWConstr[0] + workingset->nWConstr[1];
        if (i1 != 0) {
          idx_col = memspace->workspace_int.data[idx - 1];
          if (idx_col <= i1) {
            if ((workingset->nActiveConstr == i1) || (idx_col == i1)) {
              workingset->mEqRemoved++;
              workingset->indexEqRemoved[workingset->mEqRemoved - 1] =
                workingset->Wlocalidx.data[idx_col - 1];
              removeConstr(workingset, idx_col);
            } else {
              workingset->mEqRemoved++;
              ix0 = workingset->Wid.data[idx_col - 1] - 1;
              idxDiag = workingset->Wlocalidx.data[idx_col - 1];
              workingset->indexEqRemoved[workingset->mEqRemoved - 1] = idxDiag;
              workingset->isActiveConstr.data[(workingset->isActiveIdx[ix0] +
                idxDiag) - 2] = false;
              workingset->Wid.data[idx_col - 1] = workingset->Wid.data[i1 - 1];
              workingset->Wlocalidx.data[idx_col - 1] =
                workingset->Wlocalidx.data[i1 - 1];
              for (idxDiag = 0; idxDiag < i; idxDiag++) {
                workingset->ATwset.data[idxDiag + workingset->ldA * (idx_col - 1)]
                  = workingset->ATwset.data[idxDiag + workingset->ldA * (i1 - 1)];
              }

              workingset->bwset.data[idx_col - 1] = workingset->bwset.data[i1 -
                1];
              idx_col = workingset->nActiveConstr - 1;
              workingset->Wid.data[i1 - 1] = workingset->Wid.data[idx_col];
              workingset->Wlocalidx.data[i1 - 1] = workingset->
                Wlocalidx.data[idx_col];
              for (idxDiag = 0; idxDiag < i; idxDiag++) {
                workingset->ATwset.data[idxDiag + workingset->ldA * (i1 - 1)] =
                  workingset->ATwset.data[idxDiag + workingset->ldA * idx_col];
              }

              workingset->bwset.data[i1 - 1] = workingset->bwset.data[idx_col];
              workingset->nActiveConstr = idx_col;
              workingset->nWConstr[ix0]--;
            }
          }
        }
      }
    }
  }

  return nDepInd;
}

/*
 * Arguments    : h_struct_T *obj
 *                int idx_local
 * Return Type  : void
 */
static void addAeqConstr(h_struct_T *obj, int idx_local)
{
  int idx;
  int totalEq;
  totalEq = obj->nWConstr[0] + obj->nWConstr[1];
  if ((obj->nActiveConstr == totalEq) && (idx_local > obj->nWConstr[1])) {
    int i;
    int i1;
    int iAeq0;
    int iAw0;
    obj->nWConstr[1]++;
    obj->isActiveConstr.data[(obj->isActiveIdx[1] + idx_local) - 2] = true;
    obj->nActiveConstr++;
    i = obj->nActiveConstr - 1;
    obj->Wid.data[i] = 2;
    obj->Wlocalidx.data[i] = idx_local;
    iAeq0 = obj->ldA * (idx_local - 1);
    iAw0 = obj->ldA * (obj->nActiveConstr - 1);
    i1 = (unsigned char)obj->nVar;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset.data[iAw0 + idx] = obj->Aeq.data[iAeq0 + idx];
    }

    obj->bwset.data[i] = obj->beq[idx_local - 1];
  } else {
    int i;
    int i1;
    int iAeq0;
    int iAw0;
    obj->nActiveConstr++;
    i = obj->nActiveConstr - 1;
    obj->Wid.data[i] = obj->Wid.data[totalEq];
    obj->Wlocalidx.data[i] = obj->Wlocalidx.data[totalEq];
    i1 = (unsigned char)obj->nVar;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset.data[idx + obj->ldA * i] = obj->ATwset.data[idx + obj->ldA *
        totalEq];
    }

    obj->bwset.data[i] = obj->bwset.data[totalEq];
    obj->nWConstr[1]++;
    obj->isActiveConstr.data[(obj->isActiveIdx[1] + idx_local) - 2] = true;
    obj->Wid.data[totalEq] = 2;
    obj->Wlocalidx.data[totalEq] = idx_local;
    iAeq0 = obj->ldA * (idx_local - 1);
    iAw0 = obj->ldA * totalEq;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset.data[iAw0 + idx] = obj->Aeq.data[iAeq0 + idx];
    }

    obj->bwset.data[totalEq] = obj->beq[idx_local - 1];
  }
}

/*
 * Arguments    : h_struct_T *obj
 *                int TYPE
 *                int idx_local
 * Return Type  : void
 */
static void addBoundToActiveSetMatrix_(h_struct_T *obj, int TYPE, int idx_local)
{
  int colOffset;
  int i;
  int idx_bnd_local;
  obj->nWConstr[TYPE - 1]++;
  obj->isActiveConstr.data[(obj->isActiveIdx[TYPE - 1] + idx_local) - 2] = true;
  obj->nActiveConstr++;
  i = obj->nActiveConstr - 1;
  obj->Wid.data[i] = TYPE;
  obj->Wlocalidx.data[i] = idx_local;
  colOffset = obj->ldA * i - 1;
  if (TYPE == 5) {
    idx_bnd_local = obj->indexUB.data[idx_local - 1];
    obj->bwset.data[i] = obj->ub.data[idx_bnd_local - 1];
  } else {
    idx_bnd_local = obj->indexLB.data[idx_local - 1];
    obj->bwset.data[i] = obj->lb.data[idx_bnd_local - 1];
  }

  i = (unsigned char)(idx_bnd_local - 1);
  if (i - 1 >= 0) {
    memset(&obj->ATwset.data[colOffset + 1], 0, (unsigned int)i * sizeof(double));
  }

  obj->ATwset.data[idx_bnd_local + colOffset] = 2.0 * (double)(TYPE == 5) - 1.0;
  i = idx_bnd_local + 1;
  idx_bnd_local = obj->nVar;
  if (i <= idx_bnd_local) {
    memset(&obj->ATwset.data[i + colOffset], 0, (unsigned int)((((idx_bnd_local
               + colOffset) - i) - colOffset) + 1) * sizeof(double));
  }

  switch (obj->probType) {
   case 3:
   case 2:
    break;

   default:
    obj->ATwset.data[obj->nVar + colOffset] = -1.0;
    break;
  }
}

/*
 * Arguments    : const bool x[3]
 * Return Type  : bool
 */
static bool all(const bool x[3])
{
  int k;
  bool exitg1;
  bool y;
  y = true;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= 2)) {
    if (!x[k]) {
      y = false;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return y;
}

/*
 * Arguments    : const double H[81]
 *                const double f_data[]
 *                c_struct_T *solution
 *                g_struct_T *memspace
 *                h_struct_T *workingset
 *                d_struct_T *qrmanager
 *                e_struct_T *cholmanager
 *                f_struct_T *objective
 *                const k_struct_T *options
 *                const k_struct_T *runTimeOptions
 * Return Type  : void
 */
static void b_driver(const double H[81], const double f_data[], c_struct_T
                     *solution, g_struct_T *memspace, h_struct_T *workingset,
                     d_struct_T *qrmanager, e_struct_T *cholmanager, f_struct_T *
                     objective, const k_struct_T *options, const k_struct_T
                     *runTimeOptions)
{
  static const char b[7] = { 'f', 'm', 'i', 'n', 'c', 'o', 'n' };

  double d;
  double minLambda;
  int i;
  int i1;
  int idx;
  int idx_local;
  int ret;
  bool guard1;
  bool guard2;
  bool updateFval;
  solution->iterations = 0;
  ret = workingset->nVar;
  guard1 = false;
  guard2 = false;
  if (workingset->probType == 3) {
    i = (unsigned char)workingset->sizes[0];
    for (idx = 0; idx < i; idx++) {
      solution->xstar.data[workingset->indexFixed.data[idx] - 1] =
        workingset->ub.data[workingset->indexFixed.data[idx] - 1];
    }

    i = (unsigned char)workingset->sizes[3];
    for (idx = 0; idx < i; idx++) {
      if (workingset->isActiveConstr.data[(workingset->isActiveIdx[3] + idx) - 1])
      {
        solution->xstar.data[workingset->indexLB.data[idx] - 1] =
          -workingset->lb.data[workingset->indexLB.data[idx] - 1];
      }
    }

    i = (unsigned char)workingset->sizes[4];
    for (idx = 0; idx < i; idx++) {
      if (workingset->isActiveConstr.data[(workingset->isActiveIdx[4] + idx) - 1])
      {
        solution->xstar.data[workingset->indexUB.data[idx] - 1] =
          workingset->ub.data[workingset->indexUB.data[idx] - 1];
      }
    }

    PresolveWorkingSet(solution, memspace, workingset, qrmanager);
    if (solution->state >= 0) {
      guard2 = true;
    }
  } else {
    solution->state = 82;
    guard2 = true;
  }

  if (guard2) {
    solution->iterations = 0;
    solution->maxConstr = b_maxConstraintViolation(workingset,
      solution->xstar.data);
    if (solution->maxConstr > 1.0E-6) {
      phaseone(H, f_data, solution, memspace, workingset, qrmanager, cholmanager,
               objective, options->SolverName, runTimeOptions);
      if (solution->state != 0) {
        solution->maxConstr = b_maxConstraintViolation(workingset,
          solution->xstar.data);
        if (solution->maxConstr > 1.0E-6) {
          ret = workingset->mConstrMax;
          if (ret - 1 >= 0) {
            memset(&solution->lambda.data[0], 0, (unsigned int)ret * sizeof
                   (double));
          }

          d = 0.0;
          switch (objective->objtype) {
           case 5:
            d = objective->gammaScalar * solution->xstar.data[objective->nvar -
              1];
            break;

           case 3:
            linearForm_(objective->hasLinear, objective->nvar,
                        memspace->workspace_float.data, H, f_data,
                        solution->xstar.data);
            if (objective->nvar >= 1) {
              ret = objective->nvar;
              for (idx_local = 0; idx_local < ret; idx_local++) {
                d += solution->xstar.data[idx_local] *
                  memspace->workspace_float.data[idx_local];
              }
            }
            break;

           case 4:
            linearForm_(objective->hasLinear, objective->nvar,
                        memspace->workspace_float.data, H, f_data,
                        solution->xstar.data);
            i = objective->nvar + 1;
            i1 = objective->maxVar;
            for (idx = i; idx < i1; idx++) {
              memspace->workspace_float.data[idx - 1] = 0.5 * objective->beta *
                solution->xstar.data[idx - 1] + objective->rho;
            }

            if (objective->maxVar - 1 >= 1) {
              for (idx_local = 0; idx_local <= i1 - 2; idx_local++) {
                d += solution->xstar.data[idx_local] *
                  memspace->workspace_float.data[idx_local];
              }
            }
            break;
          }

          solution->fstar = d;
          solution->state = -2;
        } else {
          if (solution->maxConstr > 0.0) {
            if (ret - 1 >= 0) {
              memcpy(&solution->searchDir.data[0], &solution->xstar.data[0],
                     (unsigned int)ret * sizeof(double));
            }

            PresolveWorkingSet(solution, memspace, workingset, qrmanager);
            minLambda = b_maxConstraintViolation(workingset,
              solution->xstar.data);
            if (minLambda >= solution->maxConstr) {
              solution->maxConstr = minLambda;
              if (ret - 1 >= 0) {
                memcpy(&solution->xstar.data[0], &solution->searchDir.data[0],
                       (unsigned int)ret * sizeof(double));
              }
            }
          }

          guard1 = true;
        }
      }
    } else {
      guard1 = true;
    }
  }

  if (guard1) {
    int TYPE;
    int activeSetChangeID;
    int globalActiveConstrIdx;
    bool subProblemChanged;
    subProblemChanged = true;
    updateFval = true;
    activeSetChangeID = 0;
    TYPE = objective->objtype;
    i = workingset->nVar;
    globalActiveConstrIdx = 0;
    computeGrad_StoreHx(objective, H, f_data, solution->xstar.data);
    solution->fstar = computeFval_ReuseHx(objective,
      memspace->workspace_float.data, f_data, solution->xstar.data);
    if (solution->iterations < runTimeOptions->MaxIterations) {
      solution->state = -5;
    } else {
      solution->state = 0;
    }

    ret = workingset->mConstrMax;
    if (ret - 1 >= 0) {
      memset(&solution->lambda.data[0], 0, (unsigned int)ret * sizeof(double));
    }

    int exitg1;
    do {
      exitg1 = 0;
      if (solution->state == -5) {
        int i2;
        bool guard3;
        bool guard4;
        guard3 = false;
        guard4 = false;
        if (subProblemChanged) {
          switch (activeSetChangeID) {
           case 1:
            squareQ_appendCol(qrmanager, workingset->ATwset.data,
                              workingset->ldA * (workingset->nActiveConstr - 1)
                              + 1);
            break;

           case -1:
            deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
            break;

           default:
            factorQR(qrmanager, workingset->ATwset.data, i,
                     workingset->nActiveConstr, workingset->ldA);
            computeQ_(qrmanager, qrmanager->mrows);
            break;
          }

          ret = memcmp(&options->SolverName[0], &b[0], 7);
          compute_deltax(H, solution, memspace, qrmanager, cholmanager,
                         objective, ret == 0);
          if (solution->state != -5) {
            exitg1 = 1;
          } else if ((c_xnrm2(i, solution->searchDir.data) < 1.0E-6) ||
                     (workingset->nActiveConstr >= i)) {
            guard4 = true;
          } else {
            minLambda = feasibleratiotest(solution->xstar.data,
              solution->searchDir.data, memspace->workspace_float.data,
              memspace->workspace_float.size, workingset->nVar, workingset->ldA,
              workingset->Aineq.data, workingset->bineq.data,
              workingset->lb.data, workingset->ub.data, workingset->indexLB.data,
              workingset->indexUB.data, workingset->sizes,
              workingset->isActiveIdx, workingset->isActiveConstr.data,
              workingset->nWConstr, TYPE == 5, &updateFval, &i1, &idx_local);
            if (updateFval) {
              switch (i1) {
               case 3:
                workingset->nWConstr[2]++;
                workingset->isActiveConstr.data[(workingset->isActiveIdx[2] +
                  idx_local) - 2] = true;
                workingset->nActiveConstr++;
                i1 = workingset->nActiveConstr - 1;
                workingset->Wid.data[i1] = 3;
                workingset->Wlocalidx.data[i1] = idx_local;
                ret = workingset->ldA * (idx_local - 1);
                activeSetChangeID = workingset->ldA * i1;
                i2 = workingset->nVar;
                for (idx = 0; idx < i2; idx++) {
                  workingset->ATwset.data[activeSetChangeID + idx] =
                    workingset->Aineq.data[ret + idx];
                }

                workingset->bwset.data[i1] = workingset->bineq.data[idx_local -
                  1];
                break;

               case 4:
                addBoundToActiveSetMatrix_(workingset, 4, idx_local);
                break;

               default:
                addBoundToActiveSetMatrix_(workingset, 5, idx_local);
                break;
              }

              activeSetChangeID = 1;
            } else {
              if (objective->objtype == 5) {
                if (c_xnrm2(objective->nvar, solution->searchDir.data) > 100.0 *
                    (double)objective->nvar * 1.4901161193847656E-8) {
                  solution->state = 3;
                } else {
                  solution->state = 4;
                }
              }

              subProblemChanged = false;
              if (workingset->nActiveConstr == 0) {
                solution->state = 1;
              }
            }

            if ((i >= 1) && (!(minLambda == 0.0))) {
              for (idx_local = 0; idx_local < i; idx_local++) {
                solution->xstar.data[idx_local] += minLambda *
                  solution->searchDir.data[idx_local];
              }
            }

            computeGrad_StoreHx(objective, H, f_data, solution->xstar.data);
            guard3 = true;
          }
        } else {
          if (i - 1 >= 0) {
            memset(&solution->searchDir.data[0], 0, (unsigned int)i * sizeof
                   (double));
          }

          guard4 = true;
        }

        if (guard4) {
          compute_lambda(memspace->workspace_float.data, solution, objective,
                         qrmanager);
          if ((solution->state != -7) || (workingset->nActiveConstr > i)) {
            ret = 0;
            minLambda = 0.0;
            i1 = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
            i2 = workingset->nActiveConstr;
            for (idx = i1; idx <= i2; idx++) {
              d = solution->lambda.data[idx - 1];
              if (d < minLambda) {
                minLambda = d;
                ret = idx;
              }
            }

            if (ret == 0) {
              solution->state = 1;
            } else {
              activeSetChangeID = -1;
              globalActiveConstrIdx = ret;
              subProblemChanged = true;
              removeConstr(workingset, ret);
              if (ret < workingset->nActiveConstr + 1) {
                solution->lambda.data[ret - 1] = solution->
                  lambda.data[workingset->nActiveConstr];
              }

              solution->lambda.data[workingset->nActiveConstr] = 0.0;
            }
          } else {
            ret = workingset->nActiveConstr;
            activeSetChangeID = 0;
            globalActiveConstrIdx = workingset->nActiveConstr;
            subProblemChanged = true;
            removeConstr(workingset, workingset->nActiveConstr);
            solution->lambda.data[ret - 1] = 0.0;
          }

          guard3 = true;
        }

        if (guard3) {
          solution->iterations++;
          ret = objective->nvar;
          if ((solution->iterations >= runTimeOptions->MaxIterations) &&
              ((solution->state != 1) || (objective->objtype == 5))) {
            solution->state = 0;
          }

          if (solution->iterations - solution->iterations / 50 * 50 == 0) {
            solution->maxConstr = b_maxConstraintViolation(workingset,
              solution->xstar.data);
            minLambda = solution->maxConstr;
            if (objective->objtype == 5) {
              minLambda = solution->maxConstr - solution->xstar.data
                [objective->nvar - 1];
            }

            if (minLambda > 1.0E-6) {
              if (ret - 1 >= 0) {
                memcpy(&solution->searchDir.data[0], &solution->xstar.data[0],
                       (unsigned int)ret * sizeof(double));
              }

              updateFval = feasibleX0ForWorkingSet
                (memspace->workspace_float.data, memspace->workspace_float.size,
                 solution->searchDir.data, workingset, qrmanager);
              if ((!updateFval) && (solution->state != 0)) {
                solution->state = -2;
              }

              activeSetChangeID = 0;
              minLambda = b_maxConstraintViolation(workingset,
                solution->searchDir.data);
              if (minLambda < solution->maxConstr) {
                if (ret - 1 >= 0) {
                  memcpy(&solution->xstar.data[0], &solution->searchDir.data[0],
                         (unsigned int)ret * sizeof(double));
                }

                solution->maxConstr = minLambda;
              }
            }
          }

          updateFval = false;
        }
      } else {
        if (!updateFval) {
          solution->fstar = computeFval_ReuseHx(objective,
            memspace->workspace_float.data, f_data, solution->xstar.data);
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
}

/*
 * Arguments    : int mIneqMax
 *                int nVarMax
 *                int mConstrMax
 *                h_struct_T *obj
 * Return Type  : void
 */
static void b_factoryConstruct(int mIneqMax, int nVarMax, int mConstrMax,
  h_struct_T *obj)
{
  int i;
  obj->mConstr = 0;
  obj->mConstrOrig = 0;
  obj->mConstrMax = mConstrMax;
  obj->nVar = 9;
  obj->nVarMax = nVarMax;
  obj->ldA = nVarMax;
  obj->Aineq.size[0] = mIneqMax * nVarMax;
  obj->bineq.size[0] = mIneqMax;
  obj->lb.size[0] = nVarMax;
  obj->ub.size[0] = nVarMax;
  obj->mEqRemoved = 0;
  obj->bwset.size[0] = mConstrMax;
  obj->nActiveConstr = 0;
  obj->maxConstrWorkspace.size[0] = mConstrMax;
  for (i = 0; i < 5; i++) {
    obj->sizes[i] = 0;
    obj->sizesNormal[i] = 0;
    obj->sizesPhaseOne[i] = 0;
    obj->sizesRegularized[i] = 0;
    obj->sizesRegPhaseOne[i] = 0;
  }

  for (i = 0; i < 6; i++) {
    obj->isActiveIdx[i] = 0;
    obj->isActiveIdxNormal[i] = 0;
    obj->isActiveIdxPhaseOne[i] = 0;
    obj->isActiveIdxRegularized[i] = 0;
    obj->isActiveIdxRegPhaseOne[i] = 0;
  }

  for (i = 0; i < 5; i++) {
    obj->nWConstr[i] = 0;
  }

  obj->probType = 3;
}

/*
 * Arguments    : h_struct_T *obj
 *                const double x_data[]
 * Return Type  : double
 */
static double b_maxConstraintViolation(h_struct_T *obj, const double x_data[])
{
  double v;
  int i;
  int idx;
  int mFixed;
  int mLB;
  int mUB;
  mLB = obj->sizes[3];
  mUB = obj->sizes[4];
  mFixed = obj->sizes[0];
  if (obj->probType == 2) {
    v = c_maxConstraintViolation_AMats_(obj, x_data);
  } else {
    v = d_maxConstraintViolation_AMats_(obj, x_data);
  }

  if (mLB > 0) {
    i = (unsigned char)mLB;
    for (idx = 0; idx < i; idx++) {
      mLB = obj->indexLB.data[idx] - 1;
      v = fmax(v, -x_data[mLB] - obj->lb.data[mLB]);
    }
  }

  if (mUB > 0) {
    i = (unsigned char)mUB;
    for (idx = 0; idx < i; idx++) {
      mLB = obj->indexUB.data[idx] - 1;
      v = fmax(v, x_data[mLB] - obj->ub.data[mLB]);
    }
  }

  if (mFixed > 0) {
    i = (unsigned char)mFixed;
    for (idx = 0; idx < i; idx++) {
      v = fmax(v, fabs(x_data[obj->indexFixed.data[idx] - 1] - obj->ub.data
                       [obj->indexFixed.data[idx] - 1]));
    }
  }

  return v;
}

/*
 * Arguments    : struct_T *Flags
 *                g_struct_T *memspace
 *                j_struct_T *MeritFunction
 *                int fscales_lineq_constraint_size
 *                h_struct_T *WorkingSet
 *                c_struct_T *TrialState
 *                d_struct_T *QRManager
 * Return Type  : void
 */
static void b_test_exit(struct_T *Flags, g_struct_T *memspace, j_struct_T
  *MeritFunction, int fscales_lineq_constraint_size, h_struct_T *WorkingSet,
  c_struct_T *TrialState, d_struct_T *QRManager)
{
  double optimRelativeFactor;
  double s;
  double smax;
  int ia;
  int iac;
  int idx_max;
  int mFixed;
  int mIneq;
  int mLB;
  int mLambda;
  int mUB;
  int nVar;
  int rankR;
  bool isFeasible;
  nVar = WorkingSet->nVar;
  mFixed = WorkingSet->sizes[0];
  mIneq = WorkingSet->sizes[2];
  mLB = WorkingSet->sizes[3];
  mUB = WorkingSet->sizes[4];
  mLambda = (((WorkingSet->sizes[0] + WorkingSet->sizes[2]) + WorkingSet->sizes
              [3]) + WorkingSet->sizes[4]) + 2;
  if (mLambda >= 0) {
    memcpy(&TrialState->lambdaStopTest.data[0], &TrialState->lambdasqp.data[0],
           (unsigned int)(mLambda + 1) * sizeof(double));
  }

  computeGradLag(TrialState->gradLag.data, WorkingSet->ldA, WorkingSet->nVar,
                 TrialState->grad.data, WorkingSet->sizes[2],
                 WorkingSet->Aineq.data, WorkingSet->Aeq.data,
                 WorkingSet->indexFixed.data, WorkingSet->sizes[0],
                 WorkingSet->indexLB.data, WorkingSet->sizes[3],
                 WorkingSet->indexUB.data, WorkingSet->sizes[4],
                 TrialState->lambdaStopTest.data);
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = fabs(TrialState->grad.data[0]);
      for (rankR = 2; rankR <= nVar; rankR++) {
        s = fabs(TrialState->grad.data[rankR - 1]);
        if (s > smax) {
          idx_max = rankR;
          smax = s;
        }
      }
    }
  }

  optimRelativeFactor = fmax(1.0, fabs(TrialState->grad.data[idx_max - 1]));
  if (rtIsInf(optimRelativeFactor)) {
    optimRelativeFactor = 1.0;
  }

  MeritFunction->nlpPrimalFeasError = computePrimalFeasError
    (TrialState->xstarsqp, WorkingSet->sizes[2], TrialState->cIneq.data,
     TrialState->cEq, WorkingSet->indexLB.data, WorkingSet->sizes[3],
     WorkingSet->indexUB.data, WorkingSet->sizes[4]);
  if (TrialState->sqpIterations == 0) {
    MeritFunction->feasRelativeFactor = fmax(1.0,
      MeritFunction->nlpPrimalFeasError);
  }

  isFeasible = (MeritFunction->nlpPrimalFeasError <= 1.0E-6 *
                MeritFunction->feasRelativeFactor);
  Flags->gradOK = computeDualFeasError(WorkingSet->nVar,
    TrialState->gradLag.data, &MeritFunction->nlpDualFeasError);
  if (!Flags->gradOK) {
    Flags->done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    int fullRank_R;
    MeritFunction->nlpComplError = computeComplError
      (fscales_lineq_constraint_size, TrialState->xstarsqp, WorkingSet->sizes[2],
       TrialState->cIneq.data, WorkingSet->indexLB.data, WorkingSet->sizes[3],
       WorkingSet->indexUB.data, WorkingSet->sizes[4],
       TrialState->lambdaStopTest.data, WorkingSet->sizes[0] + 4);
    MeritFunction->firstOrderOpt = fmax(MeritFunction->nlpDualFeasError,
      MeritFunction->nlpComplError);
    if (TrialState->sqpIterations > 1) {
      computeGradLag(memspace->workspace_float.data, WorkingSet->ldA,
                     WorkingSet->nVar, TrialState->grad.data, WorkingSet->sizes
                     [2], WorkingSet->Aineq.data, WorkingSet->Aeq.data,
                     WorkingSet->indexFixed.data, WorkingSet->sizes[0],
                     WorkingSet->indexLB.data, WorkingSet->sizes[3],
                     WorkingSet->indexUB.data, WorkingSet->sizes[4],
                     TrialState->lambdaStopTestPrev.data);
      smax = 0.0;
      fullRank_R = (unsigned char)WorkingSet->nVar;
      idx_max = 0;
      while ((idx_max <= fullRank_R - 1) && ((!rtIsInf
               (memspace->workspace_float.data[idx_max])) && (!rtIsNaN
               (memspace->workspace_float.data[idx_max])))) {
        smax = fmax(smax, fabs(memspace->workspace_float.data[idx_max]));
        idx_max++;
      }

      s = computeComplError(fscales_lineq_constraint_size, TrialState->xstarsqp,
                            WorkingSet->sizes[2], TrialState->cIneq.data,
                            WorkingSet->indexLB.data, WorkingSet->sizes[3],
                            WorkingSet->indexUB.data, WorkingSet->sizes[4],
                            TrialState->lambdaStopTestPrev.data,
                            WorkingSet->sizes[0] + 4);
      if ((smax < MeritFunction->nlpDualFeasError) && (s <
           MeritFunction->nlpComplError)) {
        MeritFunction->nlpDualFeasError = smax;
        MeritFunction->nlpComplError = s;
        MeritFunction->firstOrderOpt = fmax(smax, s);
        if (mLambda >= 0) {
          memcpy(&TrialState->lambdaStopTest.data[0],
                 &TrialState->lambdaStopTestPrev.data[0], (unsigned int)(mLambda
                  + 1) * sizeof(double));
        }
      } else if (mLambda >= 0) {
        memcpy(&TrialState->lambdaStopTestPrev.data[0],
               &TrialState->lambdaStopTest.data[0], (unsigned int)(mLambda + 1) *
               sizeof(double));
      }
    } else if (mLambda >= 0) {
      memcpy(&TrialState->lambdaStopTestPrev.data[0],
             &TrialState->lambdaStopTest.data[0], (unsigned int)(mLambda + 1) *
             sizeof(double));
    }

    if (isFeasible && (MeritFunction->nlpDualFeasError <= 1.0E-6 *
                       optimRelativeFactor) && (MeritFunction->nlpComplError <=
         1.0E-6 * optimRelativeFactor)) {
      Flags->done = true;
      TrialState->sqpExitFlag = 1;
    } else {
      Flags->done = false;
      if (isFeasible && (TrialState->sqpFval < -1.0E+20)) {
        Flags->done = true;
        TrialState->sqpExitFlag = -3;
      } else {
        bool guard1;
        guard1 = false;
        if (TrialState->sqpIterations > 0) {
          bool dxTooSmall;
          bool exitg1;
          dxTooSmall = true;
          fullRank_R = (unsigned char)WorkingSet->nVar;
          idx_max = 0;
          exitg1 = false;
          while ((!exitg1) && (idx_max <= fullRank_R - 1)) {
            if (1.0E-6 * fmax(1.0, fabs(TrialState->xstarsqp[idx_max])) <= fabs
                (TrialState->delta_x.data[idx_max])) {
              dxTooSmall = false;
              exitg1 = true;
            } else {
              idx_max++;
            }
          }

          if (dxTooSmall) {
            if (!isFeasible) {
              if (Flags->stepType == 2) {
                Flags->done = true;
                TrialState->sqpExitFlag = -2;
              } else {
                Flags->stepType = 2;
                Flags->failedLineSearch = false;
                Flags->stepAccepted = false;
                guard1 = true;
              }
            } else {
              int nActiveConstr;
              nActiveConstr = WorkingSet->nActiveConstr;
              if (WorkingSet->nActiveConstr == 0) {
                Flags->done = true;
                TrialState->sqpExitFlag = 2;
              } else {
                double d;
                int ix;
                updateWorkingSetForNewQP(TrialState->xstarsqp, WorkingSet,
                  WorkingSet->sizes[2], TrialState->cIneq.data, TrialState->cEq,
                  WorkingSet->sizes[3], WorkingSet->sizes[4], WorkingSet->sizes
                  [0]);
                if (nActiveConstr - 1 >= 0) {
                  memset(&TrialState->lambda.data[0], 0, (unsigned int)
                         nActiveConstr * sizeof(double));
                }

                factorQRE(QRManager, WorkingSet->ATwset.data, nVar,
                          nActiveConstr, WorkingSet->ldA);
                computeQ_(QRManager, QRManager->mrows);
                idx_max = QRManager->ldq;
                fullRank_R = (unsigned char)nVar;
                memset(&memspace->workspace_float.data[0], 0, (unsigned int)
                       fullRank_R * sizeof(double));
                ix = 0;
                fullRank_R = QRManager->ldq * (nVar - 1) + 1;
                for (iac = 1; idx_max < 0 ? iac >= fullRank_R : iac <=
                     fullRank_R; iac += idx_max) {
                  smax = 0.0;
                  rankR = iac + nVar;
                  for (ia = iac; ia < rankR; ia++) {
                    smax += QRManager->Q.data[ia - 1] * TrialState->grad.data[ia
                      - iac];
                  }

                  memspace->workspace_float.data[ix] -= smax;
                  ix++;
                }

                if (nVar >= nActiveConstr) {
                  idx_max = nVar;
                } else {
                  idx_max = nActiveConstr;
                }

                smax = fabs(QRManager->QR.data[0]) * fmin(1.4901161193847656E-8,
                  (double)idx_max * 2.2204460492503131E-16);
                if (nVar <= nActiveConstr) {
                  fullRank_R = nVar;
                } else {
                  fullRank_R = nActiveConstr;
                }

                rankR = 0;
                idx_max = 0;
                while ((rankR < fullRank_R) && (fabs(QRManager->QR.data[idx_max])
                        > smax)) {
                  rankR++;
                  idx_max = (idx_max + QRManager->ldq) + 1;
                }

                if (rankR != 0) {
                  for (iac = rankR; iac >= 1; iac--) {
                    idx_max = (iac + (iac - 1) * QRManager->ldq) - 1;
                    memspace->workspace_float.data[iac - 1] /=
                      QRManager->QR.data[idx_max];
                    for (ia = 0; ia <= iac - 2; ia++) {
                      ix = (iac - ia) - 2;
                      memspace->workspace_float.data[ix] -=
                        memspace->workspace_float.data[iac - 1] *
                        QRManager->QR.data[(idx_max - ia) - 1];
                    }
                  }
                }

                if (nActiveConstr <= fullRank_R) {
                  fullRank_R = nActiveConstr;
                }

                for (idx_max = 0; idx_max < fullRank_R; idx_max++) {
                  TrialState->lambda.data[QRManager->jpvt.data[idx_max] - 1] =
                    memspace->workspace_float.data[idx_max];
                }

                sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
                             WorkingSet->sizes, WorkingSet->isActiveIdx,
                             WorkingSet->Wid.data, WorkingSet->Wlocalidx.data,
                             memspace->workspace_float.data);
                computeGradLag(memspace->workspace_float.data, WorkingSet->ldA,
                               nVar, TrialState->grad.data, mIneq,
                               WorkingSet->Aineq.data, WorkingSet->Aeq.data,
                               WorkingSet->indexFixed.data, mFixed,
                               WorkingSet->indexLB.data, mLB,
                               WorkingSet->indexUB.data, mUB,
                               TrialState->lambda.data);
                smax = 0.0;
                idx_max = 0;
                while ((idx_max <= (unsigned char)nVar - 1) && ((!rtIsInf
                         (memspace->workspace_float.data[idx_max])) && (!rtIsNaN
                         (memspace->workspace_float.data[idx_max])))) {
                  smax = fmax(smax, fabs(memspace->workspace_float.data[idx_max]));
                  idx_max++;
                }

                s = computeComplError(fscales_lineq_constraint_size,
                                      TrialState->xstarsqp, mIneq,
                                      TrialState->cIneq.data,
                                      WorkingSet->indexLB.data, mLB,
                                      WorkingSet->indexUB.data, mUB,
                                      TrialState->lambda.data, mFixed + 4);
                d = fmax(smax, s);
                if (d <= fmax(MeritFunction->nlpDualFeasError,
                              MeritFunction->nlpComplError)) {
                  MeritFunction->nlpDualFeasError = smax;
                  MeritFunction->nlpComplError = s;
                  MeritFunction->firstOrderOpt = d;
                  if (mLambda >= 0) {
                    memcpy(&TrialState->lambdaStopTest.data[0],
                           &TrialState->lambda.data[0], (unsigned int)(mLambda +
                            1) * sizeof(double));
                  }
                }

                if ((MeritFunction->nlpDualFeasError <= 1.0E-6 *
                     optimRelativeFactor) && (MeritFunction->nlpComplError <=
                     1.0E-6 * optimRelativeFactor)) {
                  TrialState->sqpExitFlag = 1;
                } else {
                  TrialState->sqpExitFlag = 2;
                }

                Flags->done = true;
                guard1 = true;
              }
            }
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          if (TrialState->sqpIterations >= 400) {
            Flags->done = true;
            TrialState->sqpExitFlag = 0;
          } else if (TrialState->FunctionEvaluations >= 900) {
            Flags->done = true;
            TrialState->sqpExitFlag = 0;
          }
        }
      }
    }
  }
}

/*
 * Arguments    : const double xk[9]
 *                h_struct_T *WorkingSet
 *                int mIneq
 *                const double cIneq_data[]
 *                const double cEq[3]
 * Return Type  : void
 */
static void b_updateWorkingSetForNewQP(const double xk[9], h_struct_T
  *WorkingSet, int mIneq, const double cIneq_data[], const double cEq[3])
{
  int iEq0;
  int idx;
  int iw0;
  iw0 = 0;
  iEq0 = 0;
  for (idx = 0; idx < 3; idx++) {
    double d;
    d = cEq[idx];
    WorkingSet->beq[idx] = -d;
    WorkingSet->bwset.data[idx] = -d;
    memcpy(&WorkingSet->ATwset.data[iw0], &WorkingSet->Aeq.data[iEq0], 9U *
           sizeof(double));
    iw0 += WorkingSet->ldA;
    iEq0 += WorkingSet->ldA;
  }

  iw0 = (unsigned char)mIneq;
  for (idx = 0; idx < iw0; idx++) {
    WorkingSet->bineq.data[idx] = -cIneq_data[idx];
  }

  for (idx = 0; idx < 9; idx++) {
    WorkingSet->lb.data[WorkingSet->indexLB.data[idx] - 1] = xk
      [WorkingSet->indexLB.data[idx] - 1] + 1.0;
  }

  for (idx = 0; idx < 9; idx++) {
    WorkingSet->ub.data[WorkingSet->indexUB.data[idx] - 1] = (double)
      iv[WorkingSet->indexUB.data[idx] - 1] - xk[WorkingSet->indexUB.data[idx] -
      1];
  }
}

/*
 * Arguments    : const double x[18]
 *                double y[6]
 * Return Type  : void
 */
static void b_vecnorm(const double x[18], double y[6])
{
  int j;
  for (j = 0; j < 6; j++) {
    double absxk;
    double b_y;
    double scale;
    double t;
    scale = 3.3121686421112381E-170;
    absxk = fabs(x[j]);
    if (absxk > 3.3121686421112381E-170) {
      b_y = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      b_y = t * t;
    }

    absxk = fabs(x[j + 6]);
    if (absxk > scale) {
      t = scale / absxk;
      b_y = b_y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      b_y += t * t;
    }

    absxk = fabs(x[j + 12]);
    if (absxk > scale) {
      t = scale / absxk;
      b_y = b_y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      b_y += t * t;
    }

    y[j] = scale * sqrt(b_y);
  }
}

/*
 * Arguments    : double a
 *                const double x[9]
 *                int ix0
 *                double y[3]
 * Return Type  : void
 */
static void b_xaxpy(double a, const double x[9], int ix0, double y[3])
{
  int k;
  if (!(a == 0.0)) {
    for (k = 0; k < 2; k++) {
      y[k + 1] += a * x[(ix0 + k) - 1];
    }
  }
}

/*
 * Arguments    : int m
 *                int n
 *                int k
 *                const double A_data[]
 *                int ia0
 *                int lda
 *                const double B_data[]
 *                int ldb
 *                double C_data[]
 *                int ldc
 * Return Type  : void
 */
static void b_xgemm(int m, int n, int k, const double A_data[], int ia0, int lda,
                    const double B_data[], int ldb, double C_data[], int ldc)
{
  int cr;
  int ic;
  int w;
  if ((m != 0) && (n != 0)) {
    int br;
    int i;
    int i1;
    int lastColC;
    lastColC = ldc * (n - 1);
    for (cr = 0; ldc < 0 ? cr >= lastColC : cr <= lastColC; cr += ldc) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        memset(&C_data[i + -1], 0, (unsigned int)((i1 - i) + 1) * sizeof(double));
      }
    }

    br = -1;
    for (cr = 0; ldc < 0 ? cr >= lastColC : cr <= lastColC; cr += ldc) {
      int ar;
      ar = ia0;
      i = cr + 1;
      i1 = cr + m;
      for (ic = i; ic <= i1; ic++) {
        double temp;
        temp = 0.0;
        for (w = 0; w < k; w++) {
          temp += A_data[(w + ar) - 1] * B_data[(w + br) + 1];
        }

        C_data[ic - 1] += temp;
        ar += lda;
      }

      br += ldb;
    }
  }
}

/*
 * Arguments    : int m
 *                int n
 *                const double A_data[]
 *                int lda
 *                const double x_data[]
 *                int ix0
 *                double y_data[]
 * Return Type  : void
 */
static void b_xgemv(int m, int n, const double A_data[], int lda, const double
                    x_data[], int ix0, double y_data[])
{
  int ia;
  int iac;
  int iy;
  if (n != 0) {
    int i;
    i = (unsigned char)n;
    for (iy = 0; iy < i; iy++) {
      y_data[iy] = -y_data[iy];
    }

    iy = 0;
    i = lda * (n - 1) + 1;
    for (iac = 1; lda < 0 ? iac >= i : iac <= i; iac += lda) {
      double c;
      int i1;
      c = 0.0;
      i1 = iac + m;
      for (ia = iac; ia < i1; ia++) {
        c += A_data[ia - 1] * x_data[((ix0 + ia) - iac) - 1];
      }

      y_data[iy] += c;
      iy++;
    }
  }
}

/*
 * Arguments    : int n
 *                const double x_data[]
 *                int ix0
 * Return Type  : double
 */
static double b_xnrm2(int n, const double x_data[], int ix0)
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x_data[ix0 - 1]);
    } else {
      double scale;
      int kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (k = ix0; k < kend; k++) {
        double absxk;
        absxk = fabs(x_data[k - 1]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : double cfrom
 *                double cto
 *                double A[3]
 * Return Type  : void
 */
static void b_xzlascl(double cfrom, double cto, double A[3])
{
  double cfromc;
  double ctoc;
  bool notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }

    A[0] *= mul;
    A[1] *= mul;
    A[2] *= mul;
  }
}

/*
 * Arguments    : h_struct_T *obj
 *                const double x_data[]
 * Return Type  : double
 */
static double c_maxConstraintViolation_AMats_(h_struct_T *obj, const double
  x_data[])
{
  double v;
  int idx;
  int mIneq;
  v = 0.0;
  mIneq = obj->sizes[2];
  if (obj->Aineq.size[0] != 0) {
    int i;
    if (mIneq - 1 >= 0) {
      memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0], (unsigned
              int)mIneq * sizeof(double));
    }

    xgemv(9, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data,
          obj->maxConstrWorkspace.data);
    i = (unsigned char)obj->sizes[2];
    for (idx = 0; idx < i; idx++) {
      double d;
      d = obj->maxConstrWorkspace.data[idx] - x_data[idx + 9];
      obj->maxConstrWorkspace.data[idx] = d;
      v = fmax(v, d);
    }
  }

  obj->maxConstrWorkspace.data[0] = obj->beq[0];
  obj->maxConstrWorkspace.data[1] = obj->beq[1];
  obj->maxConstrWorkspace.data[2] = obj->beq[2];
  xgemv(9, 3, obj->Aeq.data, obj->ldA, x_data, obj->maxConstrWorkspace.data);
  obj->maxConstrWorkspace.data[0] = (obj->maxConstrWorkspace.data[0] -
    x_data[mIneq + 9]) + x_data[obj->sizes[2] + 12];
  v = fmax(v, fabs(obj->maxConstrWorkspace.data[0]));
  obj->maxConstrWorkspace.data[1] = (obj->maxConstrWorkspace.data[1] -
    x_data[mIneq + 10]) + x_data[obj->sizes[2] + 13];
  v = fmax(v, fabs(obj->maxConstrWorkspace.data[1]));
  obj->maxConstrWorkspace.data[2] = (obj->maxConstrWorkspace.data[2] -
    x_data[mIneq + 11]) + x_data[obj->sizes[2] + 14];
  return fmax(v, fabs(obj->maxConstrWorkspace.data[2]));
}

/*
 * Arguments    : double a
 *                const double x[3]
 *                double y[9]
 *                int iy0
 * Return Type  : void
 */
static void c_xaxpy(double a, const double x[3], double y[9], int iy0)
{
  int k;
  if (!(a == 0.0)) {
    for (k = 0; k < 2; k++) {
      int i;
      i = (iy0 + k) - 1;
      y[i] += a * x[k + 1];
    }
  }
}

/*
 * Arguments    : int n
 *                const double x_data[]
 * Return Type  : double
 */
static double c_xnrm2(int n, const double x_data[])
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x_data[0]);
    } else {
      double scale;
      scale = 3.3121686421112381E-170;
      for (k = 0; k < n; k++) {
        double absxk;
        absxk = fabs(x_data[k]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

/*
 * Arguments    : int fscales_lineq_constraint_size
 *                const double xCurrent[9]
 *                int mIneq
 *                const double cIneq_data[]
 *                const int finiteLB_data[]
 *                int mLB
 *                const int finiteUB_data[]
 *                int mUB
 *                const double lambda_data[]
 *                int iL0
 * Return Type  : double
 */
static double computeComplError(int fscales_lineq_constraint_size, const double
  xCurrent[9], int mIneq, const double cIneq_data[], const int finiteLB_data[],
  int mLB, const int finiteUB_data[], int mUB, const double lambda_data[], int
  iL0)
{
  double nlpComplError;
  int idx;
  nlpComplError = 0.0;
  if ((mIneq + mLB) + mUB > 0) {
    double lbLambda;
    double ubLambda;
    int i;
    int lbOffset;
    int ubOffset;
    for (idx = 0; idx < fscales_lineq_constraint_size; idx++) {
      lbLambda = lambda_data[(iL0 + idx) - 1];
      ubLambda = cIneq_data[idx];
      nlpComplError = fmax(nlpComplError, fmin(fabs(ubLambda * lbLambda), fmin
        (fabs(ubLambda), lbLambda)));
    }

    lbOffset = (iL0 + mIneq) - 1;
    ubOffset = lbOffset + mLB;
    i = (unsigned char)mLB;
    for (idx = 0; idx < i; idx++) {
      lbLambda = lambda_data[lbOffset + idx];
      ubLambda = xCurrent[finiteLB_data[idx] - 1] - -1.0;
      nlpComplError = fmax(nlpComplError, fmin(fabs(ubLambda * lbLambda), fmin
        (fabs(ubLambda), lbLambda)));
    }

    i = (unsigned char)mUB;
    for (idx = 0; idx < i; idx++) {
      lbLambda = (double)iv[finiteUB_data[idx] - 1] - xCurrent[finiteUB_data[idx]
        - 1];
      ubLambda = lambda_data[ubOffset + idx];
      nlpComplError = fmax(nlpComplError, fmin(fabs(lbLambda * ubLambda), fmin
        (fabs(lbLambda), ubLambda)));
    }
  }

  return nlpComplError;
}

/*
 * Arguments    : const double x[9]
 *                const double Cineq_workspace_data[]
 *                int ineq0
 *                double Ceq_workspace[3]
 * Return Type  : int
 */
static int computeConstraints_(const double x[9], const double
  Cineq_workspace_data[], int ineq0, double Ceq_workspace[3])
{
  double scale;
  int idx;
  int idx_current;
  int status;
  bool allFinite;

  /* 'optimize_forces:71' @(W) constraint(W) */
  /*  优化约束：单位向量约束 */
  /* 'optimize_forces:31' ceq = zeros([length(W) / 3, 1]); */
  /* 'optimize_forces:32' for idx = 1:length(W) / 3 */
  for (idx = 0; idx < 3; idx++) {
    double absxk;
    double t;
    double y;

    /* 'optimize_forces:33' ceq(idx) = norm(W(3*idx-2:3*idx)) - 1; */
    scale = 3.3121686421112381E-170;
    idx_current = 3 * (idx + 1) - 3;
    absxk = fabs(x[idx_current]);
    if (absxk > 3.3121686421112381E-170) {
      y = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      y = t * t;
    }

    absxk = fabs(x[idx_current + 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }

    absxk = fabs(x[idx_current + 2]);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }

    y = scale * sqrt(y);
    Ceq_workspace[idx] = y - 1.0;
  }

  /* 'optimize_forces:35' c = []; */
  status = 1;
  allFinite = true;
  idx_current = ineq0;
  while (allFinite && (idx_current <= ineq0 - 1)) {
    scale = Cineq_workspace_data[idx_current - 1];
    allFinite = ((!rtIsInf(scale)) && (!rtIsNaN(scale)));
    idx_current++;
  }

  if (!allFinite) {
    idx_current -= 2;
    if (rtIsNaN(Cineq_workspace_data[idx_current])) {
      status = -3;
    } else if (Cineq_workspace_data[idx_current] < 0.0) {
      status = -1;
    } else {
      status = -2;
    }
  } else {
    allFinite = true;
    idx_current = 0;
    while (allFinite && (idx_current + 1 <= 3)) {
      allFinite = ((!rtIsInf(Ceq_workspace[idx_current])) && (!rtIsNaN
        (Ceq_workspace[idx_current])));
      idx_current++;
    }

    if (!allFinite) {
      idx_current--;
      if (rtIsNaN(Ceq_workspace[idx_current])) {
        status = -3;
      } else if (Ceq_workspace[idx_current] < 0.0) {
        status = -1;
      } else {
        status = -2;
      }
    }
  }

  return status;
}

/*
 * Arguments    : int nVar
 *                const double gradLag_data[]
 *                double *val
 * Return Type  : bool
 */
static bool computeDualFeasError(int nVar, const double gradLag_data[], double
  *val)
{
  int idx;
  bool exitg1;
  bool gradOK;
  gradOK = true;
  *val = 0.0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx <= (unsigned char)nVar - 1)) {
    gradOK = ((!rtIsInf(gradLag_data[idx])) && (!rtIsNaN(gradLag_data[idx])));
    if (!gradOK) {
      exitg1 = true;
    } else {
      *val = fmax(*val, fabs(gradLag_data[idx]));
      idx++;
    }
  }

  return gradOK;
}

/*
 * Arguments    : i_struct_T *obj
 *                double fCurrent
 *                const double cEqCurrent[3]
 *                double xk[9]
 *                double gradf_data[]
 *                double JacCeqTrans_data[]
 *                int ldJE
 * Return Type  : bool
 */
static bool computeFiniteDifferences(i_struct_T *obj, double fCurrent, const
  double cEqCurrent[3], double xk[9], double gradf_data[], double
  JacCeqTrans_data[], int ldJE)
{
  int idx;
  bool evalOK;
  bool exitg1;
  evalOK = true;
  obj->numEvals = 0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx < 9)) {
    double d;
    double deltaX;
    double ubDiff;
    bool guard1;
    bool modifiedStep;
    ubDiff = 1.4901161193847656E-8 * (1.0 - 2.0 * (double)(xk[idx] < 0.0)) *
      fmax(fabs(xk[idx]), 1.0);
    modifiedStep = false;
    if (xk[idx] >= -1.0) {
      int i;
      i = iv[idx];
      if (xk[idx] <= i) {
        d = xk[idx] + ubDiff;
        if ((d > i) || (d < -1.0)) {
          ubDiff = -ubDiff;
          modifiedStep = true;
          d = xk[idx] + ubDiff;
          if ((d > i) || (d < -1.0)) {
            ubDiff = (double)i - xk[idx];
            if (xk[idx] - -1.0 <= ubDiff) {
              ubDiff = -(xk[idx] - -1.0);
            }
          }
        }
      }
    }

    deltaX = ubDiff;
    evalOK = finDiffEvalAndChkErr(obj->objfun.workspace.center,
      obj->objfun.workspace.A_alpha, obj->objfun.workspace.b_alpha,
      obj->objfun.workspace.deltaT, obj->objfun.workspace.T_min, obj->cEq_1, idx
      + 1, ubDiff, xk, &obj->f_1);
    obj->numEvals++;
    guard1 = false;
    if (!evalOK) {
      if (!modifiedStep) {
        deltaX = -ubDiff;
        d = xk[idx] - ubDiff;
        if ((d >= -1.0) && (d <= iv[idx])) {
          modifiedStep = true;
        } else {
          modifiedStep = false;
        }

        if ((!obj->hasBounds) || modifiedStep) {
          evalOK = finDiffEvalAndChkErr(obj->objfun.workspace.center,
            obj->objfun.workspace.A_alpha, obj->objfun.workspace.b_alpha,
            obj->objfun.workspace.deltaT, obj->objfun.workspace.T_min,
            obj->cEq_1, idx + 1, -ubDiff, xk, &obj->f_1);
          obj->numEvals++;
        }
      }

      if (!evalOK) {
        exitg1 = true;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      gradf_data[idx] = (obj->f_1 - fCurrent) / deltaX;
      JacCeqTrans_data[idx] = (obj->cEq_1[0] - cEqCurrent[0]) / deltaX;
      JacCeqTrans_data[idx + ldJE] = (obj->cEq_1[1] - cEqCurrent[1]) / deltaX;
      JacCeqTrans_data[idx + (ldJE << 1)] = (obj->cEq_1[2] - cEqCurrent[2]) /
        deltaX;
      idx++;
    }
  }

  return evalOK;
}

/*
 * Arguments    : const f_struct_T *obj
 *                double workspace_data[]
 *                const double f_data[]
 *                const double x_data[]
 * Return Type  : double
 */
static double computeFval_ReuseHx(const f_struct_T *obj, double workspace_data[],
  const double f_data[], const double x_data[])
{
  double val;
  int k;
  val = 0.0;
  switch (obj->objtype) {
   case 5:
    val = obj->gammaScalar * x_data[obj->nvar - 1];
    break;

   case 3:
    {
      if (obj->hasLinear) {
        int ixlast;
        ixlast = obj->nvar;
        for (k = 0; k < ixlast; k++) {
          workspace_data[k] = 0.5 * obj->Hx.data[k] + f_data[k];
        }

        if (obj->nvar >= 1) {
          for (k = 0; k < ixlast; k++) {
            val += x_data[k] * workspace_data[k];
          }
        }
      } else {
        if (obj->nvar >= 1) {
          int ixlast;
          ixlast = obj->nvar;
          for (k = 0; k < ixlast; k++) {
            val += x_data[k] * obj->Hx.data[k];
          }
        }

        val *= 0.5;
      }
    }
    break;

   case 4:
    {
      int maxRegVar_tmp;
      maxRegVar_tmp = obj->maxVar;
      if (obj->hasLinear) {
        int ixlast;
        ixlast = obj->nvar;
        if (ixlast - 1 >= 0) {
          memcpy(&workspace_data[0], &f_data[0], (unsigned int)ixlast * sizeof
                 (double));
        }

        ixlast = obj->maxVar - obj->nvar;
        for (k = 0; k <= ixlast - 2; k++) {
          workspace_data[obj->nvar + k] = obj->rho;
        }

        ixlast = (unsigned char)(obj->maxVar - 1);
        for (k = 0; k < ixlast; k++) {
          workspace_data[k] += 0.5 * obj->Hx.data[k];
        }

        if (obj->maxVar - 1 >= 1) {
          for (k = 0; k <= maxRegVar_tmp - 2; k++) {
            val += x_data[k] * workspace_data[k];
          }
        }
      } else {
        int ixlast;
        if (obj->maxVar - 1 >= 1) {
          for (k = 0; k <= maxRegVar_tmp - 2; k++) {
            val += x_data[k] * obj->Hx.data[k];
          }
        }

        val *= 0.5;
        ixlast = obj->nvar + 1;
        for (k = ixlast; k < maxRegVar_tmp; k++) {
          val += x_data[k - 1] * obj->rho;
        }
      }
    }
    break;
  }

  return val;
}

/*
 * Arguments    : double workspace_data[]
 *                int ldA
 *                int nVar
 *                const double grad_data[]
 *                int mIneq
 *                const double AineqTrans_data[]
 *                const double AeqTrans_data[]
 *                const int finiteFixed_data[]
 *                int mFixed
 *                const int finiteLB_data[]
 *                int mLB
 *                const int finiteUB_data[]
 *                int mUB
 *                const double lambda_data[]
 * Return Type  : void
 */
static void computeGradLag(double workspace_data[], int ldA, int nVar, const
  double grad_data[], int mIneq, const double AineqTrans_data[], const double
  AeqTrans_data[], const int finiteFixed_data[], int mFixed, const int
  finiteLB_data[], int mLB, const int finiteUB_data[], int mUB, const double
  lambda_data[])
{
  int i;
  int i1;
  int ia;
  int iac;
  int idx;
  int ix;
  i = (unsigned char)nVar;
  memcpy(&workspace_data[0], &grad_data[0], (unsigned int)i * sizeof(double));
  i = (unsigned char)mFixed;
  for (idx = 0; idx < i; idx++) {
    workspace_data[finiteFixed_data[idx] - 1] += lambda_data[idx];
  }

  ix = mFixed;
  i = (ldA << 1) + 1;
  for (iac = 1; ldA < 0 ? iac >= i : iac <= i; iac += ldA) {
    idx = iac + nVar;
    for (ia = iac; ia < idx; ia++) {
      i1 = ia - iac;
      workspace_data[i1] += AeqTrans_data[ia - 1] * lambda_data[ix];
    }

    ix++;
  }

  if (mIneq != 0) {
    ix = mFixed + 3;
    i = ldA * (mIneq - 1) + 1;
    for (iac = 1; ldA < 0 ? iac >= i : iac <= i; iac += ldA) {
      idx = iac + nVar;
      for (ia = iac; ia < idx; ia++) {
        i1 = ia - iac;
        workspace_data[i1] += AineqTrans_data[ia - 1] * lambda_data[ix];
      }

      ix++;
    }
  }

  ix = (mFixed + mIneq) + 3;
  i = (unsigned char)mLB;
  for (idx = 0; idx < i; idx++) {
    workspace_data[finiteLB_data[idx] - 1] -= lambda_data[ix + idx];
  }

  if ((unsigned char)mLB - 1 >= 0) {
    ix += (unsigned char)mLB;
  }

  i = (unsigned char)mUB;
  for (idx = 0; idx < i; idx++) {
    workspace_data[finiteUB_data[idx] - 1] += lambda_data[ix + idx];
  }
}

/*
 * Arguments    : f_struct_T *obj
 *                const double H[81]
 *                const double f_data[]
 *                const double x_data[]
 * Return Type  : void
 */
static void computeGrad_StoreHx(f_struct_T *obj, const double H[81], const
  double f_data[], const double x_data[])
{
  double y_data[40];
  int ia;
  int iac;
  int ix;
  int ixlast;
  switch (obj->objtype) {
   case 5:
    {
      int i;
      i = obj->nvar;
      if (i - 2 >= 0) {
        memset(&obj->grad.data[0], 0, (unsigned int)(i - 1) * sizeof(double));
      }

      obj->grad.data[obj->nvar - 1] = obj->gammaScalar;
    }
    break;

   case 3:
    {
      int i;
      int iy;
      ixlast = obj->nvar - 1;
      iy = obj->nvar;
      if (obj->nvar != 0) {
        if (ixlast >= 0) {
          memset(&obj->Hx.data[0], 0, (unsigned int)(ixlast + 1) * sizeof(double));
        }

        ix = 0;
        i = obj->nvar * ixlast + 1;
        for (iac = 1; iy < 0 ? iac >= i : iac <= i; iac += iy) {
          int i1;
          i1 = iac + ixlast;
          for (ia = iac; ia <= i1; ia++) {
            int i2;
            i2 = ia - iac;
            obj->Hx.data[i2] += H[ia - 1] * x_data[ix];
          }

          ix++;
        }
      }

      i = obj->nvar;
      if (i - 1 >= 0) {
        memcpy(&obj->grad.data[0], &obj->Hx.data[0], (unsigned int)i * sizeof
               (double));
      }

      if (obj->hasLinear) {
        iy = obj->grad.size[0];
        ixlast = obj->grad.size[0];
        if (ixlast - 1 >= 0) {
          memcpy(&y_data[0], &obj->grad.data[0], (unsigned int)ixlast * sizeof
                 (double));
        }

        if (obj->nvar >= 1) {
          ixlast = obj->nvar;
          for (ix = 0; ix < ixlast; ix++) {
            y_data[ix] += f_data[ix];
          }
        }

        if (iy - 1 >= 0) {
          memcpy(&obj->grad.data[0], &y_data[0], (unsigned int)iy * sizeof
                 (double));
        }
      }
    }
    break;

   case 4:
    {
      int i;
      int iy;
      int maxRegVar;
      maxRegVar = obj->maxVar;
      ixlast = obj->nvar - 1;
      iy = obj->nvar;
      if (obj->nvar != 0) {
        if (ixlast >= 0) {
          memset(&obj->Hx.data[0], 0, (unsigned int)(ixlast + 1) * sizeof(double));
        }

        ix = 0;
        i = obj->nvar * (obj->nvar - 1) + 1;
        for (iac = 1; iy < 0 ? iac >= i : iac <= i; iac += iy) {
          int i1;
          i1 = iac + ixlast;
          for (ia = iac; ia <= i1; ia++) {
            int i2;
            i2 = ia - iac;
            obj->Hx.data[i2] += H[ia - 1] * x_data[ix];
          }

          ix++;
        }
      }

      i = obj->nvar + 1;
      for (ixlast = i; ixlast < maxRegVar; ixlast++) {
        obj->Hx.data[ixlast - 1] = obj->beta * x_data[ixlast - 1];
      }

      i = (unsigned char)(obj->maxVar - 1);
      memcpy(&obj->grad.data[0], &obj->Hx.data[0], (unsigned int)i * sizeof
             (double));
      if (obj->hasLinear) {
        iy = obj->grad.size[0];
        ixlast = obj->grad.size[0];
        if (ixlast - 1 >= 0) {
          memcpy(&y_data[0], &obj->grad.data[0], (unsigned int)ixlast * sizeof
                 (double));
        }

        if (obj->nvar >= 1) {
          ixlast = obj->nvar;
          for (ix = 0; ix < ixlast; ix++) {
            y_data[ix] += f_data[ix];
          }
        }

        if (iy - 1 >= 0) {
          memcpy(&obj->grad.data[0], &y_data[0], (unsigned int)iy * sizeof
                 (double));
        }
      }

      ixlast = (obj->maxVar - obj->nvar) - 1;
      if (ixlast >= 1) {
        iy = obj->nvar;
        for (ix = 0; ix < ixlast; ix++) {
          i = iy + ix;
          obj->grad.data[i] += obj->rho;
        }
      }
    }
    break;
  }
}

/*
 * Arguments    : const double x[9]
 *                int nVar
 *                double workspaceIneq_data[]
 *                const int *workspaceIneq_size
 *                int mLinIneq
 *                const double AineqT_data[]
 *                const double bineq_data[]
 *                int ldAi
 * Return Type  : void
 */
static void computeLinearResiduals(const double x[9], int nVar, double
  workspaceIneq_data[], const int *workspaceIneq_size, int mLinIneq, const
  double AineqT_data[], const double bineq_data[], int ldAi)
{
  double y_data[40];
  int ia;
  int iac;
  int iy;
  if (mLinIneq > 0) {
    int i;
    int loop_ub;
    loop_ub = *workspaceIneq_size;
    if (loop_ub - 1 >= 0) {
      memcpy(&y_data[0], &workspaceIneq_data[0], (unsigned int)loop_ub * sizeof
             (double));
    }

    i = (unsigned char)mLinIneq;
    memcpy(&y_data[0], &bineq_data[0], (unsigned int)i * sizeof(double));
    if (loop_ub - 1 >= 0) {
      memcpy(&workspaceIneq_data[0], &y_data[0], (unsigned int)loop_ub * sizeof
             (double));
    }

    for (iy = 0; iy < i; iy++) {
      workspaceIneq_data[iy] = -workspaceIneq_data[iy];
    }

    iy = 0;
    i = ldAi * (mLinIneq - 1) + 1;
    for (iac = 1; ldAi < 0 ? iac >= i : iac <= i; iac += ldAi) {
      double c;
      c = 0.0;
      loop_ub = iac + nVar;
      for (ia = iac; ia < loop_ub; ia++) {
        c += AineqT_data[ia - 1] * x[ia - iac];
      }

      workspaceIneq_data[iy] += c;
      iy++;
    }
  }
}

/*
 * Arguments    : double obj_penaltyParam
 *                double fval
 *                const double Cineq_workspace_data[]
 *                int mIneq
 *                const double Ceq_workspace[3]
 *                bool evalWellDefined
 * Return Type  : double
 */
static double computeMeritFcn(double obj_penaltyParam, double fval, const double
  Cineq_workspace_data[], int mIneq, const double Ceq_workspace[3], bool
  evalWellDefined)
{
  double val;
  int idx;
  if (evalWellDefined) {
    int i;
    val = 0.0;
    i = (unsigned char)mIneq;
    for (idx = 0; idx < i; idx++) {
      double d;
      d = Cineq_workspace_data[idx];
      if (d > 0.0) {
        val += d;
      }
    }

    val = fval + obj_penaltyParam * (((fabs(Ceq_workspace[0]) + fabs
      (Ceq_workspace[1])) + fabs(Ceq_workspace[2])) + val);
  } else {
    val = rtInf;
  }

  return val;
}

/*
 * Arguments    : const double x[9]
 *                int mLinIneq
 *                const double cIneq_data[]
 *                const double cEq[3]
 *                const int finiteLB_data[]
 *                int mLB
 *                const int finiteUB_data[]
 *                int mUB
 * Return Type  : double
 */
static double computePrimalFeasError(const double x[9], int mLinIneq, const
  double cIneq_data[], const double cEq[3], const int finiteLB_data[], int mLB,
  const int finiteUB_data[], int mUB)
{
  double feasError;
  int i;
  int idx;
  feasError = fmax(fmax(fmax(0.0, fabs(cEq[0])), fabs(cEq[1])), fabs(cEq[2]));
  i = (unsigned char)mLinIneq;
  for (idx = 0; idx < i; idx++) {
    feasError = fmax(feasError, cIneq_data[idx]);
  }

  i = (unsigned char)mLB;
  for (idx = 0; idx < i; idx++) {
    feasError = fmax(feasError, -1.0 - x[finiteLB_data[idx] - 1]);
  }

  i = (unsigned char)mUB;
  for (idx = 0; idx < i; idx++) {
    feasError = fmax(feasError, x[finiteUB_data[idx] - 1] - (double)
                     iv[finiteUB_data[idx] - 1]);
  }

  return feasError;
}

/*
 * Arguments    : d_struct_T *obj
 *                int nrows
 * Return Type  : void
 */
static void computeQ_(d_struct_T *obj, int nrows)
{
  double work_data[40];
  int i;
  int iQR0;
  int ia;
  int k;
  int lastc;
  int lda;
  int m;
  int n;
  k = obj->minRowCol;
  for (lastc = 0; lastc < k; lastc++) {
    iQR0 = obj->ldq * lastc + lastc;
    n = obj->mrows - lastc;
    if (n - 2 >= 0) {
      memcpy(&obj->Q.data[iQR0 + 1], &obj->QR.data[iQR0 + 1], (unsigned int)(n -
              1) * sizeof(double));
    }
  }

  m = obj->mrows;
  lda = obj->ldq;
  if (nrows >= 1) {
    int itau;
    for (iQR0 = k; iQR0 < nrows; iQR0++) {
      ia = iQR0 * lda;
      memset(&obj->Q.data[ia], 0, (unsigned int)m * sizeof(double));
      obj->Q.data[ia + iQR0] = 1.0;
    }

    itau = obj->minRowCol - 1;
    iQR0 = obj->Q.size[1];
    if (iQR0 - 1 >= 0) {
      memset(&work_data[0], 0, (unsigned int)iQR0 * sizeof(double));
    }

    for (i = obj->minRowCol; i >= 1; i--) {
      int b_i;
      int iaii;
      iaii = i + (i - 1) * lda;
      if (i < nrows) {
        int ic0;
        int lastv;
        obj->Q.data[iaii - 1] = 1.0;
        ic0 = iaii + lda;
        if (obj->tau.data[itau] != 0.0) {
          bool exitg2;
          lastv = (m - i) + 1;
          iQR0 = (iaii + m) - i;
          while ((lastv > 0) && (obj->Q.data[iQR0 - 1] == 0.0)) {
            lastv--;
            iQR0--;
          }

          lastc = nrows - i;
          exitg2 = false;
          while ((!exitg2) && (lastc > 0)) {
            int exitg1;
            iQR0 = ic0 + (lastc - 1) * lda;
            ia = iQR0;
            do {
              exitg1 = 0;
              if (ia <= (iQR0 + lastv) - 1) {
                if (obj->Q.data[ia - 1] != 0.0) {
                  exitg1 = 1;
                } else {
                  ia++;
                }
              } else {
                lastc--;
                exitg1 = 2;
              }
            } while (exitg1 == 0);

            if (exitg1 == 1) {
              exitg2 = true;
            }
          }
        } else {
          lastv = 0;
          lastc = 0;
        }

        if (lastv > 0) {
          if (lastc != 0) {
            if (lastc - 1 >= 0) {
              memset(&work_data[0], 0, (unsigned int)lastc * sizeof(double));
            }

            iQR0 = 0;
            b_i = ic0 + lda * (lastc - 1);
            for (k = ic0; lda < 0 ? k >= b_i : k <= b_i; k += lda) {
              double c;
              c = 0.0;
              n = k + lastv;
              for (ia = k; ia < n; ia++) {
                c += obj->Q.data[ia - 1] * obj->Q.data[((iaii + ia) - k) - 1];
              }

              work_data[iQR0] += c;
              iQR0++;
            }
          }

          xgerc(lastv, lastc, -obj->tau.data[itau], iaii, work_data, obj->Q.data,
                ic0, lda);
        }
      }

      if (i < m) {
        iQR0 = iaii + 1;
        b_i = (iaii + m) - i;
        for (k = iQR0; k <= b_i; k++) {
          obj->Q.data[k - 1] *= -obj->tau.data[itau];
        }
      }

      obj->Q.data[iaii - 1] = 1.0 - obj->tau.data[itau];
      b_i = (unsigned char)(i - 1);
      for (iQR0 = 0; iQR0 < b_i; iQR0++) {
        obj->Q.data[(iaii - iQR0) - 2] = 0.0;
      }

      itau--;
    }
  }
}

/*
 * Arguments    : const double H[81]
 *                c_struct_T *solution
 *                g_struct_T *memspace
 *                const d_struct_T *qrmanager
 *                e_struct_T *cholmanager
 *                const f_struct_T *objective
 *                bool alwaysPositiveDef
 * Return Type  : void
 */
static void compute_deltax(const double H[81], c_struct_T *solution, g_struct_T *
  memspace, const d_struct_T *qrmanager, e_struct_T *cholmanager, const
  f_struct_T *objective, bool alwaysPositiveDef)
{
  int ia;
  int iac;
  int idx;
  int idx_col;
  int idx_row;
  int k;
  int mNull_tmp;
  int nVar_tmp;
  nVar_tmp = qrmanager->mrows - 1;
  mNull_tmp = qrmanager->mrows - qrmanager->ncols;
  if (mNull_tmp <= 0) {
    if (nVar_tmp >= 0) {
      memset(&solution->searchDir.data[0], 0, (unsigned int)(nVar_tmp + 1) *
             sizeof(double));
    }
  } else {
    for (idx = 0; idx <= nVar_tmp; idx++) {
      solution->searchDir.data[idx] = -objective->grad.data[idx];
    }

    if (qrmanager->ncols <= 0) {
      switch (objective->objtype) {
       case 5:
        break;

       case 3:
        {
          if (alwaysPositiveDef) {
            cholmanager->ndims = qrmanager->mrows;
            for (idx = 0; idx <= nVar_tmp; idx++) {
              idx_row = (nVar_tmp + 1) * idx;
              idx_col = cholmanager->ldm * idx;
              for (k = 0; k <= nVar_tmp; k++) {
                cholmanager->FMat.data[idx_col + k] = H[idx_row + k];
              }
            }

            cholmanager->info = xpotrf(qrmanager->mrows, cholmanager->FMat.data,
              cholmanager->ldm);
          } else {
            cholmanager->ndims = qrmanager->mrows;
            for (idx = 0; idx <= nVar_tmp; idx++) {
              idx_row = qrmanager->mrows * idx;
              idx_col = cholmanager->ldm * idx;
              for (k = 0; k <= nVar_tmp; k++) {
                cholmanager->FMat.data[idx_col + k] = H[idx_row + k];
              }
            }

            fullColLDL2_(cholmanager, qrmanager->mrows);
            if (cholmanager->ConvexCheck) {
              idx = 0;
              int exitg1;
              do {
                exitg1 = 0;
                if (idx <= nVar_tmp) {
                  if (cholmanager->FMat.data[idx + cholmanager->ldm * idx] <=
                      0.0) {
                    cholmanager->info = -idx - 1;
                    exitg1 = 1;
                  } else {
                    idx++;
                  }
                } else {
                  cholmanager->ConvexCheck = false;
                  exitg1 = 1;
                }
              } while (exitg1 == 0);
            }
          }

          if (cholmanager->info != 0) {
            solution->state = -6;
          } else if (alwaysPositiveDef) {
            solve(cholmanager, solution->searchDir.data);
          } else {
            int i;
            int nVars;
            idx_col = cholmanager->ndims - 2;
            if (cholmanager->ndims != 0) {
              for (idx_row = 0; idx_row <= idx_col + 1; idx_row++) {
                nVars = idx_row + idx_row * cholmanager->ldm;
                i = idx_col - idx_row;
                for (idx = 0; idx <= i; idx++) {
                  int ix;
                  ix = (idx_row + idx) + 1;
                  solution->searchDir.data[ix] -= solution->
                    searchDir.data[idx_row] * cholmanager->FMat.data[(nVars +
                    idx) + 1];
                }
              }
            }

            i = cholmanager->ndims;
            for (idx = 0; idx < i; idx++) {
              solution->searchDir.data[idx] /= cholmanager->FMat.data[idx +
                cholmanager->ldm * idx];
            }

            idx_col = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx_row = idx_col; idx_row >= 1; idx_row--) {
                double smax;
                nVars = (idx_row - 1) * cholmanager->ldm;
                smax = solution->searchDir.data[idx_row - 1];
                i = idx_row + 1;
                for (idx = idx_col; idx >= i; idx--) {
                  smax -= cholmanager->FMat.data[(nVars + idx) - 1] *
                    solution->searchDir.data[idx - 1];
                }

                solution->searchDir.data[idx_row - 1] = smax;
              }
            }
          }
        }
        break;

       case 4:
        {
          if (alwaysPositiveDef) {
            int nVars;
            nVars = objective->nvar;
            cholmanager->ndims = objective->nvar;
            for (idx = 0; idx < nVars; idx++) {
              idx_row = nVars * idx;
              idx_col = cholmanager->ldm * idx;
              for (k = 0; k < nVars; k++) {
                cholmanager->FMat.data[idx_col + k] = H[idx_row + k];
              }
            }

            cholmanager->info = xpotrf(objective->nvar, cholmanager->FMat.data,
              cholmanager->ldm);
            if (cholmanager->info != 0) {
              solution->state = -6;
            } else {
              double smax;
              int i;
              solve(cholmanager, solution->searchDir.data);
              smax = 1.0 / objective->beta;
              idx_row = objective->nvar + 1;
              i = qrmanager->mrows;
              for (k = idx_row; k <= i; k++) {
                solution->searchDir.data[k - 1] *= smax;
              }
            }
          }
        }
        break;
      }
    } else {
      int nullStartIdx_tmp;
      nullStartIdx_tmp = qrmanager->ldq * qrmanager->ncols + 1;
      if (objective->objtype == 5) {
        for (idx = 0; idx < mNull_tmp; idx++) {
          memspace->workspace_float.data[idx] = -qrmanager->Q.data[nVar_tmp +
            qrmanager->ldq * (qrmanager->ncols + idx)];
        }

        idx_row = qrmanager->ldq;
        if (qrmanager->mrows != 0) {
          int i;
          int ix;
          memset(&solution->searchDir.data[0], 0, (unsigned int)(nVar_tmp + 1) *
                 sizeof(double));
          ix = 0;
          i = nullStartIdx_tmp + qrmanager->ldq * (mNull_tmp - 1);
          for (iac = nullStartIdx_tmp; idx_row < 0 ? iac >= i : iac <= i; iac +=
               idx_row) {
            idx = iac + nVar_tmp;
            for (ia = iac; ia <= idx; ia++) {
              int nVars;
              nVars = ia - iac;
              solution->searchDir.data[nVars] += qrmanager->Q.data[ia - 1] *
                memspace->workspace_float.data[ix];
            }

            ix++;
          }
        }
      } else {
        int i;
        int nVars;
        if (objective->objtype == 3) {
          xgemm(qrmanager->mrows, mNull_tmp, qrmanager->mrows, H,
                qrmanager->mrows, qrmanager->Q.data, nullStartIdx_tmp,
                qrmanager->ldq, memspace->workspace_float.data,
                memspace->workspace_float.size[0]);
          b_xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows, qrmanager->Q.data,
                  nullStartIdx_tmp, qrmanager->ldq,
                  memspace->workspace_float.data, memspace->
                  workspace_float.size[0], cholmanager->FMat.data,
                  cholmanager->ldm);
        } else if (alwaysPositiveDef) {
          nVars = qrmanager->mrows;
          xgemm(objective->nvar, mNull_tmp, objective->nvar, H, objective->nvar,
                qrmanager->Q.data, nullStartIdx_tmp, qrmanager->ldq,
                memspace->workspace_float.data, memspace->workspace_float.size[0]);
          i = objective->nvar + 1;
          for (idx_col = 0; idx_col < mNull_tmp; idx_col++) {
            for (idx_row = i; idx_row <= nVars; idx_row++) {
              memspace->workspace_float.data[(idx_row +
                memspace->workspace_float.size[0] * idx_col) - 1] =
                objective->beta * qrmanager->Q.data[(idx_row + qrmanager->
                Q.size[0] * (idx_col + qrmanager->ncols)) - 1];
            }
          }

          b_xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows, qrmanager->Q.data,
                  nullStartIdx_tmp, qrmanager->ldq,
                  memspace->workspace_float.data, memspace->
                  workspace_float.size[0], cholmanager->FMat.data,
                  cholmanager->ldm);
        }

        if (alwaysPositiveDef) {
          cholmanager->ndims = mNull_tmp;
          cholmanager->info = xpotrf(mNull_tmp, cholmanager->FMat.data,
            cholmanager->ldm);
        } else {
          cholmanager->ndims = mNull_tmp;
          fullColLDL2_(cholmanager, mNull_tmp);
          if (cholmanager->ConvexCheck) {
            idx = 0;
            int exitg1;
            do {
              exitg1 = 0;
              if (idx <= mNull_tmp - 1) {
                if (cholmanager->FMat.data[idx + cholmanager->ldm * idx] <= 0.0)
                {
                  cholmanager->info = -idx - 1;
                  exitg1 = 1;
                } else {
                  idx++;
                }
              } else {
                cholmanager->ConvexCheck = false;
                exitg1 = 1;
              }
            } while (exitg1 == 0);
          }
        }

        if (cholmanager->info != 0) {
          solution->state = -6;
        } else {
          double smax;
          int ix;
          k = qrmanager->ldq;
          if (qrmanager->mrows != 0) {
            memset(&memspace->workspace_float.data[0], 0, (unsigned int)
                   mNull_tmp * sizeof(double));
            idx_col = 0;
            i = nullStartIdx_tmp + qrmanager->ldq * (mNull_tmp - 1);
            for (iac = nullStartIdx_tmp; k < 0 ? iac >= i : iac <= i; iac += k)
            {
              smax = 0.0;
              idx = iac + nVar_tmp;
              for (ia = iac; ia <= idx; ia++) {
                smax += qrmanager->Q.data[ia - 1] * objective->grad.data[ia -
                  iac];
              }

              memspace->workspace_float.data[idx_col] -= smax;
              idx_col++;
            }
          }

          if (alwaysPositiveDef) {
            idx_col = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx_row = 0; idx_row < idx_col; idx_row++) {
                nVars = idx_row * cholmanager->ldm;
                smax = memspace->workspace_float.data[idx_row];
                for (idx = 0; idx < idx_row; idx++) {
                  smax -= cholmanager->FMat.data[nVars + idx] *
                    memspace->workspace_float.data[idx];
                }

                memspace->workspace_float.data[idx_row] = smax /
                  cholmanager->FMat.data[nVars + idx_row];
              }
            }

            if (cholmanager->ndims != 0) {
              for (idx_row = idx_col; idx_row >= 1; idx_row--) {
                nVars = (idx_row + (idx_row - 1) * cholmanager->ldm) - 1;
                memspace->workspace_float.data[idx_row - 1] /=
                  cholmanager->FMat.data[nVars];
                for (idx = 0; idx <= idx_row - 2; idx++) {
                  ix = (idx_row - idx) - 2;
                  memspace->workspace_float.data[ix] -=
                    memspace->workspace_float.data[idx_row - 1] *
                    cholmanager->FMat.data[(nVars - idx) - 1];
                }
              }
            }
          } else {
            idx_col = cholmanager->ndims - 2;
            if (cholmanager->ndims != 0) {
              for (idx_row = 0; idx_row <= idx_col + 1; idx_row++) {
                nVars = idx_row + idx_row * cholmanager->ldm;
                i = idx_col - idx_row;
                for (idx = 0; idx <= i; idx++) {
                  ix = (idx_row + idx) + 1;
                  memspace->workspace_float.data[ix] -=
                    memspace->workspace_float.data[idx_row] *
                    cholmanager->FMat.data[(nVars + idx) + 1];
                }
              }
            }

            i = cholmanager->ndims;
            for (idx = 0; idx < i; idx++) {
              memspace->workspace_float.data[idx] /= cholmanager->FMat.data[idx
                + cholmanager->ldm * idx];
            }

            idx_col = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx_row = idx_col; idx_row >= 1; idx_row--) {
                nVars = (idx_row - 1) * cholmanager->ldm;
                smax = memspace->workspace_float.data[idx_row - 1];
                i = idx_row + 1;
                for (idx = idx_col; idx >= i; idx--) {
                  smax -= cholmanager->FMat.data[(nVars + idx) - 1] *
                    memspace->workspace_float.data[idx - 1];
                }

                memspace->workspace_float.data[idx_row - 1] = smax;
              }
            }
          }

          if (qrmanager->mrows != 0) {
            memset(&solution->searchDir.data[0], 0, (unsigned int)(nVar_tmp + 1)
                   * sizeof(double));
            ix = 0;
            i = nullStartIdx_tmp + qrmanager->ldq * (mNull_tmp - 1);
            for (iac = nullStartIdx_tmp; k < 0 ? iac >= i : iac <= i; iac += k)
            {
              idx = iac + nVar_tmp;
              for (ia = iac; ia <= idx; ia++) {
                nVars = ia - iac;
                solution->searchDir.data[nVars] += qrmanager->Q.data[ia - 1] *
                  memspace->workspace_float.data[ix];
              }

              ix++;
            }
          }
        }
      }
    }
  }
}

/*
 * Arguments    : double workspace_data[]
 *                c_struct_T *solution
 *                const f_struct_T *objective
 *                const d_struct_T *qrmanager
 * Return Type  : void
 */
static void compute_lambda(double workspace_data[], c_struct_T *solution, const
  f_struct_T *objective, const d_struct_T *qrmanager)
{
  int ia;
  int iac;
  int idx;
  int j;
  int nActiveConstr_tmp;
  nActiveConstr_tmp = qrmanager->ncols;
  if (qrmanager->ncols > 0) {
    double c;
    int idxQR;
    bool guard1;
    guard1 = false;
    if (objective->objtype != 4) {
      bool nonDegenerate;
      idxQR = qrmanager->mrows;
      idx = qrmanager->ncols;
      if (idxQR >= idx) {
        idx = idxQR;
      }

      c = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)idx);
      if ((qrmanager->mrows > 0) && (qrmanager->ncols > 0)) {
        nonDegenerate = true;
      } else {
        nonDegenerate = false;
      }

      if (nonDegenerate) {
        bool guard2;
        idx = nActiveConstr_tmp;
        guard2 = false;
        if (qrmanager->mrows < qrmanager->ncols) {
          idxQR = qrmanager->mrows + qrmanager->ldq * (qrmanager->ncols - 1);
          while ((idx > qrmanager->mrows) && (fabs(qrmanager->QR.data[idxQR - 1])
                  >= c)) {
            idx--;
            idxQR -= qrmanager->ldq;
          }

          nonDegenerate = (idx == qrmanager->mrows);
          if (nonDegenerate) {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }

        if (guard2) {
          idxQR = idx + qrmanager->ldq * (idx - 1);
          while ((idx >= 1) && (fabs(qrmanager->QR.data[idxQR - 1]) >= c)) {
            idx--;
            idxQR = (idxQR - qrmanager->ldq) - 1;
          }

          nonDegenerate = (idx == 0);
        }
      }

      if (!nonDegenerate) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      int ldq;
      ldq = qrmanager->ldq;
      if (qrmanager->mrows != 0) {
        memset(&workspace_data[0], 0, (unsigned int)nActiveConstr_tmp * sizeof
               (double));
        idxQR = 0;
        idx = qrmanager->ldq * (qrmanager->ncols - 1) + 1;
        for (iac = 1; ldq < 0 ? iac >= idx : iac <= idx; iac += ldq) {
          c = 0.0;
          j = iac + qrmanager->mrows;
          for (ia = iac; ia < j; ia++) {
            c += qrmanager->Q.data[ia - 1] * objective->grad.data[ia - iac];
          }

          workspace_data[idxQR] += c;
          idxQR++;
        }
      }

      for (j = nActiveConstr_tmp; j >= 1; j--) {
        idx = (j + (j - 1) * ldq) - 1;
        workspace_data[j - 1] /= qrmanager->QR.data[idx];
        for (iac = 0; iac <= j - 2; iac++) {
          idxQR = (j - iac) - 2;
          workspace_data[idxQR] -= workspace_data[j - 1] * qrmanager->QR.data
            [(idx - iac) - 1];
        }
      }

      for (idx = 0; idx < nActiveConstr_tmp; idx++) {
        solution->lambda.data[idx] = -workspace_data[idx];
      }
    }
  }
}

/*
 * Arguments    : int x_data[]
 *                int xLen
 *                int workspace_data[]
 *                int xMin
 *                int xMax
 * Return Type  : void
 */
static void countsort(int x_data[], int xLen, int workspace_data[], int xMin,
                      int xMax)
{
  int idx;
  int idxFill;
  if ((xLen > 1) && (xMax > xMin)) {
    int i;
    int idxEnd;
    int idxStart;
    i = xMax - xMin;
    if (i >= 0) {
      memset(&workspace_data[0], 0, (unsigned int)(i + 1) * sizeof(int));
    }

    for (idx = 0; idx < xLen; idx++) {
      idxStart = x_data[idx] - xMin;
      workspace_data[idxStart]++;
    }

    for (idx = 2; idx <= i + 1; idx++) {
      workspace_data[idx - 1] += workspace_data[idx - 2];
    }

    idxStart = 1;
    idxEnd = workspace_data[0];
    for (idx = 0; idx < i; idx++) {
      for (idxFill = idxStart; idxFill <= idxEnd; idxFill++) {
        x_data[idxFill - 1] = idx + xMin;
      }

      idxStart = workspace_data[idx] + 1;
      idxEnd = workspace_data[idx + 1];
    }

    for (idx = idxStart; idx <= idxEnd; idx++) {
      x_data[idx - 1] = xMax;
    }
  }
}

/*
 * Arguments    : h_struct_T *obj
 *                const double x_data[]
 * Return Type  : double
 */
static double d_maxConstraintViolation_AMats_(h_struct_T *obj, const double
  x_data[])
{
  double v;
  int idx;
  int mIneq;
  v = 0.0;
  mIneq = obj->sizes[2];
  if (obj->Aineq.size[0] != 0) {
    if (mIneq - 1 >= 0) {
      memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0], (unsigned
              int)mIneq * sizeof(double));
    }

    xgemv(obj->nVar, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data,
          obj->maxConstrWorkspace.data);
    mIneq = (unsigned char)obj->sizes[2];
    for (idx = 0; idx < mIneq; idx++) {
      v = fmax(v, obj->maxConstrWorkspace.data[idx]);
    }
  }

  obj->maxConstrWorkspace.data[0] = obj->beq[0];
  obj->maxConstrWorkspace.data[1] = obj->beq[1];
  obj->maxConstrWorkspace.data[2] = obj->beq[2];
  xgemv(obj->nVar, 3, obj->Aeq.data, obj->ldA, x_data,
        obj->maxConstrWorkspace.data);
  v = fmax(v, fabs(obj->maxConstrWorkspace.data[0]));
  v = fmax(v, fabs(obj->maxConstrWorkspace.data[1]));
  return fmax(v, fabs(obj->maxConstrWorkspace.data[2]));
}

/*
 * Arguments    : const double x[3]
 * Return Type  : double
 */
static double d_xnrm2(const double x[3])
{
  double scale;
  double y;
  int k;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 2; k < 4; k++) {
    double absxk;
    absxk = fabs(x[k - 1]);
    if (absxk > scale) {
      double t;
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      double t;
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * sqrt(y);
}

/*
 * Arguments    : d_struct_T *obj
 *                int idx
 * Return Type  : void
 */
static void deleteColMoveEnd(d_struct_T *obj, int idx)
{
  double s;
  double temp;
  int b_k;
  int i;
  int k;
  if (obj->usedPivoting) {
    i = 1;
    while ((i <= obj->ncols) && (obj->jpvt.data[i - 1] != idx)) {
      i++;
    }

    idx = i;
  }

  if (idx >= obj->ncols) {
    obj->ncols--;
  } else {
    int b_i;
    int idxRotGCol;
    int ix;
    b_i = obj->ncols - 1;
    obj->jpvt.data[idx - 1] = obj->jpvt.data[b_i];
    idxRotGCol = obj->minRowCol;
    for (k = 0; k < idxRotGCol; k++) {
      obj->QR.data[k + obj->ldq * (idx - 1)] = obj->QR.data[k + obj->ldq * b_i];
    }

    obj->ncols = b_i;
    ix = obj->mrows;
    i = obj->ncols;
    if (ix <= i) {
      i = ix;
    }

    obj->minRowCol = i;
    if (idx < obj->mrows) {
      double c;
      int endIdx;
      int n;
      int temp_tmp;
      ix = obj->mrows - 1;
      endIdx = obj->ncols;
      if (ix <= endIdx) {
        endIdx = ix;
      }

      k = endIdx;
      idxRotGCol = obj->ldq * (idx - 1);
      while (k >= idx) {
        b_i = k + idxRotGCol;
        temp = obj->QR.data[b_i];
        c = xrotg(&obj->QR.data[b_i - 1], &temp, &s);
        obj->QR.data[b_i] = temp;
        b_i = obj->ldq * (k - 1);
        obj->QR.data[k + b_i] = 0.0;
        i = k + obj->ldq * idx;
        n = obj->ncols - idx;
        if (n >= 1) {
          ix = i - 1;
          for (b_k = 0; b_k < n; b_k++) {
            temp = c * obj->QR.data[ix] + s * obj->QR.data[i];
            obj->QR.data[i] = c * obj->QR.data[i] - s * obj->QR.data[ix];
            obj->QR.data[ix] = temp;
            i += obj->ldq;
            ix += obj->ldq;
          }
        }

        i = obj->ldq + b_i;
        n = obj->mrows;
        for (b_k = 0; b_k < n; b_k++) {
          ix = i + b_k;
          temp_tmp = b_i + b_k;
          temp = c * obj->Q.data[temp_tmp] + s * obj->Q.data[ix];
          obj->Q.data[ix] = c * obj->Q.data[ix] - s * obj->Q.data[temp_tmp];
          obj->Q.data[temp_tmp] = temp;
        }

        k--;
      }

      b_i = idx + 1;
      for (k = b_i; k <= endIdx; k++) {
        idxRotGCol = obj->ldq * (k - 1);
        i = k + idxRotGCol;
        temp = obj->QR.data[i];
        c = xrotg(&obj->QR.data[i - 1], &temp, &s);
        obj->QR.data[i] = temp;
        i = k * (obj->ldq + 1);
        n = obj->ncols - k;
        if (n >= 1) {
          ix = i - 1;
          for (b_k = 0; b_k < n; b_k++) {
            temp = c * obj->QR.data[ix] + s * obj->QR.data[i];
            obj->QR.data[i] = c * obj->QR.data[i] - s * obj->QR.data[ix];
            obj->QR.data[ix] = temp;
            i += obj->ldq;
            ix += obj->ldq;
          }
        }

        i = obj->ldq + idxRotGCol;
        n = obj->mrows;
        for (b_k = 0; b_k < n; b_k++) {
          ix = i + b_k;
          temp_tmp = idxRotGCol + b_k;
          temp = c * obj->Q.data[temp_tmp] + s * obj->Q.data[ix];
          obj->Q.data[ix] = c * obj->Q.data[ix] - s * obj->Q.data[temp_tmp];
          obj->Q.data[temp_tmp] = temp;
        }
      }
    }
  }
}

/*
 * Arguments    : const double v[3]
 *                double d[9]
 * Return Type  : void
 */
static void diag(const double v[3], double d[9])
{
  memset(&d[0], 0, 9U * sizeof(double));
  d[0] = v[0];
  d[4] = v[1];
  d[8] = v[2];
}

/*
 * Arguments    : const double bineq_data[]
 *                c_struct_T *TrialState
 *                j_struct_T *MeritFunction
 *                const i_coder_internal_stickyStruct *FcnEvaluator
 *                i_struct_T *FiniteDifferences
 *                g_struct_T *memspace
 *                h_struct_T *WorkingSet
 *                d_struct_T *QRManager
 *                e_struct_T *CholManager
 *                f_struct_T *QPObjective
 *                int fscales_lineq_constraint_size
 *                double Hessian[81]
 * Return Type  : void
 */
static void driver(const double bineq_data[], c_struct_T *TrialState, j_struct_T
                   *MeritFunction, const i_coder_internal_stickyStruct
                   *FcnEvaluator, i_struct_T *FiniteDifferences, g_struct_T
                   *memspace, h_struct_T *WorkingSet, d_struct_T *QRManager,
                   e_struct_T *CholManager, f_struct_T *QPObjective, int
                   fscales_lineq_constraint_size, double Hessian[81])
{
  static const signed char b_Hessian[81] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const char qpoptions_SolverName[7] = { 'f', 'm', 'i', 'n', 'c', 'o',
    'n' };

  k_struct_T expl_temp;
  struct_T Flags;
  double y_data[40];
  int i;
  int i1;
  int iCol;
  int ia;
  int ldJ_tmp;
  int mConstr;
  int mFixed;
  int mIneq;
  int mLB;
  int mUB;
  int nVar_tmp;
  int n_tmp;
  int qpoptions_MaxIterations;
  int u1;
  for (i = 0; i < 81; i++) {
    Hessian[i] = b_Hessian[i];
  }

  nVar_tmp = WorkingSet->nVar;
  mFixed = WorkingSet->sizes[0];
  mIneq = WorkingSet->sizes[2];
  mLB = WorkingSet->sizes[3];
  mUB = WorkingSet->sizes[4];
  mConstr = (((WorkingSet->sizes[0] + WorkingSet->sizes[2]) + WorkingSet->sizes
              [3]) + WorkingSet->sizes[4]) + 3;
  iCol = WorkingSet->nVar;
  u1 = ((WorkingSet->sizes[2] + WorkingSet->sizes[3]) + WorkingSet->sizes[4]) +
    (WorkingSet->sizes[0] << 1);
  if (iCol >= u1) {
    u1 = iCol;
  }

  qpoptions_MaxIterations = 10 * u1;
  TrialState->steplength = 1.0;
  Flags.gradOK = test_exit(MeritFunction, WorkingSet, TrialState, &Flags.fevalOK,
    &Flags.done, &Flags.stepAccepted, &Flags.failedLineSearch, &Flags.stepType);
  ldJ_tmp = WorkingSet->ldA;
  iCol = -1;
  i = (unsigned char)nVar_tmp;
  for (u1 = 0; u1 < 3; u1++) {
    memcpy(&TrialState->JacCeqTrans_old.data[iCol + 1], &WorkingSet->
           Aeq.data[iCol + 1], (unsigned int)i * sizeof(double));
    iCol += ldJ_tmp;
  }

  TrialState->sqpFval_old = TrialState->sqpFval;
  for (u1 = 0; u1 < 9; u1++) {
    TrialState->xstarsqp_old[u1] = TrialState->xstarsqp[u1];
    TrialState->grad_old.data[u1] = TrialState->grad.data[u1];
  }

  n_tmp = TrialState->mIneq;
  iCol = TrialState->cIneq_old.size[0];
  u1 = TrialState->cIneq_old.size[0];
  if (u1 - 1 >= 0) {
    memcpy(&y_data[0], &TrialState->cIneq_old.data[0], (unsigned int)u1 * sizeof
           (double));
  }

  if (n_tmp - 1 >= 0) {
    memcpy(&y_data[0], &TrialState->cIneq.data[0], (unsigned int)n_tmp * sizeof
           (double));
  }

  if (iCol - 1 >= 0) {
    memcpy(&TrialState->cIneq_old.data[0], &y_data[0], (unsigned int)iCol *
           sizeof(double));
  }

  TrialState->cEq_old[0] = TrialState->cEq[0];
  TrialState->cEq_old[1] = TrialState->cEq[1];
  TrialState->cEq_old[2] = TrialState->cEq[2];
  if (!Flags.done) {
    TrialState->sqpIterations = 1;
  }

  while (!Flags.done) {
    double phi_alpha;
    int i2;
    int ix;
    while (!(Flags.stepAccepted || Flags.failedLineSearch)) {
      if (Flags.stepType != 3) {
        updateWorkingSetForNewQP(TrialState->xstarsqp, WorkingSet, mIneq,
          TrialState->cIneq.data, TrialState->cEq, mLB, mUB, mFixed);
      }

      expl_temp.MaxIterations = qpoptions_MaxIterations;
      for (i1 = 0; i1 < 7; i1++) {
        expl_temp.SolverName[i1] = qpoptions_SolverName[i1];
      }

      step(&Flags, Hessian, TrialState, MeritFunction, memspace, WorkingSet,
           QRManager, CholManager, QPObjective, &expl_temp);
      if (Flags.stepAccepted) {
        for (u1 = 0; u1 < i; u1++) {
          TrialState->xstarsqp[u1] += TrialState->delta_x.data[u1];
        }

        TrialState->sqpFval = evalObjAndConstr
          (&FcnEvaluator->next.next.next.next.next.next.next.next.value.workspace,
           TrialState->xstarsqp, TrialState->cIneq.data, TrialState->iNonIneq0,
           TrialState->cEq, &iCol);
        Flags.fevalOK = (iCol == 1);
        TrialState->FunctionEvaluations++;
        computeLinearResiduals(TrialState->xstarsqp, nVar_tmp,
          TrialState->cIneq.data, &TrialState->cIneq.size[0], mIneq,
          WorkingSet->Aineq.data, bineq_data, WorkingSet->ldA);
        MeritFunction->phiFullStep = computeMeritFcn(MeritFunction->penaltyParam,
          TrialState->sqpFval, TrialState->cIneq.data, mIneq, TrialState->cEq,
          Flags.fevalOK);
      }

      if ((Flags.stepType == 1) && Flags.stepAccepted && Flags.fevalOK &&
          (MeritFunction->phi < MeritFunction->phiFullStep) &&
          (TrialState->sqpFval < TrialState->sqpFval_old)) {
        Flags.stepType = 3;
        Flags.stepAccepted = false;
      } else {
        double alpha;
        bool b;
        bool socTaken;
        if ((Flags.stepType == 3) && Flags.stepAccepted) {
          socTaken = true;
        } else {
          socTaken = false;
        }

        b = Flags.fevalOK;
        i1 = WorkingSet->nVar;
        alpha = 1.0;
        ix = 1;
        phi_alpha = MeritFunction->phiFullStep;
        if (i1 - 1 >= 0) {
          memcpy(&TrialState->searchDir.data[0], &TrialState->delta_x.data[0],
                 (unsigned int)i1 * sizeof(double));
        }

        int exitg1;
        do {
          exitg1 = 0;
          if (TrialState->FunctionEvaluations < 900) {
            if (b && (phi_alpha <= MeritFunction->phi + alpha * 0.0001 *
                      MeritFunction->phiPrimePlus)) {
              exitg1 = 1;
            } else {
              bool exitg2;
              bool tooSmallX;
              alpha *= 0.7;
              i2 = (unsigned char)i1;
              for (u1 = 0; u1 < i2; u1++) {
                TrialState->delta_x.data[u1] = alpha * TrialState->xstar.data[u1];
              }

              if (socTaken) {
                phi_alpha = alpha * alpha;
                iCol = TrialState->delta_x.size[0];
                u1 = TrialState->delta_x.size[0];
                if (u1 - 1 >= 0) {
                  memcpy(&y_data[0], &TrialState->delta_x.data[0], (unsigned int)
                         u1 * sizeof(double));
                }

                if ((i1 >= 1) && (!(phi_alpha == 0.0))) {
                  for (u1 = 0; u1 < i1; u1++) {
                    y_data[u1] += phi_alpha * TrialState->socDirection.data[u1];
                  }
                }

                if (iCol - 1 >= 0) {
                  memcpy(&TrialState->delta_x.data[0], &y_data[0], (unsigned int)
                         iCol * sizeof(double));
                }
              }

              tooSmallX = true;
              u1 = 0;
              exitg2 = false;
              while ((!exitg2) && (u1 <= (unsigned char)i1 - 1)) {
                if (1.0E-6 * fmax(1.0, fabs(TrialState->xstarsqp[u1])) <= fabs
                    (TrialState->delta_x.data[u1])) {
                  tooSmallX = false;
                  exitg2 = true;
                } else {
                  u1++;
                }
              }

              if (tooSmallX) {
                ix = -2;
                exitg1 = 1;
              } else {
                for (u1 = 0; u1 < i2; u1++) {
                  TrialState->xstarsqp[u1] = TrialState->xstarsqp_old[u1] +
                    TrialState->delta_x.data[u1];
                }

                TrialState->sqpFval = evalObjAndConstr
                  (&FcnEvaluator->next.next.next.next.next.next.next.next.value.workspace,
                   TrialState->xstarsqp, TrialState->cIneq.data,
                   TrialState->iNonIneq0, TrialState->cEq, &u1);
                computeLinearResiduals(TrialState->xstarsqp, i1,
                  TrialState->cIneq.data, &TrialState->cIneq.size[0], n_tmp,
                  WorkingSet->Aineq.data, bineq_data, WorkingSet->ldA);
                TrialState->FunctionEvaluations++;
                b = (u1 == 1);
                phi_alpha = computeMeritFcn(MeritFunction->penaltyParam,
                  TrialState->sqpFval, TrialState->cIneq.data, n_tmp,
                  TrialState->cEq, b);
              }
            }
          } else {
            ix = 0;
            exitg1 = 1;
          }
        } while (exitg1 == 0);

        Flags.fevalOK = b;
        TrialState->steplength = alpha;
        if (ix > 0) {
          Flags.stepAccepted = true;
        } else {
          Flags.failedLineSearch = true;
        }
      }
    }

    if (Flags.stepAccepted && (!Flags.failedLineSearch)) {
      for (u1 = 0; u1 < i; u1++) {
        TrialState->xstarsqp[u1] = TrialState->xstarsqp_old[u1] +
          TrialState->delta_x.data[u1];
      }

      i1 = (unsigned char)mConstr;
      for (u1 = 0; u1 < i1; u1++) {
        phi_alpha = TrialState->lambdasqp.data[u1];
        phi_alpha += TrialState->steplength * (TrialState->lambda.data[u1] -
          phi_alpha);
        TrialState->lambdasqp.data[u1] = phi_alpha;
      }

      TrialState->sqpFval_old = TrialState->sqpFval;
      for (u1 = 0; u1 < 9; u1++) {
        TrialState->xstarsqp_old[u1] = TrialState->xstarsqp[u1];
        TrialState->grad_old.data[u1] = TrialState->grad.data[u1];
      }

      iCol = TrialState->cIneq_old.size[0];
      u1 = TrialState->cIneq_old.size[0];
      if (u1 - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->cIneq_old.data[0], (unsigned int)u1 *
               sizeof(double));
      }

      if (n_tmp - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->cIneq.data[0], (unsigned int)n_tmp *
               sizeof(double));
      }

      if (iCol - 1 >= 0) {
        memcpy(&TrialState->cIneq_old.data[0], &y_data[0], (unsigned int)iCol *
               sizeof(double));
      }

      TrialState->cEq_old[0] = TrialState->cEq[0];
      TrialState->cEq_old[1] = TrialState->cEq[1];
      TrialState->cEq_old[2] = TrialState->cEq[2];
      Flags.gradOK = computeFiniteDifferences(FiniteDifferences,
        TrialState->sqpFval, TrialState->cEq, TrialState->xstarsqp,
        TrialState->grad.data, WorkingSet->Aeq.data, WorkingSet->ldA);
      TrialState->FunctionEvaluations += FiniteDifferences->numEvals;
    } else {
      TrialState->sqpFval = TrialState->sqpFval_old;
      memcpy(&TrialState->xstarsqp[0], &TrialState->xstarsqp_old[0], 9U * sizeof
             (double));
      iCol = TrialState->cIneq.size[0];
      u1 = TrialState->cIneq.size[0];
      if (u1 - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->cIneq.data[0], (unsigned int)u1 * sizeof
               (double));
      }

      if (n_tmp - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->cIneq_old.data[0], (unsigned int)n_tmp *
               sizeof(double));
      }

      if (iCol - 1 >= 0) {
        memcpy(&TrialState->cIneq.data[0], &y_data[0], (unsigned int)iCol *
               sizeof(double));
      }

      TrialState->cEq[0] = TrialState->cEq_old[0];
      TrialState->cEq[1] = TrialState->cEq_old[1];
      TrialState->cEq[2] = TrialState->cEq_old[2];
    }

    b_test_exit(&Flags, memspace, MeritFunction, fscales_lineq_constraint_size,
                WorkingSet, TrialState, QRManager);
    if ((!Flags.done) && Flags.stepAccepted) {
      Flags.stepAccepted = false;
      Flags.stepType = 1;
      Flags.failedLineSearch = false;
      memcpy(&TrialState->delta_gradLag.data[0], &TrialState->grad.data[0],
             (unsigned int)i * sizeof(double));
      iCol = TrialState->delta_gradLag.size[0];
      u1 = TrialState->delta_gradLag.size[0];
      if (u1 - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->delta_gradLag.data[0], (unsigned int)u1 *
               sizeof(double));
      }

      if (nVar_tmp >= 1) {
        for (u1 = 0; u1 < nVar_tmp; u1++) {
          y_data[u1] -= TrialState->grad_old.data[u1];
        }
      }

      ix = mFixed;
      if (iCol - 1 >= 0) {
        memcpy(&TrialState->delta_gradLag.data[0], &y_data[0], (unsigned int)
               iCol * sizeof(double));
      }

      i1 = (WorkingSet->ldA << 1) + 1;
      for (iCol = 1; ldJ_tmp < 0 ? iCol >= i1 : iCol <= i1; iCol += ldJ_tmp) {
        i2 = iCol + nVar_tmp;
        for (ia = iCol; ia < i2; ia++) {
          u1 = ia - iCol;
          TrialState->delta_gradLag.data[u1] += WorkingSet->Aeq.data[ia - 1] *
            TrialState->lambdasqp.data[ix];
        }

        ix++;
      }

      ix = mFixed;
      for (iCol = 1; ldJ_tmp < 0 ? iCol >= i1 : iCol <= i1; iCol += ldJ_tmp) {
        i2 = iCol + nVar_tmp;
        for (ia = iCol; ia < i2; ia++) {
          u1 = ia - iCol;
          TrialState->delta_gradLag.data[u1] += TrialState->
            JacCeqTrans_old.data[ia - 1] * -TrialState->lambdasqp.data[ix];
        }

        ix++;
      }

      iCol = -1;
      for (u1 = 0; u1 < 3; u1++) {
        memcpy(&TrialState->JacCeqTrans_old.data[iCol + 1],
               &WorkingSet->Aeq.data[iCol + 1], (unsigned int)i * sizeof(double));
        iCol += ldJ_tmp;
      }

      BFGSUpdate(nVar_tmp, Hessian, TrialState->delta_x.data,
                 TrialState->delta_gradLag.data, memspace->workspace_float.data);
      TrialState->sqpIterations++;
    }
  }
}

/*
 * Arguments    : const b_struct_T *c_obj_next_next_next_next_next_
 *                const double x[9]
 *                const double Cineq_workspace_data[]
 *                int ineq0
 *                double Ceq_workspace[3]
 *                int *status
 * Return Type  : double
 */
static double evalObjAndConstr(const b_struct_T *c_obj_next_next_next_next_next_,
  const double x[9], const double Cineq_workspace_data[], int ineq0, double
  Ceq_workspace[3], int *status)
{
  double fval;
  bool b;
  fval = optimize_forces_anonFcn1(c_obj_next_next_next_next_next_->center,
    c_obj_next_next_next_next_next_->A_alpha,
    c_obj_next_next_next_next_next_->b_alpha,
    c_obj_next_next_next_next_next_->deltaT,
    c_obj_next_next_next_next_next_->T_min, x);
  *status = 1;
  b = rtIsNaN(fval);
  if (rtIsInf(fval) || b) {
    if (b) {
      *status = -3;
    } else if (fval < 0.0) {
      *status = -1;
    } else {
      *status = -2;
    }
  }

  if (*status == 1) {
    *status = computeConstraints_(x, Cineq_workspace_data, ineq0, Ceq_workspace);
  }

  return fval;
}

/*
 * Arguments    : d_struct_T *obj
 *                const double A_data[]
 *                int mrows
 *                int ncols
 *                int ldA
 * Return Type  : void
 */
static void factorQR(d_struct_T *obj, const double A_data[], int mrows, int
                     ncols, int ldA)
{
  int i;
  int idx;
  int ix0;
  int k;
  int minmana;
  bool guard1;
  i = mrows * ncols;
  guard1 = false;
  if (i > 0) {
    for (idx = 0; idx < ncols; idx++) {
      ix0 = ldA * idx;
      minmana = obj->ldq * idx;
      i = (unsigned char)mrows;
      for (k = 0; k < i; k++) {
        obj->QR.data[minmana + k] = A_data[ix0 + k];
      }
    }

    guard1 = true;
  } else if (i == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    guard1 = true;
  }

  if (guard1) {
    obj->usedPivoting = false;
    obj->mrows = mrows;
    obj->ncols = ncols;
    for (idx = 0; idx < ncols; idx++) {
      obj->jpvt.data[idx] = idx + 1;
    }

    if (mrows <= ncols) {
      i = mrows;
    } else {
      i = ncols;
    }

    obj->minRowCol = i;
    ix0 = obj->QR.size[0];
    minmana = obj->QR.size[1];
    if (ix0 <= minmana) {
      minmana = ix0;
    }

    obj->tau.size[0] = minmana;
    if (minmana - 1 >= 0) {
      memset(&obj->tau.data[0], 0, (unsigned int)minmana * sizeof(double));
    }

    if (i >= 1) {
      qrf(obj->QR.data, obj->QR.size, mrows, ncols, i, obj->tau.data);
    }
  }
}

/*
 * Arguments    : d_struct_T *obj
 *                const double A_data[]
 *                int mrows
 *                int ncols
 *                int ldA
 * Return Type  : void
 */
static void factorQRE(d_struct_T *obj, const double A_data[], int mrows, int
                      ncols, int ldA)
{
  int idx;
  int k;
  int y;
  bool guard1;
  y = mrows * ncols;
  guard1 = false;
  if (y > 0) {
    for (idx = 0; idx < ncols; idx++) {
      int ix0;
      int iy0;
      ix0 = ldA * idx;
      iy0 = obj->ldq * idx;
      y = (unsigned char)mrows;
      for (k = 0; k < y; k++) {
        obj->QR.data[iy0 + k] = A_data[ix0 + k];
      }
    }

    guard1 = true;
  } else if (y == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    guard1 = true;
  }

  if (guard1) {
    obj->usedPivoting = true;
    obj->mrows = mrows;
    obj->ncols = ncols;
    if (mrows <= ncols) {
      y = mrows;
    } else {
      y = ncols;
    }

    obj->minRowCol = y;
    obj->tau.size[0] = xgeqp3(obj->QR.data, obj->QR.size, mrows, ncols,
      obj->jpvt.data, obj->tau.data);
  }
}

/*
 * Arguments    : int nVarMax
 *                int mConstrMax
 *                int mIneq
 *                c_struct_T *obj
 * Return Type  : void
 */
static void factoryConstruct(int nVarMax, int mConstrMax, int mIneq, c_struct_T *
  obj)
{
  obj->mIneq = mIneq;
  obj->iNonIneq0 = mIneq + 1;
  obj->sqpFval = 0.0;
  obj->sqpFval_old = 0.0;
  obj->cIneq.size[0] = mIneq;
  obj->cIneq_old.size[0] = mIneq;
  obj->grad.size[0] = nVarMax;
  obj->grad_old.size[0] = nVarMax;
  obj->FunctionEvaluations = 0;
  obj->sqpIterations = 0;
  obj->sqpExitFlag = 0;
  obj->lambdasqp.size[0] = mConstrMax;
  if (mConstrMax - 1 >= 0) {
    memset(&obj->lambdasqp.data[0], 0, (unsigned int)mConstrMax * sizeof(double));
  }

  obj->lambdaStopTest.size[0] = mConstrMax;
  obj->lambdaStopTestPrev.size[0] = mConstrMax;
  obj->steplength = 1.0;
  obj->delta_x.size[0] = nVarMax;
  if (nVarMax - 1 >= 0) {
    memset(&obj->delta_x.data[0], 0, (unsigned int)nVarMax * sizeof(double));
  }

  obj->socDirection.size[0] = nVarMax;
  obj->gradLag.size[0] = nVarMax;
  obj->delta_gradLag.size[0] = nVarMax;
  obj->xstar.size[0] = nVarMax;
  obj->fstar = 0.0;
  obj->lambda.size[0] = mConstrMax;
  if (mConstrMax - 1 >= 0) {
    memset(&obj->lambda.data[0], 0, (unsigned int)mConstrMax * sizeof(double));
  }

  obj->state = 0;
  obj->maxConstr = 0.0;
  obj->iterations = 0;
  obj->searchDir.size[0] = nVarMax;
}

/*
 * Arguments    : double workspace_data[]
 *                const int workspace_size[2]
 *                double xCurrent_data[]
 *                h_struct_T *workingset
 *                d_struct_T *qrmanager
 * Return Type  : bool
 */
static bool feasibleX0ForWorkingSet(double workspace_data[], const int
  workspace_size[2], double xCurrent_data[], h_struct_T *workingset, d_struct_T *
  qrmanager)
{
  double B_data[880];
  int b_i;
  int ix;
  int j;
  int jBcol;
  int k;
  int mWConstr;
  int minmn;
  int nVar;
  bool nonDegenerateWset;
  mWConstr = workingset->nActiveConstr;
  nVar = workingset->nVar;
  nonDegenerateWset = true;
  if (mWConstr != 0) {
    double d;
    double tol;
    int i;
    if (mWConstr >= nVar) {
      int i1;
      int i2;
      int ldq;
      int ldw_tmp;
      int rankQR;
      i = (unsigned char)nVar;
      for (minmn = 0; minmn < i; minmn++) {
        ix = qrmanager->ldq * minmn;
        for (jBcol = 0; jBcol < mWConstr; jBcol++) {
          qrmanager->QR.data[jBcol + ix] = workingset->ATwset.data[minmn +
            workingset->ldA * jBcol];
        }

        qrmanager->jpvt.data[minmn] = 0;
      }

      if (mWConstr * nVar == 0) {
        qrmanager->mrows = mWConstr;
        qrmanager->ncols = nVar;
        qrmanager->minRowCol = 0;
      } else {
        qrmanager->usedPivoting = true;
        qrmanager->mrows = mWConstr;
        qrmanager->ncols = nVar;
        if (mWConstr <= nVar) {
          ix = mWConstr;
        } else {
          ix = nVar;
        }

        qrmanager->minRowCol = ix;
        qrmanager->tau.size[0] = xgeqp3(qrmanager->QR.data, qrmanager->QR.size,
          mWConstr, nVar, qrmanager->jpvt.data, qrmanager->tau.data);
      }

      computeQ_(qrmanager, qrmanager->mrows);
      rankQR = 0;
      ix = qrmanager->mrows;
      minmn = qrmanager->ncols;
      if (ix <= minmn) {
        minmn = ix;
      }

      if (minmn > 0) {
        jBcol = qrmanager->ncols;
        if (ix >= jBcol) {
          jBcol = ix;
        }

        tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)jBcol)
          * fabs(qrmanager->QR.data[0]);
        while ((rankQR < minmn) && (!(fabs(qrmanager->QR.data[rankQR +
                  qrmanager->QR.size[0] * rankQR]) <= tol))) {
          rankQR++;
        }
      }

      for (minmn = 0; minmn < mWConstr; minmn++) {
        d = workingset->bwset.data[minmn];
        workspace_data[minmn] = d;
        workspace_data[minmn + workspace_size[0]] = d;
      }

      minmn = workingset->ldA;
      jBcol = 0;
      i1 = workingset->ldA * (mWConstr - 1) + 1;
      for (ix = 1; minmn < 0 ? ix >= i1 : ix <= i1; ix += minmn) {
        tol = 0.0;
        i2 = ix + nVar;
        for (j = ix; j < i2; j++) {
          tol += workingset->ATwset.data[j - 1] * xCurrent_data[j - ix];
        }

        workspace_data[jBcol] -= tol;
        jBcol++;
      }

      ldq = qrmanager->ldq;
      ldw_tmp = workspace_size[0];
      minmn = workspace_size[0] * workspace_size[1];
      if (minmn - 1 >= 0) {
        memcpy(&B_data[0], &workspace_data[0], (unsigned int)minmn * sizeof
               (double));
      }

      for (j = 0; ldw_tmp < 0 ? j >= ldw_tmp : j <= ldw_tmp; j += ldw_tmp) {
        i1 = j + 1;
        i2 = j + nVar;
        if (i1 <= i2) {
          memset(&workspace_data[i1 + -1], 0, (unsigned int)((i2 - i1) + 1) *
                 sizeof(double));
        }
      }

      minmn = -1;
      for (j = 0; ldw_tmp < 0 ? j >= ldw_tmp : j <= ldw_tmp; j += ldw_tmp) {
        ix = -1;
        i1 = j + 1;
        i2 = j + nVar;
        for (b_i = i1; b_i <= i2; b_i++) {
          tol = 0.0;
          for (jBcol = 0; jBcol < mWConstr; jBcol++) {
            tol += qrmanager->Q.data[(jBcol + ix) + 1] * B_data[(jBcol + minmn)
              + 1];
          }

          workspace_data[b_i - 1] += tol;
          ix += ldq;
        }

        minmn += ldw_tmp;
      }

      for (j = 0; j < 2; j++) {
        jBcol = ldw_tmp * j - 1;
        for (k = rankQR; k >= 1; k--) {
          minmn = ldq * (k - 1) - 1;
          i1 = k + jBcol;
          d = workspace_data[i1];
          if (d != 0.0) {
            workspace_data[i1] = d / qrmanager->QR.data[k + minmn];
            i2 = (unsigned char)(k - 1);
            for (b_i = 0; b_i < i2; b_i++) {
              ix = (b_i + jBcol) + 1;
              workspace_data[ix] -= workspace_data[i1] * qrmanager->QR.data[(b_i
                + minmn) + 1];
            }
          }
        }
      }

      i1 = rankQR + 1;
      for (b_i = i1; b_i <= nVar; b_i++) {
        workspace_data[b_i - 1] = 0.0;
        workspace_data[(b_i + workspace_size[0]) - 1] = 0.0;
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace_data[(qrmanager->jpvt.data[b_i] + workspace_size[0] * 2) - 1] =
          workspace_data[b_i];
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace_data[b_i] = workspace_data[b_i + workspace_size[0] * 2];
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace_data[(qrmanager->jpvt.data[b_i] + workspace_size[0] * 2) - 1] =
          workspace_data[b_i + workspace_size[0]];
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace_data[b_i + workspace_size[0]] = workspace_data[b_i +
          workspace_size[0] * 2];
      }
    } else {
      int i1;
      int ldq;
      int ldw_tmp;
      int rankQR;
      if (mWConstr - 1 >= 0) {
        memset(&qrmanager->jpvt.data[0], 0, (unsigned int)mWConstr * sizeof(int));
      }

      factorQRE(qrmanager, workingset->ATwset.data, nVar, mWConstr,
                workingset->ldA);
      computeQ_(qrmanager, qrmanager->minRowCol);
      rankQR = 0;
      ix = qrmanager->mrows;
      minmn = qrmanager->ncols;
      if (ix <= minmn) {
        minmn = ix;
      }

      if (minmn > 0) {
        jBcol = qrmanager->ncols;
        if (ix >= jBcol) {
          jBcol = ix;
        }

        tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)jBcol)
          * fabs(qrmanager->QR.data[0]);
        while ((rankQR < minmn) && (!(fabs(qrmanager->QR.data[rankQR +
                  qrmanager->QR.size[0] * rankQR]) <= tol))) {
          rankQR++;
        }
      }

      for (minmn = 0; minmn < mWConstr; minmn++) {
        ix = (qrmanager->jpvt.data[minmn] - 1) * workingset->ldA;
        tol = 0.0;
        i = (unsigned char)nVar;
        for (k = 0; k < i; k++) {
          tol += workingset->ATwset.data[ix + k] * xCurrent_data[k];
        }

        d = workingset->bwset.data[qrmanager->jpvt.data[minmn] - 1];
        workspace_data[minmn] = d - tol;
        workspace_data[minmn + workspace_size[0]] = d;
      }

      ldq = qrmanager->ldq;
      ldw_tmp = workspace_size[0];
      i = (unsigned char)rankQR;
      for (j = 0; j < 2; j++) {
        jBcol = ldw_tmp * j;
        for (b_i = 0; b_i < i; b_i++) {
          minmn = ldq * b_i;
          ix = b_i + jBcol;
          tol = workspace_data[ix];
          for (k = 0; k < b_i; k++) {
            tol -= qrmanager->QR.data[k + minmn] * workspace_data[k + jBcol];
          }

          workspace_data[ix] = tol / qrmanager->QR.data[b_i + minmn];
        }
      }

      minmn = workspace_size[0] * workspace_size[1];
      if (minmn - 1 >= 0) {
        memcpy(&B_data[0], &workspace_data[0], (unsigned int)minmn * sizeof
               (double));
      }

      for (j = 0; ldw_tmp < 0 ? j >= ldw_tmp : j <= ldw_tmp; j += ldw_tmp) {
        i = j + 1;
        i1 = j + nVar;
        if (i <= i1) {
          memset(&workspace_data[i + -1], 0, (unsigned int)((i1 - i) + 1) *
                 sizeof(double));
        }
      }

      minmn = 1;
      for (j = 0; ldw_tmp < 0 ? j >= ldw_tmp : j <= ldw_tmp; j += ldw_tmp) {
        ix = -1;
        i = minmn + rankQR;
        for (jBcol = minmn; jBcol < i; jBcol++) {
          int i2;
          i1 = j + 1;
          i2 = j + nVar;
          for (b_i = i1; b_i <= i2; b_i++) {
            workspace_data[b_i - 1] += B_data[jBcol - 1] * qrmanager->Q.data[(ix
              + b_i) - j];
          }

          ix += ldq;
        }

        minmn += ldw_tmp;
      }
    }

    minmn = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (minmn <= (unsigned char)nVar - 1) {
        tol = workspace_data[minmn];
        if (rtIsInf(tol) || rtIsNaN(tol)) {
          nonDegenerateWset = false;
          exitg1 = 1;
        } else {
          tol = workspace_data[minmn + workspace_size[0]];
          if (rtIsInf(tol) || rtIsNaN(tol)) {
            nonDegenerateWset = false;
            exitg1 = 1;
          } else {
            minmn++;
          }
        }
      } else {
        double v;
        for (k = 0; k < nVar; k++) {
          workspace_data[k] += xCurrent_data[k];
        }

        if (workingset->probType == 2) {
          v = 0.0;
          ix = workingset->sizes[2];
          if (workingset->Aineq.size[0] != 0) {
            if (ix - 1 >= 0) {
              memcpy(&workingset->maxConstrWorkspace.data[0],
                     &workingset->bineq.data[0], (unsigned int)ix * sizeof
                     (double));
            }

            xgemv(9, workingset->sizes[2], workingset->Aineq.data,
                  workingset->ldA, workspace_data,
                  workingset->maxConstrWorkspace.data);
            i = (unsigned char)workingset->sizes[2];
            for (minmn = 0; minmn < i; minmn++) {
              d = workingset->maxConstrWorkspace.data[minmn] -
                workspace_data[minmn + 9];
              workingset->maxConstrWorkspace.data[minmn] = d;
              v = fmax(v, d);
            }
          }

          workingset->maxConstrWorkspace.data[0] = workingset->beq[0];
          workingset->maxConstrWorkspace.data[1] = workingset->beq[1];
          workingset->maxConstrWorkspace.data[2] = workingset->beq[2];
          xgemv(9, 3, workingset->Aeq.data, workingset->ldA, workspace_data,
                workingset->maxConstrWorkspace.data);
          workingset->maxConstrWorkspace.data[0] =
            (workingset->maxConstrWorkspace.data[0] - workspace_data[ix + 9]) +
            workspace_data[workingset->sizes[2] + 12];
          v = fmax(v, fabs(workingset->maxConstrWorkspace.data[0]));
          workingset->maxConstrWorkspace.data[1] =
            (workingset->maxConstrWorkspace.data[1] - workspace_data[ix + 10]) +
            workspace_data[workingset->sizes[2] + 13];
          v = fmax(v, fabs(workingset->maxConstrWorkspace.data[1]));
          workingset->maxConstrWorkspace.data[2] =
            (workingset->maxConstrWorkspace.data[2] - workspace_data[ix + 11]) +
            workspace_data[workingset->sizes[2] + 14];
          v = fmax(v, fabs(workingset->maxConstrWorkspace.data[2]));
        } else {
          v = 0.0;
          ix = workingset->sizes[2];
          if (workingset->Aineq.size[0] != 0) {
            if (ix - 1 >= 0) {
              memcpy(&workingset->maxConstrWorkspace.data[0],
                     &workingset->bineq.data[0], (unsigned int)ix * sizeof
                     (double));
            }

            xgemv(workingset->nVar, workingset->sizes[2], workingset->Aineq.data,
                  workingset->ldA, workspace_data,
                  workingset->maxConstrWorkspace.data);
            i = (unsigned char)workingset->sizes[2];
            for (minmn = 0; minmn < i; minmn++) {
              v = fmax(v, workingset->maxConstrWorkspace.data[minmn]);
            }
          }

          workingset->maxConstrWorkspace.data[0] = workingset->beq[0];
          workingset->maxConstrWorkspace.data[1] = workingset->beq[1];
          workingset->maxConstrWorkspace.data[2] = workingset->beq[2];
          xgemv(workingset->nVar, 3, workingset->Aeq.data, workingset->ldA,
                workspace_data, workingset->maxConstrWorkspace.data);
          v = fmax(v, fabs(workingset->maxConstrWorkspace.data[0]));
          v = fmax(v, fabs(workingset->maxConstrWorkspace.data[1]));
          v = fmax(v, fabs(workingset->maxConstrWorkspace.data[2]));
        }

        if (workingset->sizes[3] > 0) {
          i = (unsigned char)workingset->sizes[3];
          for (minmn = 0; minmn < i; minmn++) {
            ix = workingset->indexLB.data[minmn] - 1;
            v = fmax(v, -workspace_data[ix] - workingset->lb.data[ix]);
          }
        }

        if (workingset->sizes[4] > 0) {
          i = (unsigned char)workingset->sizes[4];
          for (minmn = 0; minmn < i; minmn++) {
            ix = workingset->indexUB.data[minmn] - 1;
            v = fmax(v, workspace_data[ix] - workingset->ub.data[ix]);
          }
        }

        if (workingset->sizes[0] > 0) {
          i = (unsigned char)workingset->sizes[0];
          for (minmn = 0; minmn < i; minmn++) {
            v = fmax(v, fabs(workspace_data[workingset->indexFixed.data[minmn] -
                             1] - workingset->ub.data
                             [workingset->indexFixed.data[minmn] - 1]));
          }
        }

        tol = maxConstraintViolation(workingset, workspace_data, workspace_size
          [0] + 1);
        if ((v <= 2.2204460492503131E-16) || (v < tol)) {
          i = (unsigned char)nVar;
          memcpy(&xCurrent_data[0], &workspace_data[0], (unsigned int)i * sizeof
                 (double));
        } else {
          i = (unsigned char)nVar;
          for (k = 0; k < i; k++) {
            xCurrent_data[k] = workspace_data[workspace_size[0] + k];
          }
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return nonDegenerateWset;
}

/*
 * Arguments    : const double solution_xstar_data[]
 *                const double solution_searchDir_data[]
 *                double workspace_data[]
 *                const int workspace_size[2]
 *                int workingset_nVar
 *                int workingset_ldA
 *                const double workingset_Aineq_data[]
 *                const double workingset_bineq_data[]
 *                const double workingset_lb_data[]
 *                const double workingset_ub_data[]
 *                const int workingset_indexLB_data[]
 *                const int workingset_indexUB_data[]
 *                const int workingset_sizes[5]
 *                const int workingset_isActiveIdx[6]
 *                const bool workingset_isActiveConstr_data[]
 *                const int workingset_nWConstr[5]
 *                bool isPhaseOne
 *                bool *newBlocking
 *                int *constrType
 *                int *constrIdx
 * Return Type  : double
 */
static double feasibleratiotest(const double solution_xstar_data[], const double
  solution_searchDir_data[], double workspace_data[], const int workspace_size[2],
  int workingset_nVar, int workingset_ldA, const double workingset_Aineq_data[],
  const double workingset_bineq_data[], const double workingset_lb_data[], const
  double workingset_ub_data[], const int workingset_indexLB_data[], const int
  workingset_indexUB_data[], const int workingset_sizes[5], const int
  workingset_isActiveIdx[6], const bool workingset_isActiveConstr_data[], const
  int workingset_nWConstr[5], bool isPhaseOne, bool *newBlocking, int
  *constrType, int *constrIdx)
{
  double alpha;
  double c;
  double denomTol;
  double phaseOneCorrectionP;
  double phaseOneCorrectionX;
  double ratio;
  int i;
  int ia;
  int iac;
  int iyend;
  alpha = 1.0E+30;
  *newBlocking = false;
  *constrType = 0;
  *constrIdx = 0;
  denomTol = 2.2204460492503131E-13 * c_xnrm2(workingset_nVar,
    solution_searchDir_data);
  if (workingset_nWConstr[2] < workingset_sizes[2]) {
    int iy;
    int ldw_tmp;
    i = (unsigned char)workingset_sizes[2];
    if (i - 1 >= 0) {
      memcpy(&workspace_data[0], &workingset_bineq_data[0], (unsigned int)i *
             sizeof(double));
    }

    xgemv(workingset_nVar, workingset_sizes[2], workingset_Aineq_data,
          workingset_ldA, solution_xstar_data, workspace_data);
    ldw_tmp = workspace_size[0];
    iy = workspace_size[0] + 1;
    if (workingset_sizes[2] != 0) {
      iyend = workspace_size[0] + workingset_sizes[2];
      if (iy <= iyend) {
        memset(&workspace_data[iy + -1], 0, (unsigned int)((iyend - iy) + 1) *
               sizeof(double));
      }

      iy = ldw_tmp;
      iyend = workingset_ldA * (workingset_sizes[2] - 1) + 1;
      for (iac = 1; workingset_ldA < 0 ? iac >= iyend : iac <= iyend; iac +=
           workingset_ldA) {
        int i1;
        c = 0.0;
        i1 = iac + workingset_nVar;
        for (ia = iac; ia < i1; ia++) {
          c += workingset_Aineq_data[ia - 1] * solution_searchDir_data[ia - iac];
        }

        workspace_data[iy] += c;
        iy++;
      }
    }

    for (iyend = 0; iyend < i; iyend++) {
      c = workspace_data[ldw_tmp + iyend];
      if ((c > denomTol) && (!workingset_isActiveConstr_data
                             [(workingset_isActiveIdx[2] + iyend) - 1])) {
        c = fmin(fabs(workspace_data[iyend]), 1.0E-6 - workspace_data[iyend]) /
          c;
        if (c < alpha) {
          alpha = c;
          *constrType = 3;
          *constrIdx = iyend + 1;
          *newBlocking = true;
        }
      }
    }
  }

  if (workingset_nWConstr[3] < workingset_sizes[3]) {
    phaseOneCorrectionX = (double)isPhaseOne *
      solution_xstar_data[workingset_nVar - 1];
    phaseOneCorrectionP = (double)isPhaseOne *
      solution_searchDir_data[workingset_nVar - 1];
    i = workingset_sizes[3];
    for (iyend = 0; iyend <= i - 2; iyend++) {
      c = -solution_searchDir_data[workingset_indexLB_data[iyend] - 1] -
        phaseOneCorrectionP;
      if ((c > denomTol) && (!workingset_isActiveConstr_data
                             [(workingset_isActiveIdx[3] + iyend) - 1])) {
        ratio = (-solution_xstar_data[workingset_indexLB_data[iyend] - 1] -
                 workingset_lb_data[workingset_indexLB_data[iyend] - 1]) -
          phaseOneCorrectionX;
        c = fmin(fabs(ratio), 1.0E-6 - ratio) / c;
        if (c < alpha) {
          alpha = c;
          *constrType = 4;
          *constrIdx = iyend + 1;
          *newBlocking = true;
        }
      }
    }

    i = workingset_indexLB_data[workingset_sizes[3] - 1] - 1;
    c = -solution_searchDir_data[i];
    if ((c > denomTol) && (!workingset_isActiveConstr_data
                           [(workingset_isActiveIdx[3] + workingset_sizes[3]) -
                           2])) {
      ratio = -solution_xstar_data[i] - workingset_lb_data[i];
      c = fmin(fabs(ratio), 1.0E-6 - ratio) / c;
      if (c < alpha) {
        alpha = c;
        *constrType = 4;
        *constrIdx = workingset_sizes[3];
        *newBlocking = true;
      }
    }
  }

  if (workingset_nWConstr[4] < workingset_sizes[4]) {
    phaseOneCorrectionX = (double)isPhaseOne *
      solution_xstar_data[workingset_nVar - 1];
    phaseOneCorrectionP = (double)isPhaseOne *
      solution_searchDir_data[workingset_nVar - 1];
    i = (unsigned char)workingset_sizes[4];
    for (iyend = 0; iyend < i; iyend++) {
      c = solution_searchDir_data[workingset_indexUB_data[iyend] - 1] -
        phaseOneCorrectionP;
      if ((c > denomTol) && (!workingset_isActiveConstr_data
                             [(workingset_isActiveIdx[4] + iyend) - 1])) {
        ratio = (solution_xstar_data[workingset_indexUB_data[iyend] - 1] -
                 workingset_ub_data[workingset_indexUB_data[iyend] - 1]) -
          phaseOneCorrectionX;
        c = fmin(fabs(ratio), 1.0E-6 - ratio) / c;
        if (c < alpha) {
          alpha = c;
          *constrType = 5;
          *constrIdx = iyend + 1;
          *newBlocking = true;
        }
      }
    }
  }

  if (!isPhaseOne) {
    if ((*newBlocking) && (alpha > 1.0)) {
      *newBlocking = false;
    }

    alpha = fmin(alpha, 1.0);
  }

  return alpha;
}

/*
 * Arguments    : const double obj_objfun_workspace_center[3]
 *                const double obj_objfun_workspace_A_alpha[18]
 *                const double obj_objfun_workspace_b_alpha[6]
 *                const double obj_objfun_workspace_deltaT[3]
 *                const double obj_objfun_workspace_T_min[3]
 *                double cEqPlus[3]
 *                int dim
 *                double delta
 *                double xk[9]
 *                double *fplus
 * Return Type  : bool
 */
static bool finDiffEvalAndChkErr(const double obj_objfun_workspace_center[3],
  const double obj_objfun_workspace_A_alpha[18], const double
  obj_objfun_workspace_b_alpha[6], const double obj_objfun_workspace_deltaT[3],
  const double obj_objfun_workspace_T_min[3], double cEqPlus[3], int dim, double
  delta, double xk[9], double *fplus)
{
  double temp_tmp;
  int idx;
  bool evalOK;
  temp_tmp = xk[dim - 1];
  xk[dim - 1] = temp_tmp + delta;
  *fplus = optimize_forces_anonFcn1(obj_objfun_workspace_center,
    obj_objfun_workspace_A_alpha, obj_objfun_workspace_b_alpha,
    obj_objfun_workspace_deltaT, obj_objfun_workspace_T_min, xk);
  evalOK = ((!rtIsInf(*fplus)) && (!rtIsNaN(*fplus)));
  if (evalOK) {
    /* 'optimize_forces:71' @(W) constraint(W) */
    /*  优化约束：单位向量约束 */
    /* 'optimize_forces:31' ceq = zeros([length(W) / 3, 1]); */
    /* 'optimize_forces:32' for idx = 1:length(W) / 3 */
    for (idx = 0; idx < 3; idx++) {
      double absxk;
      double scale;
      double t;
      double y;
      int i;

      /* 'optimize_forces:33' ceq(idx) = norm(W(3*idx-2:3*idx)) - 1; */
      scale = 3.3121686421112381E-170;
      i = 3 * (idx + 1) - 3;
      absxk = fabs(xk[i]);
      if (absxk > 3.3121686421112381E-170) {
        y = 1.0;
        scale = absxk;
      } else {
        t = absxk / 3.3121686421112381E-170;
        y = t * t;
      }

      absxk = fabs(xk[i + 1]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      absxk = fabs(xk[i + 2]);
      if (absxk > scale) {
        t = scale / absxk;
        y = y * t * t + 1.0;
        scale = absxk;
      } else {
        t = absxk / scale;
        y += t * t;
      }

      y = scale * sqrt(y);
      cEqPlus[idx] = y - 1.0;
    }

    /* 'optimize_forces:35' c = []; */
    idx = 0;
    while (evalOK && (idx + 1 <= 3)) {
      evalOK = ((!rtIsInf(cEqPlus[idx])) && (!rtIsNaN(cEqPlus[idx])));
      idx++;
    }

    xk[dim - 1] = temp_tmp;
  }

  return evalOK;
}

/*
 * Arguments    : const double fun_workspace_center[3]
 *                const double fun_workspace_A_alpha[18]
 *                const double fun_workspace_b_alpha[6]
 *                const double fun_workspace_deltaT[3]
 *                const double fun_workspace_T_min[3]
 *                const double x0[9]
 *                const double Aineq_data[]
 *                const double bineq_data[]
 *                const int bineq_size[2]
 *                double x[9]
 *                double *exitflag
 *                double *output_iterations
 *                double *output_funcCount
 *                char output_algorithm[3]
 *                double *output_constrviolation
 *                double *output_stepsize
 *                double *output_lssteplength
 *                double *output_firstorderopt
 * Return Type  : double
 */
static double fmincon(const double fun_workspace_center[3], const double
                      fun_workspace_A_alpha[18], const double
                      fun_workspace_b_alpha[6], const double
                      fun_workspace_deltaT[3], const double fun_workspace_T_min
                      [3], const double x0[9], const double Aineq_data[], const
                      double bineq_data[], const int bineq_size[2], double x[9],
                      double *exitflag, double *output_iterations, double
                      *output_funcCount, char output_algorithm[3], double
                      *output_constrviolation, double *output_stepsize, double
                      *output_lssteplength, double *output_firstorderopt)
{
  c_struct_T TrialState;
  d_struct_T QRManager;
  e_struct_T CholManager;
  f_struct_T QPObjective;
  g_struct_T memspace;
  h_struct_T WorkingSet;
  i_coder_internal_stickyStruct FcnEvaluator;
  i_struct_T obj;
  j_struct_T MeritFunction;
  double unusedExpr[81];
  double y_data[40];
  double c;
  double fval;
  double scale;
  int b_i;
  int i;
  int ia;
  int iac;
  int ldAi;
  int mConstrMax;
  int mLinIneq_tmp;
  signed char obj_tmp[5];
  bool b;
  output_algorithm[0] = 's';
  output_algorithm[1] = 'q';
  output_algorithm[2] = 'p';
  mLinIneq_tmp = bineq_size[0] * bineq_size[1];
  mConstrMax = (mLinIneq_tmp + mLinIneq_tmp) + 28;
  factoryConstruct(mLinIneq_tmp + 16, mConstrMax, mLinIneq_tmp, &TrialState);
  memcpy(&TrialState.xstarsqp[0], &x0[0], 9U * sizeof(double));
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.center[0]
    = fun_workspace_center[0];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.center[1]
    = fun_workspace_center[1];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.center[2]
    = fun_workspace_center[2];
  memcpy
    (&FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.A_alpha
     [0], &fun_workspace_A_alpha[0], 18U * sizeof(double));
  for (i = 0; i < 6; i++) {
    FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.b_alpha[
      i] = fun_workspace_b_alpha[i];
  }

  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.deltaT[0]
    = fun_workspace_deltaT[0];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_min[0] =
    fun_workspace_T_min[0];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.deltaT[1]
    = fun_workspace_deltaT[1];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_min[1] =
    fun_workspace_T_min[1];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.deltaT[2]
    = fun_workspace_deltaT[2];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_min[2] =
    fun_workspace_T_min[2];
  QRManager.ldq = mConstrMax;
  QRManager.QR.size[0] = mConstrMax;
  QRManager.QR.size[1] = mConstrMax;
  QRManager.Q.size[0] = mConstrMax;
  QRManager.Q.size[1] = mConstrMax;
  i = mConstrMax * mConstrMax;
  memset(&QRManager.Q.data[0], 0, (unsigned int)i * sizeof(double));
  memset(&QRManager.jpvt.data[0], 0, (unsigned int)mConstrMax * sizeof(int));
  QRManager.mrows = 0;
  QRManager.ncols = 0;
  QRManager.tau.size[0] = mConstrMax;
  QRManager.minRowCol = 0;
  QRManager.usedPivoting = false;
  CholManager.FMat.size[0] = mConstrMax;
  CholManager.FMat.size[1] = mConstrMax;
  CholManager.ldm = mConstrMax;
  CholManager.ndims = 0;
  CholManager.info = 0;
  CholManager.ConvexCheck = true;
  CholManager.workspace_ = rtInf;
  QPObjective.grad.size[0] = mLinIneq_tmp + 16;
  QPObjective.maxVar = mLinIneq_tmp + 16;
  QPObjective.beta = 0.0;
  QPObjective.rho = 0.0;
  QPObjective.prev_objtype = 3;
  QPObjective.prev_nvar = 0;
  QPObjective.prev_hasLinear = false;
  QPObjective.gammaScalar = 0.0;
  QPObjective.hasLinear = true;
  QPObjective.nvar = 9;
  QPObjective.objtype = 3;
  memspace.workspace_float.size[0] = mConstrMax;
  memspace.workspace_float.size[1] = mLinIneq_tmp + 16;
  b_factoryConstruct(mLinIneq_tmp, mLinIneq_tmp + 16, mConstrMax, &WorkingSet);
  for (i = 0; i < 9; i++) {
    WorkingSet.indexLB.data[i] = i + 1;
    WorkingSet.indexUB.data[i] = i + 1;
  }

  WorkingSet.mConstr = mLinIneq_tmp + 21;
  WorkingSet.mConstrOrig = mLinIneq_tmp + 21;
  WorkingSet.mConstrMax = mConstrMax;
  obj_tmp[0] = 0;
  obj_tmp[1] = 3;
  obj_tmp[2] = (signed char)mLinIneq_tmp;
  obj_tmp[3] = 9;
  obj_tmp[4] = 9;
  WorkingSet.sizesPhaseOne[0] = 0;
  WorkingSet.sizesPhaseOne[1] = 3;
  WorkingSet.sizesPhaseOne[2] = mLinIneq_tmp;
  WorkingSet.sizesPhaseOne[3] = 10;
  WorkingSet.sizesPhaseOne[4] = 9;
  WorkingSet.sizesRegularized[0] = 0;
  WorkingSet.sizesRegularized[1] = 3;
  WorkingSet.sizesRegularized[2] = mLinIneq_tmp;
  WorkingSet.sizesRegularized[3] = mLinIneq_tmp + 15;
  WorkingSet.sizesRegularized[4] = 9;
  WorkingSet.sizesRegPhaseOne[0] = 0;
  WorkingSet.sizesRegPhaseOne[1] = 3;
  WorkingSet.sizesRegPhaseOne[2] = mLinIneq_tmp;
  WorkingSet.sizesRegPhaseOne[3] = mLinIneq_tmp + 16;
  WorkingSet.sizesRegPhaseOne[4] = 9;
  WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  WorkingSet.isActiveIdxRegPhaseOne[1] = 0;
  WorkingSet.isActiveIdxRegPhaseOne[2] = 3;
  WorkingSet.isActiveIdxRegPhaseOne[3] = mLinIneq_tmp;
  WorkingSet.isActiveIdxRegPhaseOne[4] = 9;
  WorkingSet.isActiveIdxRegPhaseOne[5] = 9;
  for (i = 0; i < 5; i++) {
    signed char i1;
    i1 = obj_tmp[i];
    WorkingSet.sizes[i] = i1;
    WorkingSet.sizesNormal[i] = i1;
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (b_i = 0; b_i < 6; b_i++) {
    WorkingSet.isActiveIdx[b_i] = WorkingSet.isActiveIdxRegPhaseOne[b_i];
    WorkingSet.isActiveIdxNormal[b_i] = WorkingSet.isActiveIdxRegPhaseOne[b_i];
  }

  WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  WorkingSet.isActiveIdxRegPhaseOne[1] = 0;
  WorkingSet.isActiveIdxRegPhaseOne[2] = 3;
  WorkingSet.isActiveIdxRegPhaseOne[3] = mLinIneq_tmp;
  WorkingSet.isActiveIdxRegPhaseOne[4] = 10;
  WorkingSet.isActiveIdxRegPhaseOne[5] = 9;
  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (b_i = 0; b_i < 6; b_i++) {
    WorkingSet.isActiveIdxPhaseOne[b_i] = WorkingSet.isActiveIdxRegPhaseOne[b_i];
  }

  WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  WorkingSet.isActiveIdxRegPhaseOne[1] = 0;
  WorkingSet.isActiveIdxRegPhaseOne[2] = 3;
  WorkingSet.isActiveIdxRegPhaseOne[3] = mLinIneq_tmp;
  WorkingSet.isActiveIdxRegPhaseOne[4] = mLinIneq_tmp + 15;
  WorkingSet.isActiveIdxRegPhaseOne[5] = 9;
  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (b_i = 0; b_i < 6; b_i++) {
    WorkingSet.isActiveIdxRegularized[b_i] =
      WorkingSet.isActiveIdxRegPhaseOne[b_i];
  }

  WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  WorkingSet.isActiveIdxRegPhaseOne[1] = 0;
  WorkingSet.isActiveIdxRegPhaseOne[2] = 3;
  WorkingSet.isActiveIdxRegPhaseOne[3] = mLinIneq_tmp;
  WorkingSet.isActiveIdxRegPhaseOne[4] = mLinIneq_tmp + 16;
  WorkingSet.isActiveIdxRegPhaseOne[5] = 9;
  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  if (mLinIneq_tmp > 0) {
    for (mConstrMax = 0; mConstrMax < mLinIneq_tmp; mConstrMax++) {
      for (i = 0; i < 9; i++) {
        WorkingSet.Aineq.data[i + WorkingSet.ldA * mConstrMax] =
          Aineq_data[mConstrMax + mLinIneq_tmp * i];
      }
    }
  }

  for (i = 0; i < 9; i++) {
    b_i = WorkingSet.indexLB.data[i];
    TrialState.xstarsqp[b_i - 1] = fmax(TrialState.xstarsqp[b_i - 1], -1.0);
  }

  for (i = 0; i < 9; i++) {
    b_i = WorkingSet.indexUB.data[i];
    TrialState.xstarsqp[b_i - 1] = fmin(TrialState.xstarsqp[b_i - 1], iv[b_i - 1]);
  }

  memcpy(&x[0], &TrialState.xstarsqp[0], 9U * sizeof(double));
  fval = optimize_forces_anonFcn1(fun_workspace_center, fun_workspace_A_alpha,
    fun_workspace_b_alpha, fun_workspace_deltaT, fun_workspace_T_min,
    TrialState.xstarsqp);
  i = 1;
  b = rtIsNaN(fval);
  if (rtIsInf(fval) || b) {
    if (b) {
      i = -3;
    } else if (fval < 0.0) {
      i = -1;
    } else {
      i = -2;
    }
  }

  TrialState.sqpFval = fval;
  if (i == 1) {
    computeConstraints_(x, TrialState.cIneq.data, TrialState.iNonIneq0,
                        TrialState.cEq);
  }

  obj.objfun.workspace.center[0] = fun_workspace_center[0];
  obj.objfun.workspace.center[1] = fun_workspace_center[1];
  obj.objfun.workspace.center[2] = fun_workspace_center[2];
  memcpy(&obj.objfun.workspace.A_alpha[0], &fun_workspace_A_alpha[0], 18U *
         sizeof(double));
  for (i = 0; i < 6; i++) {
    obj.objfun.workspace.b_alpha[i] = fun_workspace_b_alpha[i];
  }

  obj.objfun.workspace.deltaT[0] = fun_workspace_deltaT[0];
  obj.objfun.workspace.T_min[0] = fun_workspace_T_min[0];
  obj.objfun.workspace.deltaT[1] = fun_workspace_deltaT[1];
  obj.objfun.workspace.T_min[1] = fun_workspace_T_min[1];
  obj.objfun.workspace.deltaT[2] = fun_workspace_deltaT[2];
  obj.objfun.workspace.T_min[2] = fun_workspace_T_min[2];
  obj.f_1 = 0.0;
  obj.numEvals = 0;
  obj.hasBounds = true;
  computeFiniteDifferences(&obj, fval, TrialState.cEq, TrialState.xstarsqp,
    TrialState.grad.data, WorkingSet.Aeq.data, WorkingSet.ldA);
  TrialState.FunctionEvaluations = obj.numEvals + 1;
  ldAi = WorkingSet.ldA;
  if (mLinIneq_tmp > 0) {
    i = TrialState.cIneq.size[0];
    if (i - 1 >= 0) {
      memcpy(&y_data[0], &TrialState.cIneq.data[0], (unsigned int)i * sizeof
             (double));
    }

    memcpy(&y_data[0], &bineq_data[0], (unsigned int)mLinIneq_tmp * sizeof
           (double));
    if (i - 1 >= 0) {
      memcpy(&TrialState.cIneq.data[0], &y_data[0], (unsigned int)i * sizeof
             (double));
    }

    for (i = 0; i < mLinIneq_tmp; i++) {
      TrialState.cIneq.data[i] = -TrialState.cIneq.data[i];
    }

    i = 0;
    b_i = WorkingSet.ldA * (mLinIneq_tmp - 1) + 1;
    for (iac = 1; ldAi < 0 ? iac >= b_i : iac <= b_i; iac += ldAi) {
      c = 0.0;
      mConstrMax = iac + 8;
      for (ia = iac; ia <= mConstrMax; ia++) {
        c += WorkingSet.Aineq.data[ia - 1] * TrialState.xstarsqp[ia - iac];
      }

      TrialState.cIneq.data[i] += c;
      i++;
    }
  }

  b_updateWorkingSetForNewQP(x0, &WorkingSet, mLinIneq_tmp,
    TrialState.cIneq.data, TrialState.cEq);
  initActiveSet(&WorkingSet);
  MeritFunction.initFval = fval;
  MeritFunction.penaltyParam = 1.0;
  MeritFunction.threshold = 0.0001;
  MeritFunction.nPenaltyDecreases = 0;
  MeritFunction.linearizedConstrViol = 0.0;
  MeritFunction.initConstrViolationEq = (fabs(TrialState.cEq[0]) + fabs
    (TrialState.cEq[1])) + fabs(TrialState.cEq[2]);
  c = 0.0;
  for (i = 0; i < mLinIneq_tmp; i++) {
    scale = TrialState.cIneq.data[i];
    if (scale > 0.0) {
      c += scale;
    }
  }

  MeritFunction.initConstrViolationIneq = c;
  MeritFunction.phi = 0.0;
  MeritFunction.phiPrimePlus = 0.0;
  MeritFunction.phiFullStep = 0.0;
  MeritFunction.feasRelativeFactor = 0.0;
  MeritFunction.nlpPrimalFeasError = 0.0;
  MeritFunction.nlpDualFeasError = 0.0;
  MeritFunction.nlpComplError = 0.0;
  MeritFunction.firstOrderOpt = 0.0;
  driver(bineq_data, &TrialState, &MeritFunction, &FcnEvaluator, &obj, &memspace,
         &WorkingSet, &QRManager, &CholManager, &QPObjective, mLinIneq_tmp,
         unusedExpr);
  fval = TrialState.sqpFval;
  *exitflag = TrialState.sqpExitFlag;
  *output_iterations = TrialState.sqpIterations;
  *output_funcCount = TrialState.FunctionEvaluations;
  *output_constrviolation = MeritFunction.nlpPrimalFeasError;
  c = 0.0;
  scale = 3.3121686421112381E-170;
  for (i = 0; i < 9; i++) {
    double absxk;
    x[i] = TrialState.xstarsqp[i];
    absxk = fabs(TrialState.delta_x.data[i]);
    if (absxk > scale) {
      double t;
      t = scale / absxk;
      c = c * t * t + 1.0;
      scale = absxk;
    } else {
      double t;
      t = absxk / scale;
      c += t * t;
    }
  }

  *output_stepsize = scale * sqrt(c);
  *output_lssteplength = TrialState.steplength;
  *output_firstorderopt = MeritFunction.firstOrderOpt;
  return fval;
}

/*
 * Arguments    : e_struct_T *obj
 *                int NColsRemain
 * Return Type  : void
 */
static void fullColLDL2_(e_struct_T *obj, int NColsRemain)
{
  int LDimSizeP1;
  int ijA;
  int j;
  int jA;
  int k;
  LDimSizeP1 = obj->ldm;
  for (k = 0; k < NColsRemain; k++) {
    double alpha1;
    double y;
    int LD_diagOffset;
    int i;
    int offset1;
    int subMatrixDim;
    LD_diagOffset = (LDimSizeP1 + 1) * k;
    alpha1 = -1.0 / obj->FMat.data[LD_diagOffset];
    subMatrixDim = NColsRemain - k;
    offset1 = LD_diagOffset + 2;
    y = obj->workspace_;
    for (jA = 0; jA <= subMatrixDim - 2; jA++) {
      y = obj->FMat.data[(LD_diagOffset + jA) + 1];
    }

    obj->workspace_ = y;
    if (!(alpha1 == 0.0)) {
      jA = LD_diagOffset + LDimSizeP1;
      for (j = 0; j <= subMatrixDim - 2; j++) {
        if (y != 0.0) {
          double temp;
          int i1;
          temp = y * alpha1;
          i = jA + 2;
          i1 = subMatrixDim + jA;
          for (ijA = i; ijA <= i1; ijA++) {
            obj->FMat.data[ijA - 1] += y * temp;
          }
        }

        jA += obj->ldm;
      }
    }

    alpha1 = 1.0 / obj->FMat.data[LD_diagOffset];
    i = LD_diagOffset + subMatrixDim;
    for (jA = offset1; jA <= i; jA++) {
      obj->FMat.data[jA - 1] *= alpha1;
    }
  }
}

/*
 * Arguments    : h_struct_T *obj
 * Return Type  : void
 */
static void initActiveSet(h_struct_T *obj)
{
  int b_i;
  int i;
  int iATw0;
  int iAeq0;
  int idxFillStart;
  int idx_local;
  setProblemType(obj, 3);
  idxFillStart = obj->isActiveIdx[2];
  i = obj->mConstrMax;
  if (idxFillStart <= i) {
    memset(&obj->isActiveConstr.data[idxFillStart + -1], 0, (unsigned int)((i -
             idxFillStart) + 1) * sizeof(bool));
  }

  obj->nWConstr[0] = obj->sizes[0];
  obj->nWConstr[1] = 3;
  obj->nWConstr[2] = 0;
  obj->nWConstr[3] = 0;
  obj->nWConstr[4] = 0;
  obj->nActiveConstr = obj->nWConstr[0] + 3;
  i = (unsigned char)obj->sizes[0];
  for (idx_local = 0; idx_local < i; idx_local++) {
    obj->Wid.data[idx_local] = 1;
    obj->Wlocalidx.data[idx_local] = idx_local + 1;
    obj->isActiveConstr.data[idx_local] = true;
    idxFillStart = obj->ldA * idx_local;
    iAeq0 = (unsigned char)(obj->indexFixed.data[idx_local] - 1);
    if (iAeq0 - 1 >= 0) {
      memset(&obj->ATwset.data[idxFillStart], 0, (unsigned int)iAeq0 * sizeof
             (double));
    }

    obj->ATwset.data[(obj->indexFixed.data[idx_local] + idxFillStart) - 1] = 1.0;
    iAeq0 = obj->indexFixed.data[idx_local] + 1;
    iATw0 = obj->nVar;
    if (iAeq0 <= iATw0) {
      memset(&obj->ATwset.data[(iAeq0 + idxFillStart) + -1], 0, (unsigned int)
             ((((iATw0 + idxFillStart) - iAeq0) - idxFillStart) + 1) * sizeof
             (double));
    }

    obj->bwset.data[idx_local] = obj->ub.data[obj->indexFixed.data[idx_local] -
      1];
  }

  for (idx_local = 0; idx_local < 3; idx_local++) {
    idxFillStart = obj->sizes[0] + idx_local;
    obj->Wid.data[idxFillStart] = 2;
    obj->Wlocalidx.data[idxFillStart] = idx_local + 1;
    obj->isActiveConstr.data[idxFillStart] = true;
    iAeq0 = obj->ldA * idx_local;
    iATw0 = obj->ldA * idxFillStart;
    i = obj->nVar;
    for (b_i = 0; b_i < i; b_i++) {
      obj->ATwset.data[iATw0 + b_i] = obj->Aeq.data[iAeq0 + b_i];
    }

    obj->bwset.data[idxFillStart] = obj->beq[idx_local];
  }
}

/*
 * Arguments    : bool obj_hasLinear
 *                int obj_nvar
 *                double workspace_data[]
 *                const double H[81]
 *                const double f_data[]
 *                const double x_data[]
 * Return Type  : void
 */
static void linearForm_(bool obj_hasLinear, int obj_nvar, double workspace_data[],
  const double H[81], const double f_data[], const double x_data[])
{
  int ia;
  int iac;
  int ix;
  ix = 0;
  if (obj_hasLinear) {
    if (obj_nvar - 1 >= 0) {
      memcpy(&workspace_data[0], &f_data[0], (unsigned int)obj_nvar * sizeof
             (double));
    }

    ix = 1;
  }

  if (obj_nvar != 0) {
    int i;
    if ((ix != 1) && (obj_nvar - 1 >= 0)) {
      memset(&workspace_data[0], 0, (unsigned int)obj_nvar * sizeof(double));
    }

    ix = 0;
    i = obj_nvar * (obj_nvar - 1) + 1;
    for (iac = 1; obj_nvar < 0 ? iac >= i : iac <= i; iac += obj_nvar) {
      double c;
      int i1;
      c = 0.5 * x_data[ix];
      i1 = iac + obj_nvar;
      for (ia = iac; ia < i1; ia++) {
        int i2;
        i2 = ia - iac;
        workspace_data[i2] += H[ia - 1] * c;
      }

      ix++;
    }
  }
}

/*
 * Arguments    : h_struct_T *obj
 *                const double x_data[]
 *                int ix0
 * Return Type  : double
 */
static double maxConstraintViolation(h_struct_T *obj, const double x_data[], int
  ix0)
{
  double v;
  int i;
  int idx;
  int mIneq;
  if (obj->probType == 2) {
    v = 0.0;
    mIneq = obj->sizes[2];
    if (obj->Aineq.size[0] != 0) {
      if (mIneq - 1 >= 0) {
        memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0], (unsigned
                int)mIneq * sizeof(double));
      }

      b_xgemv(9, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data, ix0,
              obj->maxConstrWorkspace.data);
      i = (unsigned char)obj->sizes[2];
      for (idx = 0; idx < i; idx++) {
        double d;
        d = obj->maxConstrWorkspace.data[idx] - x_data[(ix0 + idx) + 8];
        obj->maxConstrWorkspace.data[idx] = d;
        v = fmax(v, d);
      }
    }

    obj->maxConstrWorkspace.data[0] = obj->beq[0];
    obj->maxConstrWorkspace.data[1] = obj->beq[1];
    obj->maxConstrWorkspace.data[2] = obj->beq[2];
    b_xgemv(9, 3, obj->Aeq.data, obj->ldA, x_data, ix0,
            obj->maxConstrWorkspace.data);
    i = ix0 + mIneq;
    mIneq = ix0 + obj->sizes[2];
    obj->maxConstrWorkspace.data[0] = (obj->maxConstrWorkspace.data[0] -
      x_data[i + 8]) + x_data[mIneq + 11];
    v = fmax(v, fabs(obj->maxConstrWorkspace.data[0]));
    obj->maxConstrWorkspace.data[1] = (obj->maxConstrWorkspace.data[1] -
      x_data[i + 9]) + x_data[mIneq + 12];
    v = fmax(v, fabs(obj->maxConstrWorkspace.data[1]));
    obj->maxConstrWorkspace.data[2] = (obj->maxConstrWorkspace.data[2] -
      x_data[i + 10]) + x_data[mIneq + 13];
    v = fmax(v, fabs(obj->maxConstrWorkspace.data[2]));
  } else {
    v = 0.0;
    mIneq = obj->sizes[2];
    if (obj->Aineq.size[0] != 0) {
      if (mIneq - 1 >= 0) {
        memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0], (unsigned
                int)mIneq * sizeof(double));
      }

      b_xgemv(obj->nVar, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data, ix0,
              obj->maxConstrWorkspace.data);
      i = (unsigned char)obj->sizes[2];
      for (idx = 0; idx < i; idx++) {
        v = fmax(v, obj->maxConstrWorkspace.data[idx]);
      }
    }

    obj->maxConstrWorkspace.data[0] = obj->beq[0];
    obj->maxConstrWorkspace.data[1] = obj->beq[1];
    obj->maxConstrWorkspace.data[2] = obj->beq[2];
    b_xgemv(obj->nVar, 3, obj->Aeq.data, obj->ldA, x_data, ix0,
            obj->maxConstrWorkspace.data);
    v = fmax(v, fabs(obj->maxConstrWorkspace.data[0]));
    v = fmax(v, fabs(obj->maxConstrWorkspace.data[1]));
    v = fmax(v, fabs(obj->maxConstrWorkspace.data[2]));
  }

  if (obj->sizes[3] > 0) {
    i = (unsigned char)obj->sizes[3];
    for (idx = 0; idx < i; idx++) {
      mIneq = obj->indexLB.data[idx] - 1;
      v = fmax(v, -x_data[(ix0 + mIneq) - 1] - obj->lb.data[mIneq]);
    }
  }

  if (obj->sizes[4] > 0) {
    i = (unsigned char)obj->sizes[4];
    for (idx = 0; idx < i; idx++) {
      mIneq = obj->indexUB.data[idx] - 1;
      v = fmax(v, x_data[(ix0 + mIneq) - 1] - obj->ub.data[mIneq]);
    }
  }

  if (obj->sizes[0] > 0) {
    i = (unsigned char)obj->sizes[0];
    for (idx = 0; idx < i; idx++) {
      v = fmax(v, fabs(x_data[(ix0 + obj->indexFixed.data[idx]) - 2] -
                       obj->ub.data[obj->indexFixed.data[idx] - 1]));
    }
  }

  return v;
}

/*
 * Arguments    : const double x[6]
 * Return Type  : double
 */
static double minimum(const double x[6])
{
  double ex;
  int idx;
  int k;
  if (!rtIsNaN(x[0])) {
    idx = 1;
  } else {
    bool exitg1;
    idx = 0;
    k = 2;
    exitg1 = false;
    while ((!exitg1) && (k < 7)) {
      if (!rtIsNaN(x[k - 1])) {
        idx = k;
        exitg1 = true;
      } else {
        k++;
      }
    }
  }

  if (idx == 0) {
    ex = x[0];
  } else {
    ex = x[idx - 1];
    idx++;
    for (k = idx; k < 7; k++) {
      double d;
      d = x[k - 1];
      if (ex > d) {
        ex = d;
      }
    }
  }

  return ex;
}

/*
 * Arguments    : h_struct_T *obj
 * Return Type  : void
 */
static void modifyOverheadPhaseOne_(h_struct_T *obj)
{
  int i;
  int idx;
  int idxEq;
  int idxEq_tmp;
  i = (unsigned char)obj->sizes[0];
  for (idx = 0; idx < i; idx++) {
    obj->ATwset.data[(obj->nVar + obj->ldA * idx) - 1] = 0.0;
  }

  idxEq_tmp = obj->nVar - 1;
  obj->Aeq.data[idxEq_tmp] = 0.0;
  i = obj->ldA * (obj->isActiveIdx[1] - 1);
  obj->ATwset.data[idxEq_tmp + i] = 0.0;
  idxEq = (obj->nVar + obj->ldA) - 1;
  obj->Aeq.data[idxEq] = 0.0;
  obj->ATwset.data[idxEq + i] = 0.0;
  idxEq = (obj->nVar + (obj->ldA << 1)) - 1;
  obj->Aeq.data[idxEq] = 0.0;
  obj->ATwset.data[idxEq + i] = 0.0;
  i = (unsigned char)obj->sizes[2];
  for (idx = 0; idx < i; idx++) {
    obj->Aineq.data[(obj->nVar + obj->ldA * idx) - 1] = -1.0;
  }

  obj->indexLB.data[obj->sizes[3] - 1] = obj->nVar;
  obj->lb.data[idxEq_tmp] = 1.0E-5;
  idxEq = obj->isActiveIdx[2];
  i = obj->nActiveConstr;
  for (idx = idxEq; idx <= i; idx++) {
    obj->ATwset.data[(obj->nVar + obj->ldA * (idx - 1)) - 1] = -1.0;
  }

  idxEq = obj->isActiveIdx[4] - 1;
  if (obj->nWConstr[4] > 0) {
    i = obj->sizesNormal[4];
    for (idx = i; idx >= 1; idx--) {
      idxEq_tmp = idxEq + idx;
      obj->isActiveConstr.data[idxEq_tmp] = obj->isActiveConstr.data[idxEq_tmp -
        1];
    }
  } else {
    obj->isActiveConstr.data[(obj->isActiveIdx[4] + obj->sizesNormal[4]) - 1] =
      false;
  }

  obj->isActiveConstr.data[obj->isActiveIdx[4] - 1] = false;
}

/*
 * @(W)
 *
 * Arguments    : const double center[3]
 *                const double A_alpha[18]
 *                const double b_alpha[6]
 *                const double deltaT[3]
 *                const double T_min[3]
 *                const double W[9]
 * Return Type  : double
 */
static double optimize_forces_anonFcn1(const double center[3], const double
  A_alpha[18], const double b_alpha[6], const double deltaT[3], const double
  T_min[3], const double W[9])
{
  double A_wrench[18];
  double b_W[9];
  double c_W[9];
  double dv1[9];
  double b_A_wrench[6];
  double c_A_wrench[6];
  double dv2[6];
  double dv[3];
  double d;
  double d1;
  double varargout_1;
  int W_tmp;
  int i;

  /* 'optimize_forces:64' @(W) objective(W, center, A_alpha, b_alpha, deltaT, T_min) */
  /*  Solver */
  /*  优化目标函数 */
  /* 'optimize_forces:18' W = reshape(W, [3, length(W) / 3]); */
  /* 'optimize_forces:19' W = W ./ vecnorm(W, 2, 1); */
  vecnorm(W, dv);
  for (i = 0; i < 3; i++) {
    varargout_1 = dv[i];
    b_W[3 * i] = W[3 * i] / varargout_1;
    W_tmp = 3 * i + 1;
    b_W[W_tmp] = W[W_tmp] / varargout_1;
    W_tmp = 3 * i + 2;
    b_W[W_tmp] = W[W_tmp] / varargout_1;
  }

  /* 'optimize_forces:21' A_wrench = A_alpha * pinv(W * diag(deltaT)); */
  diag(deltaT, dv1);
  for (i = 0; i < 3; i++) {
    varargout_1 = b_W[i];
    d = b_W[i + 3];
    d1 = b_W[i + 6];
    for (W_tmp = 0; W_tmp < 3; W_tmp++) {
      c_W[i + 3 * W_tmp] = (varargout_1 * dv1[3 * W_tmp] + d * dv1[3 * W_tmp + 1])
        + d1 * dv1[3 * W_tmp + 2];
    }
  }

  pinv(c_W, dv1);

  /* 'optimize_forces:22' b_wrench = A_wrench * W * T_min + b_alpha; */
  /* 'optimize_forces:24' distances = (b_wrench - A_wrench * p) ./ vecnorm(A_wrench, 2, 2); */
  /*  返回最小距离（负号用于最大化） */
  /* 'optimize_forces:27' obj_val = -min(distances); */
  for (i = 0; i < 6; i++) {
    varargout_1 = A_alpha[i];
    d = A_alpha[i + 6];
    d1 = A_alpha[i + 12];
    for (W_tmp = 0; W_tmp < 3; W_tmp++) {
      A_wrench[i + 6 * W_tmp] = (varargout_1 * dv1[3 * W_tmp] + d * dv1[3 *
        W_tmp + 1]) + d1 * dv1[3 * W_tmp + 2];
    }

    varargout_1 = 0.0;
    d = 0.0;
    for (W_tmp = 0; W_tmp < 3; W_tmp++) {
      varargout_1 += ((A_wrench[i] * b_W[3 * W_tmp] + A_wrench[i + 6] * b_W[3 *
                       W_tmp + 1]) + A_wrench[i + 12] * b_W[3 * W_tmp + 2]) *
        T_min[W_tmp];
      d += A_wrench[i + 6 * W_tmp] * center[W_tmp];
    }

    c_A_wrench[i] = d;
    b_A_wrench[i] = varargout_1 + b_alpha[i];
  }

  b_vecnorm(A_wrench, dv2);
  for (i = 0; i < 6; i++) {
    b_A_wrench[i] = (b_A_wrench[i] - c_A_wrench[i]) / dv2[i];
  }

  return -minimum(b_A_wrench);
}

/*
 * Arguments    : const double H[81]
 *                const double f_data[]
 *                c_struct_T *solution
 *                g_struct_T *memspace
 *                h_struct_T *workingset
 *                d_struct_T *qrmanager
 *                e_struct_T *cholmanager
 *                f_struct_T *objective
 *                const char options_SolverName[7]
 *                const k_struct_T *runTimeOptions
 * Return Type  : void
 */
static void phaseone(const double H[81], const double f_data[], c_struct_T
                     *solution, g_struct_T *memspace, h_struct_T *workingset,
                     d_struct_T *qrmanager, e_struct_T *cholmanager, f_struct_T *
                     objective, const char options_SolverName[7], const
                     k_struct_T *runTimeOptions)
{
  static const char b[7] = { 'f', 'm', 'i', 'n', 'c', 'o', 'n' };

  int PROBTYPE_ORIG;
  int activeSetChangeID;
  int globalActiveConstrIdx;
  int i;
  int iAw0;
  int idx;
  int idxMinLambda;
  int idx_global;
  int nVar_tmp;
  int ret;
  bool subProblemChanged;
  bool updateFval;
  PROBTYPE_ORIG = workingset->probType;
  nVar_tmp = workingset->nVar;
  solution->xstar.data[workingset->nVar] = solution->maxConstr + 1.0;
  if (workingset->probType == 3) {
    ret = 1;
  } else {
    ret = 4;
  }

  setProblemType(workingset, ret);
  ret = workingset->nWConstr[0] + workingset->nWConstr[1];
  iAw0 = ret + 1;
  idxMinLambda = workingset->nActiveConstr;
  for (idx_global = iAw0; idx_global <= idxMinLambda; idx_global++) {
    workingset->isActiveConstr.data[(workingset->isActiveIdx
      [workingset->Wid.data[idx_global - 1] - 1] + workingset->
      Wlocalidx.data[idx_global - 1]) - 2] = false;
  }

  workingset->nWConstr[2] = 0;
  workingset->nWConstr[3] = 0;
  workingset->nWConstr[4] = 0;
  workingset->nActiveConstr = ret;
  objective->prev_objtype = objective->objtype;
  objective->prev_nvar = objective->nvar;
  objective->prev_hasLinear = objective->hasLinear;
  objective->objtype = 5;
  objective->nvar = nVar_tmp + 1;
  objective->gammaScalar = 1.0;
  objective->hasLinear = true;
  subProblemChanged = true;
  updateFval = true;
  activeSetChangeID = 0;
  i = workingset->nVar;
  globalActiveConstrIdx = 0;
  computeGrad_StoreHx(objective, H, f_data, solution->xstar.data);
  solution->fstar = computeFval_ReuseHx(objective,
    memspace->workspace_float.data, f_data, solution->xstar.data);
  solution->state = -5;
  ret = workingset->mConstrMax;
  if (ret - 1 >= 0) {
    memset(&solution->lambda.data[0], 0, (unsigned int)ret * sizeof(double));
  }

  int exitg1;
  do {
    exitg1 = 0;
    if (solution->state == -5) {
      double minLambda;
      int i1;
      bool guard1;
      bool guard2;
      guard1 = false;
      guard2 = false;
      if (subProblemChanged) {
        switch (activeSetChangeID) {
         case 1:
          squareQ_appendCol(qrmanager, workingset->ATwset.data, workingset->ldA *
                            (workingset->nActiveConstr - 1) + 1);
          break;

         case -1:
          deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
          break;

         default:
          factorQR(qrmanager, workingset->ATwset.data, i,
                   workingset->nActiveConstr, workingset->ldA);
          computeQ_(qrmanager, qrmanager->mrows);
          break;
        }

        ret = memcmp(&options_SolverName[0], &b[0], 7);
        compute_deltax(H, solution, memspace, qrmanager, cholmanager, objective,
                       ret == 0);
        if (solution->state != -5) {
          exitg1 = 1;
        } else if ((c_xnrm2(i, solution->searchDir.data) <
                    1.4901161193847657E-10) || (workingset->nActiveConstr >= i))
        {
          guard2 = true;
        } else {
          minLambda = feasibleratiotest(solution->xstar.data,
            solution->searchDir.data, memspace->workspace_float.data,
            memspace->workspace_float.size, workingset->nVar, workingset->ldA,
            workingset->Aineq.data, workingset->bineq.data, workingset->lb.data,
            workingset->ub.data, workingset->indexLB.data,
            workingset->indexUB.data, workingset->sizes, workingset->isActiveIdx,
            workingset->isActiveConstr.data, workingset->nWConstr, true,
            &updateFval, &idx_global, &idxMinLambda);
          if (updateFval) {
            switch (idx_global) {
             case 3:
              workingset->nWConstr[2]++;
              workingset->isActiveConstr.data[(workingset->isActiveIdx[2] +
                idxMinLambda) - 2] = true;
              workingset->nActiveConstr++;
              idx_global = workingset->nActiveConstr - 1;
              workingset->Wid.data[idx_global] = 3;
              workingset->Wlocalidx.data[idx_global] = idxMinLambda;
              ret = workingset->ldA * (idxMinLambda - 1);
              iAw0 = workingset->ldA * idx_global;
              i1 = workingset->nVar;
              for (idx = 0; idx < i1; idx++) {
                workingset->ATwset.data[iAw0 + idx] = workingset->Aineq.data[ret
                  + idx];
              }

              workingset->bwset.data[idx_global] = workingset->
                bineq.data[idxMinLambda - 1];
              break;

             case 4:
              addBoundToActiveSetMatrix_(workingset, 4, idxMinLambda);
              break;

             default:
              addBoundToActiveSetMatrix_(workingset, 5, idxMinLambda);
              break;
            }

            activeSetChangeID = 1;
          } else {
            if (objective->objtype == 5) {
              if (c_xnrm2(objective->nvar, solution->searchDir.data) > 100.0 *
                  (double)objective->nvar * 1.4901161193847656E-8) {
                solution->state = 3;
              } else {
                solution->state = 4;
              }
            }

            subProblemChanged = false;
            if (workingset->nActiveConstr == 0) {
              solution->state = 1;
            }
          }

          if ((i >= 1) && (!(minLambda == 0.0))) {
            for (ret = 0; ret < i; ret++) {
              solution->xstar.data[ret] += minLambda * solution->
                searchDir.data[ret];
            }
          }

          computeGrad_StoreHx(objective, H, f_data, solution->xstar.data);
          updateFval = true;
          guard1 = true;
        }
      } else {
        if (i - 1 >= 0) {
          memset(&solution->searchDir.data[0], 0, (unsigned int)i * sizeof
                 (double));
        }

        guard2 = true;
      }

      if (guard2) {
        compute_lambda(memspace->workspace_float.data, solution, objective,
                       qrmanager);
        if ((solution->state != -7) || (workingset->nActiveConstr > i)) {
          idxMinLambda = 0;
          minLambda = 0.0;
          idx_global = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
          i1 = workingset->nActiveConstr;
          for (idx = idx_global; idx <= i1; idx++) {
            double d;
            d = solution->lambda.data[idx - 1];
            if (d < minLambda) {
              minLambda = d;
              idxMinLambda = idx;
            }
          }

          if (idxMinLambda == 0) {
            solution->state = 1;
          } else {
            activeSetChangeID = -1;
            globalActiveConstrIdx = idxMinLambda;
            subProblemChanged = true;
            removeConstr(workingset, idxMinLambda);
            if (idxMinLambda < workingset->nActiveConstr + 1) {
              solution->lambda.data[idxMinLambda - 1] = solution->
                lambda.data[workingset->nActiveConstr];
            }

            solution->lambda.data[workingset->nActiveConstr] = 0.0;
          }
        } else {
          idx_global = workingset->nActiveConstr;
          activeSetChangeID = 0;
          globalActiveConstrIdx = workingset->nActiveConstr;
          subProblemChanged = true;
          removeConstr(workingset, workingset->nActiveConstr);
          solution->lambda.data[idx_global - 1] = 0.0;
        }

        updateFval = false;
        guard1 = true;
      }

      if (guard1) {
        solution->iterations++;
        iAw0 = objective->nvar;
        if ((solution->iterations >= runTimeOptions->MaxIterations) &&
            ((solution->state != 1) || (objective->objtype == 5))) {
          solution->state = 0;
        }

        if (solution->iterations - solution->iterations / 50 * 50 == 0) {
          solution->maxConstr = b_maxConstraintViolation(workingset,
            solution->xstar.data);
          minLambda = solution->maxConstr;
          if (objective->objtype == 5) {
            minLambda = solution->maxConstr - solution->xstar.data
              [objective->nvar - 1];
          }

          if (minLambda > 1.0E-6) {
            bool nonDegenerateWset;
            if (iAw0 - 1 >= 0) {
              memcpy(&solution->searchDir.data[0], &solution->xstar.data[0],
                     (unsigned int)iAw0 * sizeof(double));
            }

            nonDegenerateWset = feasibleX0ForWorkingSet
              (memspace->workspace_float.data, memspace->workspace_float.size,
               solution->searchDir.data, workingset, qrmanager);
            if ((!nonDegenerateWset) && (solution->state != 0)) {
              solution->state = -2;
            }

            activeSetChangeID = 0;
            minLambda = b_maxConstraintViolation(workingset,
              solution->searchDir.data);
            if ((minLambda < solution->maxConstr) && (iAw0 - 1 >= 0)) {
              memcpy(&solution->xstar.data[0], &solution->searchDir.data[0],
                     (unsigned int)iAw0 * sizeof(double));
            }
          }
        }

        if (updateFval) {
          solution->fstar = computeFval_ReuseHx(objective,
            memspace->workspace_float.data, f_data, solution->xstar.data);
          if ((solution->fstar < 1.0E-6) && ((solution->state != 0) ||
               (objective->objtype != 5))) {
            solution->state = 2;
          }
        }
      }
    } else {
      if (!updateFval) {
        solution->fstar = computeFval_ReuseHx(objective,
          memspace->workspace_float.data, f_data, solution->xstar.data);
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (workingset->isActiveConstr.data[(workingset->isActiveIdx[3] +
       workingset->sizes[3]) - 2]) {
    bool exitg2;
    idx = workingset->sizes[0] + 4;
    exitg2 = false;
    while ((!exitg2) && (idx <= workingset->nActiveConstr)) {
      if ((workingset->Wid.data[idx - 1] == 4) && (workingset->
           Wlocalidx.data[idx - 1] == workingset->sizes[3])) {
        removeConstr(workingset, idx);
        exitg2 = true;
      } else {
        idx++;
      }
    }
  }

  ret = workingset->nActiveConstr;
  iAw0 = workingset->sizes[0] + 3;
  while ((ret > iAw0) && (ret > nVar_tmp)) {
    removeConstr(workingset, ret);
    ret--;
  }

  solution->maxConstr = solution->xstar.data[nVar_tmp];
  setProblemType(workingset, PROBTYPE_ORIG);
  objective->objtype = objective->prev_objtype;
  objective->nvar = objective->prev_nvar;
  objective->hasLinear = objective->prev_hasLinear;
}

/*
 * Arguments    : const double A[9]
 *                double X[9]
 * Return Type  : void
 */
static void pinv(const double A[9], double X[9])
{
  double absx;
  int br;
  int i;
  int ib;
  int ic;
  int j;
  int vcol;
  bool p;
  p = true;
  for (br = 0; br < 9; br++) {
    X[br] = 0.0;
    if (p) {
      absx = A[br];
      if (rtIsInf(absx) || rtIsNaN(absx)) {
        p = false;
      }
    } else {
      p = false;
    }
  }

  if (!p) {
    for (i = 0; i < 9; i++) {
      X[i] = rtNaN;
    }
  } else {
    double U[9];
    double V[9];
    double s[3];
    int r;
    bool exitg1;
    svd(A, U, s, V);
    absx = fabs(s[0]);
    if (rtIsInf(absx) || rtIsNaN(absx)) {
      absx = rtNaN;
    } else if (absx < 4.4501477170144028E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &vcol);
      absx = ldexp(1.0, vcol - 53);
    }

    absx *= 3.0;
    vcol = 0;
    exitg1 = false;
    while ((!exitg1) && (vcol < 3)) {
      if (rtIsInf(s[vcol]) || rtIsNaN(s[vcol])) {
        absx = 1.7976931348623157E+308;
        exitg1 = true;
      } else {
        vcol++;
      }
    }

    r = -1;
    br = 0;
    while ((br < 3) && (s[br] > absx)) {
      r++;
      br++;
    }

    if (r + 1 > 0) {
      vcol = 1;
      for (j = 0; j <= r; j++) {
        absx = 1.0 / s[j];
        i = vcol + 2;
        for (br = vcol; br <= i; br++) {
          V[br - 1] *= absx;
        }

        vcol += 3;
      }

      for (vcol = 0; vcol <= 6; vcol += 3) {
        i = vcol + 1;
        j = vcol + 3;
        if (i <= j) {
          memset(&X[i + -1], 0, (unsigned int)((j - i) + 1) * sizeof(double));
        }
      }

      br = 0;
      for (vcol = 0; vcol <= 6; vcol += 3) {
        int ar;
        ar = -1;
        br++;
        i = br + 3 * r;
        for (ib = br; ib <= i; ib += 3) {
          int i1;
          j = vcol + 1;
          i1 = vcol + 3;
          for (ic = j; ic <= i1; ic++) {
            X[ic - 1] += U[ib - 1] * V[(ar + ic) - vcol];
          }

          ar += 3;
        }
      }
    }
  }
}

/*
 * Arguments    : double A_data[]
 *                const int A_size[2]
 *                int m
 *                int n
 *                int nfxd
 *                double tau_data[]
 * Return Type  : void
 */
static void qrf(double A_data[], const int A_size[2], int m, int n, int nfxd,
                double tau_data[])
{
  double work_data[40];
  double atmp;
  int i;
  int lda;
  int loop_ub;
  lda = A_size[0];
  loop_ub = A_size[1];
  if (loop_ub - 1 >= 0) {
    memset(&work_data[0], 0, (unsigned int)loop_ub * sizeof(double));
  }

  loop_ub = (unsigned char)nfxd;
  for (i = 0; i < loop_ub; i++) {
    double d;
    int ii;
    int mmi;
    ii = i * lda + i;
    mmi = m - i;
    if (i + 1 < m) {
      atmp = A_data[ii];
      d = xzlarfg(mmi, &atmp, A_data, ii + 2);
      tau_data[i] = d;
      A_data[ii] = atmp;
    } else {
      d = 0.0;
      tau_data[i] = 0.0;
    }

    if (i + 1 < n) {
      atmp = A_data[ii];
      A_data[ii] = 1.0;
      xzlarf(mmi, (n - i) - 1, ii + 1, d, A_data, (ii + lda) + 1, lda, work_data);
      A_data[ii] = atmp;
    }
  }
}

/*
 * Arguments    : const double Hessian[81]
 *                const double grad_data[]
 *                c_struct_T *TrialState
 *                j_struct_T *MeritFunction
 *                g_struct_T *memspace
 *                h_struct_T *WorkingSet
 *                d_struct_T *QRManager
 *                e_struct_T *CholManager
 *                f_struct_T *QPObjective
 *                k_struct_T *qpoptions
 * Return Type  : void
 */
static void relaxed(const double Hessian[81], const double grad_data[],
                    c_struct_T *TrialState, j_struct_T *MeritFunction,
                    g_struct_T *memspace, h_struct_T *WorkingSet, d_struct_T
                    *QRManager, e_struct_T *CholManager, f_struct_T *QPObjective,
                    k_struct_T *qpoptions)
{
  double beta;
  double d;
  double s;
  double smax;
  int b_mIneq;
  int i;
  int idx;
  int idx_max;
  int mIneq;
  int nActiveLBArtificial;
  int nVarOrig;
  bool b_tf;
  bool tf;
  nVarOrig = WorkingSet->nVar;
  mIneq = WorkingSet->sizes[2];
  beta = 0.0;
  i = (unsigned char)WorkingSet->nVar;
  for (idx = 0; idx < i; idx++) {
    beta += Hessian[idx + 9 * idx];
  }

  beta /= (double)WorkingSet->nVar;
  if (TrialState->sqpIterations <= 1) {
    b_mIneq = QPObjective->nvar;
    if (QPObjective->nvar < 1) {
      idx_max = 0;
    } else {
      idx_max = 1;
      if (QPObjective->nvar > 1) {
        smax = fabs(grad_data[0]);
        for (idx = 2; idx <= b_mIneq; idx++) {
          s = fabs(grad_data[idx - 1]);
          if (s > smax) {
            idx_max = idx;
            smax = s;
          }
        }
      }
    }

    smax = 100.0 * fmax(1.0, fabs(grad_data[idx_max - 1]));
  } else {
    b_mIneq = WorkingSet->mConstr;
    if (WorkingSet->mConstr < 1) {
      idx_max = 0;
    } else {
      idx_max = 1;
      if (WorkingSet->mConstr > 1) {
        smax = fabs(TrialState->lambdasqp.data[0]);
        for (idx = 2; idx <= b_mIneq; idx++) {
          s = fabs(TrialState->lambdasqp.data[idx - 1]);
          if (s > smax) {
            idx_max = idx;
            smax = s;
          }
        }
      }
    }

    smax = fabs(TrialState->lambdasqp.data[idx_max - 1]);
  }

  QPObjective->nvar = WorkingSet->nVar;
  QPObjective->beta = beta;
  QPObjective->rho = smax;
  QPObjective->hasLinear = true;
  QPObjective->objtype = 4;
  setProblemType(WorkingSet, 2);
  b_mIneq = WorkingSet->sizes[2] + 1;
  idx_max = WorkingSet->sizes[3] - WorkingSet->sizes[2];
  i = (unsigned char)WorkingSet->sizes[2];
  if (i - 1 >= 0) {
    memcpy(&memspace->workspace_float.data[0], &WorkingSet->bineq.data[0],
           (unsigned int)i * sizeof(double));
  }

  xgemv(nVarOrig, WorkingSet->sizes[2], WorkingSet->Aineq.data, WorkingSet->ldA,
        TrialState->xstar.data, memspace->workspace_float.data);
  for (idx = 0; idx < i; idx++) {
    d = memspace->workspace_float.data[idx];
    TrialState->xstar.data[nVarOrig + idx] = (double)(d > 0.0) * d;
  }

  memspace->workspace_float.data[0] = WorkingSet->beq[0];
  memspace->workspace_float.data[1] = WorkingSet->beq[1];
  memspace->workspace_float.data[2] = WorkingSet->beq[2];
  xgemv(nVarOrig, 3, WorkingSet->Aeq.data, WorkingSet->ldA,
        TrialState->xstar.data, memspace->workspace_float.data);
  if (memspace->workspace_float.data[0] <= 0.0) {
    TrialState->xstar.data[(nVarOrig + b_mIneq) - 1] = 0.0;
    TrialState->xstar.data[(nVarOrig + b_mIneq) + 2] =
      -memspace->workspace_float.data[0];
    addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 6);
    if (memspace->workspace_float.data[0] >= -1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 3);
    }
  } else {
    i = nVarOrig + b_mIneq;
    TrialState->xstar.data[i - 1] = memspace->workspace_float.data[0];
    TrialState->xstar.data[i + 2] = 0.0;
    addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 3);
    if (memspace->workspace_float.data[0] <= 1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 6);
    }
  }

  if (memspace->workspace_float.data[1] <= 0.0) {
    TrialState->xstar.data[nVarOrig + b_mIneq] = 0.0;
    TrialState->xstar.data[(nVarOrig + b_mIneq) + 3] =
      -memspace->workspace_float.data[1];
    addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 5);
    if (memspace->workspace_float.data[1] >= -1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 2);
    }
  } else {
    i = nVarOrig + b_mIneq;
    TrialState->xstar.data[i] = memspace->workspace_float.data[1];
    TrialState->xstar.data[i + 3] = 0.0;
    addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 2);
    if (memspace->workspace_float.data[1] <= 1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 5);
    }
  }

  if (memspace->workspace_float.data[2] <= 0.0) {
    TrialState->xstar.data[(nVarOrig + b_mIneq) + 1] = 0.0;
    TrialState->xstar.data[(nVarOrig + b_mIneq) + 4] =
      -memspace->workspace_float.data[2];
    addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 4);
    if (memspace->workspace_float.data[2] >= -1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 1);
    }
  } else {
    i = nVarOrig + b_mIneq;
    TrialState->xstar.data[i + 1] = memspace->workspace_float.data[2];
    TrialState->xstar.data[i + 4] = 0.0;
    addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 1);
    if (memspace->workspace_float.data[2] <= 1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, (idx_max + b_mIneq) - 4);
    }
  }

  b_mIneq = qpoptions->MaxIterations;
  qpoptions->MaxIterations = (qpoptions->MaxIterations + WorkingSet->nVar) -
    nVarOrig;
  b_driver(Hessian, grad_data, TrialState, memspace, WorkingSet, QRManager,
           CholManager, QPObjective, qpoptions, qpoptions);
  qpoptions->MaxIterations = b_mIneq;
  idx_max = (WorkingSet->isActiveIdx[3] + WorkingSet->sizes[3]) - 7;
  memspace->workspace_int.data[0] = WorkingSet->isActiveConstr.data[idx_max];
  memspace->workspace_int.data[3] = WorkingSet->isActiveConstr.data[idx_max + 3];
  tf = WorkingSet->isActiveConstr.data[idx_max + 1];
  b_tf = WorkingSet->isActiveConstr.data[idx_max + 4];
  memspace->workspace_int.data[1] = tf;
  memspace->workspace_int.data[4] = b_tf;
  nActiveLBArtificial = ((WorkingSet->isActiveConstr.data[idx_max] +
    WorkingSet->isActiveConstr.data[idx_max + 3]) + tf) + b_tf;
  tf = WorkingSet->isActiveConstr.data[idx_max + 2];
  b_tf = WorkingSet->isActiveConstr.data[idx_max + 5];
  memspace->workspace_int.data[2] = tf;
  memspace->workspace_int.data[5] = b_tf;
  nActiveLBArtificial = (nActiveLBArtificial + tf) + b_tf;
  i = (unsigned char)WorkingSet->sizes[2];
  for (idx = 0; idx < i; idx++) {
    tf = WorkingSet->isActiveConstr.data[(idx_max - WorkingSet->sizes[2]) + idx];
    memspace->workspace_int.data[idx + 6] = tf;
    nActiveLBArtificial += tf;
  }

  if (TrialState->state != -6) {
    double penaltyParamTrial;
    double qpfvalQuadExcess;
    int ix0_tmp;
    b_mIneq = (WorkingSet->nVarMax - nVarOrig) - 1;
    ix0_tmp = nVarOrig + 1;
    s = 0.0;
    if (b_mIneq >= 1) {
      idx_max = nVarOrig + b_mIneq;
      for (idx = ix0_tmp; idx <= idx_max; idx++) {
        s += fabs(TrialState->xstar.data[idx - 1]);
      }
    }

    qpfvalQuadExcess = 0.0;
    if (b_mIneq >= 1) {
      i = (unsigned char)b_mIneq;
      for (idx = 0; idx < i; idx++) {
        d = TrialState->xstar.data[nVarOrig + idx];
        qpfvalQuadExcess += d * d;
      }
    }

    beta = (TrialState->fstar - smax * s) - beta / 2.0 * qpfvalQuadExcess;
    penaltyParamTrial = MeritFunction->penaltyParam;
    smax = 0.0;
    i = (unsigned char)mIneq;
    for (idx = 0; idx < i; idx++) {
      d = TrialState->cIneq.data[idx];
      if (d > 0.0) {
        smax += d;
      }
    }

    qpfvalQuadExcess = ((fabs(TrialState->cEq[0]) + fabs(TrialState->cEq[1])) +
                        fabs(TrialState->cEq[2])) + smax;
    smax = MeritFunction->linearizedConstrViol;
    s = 0.0;
    if (b_mIneq >= 1) {
      idx_max = nVarOrig + b_mIneq;
      for (idx = ix0_tmp; idx <= idx_max; idx++) {
        s += fabs(TrialState->xstar.data[idx - 1]);
      }
    }

    MeritFunction->linearizedConstrViol = s;
    smax = (qpfvalQuadExcess + smax) - s;
    if ((smax > 2.2204460492503131E-16) && (beta > 0.0)) {
      if (TrialState->sqpFval == 0.0) {
        d = 1.0;
      } else {
        d = 1.5;
      }

      penaltyParamTrial = d * beta / smax;
    }

    if (penaltyParamTrial < MeritFunction->penaltyParam) {
      MeritFunction->phi = TrialState->sqpFval + penaltyParamTrial *
        qpfvalQuadExcess;
      if ((MeritFunction->initFval + penaltyParamTrial *
           (MeritFunction->initConstrViolationEq +
            MeritFunction->initConstrViolationIneq)) - MeritFunction->phi >
          (double)MeritFunction->nPenaltyDecreases * MeritFunction->threshold) {
        MeritFunction->nPenaltyDecreases++;
        if ((MeritFunction->nPenaltyDecreases << 1) > TrialState->sqpIterations)
        {
          MeritFunction->threshold *= 10.0;
        }

        MeritFunction->penaltyParam = fmax(penaltyParamTrial, 1.0E-10);
      } else {
        MeritFunction->phi = TrialState->sqpFval + MeritFunction->penaltyParam *
          qpfvalQuadExcess;
      }
    } else {
      MeritFunction->penaltyParam = fmax(penaltyParamTrial, 1.0E-10);
      MeritFunction->phi = TrialState->sqpFval + MeritFunction->penaltyParam *
        qpfvalQuadExcess;
    }

    MeritFunction->phiPrimePlus = fmin(beta - MeritFunction->penaltyParam *
      qpfvalQuadExcess, 0.0);
    b_mIneq = WorkingSet->isActiveIdx[1] - 1;
    for (idx = 0; idx < 3; idx++) {
      if ((memspace->workspace_int.data[idx] != 0) &&
          (memspace->workspace_int.data[idx + 3] != 0)) {
        tf = true;
      } else {
        tf = false;
      }

      i = b_mIneq + idx;
      TrialState->lambda.data[i] *= (double)tf;
    }

    b_mIneq = WorkingSet->isActiveIdx[2];
    idx_max = WorkingSet->nActiveConstr;
    for (idx = b_mIneq; idx <= idx_max; idx++) {
      if (WorkingSet->Wid.data[idx - 1] == 3) {
        TrialState->lambda.data[idx - 1] *= (double)memspace->
          workspace_int.data[WorkingSet->Wlocalidx.data[idx - 1] + 5];
      }
    }
  }

  idx_max = WorkingSet->sizes[0];
  i = WorkingSet->sizes[3] - WorkingSet->sizes[2];
  idx = WorkingSet->nActiveConstr;
  while ((idx > idx_max + 3) && (nActiveLBArtificial > 0)) {
    if ((WorkingSet->Wid.data[idx - 1] == 4) && (WorkingSet->Wlocalidx.data[idx
         - 1] > i - 6)) {
      b_mIneq = WorkingSet->nActiveConstr - 1;
      smax = TrialState->lambda.data[b_mIneq];
      TrialState->lambda.data[b_mIneq] = 0.0;
      TrialState->lambda.data[idx - 1] = smax;
      removeConstr(WorkingSet, idx);
      nActiveLBArtificial--;
    }

    idx--;
  }

  QPObjective->nvar = nVarOrig;
  QPObjective->hasLinear = true;
  QPObjective->objtype = 3;
  setProblemType(WorkingSet, 3);
  sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
               WorkingSet->sizes, WorkingSet->isActiveIdx, WorkingSet->Wid.data,
               WorkingSet->Wlocalidx.data, memspace->workspace_float.data);
}

/*
 * Arguments    : h_struct_T *obj
 *                int idx_global
 * Return Type  : void
 */
static void removeConstr(h_struct_T *obj, int idx_global)
{
  int TYPE_tmp;
  int idx;
  TYPE_tmp = obj->Wid.data[idx_global - 1] - 1;
  obj->isActiveConstr.data[(obj->isActiveIdx[TYPE_tmp] + obj->
    Wlocalidx.data[idx_global - 1]) - 2] = false;
  if (idx_global < obj->nActiveConstr) {
    int i;
    int i1;
    i = obj->nActiveConstr - 1;
    obj->Wid.data[idx_global - 1] = obj->Wid.data[i];
    obj->Wlocalidx.data[idx_global - 1] = obj->Wlocalidx.data[i];
    i1 = (unsigned char)obj->nVar;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset.data[idx + obj->ldA * (idx_global - 1)] = obj->ATwset.data[idx
        + obj->ldA * i];
    }

    obj->bwset.data[idx_global - 1] = obj->bwset.data[i];
  }

  obj->nActiveConstr--;
  obj->nWConstr[TYPE_tmp]--;
}

/*
 * Arguments    : h_struct_T *obj
 *                int PROBLEM_TYPE
 * Return Type  : void
 */
static void setProblemType(h_struct_T *obj, int PROBLEM_TYPE)
{
  int i;
  int idx;
  int idx_col;
  switch (PROBLEM_TYPE) {
   case 3:
    {
      obj->nVar = 9;
      obj->mConstr = obj->mConstrOrig;
      if (obj->nWConstr[4] > 0) {
        int idxUpperExisting;
        idxUpperExisting = obj->isActiveIdx[4] - 2;
        i = (unsigned char)obj->sizesNormal[4];
        for (idx = 0; idx < i; idx++) {
          int i1;
          i1 = (idxUpperExisting + idx) + 1;
          obj->isActiveConstr.data[(obj->isActiveIdxNormal[4] + idx) - 1] =
            obj->isActiveConstr.data[i1];
          obj->isActiveConstr.data[i1] = false;
        }
      }

      for (i = 0; i < 5; i++) {
        obj->sizes[i] = obj->sizesNormal[i];
      }

      for (i = 0; i < 6; i++) {
        obj->isActiveIdx[i] = obj->isActiveIdxNormal[i];
      }
    }
    break;

   case 1:
    obj->nVar = 10;
    obj->mConstr = obj->mConstrOrig + 1;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesPhaseOne[i];
    }

    modifyOverheadPhaseOne_(obj);
    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxPhaseOne[i];
    }
    break;

   case 2:
    {
      obj->nVar = obj->nVarMax - 1;
      obj->mConstr = obj->mConstrMax - 1;
      for (i = 0; i < 5; i++) {
        obj->sizes[i] = obj->sizesRegularized[i];
      }

      if (obj->probType != 4) {
        int colOffsetATw;
        int i1;
        int i2;
        int idxUpperExisting;
        int mIneq;
        int offsetEq1;
        int offsetEq2;
        mIneq = obj->sizes[2] + 8;
        offsetEq1 = obj->sizes[2] + 9;
        offsetEq2 = obj->sizes[2] + 12;
        i = (unsigned char)obj->sizes[0];
        for (idx_col = 0; idx_col < i; idx_col++) {
          colOffsetATw = obj->ldA * idx_col;
          i1 = obj->nVar;
          memset(&obj->ATwset.data[colOffsetATw + 9], 0, (unsigned int)(i1 - 9) *
                 sizeof(double));
        }

        i = (unsigned char)obj->sizes[2];
        for (idx_col = 0; idx_col < i; idx_col++) {
          idxUpperExisting = obj->ldA * idx_col - 1;
          i1 = idx_col + 9;
          if (i1 >= 10) {
            memset(&obj->Aineq.data[idxUpperExisting + 10], 0, (unsigned int)(i1
                    - 9) * sizeof(double));
          }

          obj->Aineq.data[(idx_col + idxUpperExisting) + 10] = -1.0;
          i1 = idx_col + 11;
          i2 = obj->nVar;
          if (i1 <= i2) {
            memset(&obj->Aineq.data[i1 + idxUpperExisting], 0, (unsigned int)
                   ((((i2 + idxUpperExisting) - i1) - idxUpperExisting) + 1) *
                   sizeof(double));
          }
        }

        i = mIneq + 2;
        i1 = mIneq + 5;
        for (idx_col = 0; idx_col < 3; idx_col++) {
          int i3;
          int i4;
          idxUpperExisting = obj->ldA * idx_col - 1;
          colOffsetATw = idxUpperExisting + obj->ldA * (obj->isActiveIdx[1] - 1);
          if (offsetEq1 >= 10) {
            memset(&obj->Aeq.data[idxUpperExisting + 10], 0, (unsigned int)
                   (offsetEq1 - 9) * sizeof(double));
            memset(&obj->ATwset.data[colOffsetATw + 10], 0, (unsigned int)
                   (offsetEq1 - 9) * sizeof(double));
          }

          i2 = mIneq + idx_col;
          idx = i2 + 1;
          if (i <= idx) {
            memset(&obj->Aeq.data[i + idxUpperExisting], 0, (unsigned int)
                   ((((idx + idxUpperExisting) - i) - idxUpperExisting) + 1) *
                   sizeof(double));
            memset(&obj->ATwset.data[i + colOffsetATw], 0, (unsigned int)((((idx
                       + colOffsetATw) - i) - colOffsetATw) + 1) * sizeof(double));
          }

          idx = i2 + idxUpperExisting;
          obj->Aeq.data[idx + 2] = -1.0;
          i3 = i2 + colOffsetATw;
          obj->ATwset.data[i3 + 2] = -1.0;
          i4 = i2 + 3;
          if (i4 <= offsetEq2) {
            memset(&obj->Aeq.data[i4 + idxUpperExisting], 0, (unsigned int)
                   ((((offsetEq2 + idxUpperExisting) - i4) - idxUpperExisting) +
                    1) * sizeof(double));
            memset(&obj->ATwset.data[i4 + colOffsetATw], 0, (unsigned int)
                   ((((offsetEq2 + colOffsetATw) - i4) - colOffsetATw) + 1) *
                   sizeof(double));
          }

          i4 = i2 + 4;
          if (i1 <= i4) {
            memset(&obj->Aeq.data[i1 + idxUpperExisting], 0, (unsigned int)
                   ((((i4 + idxUpperExisting) - i1) - idxUpperExisting) + 1) *
                   sizeof(double));
            memset(&obj->ATwset.data[i1 + colOffsetATw], 0, (unsigned int)((((i4
                       + colOffsetATw) - i1) - colOffsetATw) + 1) * sizeof
                   (double));
          }

          obj->Aeq.data[idx + 5] = 1.0;
          obj->ATwset.data[i3 + 5] = 1.0;
          i2 += 6;
          idx = obj->nVar;
          if (i2 <= idx) {
            memset(&obj->Aeq.data[i2 + idxUpperExisting], 0, (unsigned int)
                   ((((idx + idxUpperExisting) - i2) - idxUpperExisting) + 1) *
                   sizeof(double));
            memset(&obj->ATwset.data[i2 + colOffsetATw], 0, (unsigned int)
                   ((((idx + colOffsetATw) - i2) - colOffsetATw) + 1) * sizeof
                   (double));
          }
        }

        idxUpperExisting = 9;
        i = obj->sizesNormal[3] + 1;
        i1 = obj->sizesRegularized[3];
        for (idx = i; idx <= i1; idx++) {
          idxUpperExisting++;
          obj->indexLB.data[idx - 1] = idxUpperExisting;
        }

        if (obj->nWConstr[4] > 0) {
          i = (unsigned char)obj->sizesRegularized[4];
          for (idx = 0; idx < i; idx++) {
            obj->isActiveConstr.data[obj->isActiveIdxRegularized[4] + idx] =
              obj->isActiveConstr.data[(obj->isActiveIdx[4] + idx) - 1];
          }
        }

        i = obj->isActiveIdx[4];
        i1 = obj->isActiveIdxRegularized[4];
        if (i <= i1 - 1) {
          memset(&obj->isActiveConstr.data[i + -1], 0, (unsigned int)(i1 - i) *
                 sizeof(bool));
        }

        i = obj->sizes[2] + 15;
        memset(&obj->lb.data[9], 0, (unsigned int)(i - 9) * sizeof(double));
        idxUpperExisting = obj->isActiveIdx[2];
        i = obj->nActiveConstr;
        for (idx_col = idxUpperExisting; idx_col <= i; idx_col++) {
          colOffsetATw = obj->ldA * (idx_col - 1) - 1;
          if (obj->Wid.data[idx_col - 1] == 3) {
            i1 = obj->Wlocalidx.data[idx_col - 1];
            i2 = i1 + 8;
            if (i2 >= 10) {
              memset(&obj->ATwset.data[colOffsetATw + 10], 0, (unsigned int)(i2
                      - 9) * sizeof(double));
            }

            obj->ATwset.data[(i1 + colOffsetATw) + 9] = -1.0;
            i1 += 10;
            i2 = obj->nVar;
            if (i1 <= i2) {
              memset(&obj->ATwset.data[i1 + colOffsetATw], 0, (unsigned int)
                     ((((i2 + colOffsetATw) - i1) - colOffsetATw) + 1) * sizeof
                     (double));
            }
          } else {
            i1 = obj->nVar;
            memset(&obj->ATwset.data[colOffsetATw + 10], 0, (unsigned int)(i1 -
                    9) * sizeof(double));
          }
        }
      }

      for (i = 0; i < 6; i++) {
        obj->isActiveIdx[i] = obj->isActiveIdxRegularized[i];
      }
    }
    break;

   default:
    obj->nVar = obj->nVarMax;
    obj->mConstr = obj->mConstrMax;
    for (i = 0; i < 5; i++) {
      obj->sizes[i] = obj->sizesRegPhaseOne[i];
    }

    modifyOverheadPhaseOne_(obj);
    for (i = 0; i < 6; i++) {
      obj->isActiveIdx[i] = obj->isActiveIdxRegPhaseOne[i];
    }
    break;
  }

  obj->probType = PROBLEM_TYPE;
}

/*
 * Arguments    : const double Hessian[81]
 *                const double grad_data[]
 *                c_struct_T *TrialState
 *                g_struct_T *memspace
 *                h_struct_T *WorkingSet
 *                d_struct_T *QRManager
 *                e_struct_T *CholManager
 *                f_struct_T *QPObjective
 *                const k_struct_T *qpoptions
 * Return Type  : bool
 */
static bool soc(const double Hessian[81], const double grad_data[], c_struct_T
                *TrialState, g_struct_T *memspace, h_struct_T *WorkingSet,
                d_struct_T *QRManager, e_struct_T *CholManager, f_struct_T
                *QPObjective, const k_struct_T *qpoptions)
{
  double c;
  int i;
  int i1;
  int ia;
  int idxIneqOffset;
  int idx_Aineq;
  int idx_Partition;
  int idx_lower;
  int iy;
  int mConstrMax;
  int mLB;
  int nWIneq_old;
  int nWLower_old;
  int nWUpper_old;
  bool success;
  nWIneq_old = WorkingSet->nWConstr[2];
  nWLower_old = WorkingSet->nWConstr[3];
  nWUpper_old = WorkingSet->nWConstr[4];
  mLB = WorkingSet->nVar;
  mConstrMax = WorkingSet->mConstrMax;
  i = (unsigned char)WorkingSet->nVar;
  memcpy(&TrialState->xstarsqp[0], &TrialState->xstarsqp_old[0], (unsigned int)i
         * sizeof(double));
  memcpy(&TrialState->socDirection.data[0], &TrialState->xstar.data[0],
         (unsigned int)i * sizeof(double));
  if (mConstrMax - 1 >= 0) {
    memcpy(&TrialState->lambdaStopTest.data[0], &TrialState->lambda.data[0],
           (unsigned int)mConstrMax * sizeof(double));
  }

  idxIneqOffset = WorkingSet->isActiveIdx[2];
  WorkingSet->beq[0] = -TrialState->cEq[0];
  WorkingSet->beq[1] = -TrialState->cEq[1];
  WorkingSet->beq[2] = -TrialState->cEq[2];
  idx_lower = WorkingSet->ldA;
  iy = 0;
  i1 = (WorkingSet->ldA << 1) + 1;
  for (idx_Partition = 1; idx_lower < 0 ? idx_Partition >= i1 : idx_Partition <=
       i1; idx_Partition += idx_lower) {
    c = 0.0;
    idx_Aineq = idx_Partition + WorkingSet->nVar;
    for (ia = idx_Partition; ia < idx_Aineq; ia++) {
      c += WorkingSet->Aeq.data[ia - 1] * TrialState->searchDir.data[ia -
        idx_Partition];
    }

    WorkingSet->beq[iy] += c;
    iy++;
  }

  WorkingSet->bwset.data[WorkingSet->sizes[0]] = WorkingSet->beq[0];
  WorkingSet->bwset.data[WorkingSet->sizes[0] + 1] = WorkingSet->beq[1];
  WorkingSet->bwset.data[WorkingSet->sizes[0] + 2] = WorkingSet->beq[2];
  if (WorkingSet->sizes[2] > 0) {
    i1 = (unsigned char)WorkingSet->sizes[2];
    for (ia = 0; ia < i1; ia++) {
      WorkingSet->bineq.data[ia] = -TrialState->cIneq.data[ia];
    }

    iy = 0;
    i1 = WorkingSet->ldA * (WorkingSet->sizes[2] - 1) + 1;
    for (idx_Partition = 1; idx_lower < 0 ? idx_Partition >= i1 : idx_Partition <=
         i1; idx_Partition += idx_lower) {
      c = 0.0;
      idx_Aineq = idx_Partition + WorkingSet->nVar;
      for (ia = idx_Partition; ia < idx_Aineq; ia++) {
        c += WorkingSet->Aineq.data[ia - 1] * TrialState->searchDir.data[ia -
          idx_Partition];
      }

      WorkingSet->bineq.data[iy] += c;
      iy++;
    }

    idx_Aineq = 1;
    idx_lower = WorkingSet->sizes[2] + 1;
    iy = (WorkingSet->sizes[2] + WorkingSet->sizes[3]) + 1;
    i1 = WorkingSet->nActiveConstr;
    for (ia = idxIneqOffset; ia <= i1; ia++) {
      switch (WorkingSet->Wid.data[ia - 1]) {
       case 3:
        idx_Partition = idx_Aineq;
        idx_Aineq++;
        WorkingSet->bwset.data[ia - 1] = WorkingSet->bineq.data
          [WorkingSet->Wlocalidx.data[ia - 1] - 1];
        break;

       case 4:
        idx_Partition = idx_lower;
        idx_lower++;
        break;

       default:
        idx_Partition = iy;
        iy++;
        break;
      }

      TrialState->workingset_old.data[idx_Partition - 1] =
        WorkingSet->Wlocalidx.data[ia - 1];
    }
  }

  memcpy(&TrialState->xstar.data[0], &TrialState->xstarsqp[0], (unsigned int)i *
         sizeof(double));
  b_driver(Hessian, grad_data, TrialState, memspace, WorkingSet, QRManager,
           CholManager, QPObjective, qpoptions, qpoptions);
  while ((WorkingSet->mEqRemoved > 0) && (WorkingSet->indexEqRemoved
          [WorkingSet->mEqRemoved - 1] >= 1)) {
    addAeqConstr(WorkingSet, WorkingSet->indexEqRemoved[WorkingSet->mEqRemoved -
                 1]);
    WorkingSet->mEqRemoved--;
  }

  i = (unsigned char)mLB;
  for (ia = 0; ia < i; ia++) {
    double oldDirIdx;
    c = TrialState->socDirection.data[ia];
    oldDirIdx = c;
    c = TrialState->xstar.data[ia] - c;
    TrialState->socDirection.data[ia] = c;
    TrialState->xstar.data[ia] = oldDirIdx;
  }

  success = (c_xnrm2(mLB, TrialState->socDirection.data) <= 2.0 * c_xnrm2(mLB,
              TrialState->xstar.data));
  idx_Partition = WorkingSet->sizes[2];
  mLB = WorkingSet->sizes[3];
  WorkingSet->beq[0] = -TrialState->cEq[0];
  WorkingSet->beq[1] = -TrialState->cEq[1];
  WorkingSet->beq[2] = -TrialState->cEq[2];
  WorkingSet->bwset.data[WorkingSet->sizes[0]] = WorkingSet->beq[0];
  WorkingSet->bwset.data[WorkingSet->sizes[0] + 1] = WorkingSet->beq[1];
  WorkingSet->bwset.data[WorkingSet->sizes[0] + 2] = WorkingSet->beq[2];
  if (WorkingSet->sizes[2] > 0) {
    i = (unsigned char)WorkingSet->sizes[2];
    for (ia = 0; ia < i; ia++) {
      WorkingSet->bineq.data[ia] = -TrialState->cIneq.data[ia];
    }

    if (!success) {
      iy = WorkingSet->nWConstr[0] + WorkingSet->nWConstr[1];
      idx_Aineq = iy + 1;
      idx_lower = WorkingSet->nActiveConstr;
      for (idxIneqOffset = idx_Aineq; idxIneqOffset <= idx_lower; idxIneqOffset
           ++) {
        WorkingSet->isActiveConstr.data[(WorkingSet->isActiveIdx
          [WorkingSet->Wid.data[idxIneqOffset - 1] - 1] +
          WorkingSet->Wlocalidx.data[idxIneqOffset - 1]) - 2] = false;
      }

      WorkingSet->nWConstr[2] = 0;
      WorkingSet->nWConstr[3] = 0;
      WorkingSet->nWConstr[4] = 0;
      WorkingSet->nActiveConstr = iy;
      for (ia = 0; ia < nWIneq_old; ia++) {
        iy = TrialState->workingset_old.data[ia];
        WorkingSet->nWConstr[2]++;
        WorkingSet->isActiveConstr.data[(WorkingSet->isActiveIdx[2] + iy) - 2] =
          true;
        WorkingSet->nActiveConstr++;
        i = WorkingSet->nActiveConstr - 1;
        WorkingSet->Wid.data[i] = 3;
        WorkingSet->Wlocalidx.data[i] = iy;
        idx_Aineq = WorkingSet->ldA * (iy - 1);
        idx_lower = WorkingSet->ldA * i;
        i1 = WorkingSet->nVar;
        for (idxIneqOffset = 0; idxIneqOffset < i1; idxIneqOffset++) {
          WorkingSet->ATwset.data[idx_lower + idxIneqOffset] =
            WorkingSet->Aineq.data[idx_Aineq + idxIneqOffset];
        }

        WorkingSet->bwset.data[i] = WorkingSet->bineq.data[iy - 1];
      }

      for (ia = 0; ia < nWLower_old; ia++) {
        addBoundToActiveSetMatrix_(WorkingSet, 4,
          TrialState->workingset_old.data[ia + idx_Partition]);
      }

      for (ia = 0; ia < nWUpper_old; ia++) {
        addBoundToActiveSetMatrix_(WorkingSet, 5,
          TrialState->workingset_old.data[(ia + idx_Partition) + mLB]);
      }
    }
  }

  if (!success) {
    if (mConstrMax - 1 >= 0) {
      memcpy(&TrialState->lambda.data[0], &TrialState->lambdaStopTest.data[0],
             (unsigned int)mConstrMax * sizeof(double));
    }
  } else {
    sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
                 WorkingSet->sizes, WorkingSet->isActiveIdx,
                 WorkingSet->Wid.data, WorkingSet->Wlocalidx.data,
                 memspace->workspace_float.data);
  }

  return success;
}

/*
 * Arguments    : const e_struct_T *obj
 *                double rhs_data[]
 * Return Type  : void
 */
static void solve(const e_struct_T *obj, double rhs_data[])
{
  int i;
  int j;
  int jA;
  int n_tmp;
  n_tmp = obj->ndims;
  if (obj->ndims != 0) {
    for (j = 0; j < n_tmp; j++) {
      double temp;
      jA = j * obj->ldm;
      temp = rhs_data[j];
      for (i = 0; i < j; i++) {
        temp -= obj->FMat.data[jA + i] * rhs_data[i];
      }

      rhs_data[j] = temp / obj->FMat.data[jA + j];
    }
  }

  if (obj->ndims != 0) {
    for (j = n_tmp; j >= 1; j--) {
      jA = (j + (j - 1) * obj->ldm) - 1;
      rhs_data[j - 1] /= obj->FMat.data[jA];
      for (i = 0; i <= j - 2; i++) {
        int ix;
        ix = (j - i) - 2;
        rhs_data[ix] -= rhs_data[j - 1] * obj->FMat.data[(jA - i) - 1];
      }
    }
  }
}

/*
 * Arguments    : double lambda_data[]
 *                int WorkingSet_nActiveConstr
 *                const int WorkingSet_sizes[5]
 *                const int WorkingSet_isActiveIdx[6]
 *                const int WorkingSet_Wid_data[]
 *                const int WorkingSet_Wlocalidx_data[]
 *                double workspace_data[]
 * Return Type  : void
 */
static void sortLambdaQP(double lambda_data[], int WorkingSet_nActiveConstr,
  const int WorkingSet_sizes[5], const int WorkingSet_isActiveIdx[6], const int
  WorkingSet_Wid_data[], const int WorkingSet_Wlocalidx_data[], double
  workspace_data[])
{
  if (WorkingSet_nActiveConstr != 0) {
    int idx;
    int idxOffset;
    int mAll;
    mAll = (((WorkingSet_sizes[0] + WorkingSet_sizes[3]) + WorkingSet_sizes[4])
            + WorkingSet_sizes[2]) + 2;
    idx = (unsigned char)(mAll + 1);
    if (idx - 1 >= 0) {
      memcpy(&workspace_data[0], &lambda_data[0], (unsigned int)idx * sizeof
             (double));
    }

    if (mAll >= 0) {
      memset(&lambda_data[0], 0, (unsigned int)(mAll + 1) * sizeof(double));
    }

    mAll = 0;
    idx = 0;
    while ((idx + 1 <= WorkingSet_nActiveConstr) && (WorkingSet_Wid_data[idx] <=
            2)) {
      if (WorkingSet_Wid_data[idx] == 1) {
        idxOffset = 1;
      } else {
        idxOffset = WorkingSet_isActiveIdx[1];
      }

      lambda_data[(idxOffset + WorkingSet_Wlocalidx_data[idx]) - 2] =
        workspace_data[mAll];
      mAll++;
      idx++;
    }

    while (idx + 1 <= WorkingSet_nActiveConstr) {
      switch (WorkingSet_Wid_data[idx]) {
       case 3:
        idxOffset = WorkingSet_isActiveIdx[2];
        break;

       case 4:
        idxOffset = WorkingSet_isActiveIdx[3];
        break;

       default:
        idxOffset = WorkingSet_isActiveIdx[4];
        break;
      }

      lambda_data[(idxOffset + WorkingSet_Wlocalidx_data[idx]) - 2] =
        workspace_data[mAll];
      mAll++;
      idx++;
    }
  }
}

/*
 * Arguments    : d_struct_T *obj
 *                const double vec_data[]
 *                int iv0
 * Return Type  : void
 */
static void squareQ_appendCol(d_struct_T *obj, const double vec_data[], int iv0)
{
  double c;
  double s;
  double temp;
  int Qk0;
  int i;
  int ia;
  int iac;
  int idx;
  int iy;
  int iyend;
  int k;
  iyend = obj->mrows;
  Qk0 = obj->ncols + 1;
  if (iyend <= Qk0) {
    Qk0 = iyend;
  }

  obj->minRowCol = Qk0;
  iy = obj->ldq * obj->ncols;
  Qk0 = obj->ldq;
  if (obj->mrows != 0) {
    iyend = iy + obj->mrows;
    if (iy + 1 <= iyend) {
      memset(&obj->QR.data[iy], 0, (unsigned int)(iyend - iy) * sizeof(double));
    }

    i = obj->ldq * (obj->mrows - 1) + 1;
    for (iac = 1; Qk0 < 0 ? iac >= i : iac <= i; iac += Qk0) {
      c = 0.0;
      iyend = iac + obj->mrows;
      for (ia = iac; ia < iyend; ia++) {
        c += obj->Q.data[ia - 1] * vec_data[((iv0 + ia) - iac) - 1];
      }

      obj->QR.data[iy] += c;
      iy++;
    }
  }

  obj->ncols++;
  i = obj->ncols - 1;
  obj->jpvt.data[i] = obj->ncols;
  for (idx = obj->mrows - 2; idx + 2 > obj->ncols; idx--) {
    iyend = idx + obj->ldq * i;
    temp = obj->QR.data[iyend + 1];
    c = xrotg(&obj->QR.data[iyend], &temp, &s);
    obj->QR.data[iyend + 1] = temp;
    Qk0 = obj->ldq * idx;
    iyend = obj->mrows;
    if (obj->mrows >= 1) {
      iy = obj->ldq + Qk0;
      for (k = 0; k < iyend; k++) {
        iac = iy + k;
        ia = Qk0 + k;
        temp = c * obj->Q.data[ia] + s * obj->Q.data[iac];
        obj->Q.data[iac] = c * obj->Q.data[iac] - s * obj->Q.data[ia];
        obj->Q.data[ia] = temp;
      }
    }
  }
}

/*
 * Arguments    : struct_T *stepFlags
 *                double Hessian[81]
 *                c_struct_T *TrialState
 *                j_struct_T *MeritFunction
 *                g_struct_T *memspace
 *                h_struct_T *WorkingSet
 *                d_struct_T *QRManager
 *                e_struct_T *CholManager
 *                f_struct_T *QPObjective
 *                k_struct_T *qpoptions
 * Return Type  : void
 */
static void step(struct_T *stepFlags, double Hessian[81], c_struct_T *TrialState,
                 j_struct_T *MeritFunction, g_struct_T *memspace, h_struct_T
                 *WorkingSet, d_struct_T *QRManager, e_struct_T *CholManager,
                 f_struct_T *QPObjective, k_struct_T *qpoptions)
{
  h_struct_T obj;
  double y_data[40];
  double tmp_data[22];
  double constrViolationIneq;
  double linearizedConstrViolPrev;
  int iH0;
  int idxEndIneq;
  int idxStartIneq;
  int idx_global;
  int nVar;
  bool checkBoundViolation;
  stepFlags->stepAccepted = true;
  checkBoundViolation = true;
  nVar = WorkingSet->nVar - 1;
  if (stepFlags->stepType != 3) {
    idx_global = (unsigned char)WorkingSet->nVar;
    memcpy(&TrialState->xstar.data[0], &TrialState->xstarsqp[0], (unsigned int)
           idx_global * sizeof(double));
  } else if (nVar >= 0) {
    memcpy(&TrialState->searchDir.data[0], &TrialState->xstar.data[0], (unsigned
            int)(nVar + 1) * sizeof(double));
  }

  int exitg1;
  bool guard1;
  do {
    exitg1 = 0;
    guard1 = false;
    switch (stepFlags->stepType) {
     case 1:
      {
        bool nonlinEqRemoved;
        b_driver(Hessian, TrialState->grad.data, TrialState, memspace,
                 WorkingSet, QRManager, CholManager, QPObjective, qpoptions,
                 qpoptions);
        if (WorkingSet->probType == 2) {
          obj = *WorkingSet;
          constrViolationIneq = c_maxConstraintViolation_AMats_(&obj,
            TrialState->xstar.data);
        } else {
          obj = *WorkingSet;
          constrViolationIneq = d_maxConstraintViolation_AMats_(&obj,
            TrialState->xstar.data);
        }

        if (WorkingSet->sizes[3] > 0) {
          idx_global = (unsigned char)WorkingSet->sizes[3];
          for (iH0 = 0; iH0 < idx_global; iH0++) {
            constrViolationIneq = fmax(constrViolationIneq,
              -TrialState->xstar.data[obj.indexLB.data[iH0] - 1] -
              obj.lb.data[obj.indexLB.data[iH0] - 1]);
          }
        }

        if (WorkingSet->sizes[4] > 0) {
          idx_global = (unsigned char)WorkingSet->sizes[4];
          for (iH0 = 0; iH0 < idx_global; iH0++) {
            constrViolationIneq = fmax(constrViolationIneq,
              TrialState->xstar.data[obj.indexUB.data[iH0] - 1] -
              obj.ub.data[obj.indexUB.data[iH0] - 1]);
          }
        }

        if (WorkingSet->sizes[0] > 0) {
          idx_global = (unsigned char)WorkingSet->sizes[0];
          for (iH0 = 0; iH0 < idx_global; iH0++) {
            constrViolationIneq = fmax(constrViolationIneq, fabs
              (TrialState->xstar.data[obj.indexFixed.data[iH0] - 1] -
               obj.ub.data[obj.indexFixed.data[iH0] - 1]));
          }
        }

        if ((TrialState->state > 0) || ((TrialState->state == 0) &&
             (constrViolationIneq <= 1.0E-6))) {
          double constrViolation;
          double penaltyParamTrial;
          penaltyParamTrial = MeritFunction->penaltyParam;
          constrViolationIneq = 0.0;
          idx_global = (unsigned char)WorkingSet->sizes[2];
          for (iH0 = 0; iH0 < idx_global; iH0++) {
            linearizedConstrViolPrev = TrialState->cIneq.data[iH0];
            if (linearizedConstrViolPrev > 0.0) {
              constrViolationIneq += linearizedConstrViolPrev;
            }
          }

          constrViolation = ((fabs(TrialState->cEq[0]) + fabs(TrialState->cEq[1]))
                             + fabs(TrialState->cEq[2])) + constrViolationIneq;
          linearizedConstrViolPrev = MeritFunction->linearizedConstrViol;
          MeritFunction->linearizedConstrViol = 0.0;
          constrViolationIneq = constrViolation + linearizedConstrViolPrev;
          if ((constrViolationIneq > 2.2204460492503131E-16) &&
              (TrialState->fstar > 0.0)) {
            if (TrialState->sqpFval == 0.0) {
              linearizedConstrViolPrev = 1.0;
            } else {
              linearizedConstrViolPrev = 1.5;
            }

            penaltyParamTrial = linearizedConstrViolPrev * TrialState->fstar /
              constrViolationIneq;
          }

          if (penaltyParamTrial < MeritFunction->penaltyParam) {
            MeritFunction->phi = TrialState->sqpFval + penaltyParamTrial *
              constrViolation;
            if ((MeritFunction->initFval + penaltyParamTrial *
                 (MeritFunction->initConstrViolationEq +
                  MeritFunction->initConstrViolationIneq)) - MeritFunction->phi >
                (double)MeritFunction->nPenaltyDecreases *
                MeritFunction->threshold) {
              MeritFunction->nPenaltyDecreases++;
              if ((MeritFunction->nPenaltyDecreases << 1) >
                  TrialState->sqpIterations) {
                MeritFunction->threshold *= 10.0;
              }

              MeritFunction->penaltyParam = fmax(penaltyParamTrial, 1.0E-10);
            } else {
              MeritFunction->phi = TrialState->sqpFval +
                MeritFunction->penaltyParam * constrViolation;
            }
          } else {
            MeritFunction->penaltyParam = fmax(penaltyParamTrial, 1.0E-10);
            MeritFunction->phi = TrialState->sqpFval +
              MeritFunction->penaltyParam * constrViolation;
          }

          MeritFunction->phiPrimePlus = fmin(TrialState->fstar -
            MeritFunction->penaltyParam * constrViolation, 0.0);
        } else if (TrialState->state != -6) {
          stepFlags->stepType = 2;
        }

        sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
                     WorkingSet->sizes, WorkingSet->isActiveIdx,
                     WorkingSet->Wid.data, WorkingSet->Wlocalidx.data,
                     memspace->workspace_float.data);
        nonlinEqRemoved = (WorkingSet->mEqRemoved > 0);
        while ((WorkingSet->mEqRemoved > 0) && (WorkingSet->
                indexEqRemoved[WorkingSet->mEqRemoved - 1] >= 1)) {
          addAeqConstr(WorkingSet, WorkingSet->indexEqRemoved
                       [WorkingSet->mEqRemoved - 1]);
          WorkingSet->mEqRemoved--;
        }

        if (nonlinEqRemoved) {
          WorkingSet->Wlocalidx.data[WorkingSet->sizes[0]] = 1;
          WorkingSet->Wlocalidx.data[WorkingSet->sizes[0] + 1] = 2;
          WorkingSet->Wlocalidx.data[WorkingSet->sizes[0] + 2] = 3;
        }

        if (stepFlags->stepType != 2) {
          idxStartIneq = TrialState->delta_x.size[0];
          idxEndIneq = TrialState->delta_x.size[0];
          if (idxEndIneq - 1 >= 0) {
            memcpy(&y_data[0], &TrialState->delta_x.data[0], (unsigned int)
                   idxEndIneq * sizeof(double));
          }

          if (nVar >= 0) {
            memcpy(&y_data[0], &TrialState->xstar.data[0], (unsigned int)(nVar +
                    1) * sizeof(double));
          }

          if (idxStartIneq - 1 >= 0) {
            memcpy(&TrialState->delta_x.data[0], &y_data[0], (unsigned int)
                   idxStartIneq * sizeof(double));
          }

          guard1 = true;
        }
      }
      break;

     case 2:
      iH0 = WorkingSet->nWConstr[0] + WorkingSet->nWConstr[1];
      idxStartIneq = iH0 + 1;
      idxEndIneq = WorkingSet->nActiveConstr;
      for (idx_global = idxStartIneq; idx_global <= idxEndIneq; idx_global++) {
        WorkingSet->isActiveConstr.data[(WorkingSet->isActiveIdx
          [WorkingSet->Wid.data[idx_global - 1] - 1] +
          WorkingSet->Wlocalidx.data[idx_global - 1]) - 2] = false;
      }

      WorkingSet->nWConstr[2] = 0;
      WorkingSet->nWConstr[3] = 0;
      WorkingSet->nWConstr[4] = 0;
      WorkingSet->nActiveConstr = iH0;
      idxStartIneq = TrialState->xstar.size[0];
      idxEndIneq = TrialState->xstar.size[0];
      if (idxEndIneq - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->xstar.data[0], (unsigned int)idxEndIneq *
               sizeof(double));
      }

      idx_global = (unsigned char)WorkingSet->sizes[3];
      for (iH0 = 0; iH0 < idx_global; iH0++) {
        linearizedConstrViolPrev = WorkingSet->lb.data[WorkingSet->
          indexLB.data[iH0] - 1];
        if (-y_data[WorkingSet->indexLB.data[iH0] - 1] >
            linearizedConstrViolPrev) {
          y_data[WorkingSet->indexLB.data[iH0] - 1] = (WorkingSet->
            ub.data[WorkingSet->indexLB.data[iH0] - 1] -
            linearizedConstrViolPrev) / 2.0;
        }
      }

      idx_global = (unsigned char)WorkingSet->sizes[4];
      for (iH0 = 0; iH0 < idx_global; iH0++) {
        linearizedConstrViolPrev = WorkingSet->ub.data[WorkingSet->
          indexUB.data[iH0] - 1];
        if (y_data[WorkingSet->indexUB.data[iH0] - 1] > linearizedConstrViolPrev)
        {
          y_data[WorkingSet->indexUB.data[iH0] - 1] = (linearizedConstrViolPrev
            - WorkingSet->lb.data[WorkingSet->indexUB.data[iH0] - 1]) / 2.0;
        }
      }

      if (idxStartIneq - 1 >= 0) {
        memcpy(&TrialState->xstar.data[0], &y_data[0], (unsigned int)
               idxStartIneq * sizeof(double));
      }

      relaxed(Hessian, TrialState->grad.data, TrialState, MeritFunction,
              memspace, WorkingSet, QRManager, CholManager, QPObjective,
              qpoptions);
      idxStartIneq = TrialState->delta_x.size[0];
      idxEndIneq = TrialState->delta_x.size[0];
      if (idxEndIneq - 1 >= 0) {
        memcpy(&y_data[0], &TrialState->delta_x.data[0], (unsigned int)
               idxEndIneq * sizeof(double));
      }

      if (nVar >= 0) {
        memcpy(&y_data[0], &TrialState->xstar.data[0], (unsigned int)(nVar + 1) *
               sizeof(double));
      }

      if (idxStartIneq - 1 >= 0) {
        memcpy(&TrialState->delta_x.data[0], &y_data[0], (unsigned int)
               idxStartIneq * sizeof(double));
      }

      guard1 = true;
      break;

     default:
      idxStartIneq = TrialState->grad.size[0];
      if (idxStartIneq - 1 >= 0) {
        memcpy(&tmp_data[0], &TrialState->grad.data[0], (unsigned int)
               idxStartIneq * sizeof(double));
      }

      checkBoundViolation = soc(Hessian, tmp_data, TrialState, memspace,
        WorkingSet, QRManager, CholManager, QPObjective, qpoptions);
      stepFlags->stepAccepted = checkBoundViolation;
      if (stepFlags->stepAccepted && (TrialState->state != -6)) {
        idx_global = (unsigned char)(nVar + 1);
        for (iH0 = 0; iH0 < idx_global; iH0++) {
          TrialState->delta_x.data[iH0] = TrialState->xstar.data[iH0] +
            TrialState->socDirection.data[iH0];
        }
      }

      guard1 = true;
      break;
    }

    if (guard1) {
      if (TrialState->state != -6) {
        exitg1 = 1;
      } else {
        constrViolationIneq = 0.0;
        linearizedConstrViolPrev = 1.0;
        for (iH0 = 0; iH0 < 9; iH0++) {
          constrViolationIneq = fmax(constrViolationIneq, fabs
            (TrialState->grad.data[iH0]));
          linearizedConstrViolPrev = fmax(linearizedConstrViolPrev, fabs
            (TrialState->xstar.data[iH0]));
        }

        constrViolationIneq = fmax(2.2204460492503131E-16, constrViolationIneq /
          linearizedConstrViolPrev);
        for (idxEndIneq = 0; idxEndIneq < 9; idxEndIneq++) {
          iH0 = 9 * idxEndIneq;
          for (idxStartIneq = 0; idxStartIneq < idxEndIneq; idxStartIneq++) {
            Hessian[iH0 + idxStartIneq] = 0.0;
          }

          iH0 = idxEndIneq + 9 * idxEndIneq;
          Hessian[iH0] = constrViolationIneq;
          idx_global = 7 - idxEndIneq;
          if (idx_global >= 0) {
            memset(&Hessian[iH0 + 1], 0, (unsigned int)(idx_global + 1) * sizeof
                   (double));
          }
        }
      }
    }
  } while (exitg1 == 0);

  if (checkBoundViolation) {
    idxStartIneq = TrialState->delta_x.size[0];
    idxEndIneq = TrialState->delta_x.size[0];
    if (idxEndIneq - 1 >= 0) {
      memcpy(&tmp_data[0], &TrialState->delta_x.data[0], (unsigned int)
             idxEndIneq * sizeof(double));
    }

    idxEndIneq = TrialState->xstar.size[0];
    iH0 = TrialState->xstar.size[0];
    if (iH0 - 1 >= 0) {
      memcpy(&y_data[0], &TrialState->xstar.data[0], (unsigned int)iH0 * sizeof
             (double));
    }

    idx_global = (unsigned char)WorkingSet->sizes[3];
    for (iH0 = 0; iH0 < idx_global; iH0++) {
      constrViolationIneq = tmp_data[WorkingSet->indexLB.data[iH0] - 1];
      linearizedConstrViolPrev = (TrialState->xstarsqp[WorkingSet->
        indexLB.data[iH0] - 1] + constrViolationIneq) - -1.0;
      if (linearizedConstrViolPrev < 0.0) {
        tmp_data[WorkingSet->indexLB.data[iH0] - 1] = constrViolationIneq -
          linearizedConstrViolPrev;
        y_data[WorkingSet->indexLB.data[iH0] - 1] -= linearizedConstrViolPrev;
      }
    }

    idx_global = (unsigned char)WorkingSet->sizes[4];
    for (iH0 = 0; iH0 < idx_global; iH0++) {
      constrViolationIneq = tmp_data[WorkingSet->indexUB.data[iH0] - 1];
      linearizedConstrViolPrev = ((double)iv[WorkingSet->indexUB.data[iH0] - 1]
        - TrialState->xstarsqp[WorkingSet->indexUB.data[iH0] - 1]) -
        constrViolationIneq;
      if (linearizedConstrViolPrev < 0.0) {
        tmp_data[WorkingSet->indexUB.data[iH0] - 1] = constrViolationIneq +
          linearizedConstrViolPrev;
        y_data[WorkingSet->indexUB.data[iH0] - 1] += linearizedConstrViolPrev;
      }
    }

    if (idxStartIneq - 1 >= 0) {
      memcpy(&TrialState->delta_x.data[0], &tmp_data[0], (unsigned int)
             idxStartIneq * sizeof(double));
    }

    if (idxEndIneq - 1 >= 0) {
      memcpy(&TrialState->xstar.data[0], &y_data[0], (unsigned int)idxEndIneq *
             sizeof(double));
    }
  }
}

/*
 * Arguments    : const double A[9]
 *                double U[9]
 *                double s[3]
 *                double V[9]
 * Return Type  : void
 */
static void svd(const double A[9], double U[9], double s[3], double V[9])
{
  double Vf[9];
  double b_A[9];
  double e[3];
  double work[3];
  double anrm;
  double cscale;
  double nrm;
  double rt;
  double snorm;
  double sqds;
  int ii;
  int jj;
  int k;
  int m;
  int q;
  int qjj;
  int qp1;
  int qq;
  int qs;
  bool doscale;
  s[0] = 0.0;
  e[0] = 0.0;
  work[0] = 0.0;
  s[1] = 0.0;
  e[1] = 0.0;
  work[1] = 0.0;
  s[2] = 0.0;
  e[2] = 0.0;
  work[2] = 0.0;
  for (qjj = 0; qjj < 9; qjj++) {
    b_A[qjj] = A[qjj];
    U[qjj] = 0.0;
    Vf[qjj] = 0.0;
  }

  doscale = false;
  anrm = xzlangeM(A);
  cscale = anrm;
  if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
    doscale = true;
    cscale = 6.7178761075670888E-139;
    xzlascl(anrm, cscale, b_A);
  } else if (anrm > 1.4885657073574029E+138) {
    doscale = true;
    cscale = 1.4885657073574029E+138;
    xzlascl(anrm, cscale, b_A);
  }

  for (q = 0; q < 2; q++) {
    bool apply_transform;
    qp1 = q + 2;
    qs = q + 3 * q;
    qq = qs + 1;
    apply_transform = false;
    nrm = xnrm2(3 - q, b_A, qs + 1);
    if (nrm > 0.0) {
      apply_transform = true;
      if (b_A[qs] < 0.0) {
        nrm = -nrm;
      }

      s[q] = nrm;
      if (fabs(nrm) >= 1.0020841800044864E-292) {
        nrm = 1.0 / nrm;
        qjj = (qs - q) + 3;
        for (k = qq; k <= qjj; k++) {
          b_A[k - 1] *= nrm;
        }
      } else {
        qjj = (qs - q) + 3;
        for (k = qq; k <= qjj; k++) {
          b_A[k - 1] /= s[q];
        }
      }

      b_A[qs]++;
      s[q] = -s[q];
    } else {
      s[q] = 0.0;
    }

    for (jj = qp1; jj < 4; jj++) {
      qjj = q + 3 * (jj - 1);
      if (apply_transform) {
        xaxpy(3 - q, -(xdotc(3 - q, b_A, qs + 1, b_A, qjj + 1) / b_A[qs]), qs +
              1, b_A, qjj + 1);
      }

      e[jj - 1] = b_A[qjj];
    }

    for (ii = q + 1; ii < 4; ii++) {
      qjj = (ii + 3 * q) - 1;
      U[qjj] = b_A[qjj];
    }

    if (q + 1 <= 1) {
      nrm = d_xnrm2(e);
      if (nrm == 0.0) {
        e[0] = 0.0;
      } else {
        if (e[1] < 0.0) {
          e[0] = -nrm;
        } else {
          e[0] = nrm;
        }

        nrm = e[0];
        if (fabs(e[0]) >= 1.0020841800044864E-292) {
          nrm = 1.0 / e[0];
          for (k = qp1; k < 4; k++) {
            e[k - 1] *= nrm;
          }
        } else {
          for (k = qp1; k < 4; k++) {
            e[k - 1] /= nrm;
          }
        }

        e[1]++;
        e[0] = -e[0];
        for (ii = qp1; ii < 4; ii++) {
          work[ii - 1] = 0.0;
        }

        for (jj = qp1; jj < 4; jj++) {
          b_xaxpy(e[jj - 1], b_A, 3 * (jj - 1) + 2, work);
        }

        for (jj = qp1; jj < 4; jj++) {
          c_xaxpy(-e[jj - 1] / e[1], work, b_A, 3 * (jj - 1) + 2);
        }
      }

      for (ii = qp1; ii < 4; ii++) {
        Vf[ii - 1] = e[ii - 1];
      }
    }
  }

  m = 1;
  s[2] = b_A[8];
  e[1] = b_A[7];
  e[2] = 0.0;
  U[6] = 0.0;
  U[7] = 0.0;
  U[8] = 1.0;
  for (q = 1; q >= 0; q--) {
    qp1 = q + 2;
    qq = q + 3 * q;
    if (s[q] != 0.0) {
      for (jj = qp1; jj < 4; jj++) {
        qjj = (q + 3 * (jj - 1)) + 1;
        xaxpy(3 - q, -(xdotc(3 - q, U, qq + 1, U, qjj) / U[qq]), qq + 1, U, qjj);
      }

      for (ii = q + 1; ii < 4; ii++) {
        qjj = (ii + 3 * q) - 1;
        U[qjj] = -U[qjj];
      }

      U[qq]++;
      if (q - 1 >= 0) {
        U[3 * q] = 0.0;
      }
    } else {
      U[3 * q] = 0.0;
      U[3 * q + 1] = 0.0;
      U[3 * q + 2] = 0.0;
      U[qq] = 1.0;
    }
  }

  for (q = 2; q >= 0; q--) {
    if ((q + 1 <= 1) && (e[0] != 0.0)) {
      xaxpy(2, -(xdotc(2, Vf, 2, Vf, 5) / Vf[1]), 2, Vf, 5);
      xaxpy(2, -(xdotc(2, Vf, 2, Vf, 8) / Vf[1]), 2, Vf, 8);
    }

    Vf[3 * q] = 0.0;
    Vf[3 * q + 1] = 0.0;
    Vf[3 * q + 2] = 0.0;
    Vf[q + 3 * q] = 1.0;
  }

  qp1 = 0;
  snorm = 0.0;
  for (q = 0; q < 3; q++) {
    nrm = s[q];
    if (nrm != 0.0) {
      rt = fabs(nrm);
      nrm /= rt;
      s[q] = rt;
      if (q + 1 < 3) {
        e[q] /= nrm;
      }

      jj = 3 * q;
      qjj = jj + 3;
      for (k = jj + 1; k <= qjj; k++) {
        U[k - 1] *= nrm;
      }
    }

    if (q + 1 < 3) {
      nrm = e[q];
      if (nrm != 0.0) {
        rt = fabs(nrm);
        nrm = rt / nrm;
        e[q] = rt;
        s[q + 1] *= nrm;
        jj = 3 * (q + 1);
        qjj = jj + 3;
        for (k = jj + 1; k <= qjj; k++) {
          Vf[k - 1] *= nrm;
        }
      }
    }

    snorm = fmax(snorm, fmax(fabs(s[q]), fabs(e[q])));
  }

  while ((m + 2 > 0) && (qp1 < 75)) {
    bool exitg1;
    qq = m + 1;
    ii = m + 1;
    exitg1 = false;
    while (!(exitg1 || (ii == 0))) {
      nrm = fabs(e[ii - 1]);
      if ((nrm <= 2.2204460492503131E-16 * (fabs(s[ii - 1]) + fabs(s[ii]))) ||
          (nrm <= 1.0020841800044864E-292) || ((qp1 > 20) && (nrm <=
            2.2204460492503131E-16 * snorm))) {
        e[ii - 1] = 0.0;
        exitg1 = true;
      } else {
        ii--;
      }
    }

    if (ii == m + 1) {
      qjj = 4;
    } else {
      qs = m + 2;
      jj = m + 2;
      exitg1 = false;
      while ((!exitg1) && (jj >= ii)) {
        qs = jj;
        if (jj == ii) {
          exitg1 = true;
        } else {
          nrm = 0.0;
          if (jj < m + 2) {
            nrm = fabs(e[jj - 1]);
          }

          if (jj > ii + 1) {
            nrm += fabs(e[jj - 2]);
          }

          rt = fabs(s[jj - 1]);
          if ((rt <= 2.2204460492503131E-16 * nrm) || (rt <=
               1.0020841800044864E-292)) {
            s[jj - 1] = 0.0;
            exitg1 = true;
          } else {
            jj--;
          }
        }
      }

      if (qs == ii) {
        qjj = 3;
      } else if (qs == m + 2) {
        qjj = 1;
      } else {
        qjj = 2;
        ii = qs;
      }
    }

    switch (qjj) {
     case 1:
      {
        rt = e[m];
        e[m] = 0.0;
        for (k = qq; k >= ii + 1; k--) {
          double sm;
          sm = xrotg(&s[k - 1], &rt, &sqds);
          if (k > ii + 1) {
            rt = -sqds * e[0];
            e[0] *= sm;
          }

          xrot(Vf, 3 * (k - 1) + 1, 3 * (m + 1) + 1, sm, sqds);
        }
      }
      break;

     case 2:
      {
        rt = e[ii - 1];
        e[ii - 1] = 0.0;
        for (k = ii + 1; k <= m + 2; k++) {
          double b;
          double sm;
          sm = xrotg(&s[k - 1], &rt, &sqds);
          b = e[k - 1];
          rt = -sqds * b;
          e[k - 1] = b * sm;
          xrot(U, 3 * (k - 1) + 1, 3 * (ii - 1) + 1, sm, sqds);
        }
      }
      break;

     case 3:
      {
        double b;
        double scale;
        double sm;
        nrm = s[m + 1];
        scale = fmax(fmax(fmax(fmax(fabs(nrm), fabs(s[m])), fabs(e[m])), fabs
                          (s[ii])), fabs(e[ii]));
        sm = nrm / scale;
        nrm = s[m] / scale;
        rt = e[m] / scale;
        sqds = s[ii] / scale;
        b = ((nrm + sm) * (nrm - sm) + rt * rt) / 2.0;
        nrm = sm * rt;
        nrm *= nrm;
        if ((b != 0.0) || (nrm != 0.0)) {
          rt = sqrt(b * b + nrm);
          if (b < 0.0) {
            rt = -rt;
          }

          rt = nrm / (b + rt);
        } else {
          rt = 0.0;
        }

        rt += (sqds + sm) * (sqds - sm);
        nrm = sqds * (e[ii] / scale);
        for (k = ii + 1; k <= qq; k++) {
          sm = xrotg(&rt, &nrm, &sqds);
          if (k > ii + 1) {
            e[0] = rt;
          }

          nrm = e[k - 1];
          b = s[k - 1];
          e[k - 1] = sm * nrm - sqds * b;
          rt = sqds * s[k];
          s[k] *= sm;
          qjj = 3 * (k - 1) + 1;
          jj = 3 * k + 1;
          xrot(Vf, qjj, jj, sm, sqds);
          s[k - 1] = sm * b + sqds * nrm;
          sm = xrotg(&s[k - 1], &rt, &sqds);
          b = e[k - 1];
          rt = sm * b + sqds * s[k];
          s[k] = -sqds * b + sm * s[k];
          nrm = sqds * e[k];
          e[k] *= sm;
          xrot(U, qjj, jj, sm, sqds);
        }

        e[m] = rt;
        qp1++;
      }
      break;

     default:
      if (s[ii] < 0.0) {
        s[ii] = -s[ii];
        jj = 3 * ii;
        qjj = jj + 3;
        for (k = jj + 1; k <= qjj; k++) {
          Vf[k - 1] = -Vf[k - 1];
        }
      }

      qp1 = ii + 1;
      while ((ii + 1 < 3) && (s[ii] < s[qp1])) {
        rt = s[ii];
        s[ii] = s[qp1];
        s[qp1] = rt;
        qjj = 3 * ii + 1;
        jj = 3 * (ii + 1) + 1;
        xswap(Vf, qjj, jj);
        xswap(U, qjj, jj);
        ii = qp1;
        qp1++;
      }

      qp1 = 0;
      m--;
      break;
    }
  }

  if (doscale) {
    b_xzlascl(cscale, anrm, s);
  }

  for (qjj = 0; qjj < 3; qjj++) {
    V[3 * qjj] = Vf[3 * qjj];
    jj = 3 * qjj + 1;
    V[jj] = Vf[jj];
    jj = 3 * qjj + 2;
    V[jj] = Vf[jj];
  }
}

/*
 * Arguments    : j_struct_T *MeritFunction
 *                const h_struct_T *WorkingSet
 *                c_struct_T *TrialState
 *                bool *Flags_fevalOK
 *                bool *Flags_done
 *                bool *Flags_stepAccepted
 *                bool *Flags_failedLineSearch
 *                int *Flags_stepType
 * Return Type  : bool
 */
static bool test_exit(j_struct_T *MeritFunction, const h_struct_T *WorkingSet,
                      c_struct_T *TrialState, bool *Flags_fevalOK, bool
                      *Flags_done, bool *Flags_stepAccepted, bool
                      *Flags_failedLineSearch, int *Flags_stepType)
{
  double smax;
  int idx_max;
  int k;
  int mLambda;
  int nVar;
  bool Flags_gradOK;
  bool isFeasible;
  *Flags_fevalOK = true;
  *Flags_done = false;
  *Flags_stepAccepted = false;
  *Flags_failedLineSearch = false;
  *Flags_stepType = 1;
  nVar = WorkingSet->nVar;
  mLambda = (((WorkingSet->sizes[0] + WorkingSet->sizes[2]) + WorkingSet->sizes
              [3]) + WorkingSet->sizes[4]) + 2;
  if (mLambda >= 0) {
    memcpy(&TrialState->lambdaStopTest.data[0], &TrialState->lambdasqp.data[0],
           (unsigned int)(mLambda + 1) * sizeof(double));
  }

  computeGradLag(TrialState->gradLag.data, WorkingSet->ldA, WorkingSet->nVar,
                 TrialState->grad.data, WorkingSet->sizes[2],
                 WorkingSet->Aineq.data, WorkingSet->Aeq.data,
                 WorkingSet->indexFixed.data, WorkingSet->sizes[0],
                 WorkingSet->indexLB.data, WorkingSet->sizes[3],
                 WorkingSet->indexUB.data, WorkingSet->sizes[4],
                 TrialState->lambdaStopTest.data);
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = fabs(TrialState->grad.data[0]);
      for (k = 2; k <= nVar; k++) {
        double s;
        s = fabs(TrialState->grad.data[k - 1]);
        if (s > smax) {
          idx_max = k;
          smax = s;
        }
      }
    }
  }

  smax = fmax(1.0, fabs(TrialState->grad.data[idx_max - 1]));
  if (rtIsInf(smax)) {
    smax = 1.0;
  }

  MeritFunction->nlpPrimalFeasError = computePrimalFeasError
    (TrialState->xstarsqp, WorkingSet->sizes[2], TrialState->cIneq.data,
     TrialState->cEq, WorkingSet->indexLB.data, WorkingSet->sizes[3],
     WorkingSet->indexUB.data, WorkingSet->sizes[4]);
  MeritFunction->feasRelativeFactor = fmax(1.0,
    MeritFunction->nlpPrimalFeasError);
  isFeasible = (MeritFunction->nlpPrimalFeasError <= 1.0E-6 *
                MeritFunction->feasRelativeFactor);
  Flags_gradOK = computeDualFeasError(WorkingSet->nVar, TrialState->gradLag.data,
    &MeritFunction->nlpDualFeasError);
  if (!Flags_gradOK) {
    *Flags_done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    MeritFunction->nlpComplError = 0.0;
    MeritFunction->firstOrderOpt = fmax(MeritFunction->nlpDualFeasError, 0.0);
    if (mLambda >= 0) {
      memcpy(&TrialState->lambdaStopTestPrev.data[0],
             &TrialState->lambdaStopTest.data[0], (unsigned int)(mLambda + 1) *
             sizeof(double));
    }

    if (isFeasible && (MeritFunction->nlpDualFeasError <= 1.0E-6 * smax)) {
      *Flags_done = true;
      TrialState->sqpExitFlag = 1;
    } else if (isFeasible && (TrialState->sqpFval < -1.0E+20)) {
      *Flags_done = true;
      TrialState->sqpExitFlag = -3;
    } else if (TrialState->FunctionEvaluations >= 900) {
      *Flags_done = true;
      TrialState->sqpExitFlag = 0;
    }
  }

  return Flags_gradOK;
}

/*
 * Arguments    : const double xk[9]
 *                h_struct_T *WorkingSet
 *                int mIneq
 *                const double cIneq_data[]
 *                const double cEq[3]
 *                int mLB
 *                int mUB
 *                int mFixed
 * Return Type  : void
 */
static void updateWorkingSetForNewQP(const double xk[9], h_struct_T *WorkingSet,
  int mIneq, const double cIneq_data[], const double cEq[3], int mLB, int mUB,
  int mFixed)
{
  double d;
  int b_i;
  int i;
  int iEq0;
  int idx;
  int iw0;
  iw0 = WorkingSet->ldA * mFixed;
  iEq0 = 0;
  for (idx = 0; idx < 3; idx++) {
    d = cEq[idx];
    WorkingSet->beq[idx] = -d;
    WorkingSet->bwset.data[mFixed + idx] = -d;
    i = WorkingSet->nVar;
    for (b_i = 0; b_i < i; b_i++) {
      WorkingSet->ATwset.data[iw0 + b_i] = WorkingSet->Aeq.data[iEq0 + b_i];
    }

    iw0 += WorkingSet->ldA;
    iEq0 += WorkingSet->ldA;
  }

  i = (unsigned char)mIneq;
  for (idx = 0; idx < i; idx++) {
    WorkingSet->bineq.data[idx] = -cIneq_data[idx];
  }

  i = (unsigned char)mLB;
  for (idx = 0; idx < i; idx++) {
    WorkingSet->lb.data[WorkingSet->indexLB.data[idx] - 1] = xk
      [WorkingSet->indexLB.data[idx] - 1] + 1.0;
  }

  i = (unsigned char)mUB;
  for (idx = 0; idx < i; idx++) {
    WorkingSet->ub.data[WorkingSet->indexUB.data[idx] - 1] = (double)
      iv[WorkingSet->indexUB.data[idx] - 1] - xk[WorkingSet->indexUB.data[idx] -
      1];
  }

  i = (unsigned char)mFixed;
  for (idx = 0; idx < i; idx++) {
    d = (double)iv[WorkingSet->indexFixed.data[idx] - 1] - xk
      [WorkingSet->indexFixed.data[idx] - 1];
    WorkingSet->ub.data[WorkingSet->indexFixed.data[idx] - 1] = d;
    WorkingSet->bwset.data[idx] = d;
  }

  if (WorkingSet->nActiveConstr > mFixed + 3) {
    iw0 = mFixed + 4;
    i = WorkingSet->nActiveConstr;
    for (idx = iw0; idx <= i; idx++) {
      switch (WorkingSet->Wid.data[idx - 1]) {
       case 4:
        WorkingSet->bwset.data[idx - 1] = WorkingSet->lb.data
          [WorkingSet->indexLB.data[WorkingSet->Wlocalidx.data[idx - 1] - 1] - 1];
        break;

       case 5:
        WorkingSet->bwset.data[idx - 1] = WorkingSet->ub.data
          [WorkingSet->indexUB.data[WorkingSet->Wlocalidx.data[idx - 1] - 1] - 1];
        break;

       default:
        WorkingSet->bwset.data[idx - 1] = WorkingSet->bineq.data
          [WorkingSet->Wlocalidx.data[idx - 1] - 1];
        break;
      }
    }
  }
}

/*
 * Arguments    : const double x[9]
 *                double y[3]
 * Return Type  : void
 */
static void vecnorm(const double x[9], double y[3])
{
  int j;
  int k;
  for (j = 0; j < 3; j++) {
    double b_y;
    double scale;
    int ix0;
    int kend;
    ix0 = j * 3;
    b_y = 0.0;
    scale = 3.3121686421112381E-170;
    kend = ix0 + 3;
    for (k = ix0 + 1; k <= kend; k++) {
      double absxk;
      absxk = fabs(x[k - 1]);
      if (absxk > scale) {
        double t;
        t = scale / absxk;
        b_y = b_y * t * t + 1.0;
        scale = absxk;
      } else {
        double t;
        t = absxk / scale;
        b_y += t * t;
      }
    }

    y[j] = scale * sqrt(b_y);
  }
}

/*
 * Arguments    : int n
 *                double a
 *                int ix0
 *                double y[9]
 *                int iy0
 * Return Type  : void
 */
static void xaxpy(int n, double a, int ix0, double y[9], int iy0)
{
  int k;
  if (!(a == 0.0)) {
    for (k = 0; k < n; k++) {
      int i;
      i = (iy0 + k) - 1;
      y[i] += a * y[(ix0 + k) - 1];
    }
  }
}

/*
 * Arguments    : int n
 *                const double x[9]
 *                int ix0
 *                const double y[9]
 *                int iy0
 * Return Type  : double
 */
static double xdotc(int n, const double x[9], int ix0, const double y[9], int
                    iy0)
{
  double d;
  int i;
  int k;
  d = 0.0;
  i = (unsigned char)n;
  for (k = 0; k < i; k++) {
    d += x[(ix0 + k) - 1] * y[(iy0 + k) - 1];
  }

  return d;
}

/*
 * Arguments    : int m
 *                int n
 *                int k
 *                const double A[81]
 *                int lda
 *                const double B_data[]
 *                int ib0
 *                int ldb
 *                double C_data[]
 *                int ldc
 * Return Type  : void
 */
static void xgemm(int m, int n, int k, const double A[81], int lda, const double
                  B_data[], int ib0, int ldb, double C_data[], int ldc)
{
  int cr;
  int ib;
  int ic;
  if ((m != 0) && (n != 0)) {
    int br;
    int i;
    int i1;
    int lastColC;
    br = ib0;
    lastColC = ldc * (n - 1);
    for (cr = 0; ldc < 0 ? cr >= lastColC : cr <= lastColC; cr += ldc) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        memset(&C_data[i + -1], 0, (unsigned int)((i1 - i) + 1) * sizeof(double));
      }
    }

    for (cr = 0; ldc < 0 ? cr >= lastColC : cr <= lastColC; cr += ldc) {
      int ar;
      ar = -1;
      i = br + k;
      for (ib = br; ib < i; ib++) {
        int i2;
        i1 = cr + 1;
        i2 = cr + m;
        for (ic = i1; ic <= i2; ic++) {
          C_data[ic - 1] += B_data[ib - 1] * A[(ar + ic) - cr];
        }

        ar += lda;
      }

      br += ldb;
    }
  }
}

/*
 * Arguments    : int m
 *                int n
 *                const double A_data[]
 *                int lda
 *                const double x_data[]
 *                double y_data[]
 * Return Type  : void
 */
static void xgemv(int m, int n, const double A_data[], int lda, const double
                  x_data[], double y_data[])
{
  int ia;
  int iac;
  int iy;
  if (n != 0) {
    int i;
    i = (unsigned char)n;
    for (iy = 0; iy < i; iy++) {
      y_data[iy] = -y_data[iy];
    }

    iy = 0;
    i = lda * (n - 1) + 1;
    for (iac = 1; lda < 0 ? iac >= i : iac <= i; iac += lda) {
      double c;
      int i1;
      c = 0.0;
      i1 = iac + m;
      for (ia = iac; ia < i1; ia++) {
        c += A_data[ia - 1] * x_data[ia - iac];
      }

      y_data[iy] += c;
      iy++;
    }
  }
}

/*
 * Arguments    : double A_data[]
 *                const int A_size[2]
 *                int m
 *                int n
 *                int jpvt_data[]
 *                double tau_data[]
 * Return Type  : int
 */
static int xgeqp3(double A_data[], const int A_size[2], int m, int n, int
                  jpvt_data[], double tau_data[])
{
  double vn1_data[40];
  double vn2_data[40];
  double work_data[40];
  double temp;
  int b_i;
  int iy;
  int k;
  int ma;
  int minmn_tmp;
  int pvt;
  int tau_size;
  ma = A_size[0];
  iy = A_size[0];
  tau_size = A_size[1];
  if (iy <= tau_size) {
    tau_size = iy;
  }

  if (m <= n) {
    minmn_tmp = m;
  } else {
    minmn_tmp = n;
  }

  if (tau_size - 1 >= 0) {
    memset(&tau_data[0], 0, (unsigned int)tau_size * sizeof(double));
  }

  if (minmn_tmp < 1) {
    for (pvt = 0; pvt < n; pvt++) {
      jpvt_data[pvt] = pvt + 1;
    }
  } else {
    int i;
    int ix;
    int nfxd;
    int temp_tmp;
    nfxd = 0;
    for (pvt = 0; pvt < n; pvt++) {
      if (jpvt_data[pvt] != 0) {
        nfxd++;
        if (pvt + 1 != nfxd) {
          ix = pvt * ma;
          iy = (nfxd - 1) * ma;
          for (k = 0; k < m; k++) {
            temp_tmp = ix + k;
            temp = A_data[temp_tmp];
            i = iy + k;
            A_data[temp_tmp] = A_data[i];
            A_data[i] = temp;
          }

          jpvt_data[pvt] = jpvt_data[nfxd - 1];
          jpvt_data[nfxd - 1] = pvt + 1;
        } else {
          jpvt_data[pvt] = pvt + 1;
        }
      } else {
        jpvt_data[pvt] = pvt + 1;
      }
    }

    if (nfxd > minmn_tmp) {
      nfxd = minmn_tmp;
    }

    qrf(A_data, A_size, m, n, nfxd, tau_data);
    if (nfxd < minmn_tmp) {
      double d;
      ma = A_size[0];
      iy = A_size[1];
      if (iy - 1 >= 0) {
        memset(&work_data[0], 0, (unsigned int)iy * sizeof(double));
        memset(&vn1_data[0], 0, (unsigned int)iy * sizeof(double));
        memset(&vn2_data[0], 0, (unsigned int)iy * sizeof(double));
      }

      i = nfxd + 1;
      for (pvt = i; pvt <= n; pvt++) {
        d = b_xnrm2(m - nfxd, A_data, (nfxd + (pvt - 1) * ma) + 1);
        vn1_data[pvt - 1] = d;
        vn2_data[pvt - 1] = d;
      }

      for (b_i = i; b_i <= minmn_tmp; b_i++) {
        double s;
        int ii;
        int ip1;
        int mmi;
        int nmi;
        ip1 = b_i + 1;
        nfxd = (b_i - 1) * ma;
        ii = (nfxd + b_i) - 1;
        nmi = (n - b_i) + 1;
        mmi = m - b_i;
        if (nmi < 1) {
          iy = -2;
        } else {
          iy = -1;
          if (nmi > 1) {
            temp = fabs(vn1_data[b_i - 1]);
            for (k = 2; k <= nmi; k++) {
              s = fabs(vn1_data[(b_i + k) - 2]);
              if (s > temp) {
                iy = k - 2;
                temp = s;
              }
            }
          }
        }

        pvt = b_i + iy;
        if (pvt + 1 != b_i) {
          ix = pvt * ma;
          for (k = 0; k < m; k++) {
            temp_tmp = ix + k;
            temp = A_data[temp_tmp];
            iy = nfxd + k;
            A_data[temp_tmp] = A_data[iy];
            A_data[iy] = temp;
          }

          iy = jpvt_data[pvt];
          jpvt_data[pvt] = jpvt_data[b_i - 1];
          jpvt_data[b_i - 1] = iy;
          vn1_data[pvt] = vn1_data[b_i - 1];
          vn2_data[pvt] = vn2_data[b_i - 1];
        }

        if (b_i < m) {
          temp = A_data[ii];
          d = xzlarfg(mmi + 1, &temp, A_data, ii + 2);
          tau_data[b_i - 1] = d;
          A_data[ii] = temp;
        } else {
          d = 0.0;
          tau_data[b_i - 1] = 0.0;
        }

        if (b_i < n) {
          temp = A_data[ii];
          A_data[ii] = 1.0;
          xzlarf(mmi + 1, nmi - 1, ii + 1, d, A_data, (ii + ma) + 1, ma,
                 work_data);
          A_data[ii] = temp;
        }

        for (pvt = ip1; pvt <= n; pvt++) {
          iy = b_i + (pvt - 1) * ma;
          d = vn1_data[pvt - 1];
          if (d != 0.0) {
            temp = fabs(A_data[iy - 1]) / d;
            temp = 1.0 - temp * temp;
            if (temp < 0.0) {
              temp = 0.0;
            }

            s = d / vn2_data[pvt - 1];
            s = temp * (s * s);
            if (s <= 1.4901161193847656E-8) {
              if (b_i < m) {
                d = b_xnrm2(mmi, A_data, iy + 1);
                vn1_data[pvt - 1] = d;
                vn2_data[pvt - 1] = d;
              } else {
                vn1_data[pvt - 1] = 0.0;
                vn2_data[pvt - 1] = 0.0;
              }
            } else {
              vn1_data[pvt - 1] = d * sqrt(temp);
            }
          }
        }
      }
    }
  }

  return tau_size;
}

/*
 * Arguments    : int m
 *                int n
 *                double alpha1
 *                int ix0
 *                const double y_data[]
 *                double A_data[]
 *                int ia0
 *                int lda
 * Return Type  : void
 */
static void xgerc(int m, int n, double alpha1, int ix0, const double y_data[],
                  double A_data[], int ia0, int lda)
{
  int ijA;
  int j;
  if (!(alpha1 == 0.0)) {
    int jA;
    jA = ia0;
    for (j = 0; j < n; j++) {
      double temp;
      temp = y_data[j];
      if (temp != 0.0) {
        int i;
        temp *= alpha1;
        i = m + jA;
        for (ijA = jA; ijA < i; ijA++) {
          A_data[ijA - 1] += A_data[((ix0 + ijA) - jA) - 1] * temp;
        }
      }

      jA += lda;
    }
  }
}

/*
 * Arguments    : int n
 *                const double x[9]
 *                int ix0
 * Return Type  : double
 */
static double xnrm2(int n, const double x[9], int ix0)
{
  double scale;
  double y;
  int k;
  int kend;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  kend = ix0 + n;
  for (k = ix0; k < kend; k++) {
    double absxk;
    absxk = fabs(x[k - 1]);
    if (absxk > scale) {
      double t;
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      double t;
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * sqrt(y);
}

/*
 * Arguments    : int n
 *                double A_data[]
 *                int lda
 * Return Type  : int
 */
static int xpotrf(int n, double A_data[], int lda)
{
  int ia;
  int iac;
  int info;
  int j;
  int k;
  bool exitg1;
  info = 0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j <= n - 1)) {
    double c;
    double ssq;
    int idxA1j;
    int idxAjj;
    idxA1j = j * lda;
    idxAjj = idxA1j + j;
    ssq = 0.0;
    if (j >= 1) {
      for (k = 0; k < j; k++) {
        c = A_data[idxA1j + k];
        ssq += c * c;
      }
    }

    ssq = A_data[idxAjj] - ssq;
    if (ssq > 0.0) {
      ssq = sqrt(ssq);
      A_data[idxAjj] = ssq;
      if (j + 1 < n) {
        int i;
        int ia0;
        int idxAjjp1;
        int nmj;
        nmj = (n - j) - 2;
        ia0 = (idxA1j + lda) + 1;
        idxAjjp1 = idxAjj + lda;
        if ((j != 0) && (nmj + 1 != 0)) {
          idxAjj = idxAjjp1;
          i = ia0 + lda * nmj;
          for (iac = ia0; lda < 0 ? iac >= i : iac <= i; iac += lda) {
            c = 0.0;
            k = iac + j;
            for (ia = iac; ia < k; ia++) {
              c += A_data[ia - 1] * A_data[(idxA1j + ia) - iac];
            }

            A_data[idxAjj] -= c;
            idxAjj += lda;
          }
        }

        ssq = 1.0 / ssq;
        i = (idxAjjp1 + lda * nmj) + 1;
        for (k = idxAjjp1 + 1; lda < 0 ? k >= i : k <= i; k += lda) {
          A_data[k - 1] *= ssq;
        }
      }

      j++;
    } else {
      A_data[idxAjj] = ssq;
      info = j + 1;
      exitg1 = true;
    }
  }

  return info;
}

/*
 * Arguments    : double x[9]
 *                int ix0
 *                int iy0
 *                double c
 *                double s
 * Return Type  : void
 */
static void xrot(double x[9], int ix0, int iy0, double c, double s)
{
  double temp;
  double temp_tmp;
  temp = x[iy0 - 1];
  temp_tmp = x[ix0 - 1];
  x[iy0 - 1] = c * temp - s * temp_tmp;
  x[ix0 - 1] = c * temp_tmp + s * temp;
  temp = c * x[ix0] + s * x[iy0];
  x[iy0] = c * x[iy0] - s * x[ix0];
  x[ix0] = temp;
  temp = x[iy0 + 1];
  temp_tmp = x[ix0 + 1];
  x[iy0 + 1] = c * temp - s * temp_tmp;
  x[ix0 + 1] = c * temp_tmp + s * temp;
}

/*
 * Arguments    : double *a
 *                double *b
 *                double *s
 * Return Type  : double
 */
static double xrotg(double *a, double *b, double *s)
{
  double absa;
  double absb;
  double c;
  double scale;
  c = *b;
  absa = fabs(*a);
  absb = fabs(*b);
  if (absa > absb) {
    c = *a;
  }

  scale = absa + absb;
  if (scale == 0.0) {
    *s = 0.0;
    c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    double ads;
    double bds;
    ads = absa / scale;
    bds = absb / scale;
    scale *= sqrt(ads * ads + bds * bds);
    if (c < 0.0) {
      scale = -scale;
    }

    c = *a / scale;
    *s = *b / scale;
    if (absa > absb) {
      *b = *s;
    } else if (c != 0.0) {
      *b = 1.0 / c;
    } else {
      *b = 1.0;
    }

    *a = scale;
  }

  return c;
}

/*
 * Arguments    : double x[9]
 *                int ix0
 *                int iy0
 * Return Type  : void
 */
static void xswap(double x[9], int ix0, int iy0)
{
  double temp;
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

/*
 * Arguments    : const double x[9]
 * Return Type  : double
 */
static double xzlangeM(const double x[9])
{
  double y;
  int k;
  bool exitg1;
  y = 0.0;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k < 9)) {
    double absxk;
    absxk = fabs(x[k]);
    if (rtIsNaN(absxk)) {
      y = rtNaN;
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

/*
 * Arguments    : int m
 *                int n
 *                int iv0
 *                double tau
 *                double C_data[]
 *                int ic0
 *                int ldc
 *                double work_data[]
 * Return Type  : void
 */
static void xzlarf(int m, int n, int iv0, double tau, double C_data[], int ic0,
                   int ldc, double work_data[])
{
  int i;
  int ia;
  int iac;
  int lastc;
  int lastv;
  if (tau != 0.0) {
    bool exitg2;
    lastv = m;
    i = iv0 + m;
    while ((lastv > 0) && (C_data[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n;
    exitg2 = false;
    while ((!exitg2) && (lastc > 0)) {
      int exitg1;
      i = ic0 + (lastc - 1) * ldc;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C_data[ia - 1] != 0.0) {
            exitg1 = 1;
          } else {
            ia++;
          }
        } else {
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    if (lastc != 0) {
      int b_i;
      if (lastc - 1 >= 0) {
        memset(&work_data[0], 0, (unsigned int)lastc * sizeof(double));
      }

      i = 0;
      b_i = ic0 + ldc * (lastc - 1);
      for (iac = ic0; ldc < 0 ? iac >= b_i : iac <= b_i; iac += ldc) {
        double c;
        int i1;
        c = 0.0;
        i1 = iac + lastv;
        for (ia = iac; ia < i1; ia++) {
          c += C_data[ia - 1] * C_data[((iv0 + ia) - iac) - 1];
        }

        work_data[i] += c;
        i++;
      }
    }

    xgerc(lastv, lastc, -tau, iv0, work_data, C_data, ic0, ldc);
  }
}

/*
 * Arguments    : int n
 *                double *alpha1
 *                double x_data[]
 *                int ix0
 * Return Type  : double
 */
static double xzlarfg(int n, double *alpha1, double x_data[], int ix0)
{
  double tau;
  int k;
  tau = 0.0;
  if (n > 0) {
    double xnorm;
    xnorm = b_xnrm2(n - 1, x_data, ix0);
    if (xnorm != 0.0) {
      double a_tmp;
      a_tmp = fabs(*alpha1);
      xnorm = fabs(xnorm);
      if (a_tmp < xnorm) {
        a_tmp /= xnorm;
        xnorm *= sqrt(a_tmp * a_tmp + 1.0);
      } else if (a_tmp > xnorm) {
        xnorm /= a_tmp;
        xnorm = a_tmp * sqrt(xnorm * xnorm + 1.0);
      } else if (rtIsNaN(xnorm)) {
        xnorm = rtNaN;
      } else {
        xnorm = a_tmp * 1.4142135623730951;
      }

      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        int i;
        int knt;
        knt = 0;
        i = (ix0 + n) - 2;
        do {
          knt++;
          for (k = ix0; k <= i; k++) {
            x_data[k - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while ((fabs(xnorm) < 1.0020841800044864E-292) && (knt < 20));

        a_tmp = fabs(*alpha1);
        xnorm = fabs(b_xnrm2(n - 1, x_data, ix0));
        if (a_tmp < xnorm) {
          a_tmp /= xnorm;
          xnorm *= sqrt(a_tmp * a_tmp + 1.0);
        } else if (a_tmp > xnorm) {
          xnorm /= a_tmp;
          xnorm = a_tmp * sqrt(xnorm * xnorm + 1.0);
        } else if (rtIsNaN(xnorm)) {
          xnorm = rtNaN;
        } else {
          xnorm = a_tmp * 1.4142135623730951;
        }

        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        a_tmp = 1.0 / (*alpha1 - xnorm);
        for (k = ix0; k <= i; k++) {
          x_data[k - 1] *= a_tmp;
        }

        for (k = 0; k < knt; k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        int i;
        tau = (xnorm - *alpha1) / xnorm;
        a_tmp = 1.0 / (*alpha1 - xnorm);
        i = (ix0 + n) - 2;
        for (k = ix0; k <= i; k++) {
          x_data[k - 1] *= a_tmp;
        }

        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

/*
 * Arguments    : double cfrom
 *                double cto
 *                double A[9]
 * Return Type  : void
 */
static void xzlascl(double cfrom, double cto, double A[9])
{
  double cfromc;
  double ctoc;
  int j;
  bool notdone;
  cfromc = cfrom;
  ctoc = cto;
  notdone = true;
  while (notdone) {
    double cfrom1;
    double cto1;
    double mul;
    cfrom1 = cfromc * 2.0041683600089728E-292;
    cto1 = ctoc / 4.9896007738368E+291;
    if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
      mul = 2.0041683600089728E-292;
      cfromc = cfrom1;
    } else if (fabs(cto1) > fabs(cfromc)) {
      mul = 4.9896007738368E+291;
      ctoc = cto1;
    } else {
      mul = ctoc / cfromc;
      notdone = false;
    }

    for (j = 0; j < 3; j++) {
      int offset;
      offset = j * 3 - 1;
      A[offset + 1] *= mul;
      A[offset + 2] *= mul;
      A[offset + 3] *= mul;
    }
  }
}

/*
 * function [result, radius, exitflag] = optimize_forces(center, initial_guess, T_min, T_max, psi)
 *
 * Parameters
 *  center: 球心
 *  initial_guess: W 初始猜测值
 *  initial_guess = -initial_guess; % 绳向与力方向相反
 *
 * Arguments    : const double center[3]
 *                const double initial_guess[9]
 *                const double T_min[3]
 *                const double T_max[3]
 *                const double psi[3]
 *                double result[9]
 *                double *radius
 *                double *exitflag
 * Return Type  : void
 */
void optimize_forces(const double center[3], const double initial_guess[9],
                     const double T_min[3], const double T_max[3], const double
                     psi[3], double result[9], double *radius, double *exitflag)
{
  double A_data[54];
  double f_objective_workspace_A_alpha[18];
  double qxy_ub[9];
  double b_data[6];
  double f_objective_workspace_b_alpha[6];
  double f_objective_workspace_deltaT[3];
  double b_expl_temp;
  double c_expl_temp;
  double d_expl_temp;
  double e_expl_temp;
  double f_expl_temp;
  double g_expl_temp;
  double maxfval;
  int b_size[2];
  int c_f_objective_workspace_A_alpha;
  int i;
  signed char b_i;
  bool bv[3];

  /* 'optimize_forces:6' num = length(initial_guess) / 3; */
  /*  Preprocess */
  /* 'optimize_forces:8' A_alpha = zeros([num * 2, num]); */
  memset(&f_objective_workspace_A_alpha[0], 0, 18U * sizeof(double));

  /* 'optimize_forces:9' b_alpha = zeros([num * 2, 1]); */
  for (i = 0; i < 6; i++) {
    f_objective_workspace_b_alpha[i] = 0.0;
  }

  /* 'optimize_forces:10' for i = 1:num */
  /* 'optimize_forces:14' deltaT = T_max - T_min; */
  /*  不等式约束 */
  /* 'optimize_forces:39' if all(isfinite(psi)) */
  for (i = 0; i < 3; i++) {
    /* 'optimize_forces:11' A_alpha(2*i-1:2*i, i) = [1; -1]; */
    b_i = (signed char)((i + 1) << 1);

    /* 'optimize_forces:12' b_alpha(2*i-1:2*i) = [1; 0]; */
    c_f_objective_workspace_A_alpha = b_i + 6 * i;
    f_objective_workspace_A_alpha[c_f_objective_workspace_A_alpha - 2] = 1.0;
    f_objective_workspace_b_alpha[b_i - 2] = 1.0;
    f_objective_workspace_A_alpha[c_f_objective_workspace_A_alpha - 1] = -1.0;
    f_objective_workspace_b_alpha[b_i - 1] = 0.0;
    f_objective_workspace_deltaT[i] = T_max[i] - T_min[i];
    maxfval = psi[i];
    bv[i] = ((!rtIsInf(maxfval)) && (!rtIsNaN(maxfval)));
  }

  if (all(bv)) {
    double qxy_lb[9];

    /* 'optimize_forces:40' range = deg2rad(60); */
    /* 'optimize_forces:41' qxy_ub = [sin(psi + range) -cos(psi + range) zeros(num, 1)]; */
    /* 'optimize_forces:42' qxy_lb = [sin(psi - range) -cos(psi - range) zeros(num, 1)]; */
    qxy_ub[0] = sin(psi[0] + 1.0471975511965976);
    qxy_ub[3] = -cos(psi[0] + 1.0471975511965976);
    qxy_ub[6] = 0.0;
    qxy_lb[0] = sin(psi[0] - 1.0471975511965976);
    qxy_lb[3] = -cos(psi[0] - 1.0471975511965976);
    qxy_lb[6] = 0.0;
    qxy_ub[1] = sin(psi[1] + 1.0471975511965976);
    qxy_ub[4] = -cos(psi[1] + 1.0471975511965976);
    qxy_ub[7] = 0.0;
    qxy_lb[1] = sin(psi[1] - 1.0471975511965976);
    qxy_lb[4] = -cos(psi[1] - 1.0471975511965976);
    qxy_lb[7] = 0.0;
    qxy_ub[2] = sin(psi[2] + 1.0471975511965976);
    qxy_ub[5] = -cos(psi[2] + 1.0471975511965976);
    qxy_ub[8] = 0.0;
    qxy_lb[2] = sin(psi[2] - 1.0471975511965976);
    qxy_lb[5] = -cos(psi[2] - 1.0471975511965976);
    qxy_lb[8] = 0.0;

    /* 'optimize_forces:43' A = zeros([num * 2, num * 3]); */
    memset(&A_data[0], 0, 54U * sizeof(double));

    /* 'optimize_forces:44' for i = 1:num */
    for (i = 0; i < 3; i++) {
      /* 'optimize_forces:45' A(i, 3*i-2:3*i) = -qxy_ub(i, :); */
      b_i = (signed char)(3 * (i + 1));

      /* 'optimize_forces:46' A(i + num, 3*i-2:3*i) = qxy_lb(i, :); */
      c_f_objective_workspace_A_alpha = i + 6 * (b_i - 3);
      A_data[c_f_objective_workspace_A_alpha] = -qxy_ub[i];
      A_data[c_f_objective_workspace_A_alpha + 3] = qxy_lb[i];
      c_f_objective_workspace_A_alpha = i + 6 * (b_i - 2);
      A_data[c_f_objective_workspace_A_alpha] = -qxy_ub[i + 3];
      A_data[c_f_objective_workspace_A_alpha + 3] = qxy_lb[i + 3];
      c_f_objective_workspace_A_alpha = i + 6 * (b_i - 1);
      A_data[c_f_objective_workspace_A_alpha] = -0.0;
      A_data[c_f_objective_workspace_A_alpha + 3] = 0.0;
    }

    /* 'optimize_forces:48' b = zeros(num * 2, 1); */
    b_size[0] = 6;
    b_size[1] = 1;
    for (c_f_objective_workspace_A_alpha = 0; c_f_objective_workspace_A_alpha <
         6; c_f_objective_workspace_A_alpha++) {
      b_data[c_f_objective_workspace_A_alpha] = 0.0;
    }
  } else {
    /* 'optimize_forces:49' else */
    /* 'optimize_forces:50' A = []; */
    /* 'optimize_forces:51' b = []; */
    b_size[0] = 0;
    b_size[1] = 0;
  }

  /*  等式约束     */
  /* 'optimize_forces:56' Aeq = []; */
  /* 'optimize_forces:57' beq = []; */
  /*  上下限 */
  /* 'optimize_forces:60' ub = repmat([1 1 0]', num, 1); */
  /* 'optimize_forces:61' lb = -ones(size(initial_guess)); */
  /*  repmat([-1 -1 0]', num, 1); */
  /*  目标函数 */
  /* 'optimize_forces:64' f_objective = @(W) objective(W, center, A_alpha, b_alpha, deltaT, T_min); */
  /*  优化选项 */
  /* 'optimize_forces:67' options = optimoptions('fmincon', 'Algorithm', 'sqp'); */
  /*  优化求解 */
  /* 'optimize_forces:70' [result, maxfval, exitflag, ~] = fmincon(f_objective, -initial_guess, A, b, Aeq, beq, ... */
  /* 'optimize_forces:71'                     lb, ub, @(W) constraint(W), options); */
  for (c_f_objective_workspace_A_alpha = 0; c_f_objective_workspace_A_alpha < 9;
       c_f_objective_workspace_A_alpha++) {
    qxy_ub[c_f_objective_workspace_A_alpha] =
      -initial_guess[c_f_objective_workspace_A_alpha];
  }

  char expl_temp[3];
  maxfval = fmincon(center, f_objective_workspace_A_alpha,
                    f_objective_workspace_b_alpha, f_objective_workspace_deltaT,
                    T_min, qxy_ub, A_data, b_data, b_size, result, exitflag,
                    &b_expl_temp, &c_expl_temp, expl_temp, &d_expl_temp,
                    &e_expl_temp, &f_expl_temp, &g_expl_temp);

  /*  输出结果 */
  /* 'optimize_forces:75' result = -result; */
  for (c_f_objective_workspace_A_alpha = 0; c_f_objective_workspace_A_alpha < 9;
       c_f_objective_workspace_A_alpha++) {
    result[c_f_objective_workspace_A_alpha] =
      -result[c_f_objective_workspace_A_alpha];
  }

  /*  绳向与力方向相反 */
  /* 'optimize_forces:76' radius = -maxfval; */
  *radius = -maxfval;

  /* 'optimize_forces:77' disp('exitflag:'); */
  /* 'optimize_forces:78' disp(exitflag); */
  /* 'optimize_forces:79' disp('Result:'); */
  /* 'optimize_forces:80' res = reshape(-result, length(result) / num, num); */
  /* 'optimize_forces:81' res = res ./ vecnorm(res, 2, 1); */
  /* 'optimize_forces:82' angle = zeros(2, 3); */
  /* 'optimize_forces:83' for i = 1:3 */
  /* 'optimize_forces:86' disp(rad2deg(angle)); */
  /* 'optimize_forces:87' disp('Margin:') */
  /* 'optimize_forces:88' disp(radius); */
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void optimize_forces_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void optimize_forces_terminate(void)
{
}

/*
 * File trailer for optimize_forces.c
 *
 * [EOF]
 */
