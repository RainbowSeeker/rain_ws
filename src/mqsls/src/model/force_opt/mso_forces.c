/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: mso_forces.c
 *
 * MATLAB Coder version            : 24.2
 * C/C++ source code generated on  : 2025-01-02 16:21:33
 */

/* Include Files */
#include "mso_forces.h"
#include "rt_nonfinite.h"
#include "rt_nonfinite.h"
#include <math.h>
#include <string.h>

/* Type Definitions */
#ifndef typedef_struct_T
#define typedef_struct_T

typedef struct {
  double workspace_float[880];
  int workspace_int[40];
  int workspace_sort[40];
} struct_T;

#endif                                 /* typedef_struct_T */

#ifndef typedef_b_struct_T
#define typedef_b_struct_T

typedef struct {
  double center[3];
  double T_max[3];
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

#ifndef typedef_c_struct_T
#define typedef_c_struct_T

typedef struct {
  double sqpFval;
  double sqpFval_old;
  double xstarsqp[9];
  double xstarsqp_old[9];
  double cIneq[6];
  double cIneq_old[6];
  double cEq[3];
  double cEq_old[3];
  double grad[22];
  double grad_old[22];
  int FunctionEvaluations;
  int sqpIterations;
  int sqpExitFlag;
  double lambdasqp[40];
  double lambdaStopTest[40];
  double lambdaStopTestPrev[40];
  double steplength;
  double delta_x[22];
  double socDirection[22];
  int workingset_old[40];
  double JacCeqTrans_old[66];
  double gradLag[22];
  double delta_gradLag[22];
  double xstar[22];
  double fstar;
  double lambda[40];
  int state;
  double maxConstr;
  int iterations;
  double searchDir[22];
} c_struct_T;

#endif                                 /* typedef_c_struct_T */

#ifndef typedef_d_struct_T
#define typedef_d_struct_T

typedef struct {
  anonymous_function objfun;
  double f_1;
  int numEvals;
  bool hasBounds;
} d_struct_T;

#endif                                 /* typedef_d_struct_T */

#ifndef typedef_e_struct_T
#define typedef_e_struct_T

typedef struct {
  int mConstr;
  int nVar;
  double Aineq[132];
  double bineq[6];
  double Aeq[66];
  double beq[3];
  double lb[22];
  double ub[22];
  int indexLB[22];
  int indexUB[22];
  int mEqRemoved;
  int indexEqRemoved[3];
  double ATwset[880];
  double bwset[40];
  int nActiveConstr;
  double maxConstrWorkspace[40];
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
  bool isActiveConstr[40];
  int Wid[40];
  int Wlocalidx[40];
  int nWConstr[5];
  int probType;
} e_struct_T;

#endif                                 /* typedef_e_struct_T */

#ifndef typedef_f_struct_T
#define typedef_f_struct_T

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
} f_struct_T;

#endif                                 /* typedef_f_struct_T */

#ifndef typedef_g_struct_T
#define typedef_g_struct_T

typedef struct {
  double QR[1600];
  double Q[1600];
  int jpvt[40];
  int mrows;
  int ncols;
  double tau[40];
  int minRowCol;
  bool usedPivoting;
} g_struct_T;

#endif                                 /* typedef_g_struct_T */

#ifndef typedef_h_struct_T
#define typedef_h_struct_T

typedef struct {
  double FMat[1600];
  int ndims;
  int info;
  bool ConvexCheck;
  double workspace_;
} h_struct_T;

#endif                                 /* typedef_h_struct_T */

#ifndef typedef_i_struct_T
#define typedef_i_struct_T

typedef struct {
  double grad[22];
  double Hx[21];
  bool hasLinear;
  int nvar;
  double beta;
  double rho;
  int objtype;
  int prev_objtype;
  int prev_nvar;
  bool prev_hasLinear;
  double gammaScalar;
} i_struct_T;

#endif                                 /* typedef_i_struct_T */

#ifndef typedef_j_struct_T
#define typedef_j_struct_T

typedef struct {
  bool fevalOK;
  bool done;
  bool stepAccepted;
  bool failedLineSearch;
  int stepType;
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

static const signed char iv1[18] = { 1, 0, 0, -1, 0, 0, 0, 1, 0, 0, -1, 0, 0, 0,
  1, 0, 0, -1 };

/* Function Declarations */
static bool BFGSUpdate(int nvar, double Bk[81], const double sk[22], double yk
  [22], double workspace[880]);
static void PresolveWorkingSet(c_struct_T *solution, struct_T *memspace,
  e_struct_T *workingset, g_struct_T *qrmanager);
static void RemoveDependentIneq_(e_struct_T *workingset, g_struct_T *qrmanager,
  struct_T *memspace, double tolfactor);
static void addAeqConstr(e_struct_T *obj, int idx_local);
static void addBoundToActiveSetMatrix_(e_struct_T *obj, int TYPE, int idx_local);
static void b_computeGradLag(double workspace[880], int nVar, const double grad
  [22], const double AineqTrans[132], const double AeqTrans[66], const int
  finiteLB[22], int mLB, const int finiteUB[22], const double lambda[40]);
static void b_driver(c_struct_T *TrialState, f_struct_T *MeritFunction, const
                     i_coder_internal_stickyStruct *FcnEvaluator, d_struct_T
                     *FiniteDifferences, struct_T *memspace, e_struct_T
                     *WorkingSet, double Hessian[81], g_struct_T *QRManager,
                     h_struct_T *CholManager, i_struct_T *QPObjective);
static double b_maxConstraintViolation(e_struct_T *obj, const double x[22]);
static void b_test_exit(j_struct_T *Flags, struct_T *memspace, f_struct_T
  *MeritFunction, e_struct_T *WorkingSet, c_struct_T *TrialState, g_struct_T
  *QRManager);
static void b_vecnorm(const double x[18], double y[6]);
static void b_xaxpy(double a, const double x[9], int ix0, double y[3]);
static void b_xgemm(int m, int n, int k, const double A[1600], int ia0, const
                    double B[880], double C[1600]);
static double b_xnrm2(int n, const double x[1600], int ix0);
static void b_xzlascl(double cfrom, double cto, double A[3]);
static double c_maxConstraintViolation_AMats_(e_struct_T *obj, const double x[22]);
static void c_xaxpy(double a, const double x[3], double y[9], int iy0);
static double c_xnrm2(int n, const double x[22]);
static int checkVectorNonFinite(const double vec[3]);
static double computeComplError(const double xCurrent[9], const double cIneq[6],
  const int finiteLB[22], int mLB, const int finiteUB[22], const double lambda
  [40]);
static bool computeFiniteDifferences(d_struct_T *obj, double fCurrent, double
  xk[9], double gradf[22]);
static double computeFval_ReuseHx(const i_struct_T *obj, double workspace[880],
  const double f[22], const double x[22]);
static void computeGradLag(double workspace[22], int nVar, const double grad[22],
  const double AineqTrans[132], const double AeqTrans[66], const int finiteLB[22],
  int mLB, const int finiteUB[22], const double lambda[40]);
static void computeGrad_StoreHx(i_struct_T *obj, const double H[81], const
  double f[22], const double x[22]);
static void computeLinearResiduals(const double x[9], int nVar, double
  workspaceIneq[6], const double AineqT[132]);
static void computeQ_(g_struct_T *obj, int nrows);
static void compute_deltax(const double H[81], c_struct_T *solution, struct_T
  *memspace, const g_struct_T *qrmanager, h_struct_T *cholmanager, const
  i_struct_T *objective, bool alwaysPositiveDef);
static void compute_lambda(double workspace[880], c_struct_T *solution, const
  i_struct_T *objective, const g_struct_T *qrmanager);
static void countsort(int x[40], int xLen, int workspace[40], int xMin, int xMax);
static double d_maxConstraintViolation_AMats_(e_struct_T *obj, const double x[22]);
static double d_xnrm2(const double x[3]);
static void deleteColMoveEnd(g_struct_T *obj, int idx);
static int div_nde_s32_floor(int numerator, int denominator);
static void driver(const double H[81], const double f[22], c_struct_T *solution,
                   struct_T *memspace, e_struct_T *workingset, g_struct_T
                   *qrmanager, h_struct_T *cholmanager, i_struct_T *objective,
                   const k_struct_T *options, const k_struct_T *runTimeOptions);
static double evalObjAndConstr(const b_struct_T *c_obj_next_next_next_next_next_,
  const double x[9], double Ceq_workspace[3], int *status);
static double evalObjAndConstrAndDerivatives(const b_struct_T
  *c_obj_next_next_next_next_next_, const double x[9], double Ceq_workspace[3],
  double JacEqTrans_workspace[66], int *status);
static void factorQR(g_struct_T *obj, const double A[880], int mrows, int ncols);
static void factorQRE(g_struct_T *obj, const double A[880], int mrows, int ncols);
static bool feasibleX0ForWorkingSet(double workspace[880], double xCurrent[22],
  e_struct_T *workingset, g_struct_T *qrmanager);
static double feasibleratiotest(const double solution_xstar[22], const double
  solution_searchDir[22], double workspace[880], int workingset_nVar, const
  double workingset_Aineq[132], const double workingset_bineq[6], const double
  workingset_lb[22], const double workingset_ub[22], const int
  workingset_indexLB[22], const int workingset_indexUB[22], const int
  workingset_sizes[5], const int workingset_isActiveIdx[6], const bool
  workingset_isActiveConstr[40], const int workingset_nWConstr[5], bool
  isPhaseOne, bool *newBlocking, int *constrType, int *constrIdx);
static double fmincon(const double fun_workspace_center[3], const double
                      fun_workspace_T_max[3], const double fun_workspace_T_min[3],
                      const double x0[9], const double Aineq[54], double x[9],
                      double *exitflag, double *output_iterations, double
                      *output_funcCount, char output_algorithm[3], double
                      *output_constrviolation, double *output_stepsize, double
                      *output_lssteplength, double *output_firstorderopt);
static void fullColLDL2_(h_struct_T *obj, int NColsRemain);
static void linearForm_(bool obj_hasLinear, int obj_nvar, double workspace[880],
  const double H[81], const double f[22], const double x[22]);
static double maxConstraintViolation(e_struct_T *obj, const double x[880]);
static double minimum(const double x[6]);
static void modifyOverheadPhaseOne_(e_struct_T *obj);
static void mso_forces_anonFcn2(const double W[9], double varargout_2[3], double
  varargout_4[27]);
static void phaseone(const double H[81], const double f[22], c_struct_T
                     *solution, struct_T *memspace, e_struct_T *workingset,
                     g_struct_T *qrmanager, h_struct_T *cholmanager, i_struct_T *
                     objective, const char options_SolverName[7], const
                     k_struct_T *runTimeOptions);
static void pinv(const double A[9], double X[9]);
static void qrf(double A[1600], int m, int n, int nfxd, double tau[40]);
static void relaxed(const double Hessian[81], const double grad[22], c_struct_T *
                    TrialState, f_struct_T *MeritFunction, struct_T *memspace,
                    e_struct_T *WorkingSet, g_struct_T *QRManager, h_struct_T
                    *CholManager, i_struct_T *QPObjective, k_struct_T *qpoptions);
static void removeConstr(e_struct_T *obj, int idx_global);
static void setProblemType(e_struct_T *obj, int PROBLEM_TYPE);
static bool soc(const double Hessian[81], const double grad[22], c_struct_T
                *TrialState, struct_T *memspace, e_struct_T *WorkingSet,
                g_struct_T *QRManager, h_struct_T *CholManager, i_struct_T
                *QPObjective, const k_struct_T *qpoptions);
static void solve(const h_struct_T *obj, double rhs[22]);
static void sortLambdaQP(double lambda[40], int WorkingSet_nActiveConstr, const
  int WorkingSet_sizes[5], const int WorkingSet_isActiveIdx[6], const int
  WorkingSet_Wid[40], const int WorkingSet_Wlocalidx[40], double workspace[880]);
static void squareQ_appendCol(g_struct_T *obj, const double vec[880], int iv0);
static void step(j_struct_T *stepFlags, double Hessian[81], c_struct_T
                 *TrialState, f_struct_T *MeritFunction, struct_T *memspace,
                 e_struct_T *WorkingSet, g_struct_T *QRManager, h_struct_T
                 *CholManager, i_struct_T *QPObjective, k_struct_T *qpoptions);
static void sum(const double x[9], double y[3]);
static void svd(const double A[9], double U[9], double s[3], double V[9]);
static bool test_exit(f_struct_T *MeritFunction, const e_struct_T *WorkingSet,
                      c_struct_T *TrialState, bool *Flags_fevalOK, bool
                      *Flags_done, bool *Flags_stepAccepted, bool
                      *Flags_failedLineSearch, int *Flags_stepType);
static void updateWorkingSetForNewQP(const double xk[9], e_struct_T *WorkingSet,
  const double cIneq[6], const double cEq[3], int mLB);
static void vecnorm(const double x[9], double y[3]);
static void xaxpy(int n, double a, int ix0, double y[9], int iy0);
static double xdotc(int n, const double x[9], int ix0, const double y[9], int
                    iy0);
static void xgemm(int m, int n, int k, const double A[81], int lda, const double
                  B[1600], int ib0, double C[880]);
static void xgemv(int m, const double A[66], const double x[880], double y[40]);
static void xgeqp3(double A[1600], int m, int n, int jpvt[40], double tau[40]);
static double xnrm2(int n, const double x[9], int ix0);
static int xpotrf(int n, double A[1600]);
static void xrot(double x[9], int ix0, int iy0, double c, double s);
static double xrotg(double *a, double *b, double *s);
static void xswap(double x[9], int ix0, int iy0);
static double xzlangeM(const double x[9]);
static void xzlarf(int m, int n, int iv0, double tau, double C[1600], int ic0,
                   double work[40]);
static double xzlarfg(int n, double *alpha1, double x[1600], int ix0);
static void xzlascl(double cfrom, double cto, double A[9]);

/* Function Definitions */
/*
 * Arguments    : int nvar
 *                double Bk[81]
 *                const double sk[22]
 *                double yk[22]
 *                double workspace[880]
 * Return Type  : bool
 */
static bool BFGSUpdate(int nvar, double Bk[81], const double sk[22], double yk
  [22], double workspace[880])
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
  i = (unsigned char)nvar;
  for (ix = 0; ix < i; ix++) {
    dotSY += sk[ix] * yk[ix];
    workspace[ix] = 0.0;
  }

  ix = 0;
  i1 = 9 * (nvar - 1) + 1;
  for (iac = 1; iac <= i1; iac += 9) {
    j = iac + nvar;
    for (ia = iac; ia < j; ia++) {
      ijA = ia - iac;
      workspace[ijA] += Bk[ia - 1] * sk[ix];
    }

    ix++;
  }

  curvatureS = 0.0;
  if (nvar >= 1) {
    for (ix = 0; ix < nvar; ix++) {
      curvatureS += sk[ix] * workspace[ix];
    }
  }

  if (dotSY < 0.2 * curvatureS) {
    theta = 0.8 * curvatureS / (curvatureS - dotSY);
    for (ix = 0; ix < i; ix++) {
      yk[ix] *= theta;
    }

    if (!(1.0 - theta == 0.0)) {
      for (ix = 0; ix < nvar; ix++) {
        yk[ix] += (1.0 - theta) * workspace[ix];
      }
    }

    dotSY = 0.0;
    for (ix = 0; ix < i; ix++) {
      dotSY += sk[ix] * yk[ix];
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
        theta = workspace[j];
        if (theta != 0.0) {
          theta *= curvatureS;
          i1 = nvar + ix;
          for (ijA = ix; ijA < i1; ijA++) {
            Bk[ijA - 1] += workspace[ijA - ix] * theta;
          }
        }

        ix += 9;
      }
    }

    curvatureS = 1.0 / dotSY;
    if (!(curvatureS == 0.0)) {
      ix = 1;
      for (j = 0; j < i; j++) {
        theta = yk[j];
        if (theta != 0.0) {
          theta *= curvatureS;
          i1 = nvar + ix;
          for (ijA = ix; ijA < i1; ijA++) {
            Bk[ijA - 1] += yk[ijA - ix] * theta;
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
 *                struct_T *memspace
 *                e_struct_T *workingset
 *                g_struct_T *qrmanager
 * Return Type  : void
 */
static void PresolveWorkingSet(c_struct_T *solution, struct_T *memspace,
  e_struct_T *workingset, g_struct_T *qrmanager)
{
  double tol;
  int idxDiag;
  int idx_col;
  int ix;
  int ix0;
  int k;
  int mTotalWorkingEq_tmp;
  int mWorkingFixed;
  int nDepInd;
  solution->state = 82;
  mWorkingFixed = workingset->nWConstr[0];
  mTotalWorkingEq_tmp = workingset->nWConstr[0] + workingset->nWConstr[1];
  nDepInd = 0;
  if (mTotalWorkingEq_tmp > 0) {
    int i;
    int i1;
    int u0;
    i = (unsigned char)workingset->nVar;
    for (ix = 0; ix < mTotalWorkingEq_tmp; ix++) {
      for (idx_col = 0; idx_col < i; idx_col++) {
        qrmanager->QR[ix + 40 * idx_col] = workingset->ATwset[idx_col + 22 * ix];
      }
    }

    nDepInd = mTotalWorkingEq_tmp - workingset->nVar;
    if (nDepInd <= 0) {
      nDepInd = 0;
    }

    memset(&qrmanager->jpvt[0], 0, (unsigned int)i * sizeof(int));
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
      xgeqp3(qrmanager->QR, mTotalWorkingEq_tmp, workingset->nVar,
             qrmanager->jpvt, qrmanager->tau);
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

    idxDiag = u0 + 40 * (u0 - 1);
    while ((idxDiag > 0) && (fabs(qrmanager->QR[idxDiag - 1]) < tol * fabs
            (qrmanager->QR[0]))) {
      idxDiag -= 41;
      nDepInd++;
    }

    if (nDepInd > 0) {
      bool exitg1;
      computeQ_(qrmanager, qrmanager->mrows);
      idx_col = 0;
      exitg1 = false;
      while ((!exitg1) && (idx_col <= nDepInd - 1)) {
        double qtb;
        ix = 40 * ((mTotalWorkingEq_tmp - idx_col) - 1);
        qtb = 0.0;
        for (k = 0; k < mTotalWorkingEq_tmp; k++) {
          qtb += qrmanager->Q[ix + k] * workingset->bwset[k];
        }

        if (fabs(qtb) >= tol) {
          nDepInd = -1;
          exitg1 = true;
        } else {
          idx_col++;
        }
      }
    }

    if (nDepInd > 0) {
      for (idx_col = 0; idx_col < mTotalWorkingEq_tmp; idx_col++) {
        idxDiag = 40 * idx_col;
        ix0 = 22 * idx_col;
        for (k = 0; k < i; k++) {
          qrmanager->QR[idxDiag + k] = workingset->ATwset[ix0 + k];
        }
      }

      for (idx_col = 0; idx_col < mWorkingFixed; idx_col++) {
        qrmanager->jpvt[idx_col] = 1;
      }

      ix0 = workingset->nWConstr[0] + 1;
      if (ix0 <= mTotalWorkingEq_tmp) {
        memset(&qrmanager->jpvt[ix0 + -1], 0, (unsigned int)
               ((mTotalWorkingEq_tmp - ix0) + 1) * sizeof(int));
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
        xgeqp3(qrmanager->QR, workingset->nVar, mTotalWorkingEq_tmp,
               qrmanager->jpvt, qrmanager->tau);
      }

      for (idx_col = 0; idx_col < nDepInd; idx_col++) {
        memspace->workspace_int[idx_col] = qrmanager->jpvt[(mTotalWorkingEq_tmp
          - nDepInd) + idx_col];
      }

      countsort(memspace->workspace_int, nDepInd, memspace->workspace_sort, 1,
                mTotalWorkingEq_tmp);
      for (idx_col = nDepInd; idx_col >= 1; idx_col--) {
        i1 = workingset->nWConstr[0] + workingset->nWConstr[1];
        if (i1 != 0) {
          ix0 = memspace->workspace_int[idx_col - 1];
          if (ix0 <= i1) {
            if ((workingset->nActiveConstr == i1) || (ix0 == i1)) {
              workingset->mEqRemoved++;
              workingset->indexEqRemoved[workingset->mEqRemoved - 1] =
                workingset->Wlocalidx[ix0 - 1];
              removeConstr(workingset, ix0);
            } else {
              workingset->mEqRemoved++;
              idxDiag = workingset->Wid[ix0 - 1] - 1;
              ix = workingset->Wlocalidx[ix0 - 1];
              workingset->indexEqRemoved[workingset->mEqRemoved - 1] = ix;
              workingset->isActiveConstr[(workingset->isActiveIdx[idxDiag] + ix)
                - 2] = false;
              workingset->Wid[ix0 - 1] = workingset->Wid[i1 - 1];
              workingset->Wlocalidx[ix0 - 1] = workingset->Wlocalidx[i1 - 1];
              for (ix = 0; ix < i; ix++) {
                workingset->ATwset[ix + 22 * (ix0 - 1)] = workingset->ATwset[ix
                  + 22 * (i1 - 1)];
              }

              workingset->bwset[ix0 - 1] = workingset->bwset[i1 - 1];
              ix0 = workingset->nActiveConstr - 1;
              workingset->Wid[i1 - 1] = workingset->Wid[ix0];
              workingset->Wlocalidx[i1 - 1] = workingset->Wlocalidx[ix0];
              for (ix = 0; ix < i; ix++) {
                workingset->ATwset[ix + 22 * (i1 - 1)] = workingset->ATwset[ix +
                  22 * ix0];
              }

              workingset->bwset[i1 - 1] = workingset->bwset[ix0];
              workingset->nActiveConstr = ix0;
              workingset->nWConstr[idxDiag]--;
            }
          }
        }
      }
    }
  }

  if ((nDepInd != -1) && (workingset->nActiveConstr <= 40)) {
    bool guard1;
    bool okWorkingSet;
    RemoveDependentIneq_(workingset, qrmanager, memspace, 1.0);
    okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_float,
      solution->xstar, workingset, qrmanager);
    guard1 = false;
    if (!okWorkingSet) {
      RemoveDependentIneq_(workingset, qrmanager, memspace, 10.0);
      okWorkingSet = feasibleX0ForWorkingSet(memspace->workspace_float,
        solution->xstar, workingset, qrmanager);
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
      tol = b_maxConstraintViolation(workingset, solution->xstar);
      if (tol > 1.0E-6) {
        solution->state = -2;
      }
    }
  } else {
    solution->state = -3;
    idxDiag = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
    ix = workingset->nActiveConstr;
    for (ix0 = idxDiag; ix0 <= ix; ix0++) {
      workingset->isActiveConstr[(workingset->isActiveIdx[workingset->Wid[ix0 -
        1] - 1] + workingset->Wlocalidx[ix0 - 1]) - 2] = false;
    }

    workingset->nWConstr[2] = 0;
    workingset->nWConstr[3] = 0;
    workingset->nWConstr[4] = 0;
    workingset->nActiveConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  }
}

/*
 * Arguments    : e_struct_T *workingset
 *                g_struct_T *qrmanager
 *                struct_T *memspace
 *                double tolfactor
 * Return Type  : void
 */
static void RemoveDependentIneq_(e_struct_T *workingset, g_struct_T *qrmanager,
  struct_T *memspace, double tolfactor)
{
  int idx;
  int idxDiag;
  int idx_col;
  int k;
  int nActiveConstr_tmp;
  int nFixedConstr;
  int nVar;
  nActiveConstr_tmp = workingset->nActiveConstr;
  nFixedConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  nVar = workingset->nVar;
  if ((workingset->nWConstr[2] + workingset->nWConstr[3]) + workingset->
      nWConstr[4] > 0) {
    double maxDiag;
    double tol;
    int nDepIneq;
    idxDiag = workingset->nVar;
    nDepIneq = workingset->nActiveConstr;
    if (idxDiag >= nDepIneq) {
      nDepIneq = idxDiag;
    }

    tol = tolfactor * fmin(1.4901161193847656E-8, 2.2204460492503131E-15 *
      (double)nDepIneq);
    for (idx = 0; idx < nFixedConstr; idx++) {
      qrmanager->jpvt[idx] = 1;
    }

    idxDiag = nFixedConstr + 1;
    if (idxDiag <= nActiveConstr_tmp) {
      memset(&qrmanager->jpvt[idxDiag + -1], 0, (unsigned int)
             ((nActiveConstr_tmp - idxDiag) + 1) * sizeof(int));
    }

    for (idx_col = 0; idx_col < nActiveConstr_tmp; idx_col++) {
      nDepIneq = 40 * idx_col;
      idx = 22 * idx_col;
      idxDiag = (unsigned char)nVar;
      for (k = 0; k < idxDiag; k++) {
        qrmanager->QR[nDepIneq + k] = workingset->ATwset[idx + k];
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
      nDepIneq = workingset->nActiveConstr;
      if (idxDiag <= nDepIneq) {
        nDepIneq = idxDiag;
      }

      qrmanager->minRowCol = nDepIneq;
      xgeqp3(qrmanager->QR, workingset->nVar, workingset->nActiveConstr,
             qrmanager->jpvt, qrmanager->tau);
    }

    nDepIneq = 0;
    for (idx = workingset->nActiveConstr - 1; idx + 1 > nVar; idx--) {
      nDepIneq++;
      memspace->workspace_int[nDepIneq - 1] = qrmanager->jpvt[idx];
    }

    maxDiag = fabs(qrmanager->QR[0]);
    for (idxDiag = 0; idxDiag < idx; idxDiag++) {
      maxDiag = fmax(maxDiag, fabs(qrmanager->QR[(40 * (idxDiag + 1) + idxDiag)
        + 1]));
    }

    if (idx + 1 <= workingset->nVar) {
      idxDiag = idx + 40 * idx;
      while ((idx + 1 > nFixedConstr) && (fabs(qrmanager->QR[idxDiag]) < tol *
              maxDiag)) {
        nDepIneq++;
        memspace->workspace_int[nDepIneq - 1] = qrmanager->jpvt[idx];
        idx--;
        idxDiag -= 41;
      }
    }

    countsort(memspace->workspace_int, nDepIneq, memspace->workspace_sort,
              nFixedConstr + 1, workingset->nActiveConstr);
    for (idx = nDepIneq; idx >= 1; idx--) {
      removeConstr(workingset, memspace->workspace_int[idx - 1]);
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                int idx_local
 * Return Type  : void
 */
static void addAeqConstr(e_struct_T *obj, int idx_local)
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
    obj->isActiveConstr[(obj->isActiveIdx[1] + idx_local) - 2] = true;
    obj->nActiveConstr++;
    i = obj->nActiveConstr - 1;
    obj->Wid[i] = 2;
    obj->Wlocalidx[i] = idx_local;
    iAeq0 = 22 * (idx_local - 1);
    iAw0 = 22 * (obj->nActiveConstr - 1);
    i1 = (unsigned char)obj->nVar;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset[iAw0 + idx] = obj->Aeq[iAeq0 + idx];
    }

    obj->bwset[i] = obj->beq[idx_local - 1];
  } else {
    int i;
    int i1;
    int iAeq0;
    int iAw0;
    obj->nActiveConstr++;
    i = obj->nActiveConstr - 1;
    obj->Wid[i] = obj->Wid[totalEq];
    obj->Wlocalidx[i] = obj->Wlocalidx[totalEq];
    i1 = (unsigned char)obj->nVar;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset[idx + 22 * i] = obj->ATwset[idx + 22 * totalEq];
    }

    obj->bwset[i] = obj->bwset[totalEq];
    obj->nWConstr[1]++;
    obj->isActiveConstr[(obj->isActiveIdx[1] + idx_local) - 2] = true;
    obj->Wid[totalEq] = 2;
    obj->Wlocalidx[totalEq] = idx_local;
    iAeq0 = 22 * (idx_local - 1);
    iAw0 = 22 * totalEq;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset[iAw0 + idx] = obj->Aeq[iAeq0 + idx];
    }

    obj->bwset[totalEq] = obj->beq[idx_local - 1];
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                int TYPE
 *                int idx_local
 * Return Type  : void
 */
static void addBoundToActiveSetMatrix_(e_struct_T *obj, int TYPE, int idx_local)
{
  int colOffset;
  int i;
  int idx_bnd_local;
  obj->nWConstr[TYPE - 1]++;
  obj->isActiveConstr[(obj->isActiveIdx[TYPE - 1] + idx_local) - 2] = true;
  obj->nActiveConstr++;
  i = obj->nActiveConstr - 1;
  obj->Wid[i] = TYPE;
  obj->Wlocalidx[i] = idx_local;
  colOffset = 22 * i - 1;
  if (TYPE == 5) {
    idx_bnd_local = obj->indexUB[idx_local - 1];
    obj->bwset[i] = obj->ub[idx_bnd_local - 1];
  } else {
    idx_bnd_local = obj->indexLB[idx_local - 1];
    obj->bwset[i] = obj->lb[idx_bnd_local - 1];
  }

  if (idx_bnd_local - 2 >= 0) {
    memset(&obj->ATwset[colOffset + 1], 0, (unsigned int)(idx_bnd_local - 1) *
           sizeof(double));
  }

  obj->ATwset[idx_bnd_local + colOffset] = 2.0 * (double)(TYPE == 5) - 1.0;
  i = idx_bnd_local + 1;
  idx_bnd_local = obj->nVar;
  if (i <= idx_bnd_local) {
    memset(&obj->ATwset[i + colOffset], 0, (unsigned int)((((idx_bnd_local +
               colOffset) - i) - colOffset) + 1) * sizeof(double));
  }

  switch (obj->probType) {
   case 3:
   case 2:
    break;

   default:
    obj->ATwset[obj->nVar + colOffset] = -1.0;
    break;
  }
}

/*
 * Arguments    : double workspace[880]
 *                int nVar
 *                const double grad[22]
 *                const double AineqTrans[132]
 *                const double AeqTrans[66]
 *                const int finiteLB[22]
 *                int mLB
 *                const int finiteUB[22]
 *                const double lambda[40]
 * Return Type  : void
 */
static void b_computeGradLag(double workspace[880], int nVar, const double grad
  [22], const double AineqTrans[132], const double AeqTrans[66], const int
  finiteLB[22], int mLB, const int finiteUB[22], const double lambda[40])
{
  int i;
  int i1;
  int ia;
  int iac;
  int ix;
  i = (unsigned char)nVar;
  memcpy(&workspace[0], &grad[0], (unsigned int)i * sizeof(double));
  ix = 0;
  for (iac = 0; iac <= 44; iac += 22) {
    i = iac + nVar;
    for (ia = iac + 1; ia <= i; ia++) {
      i1 = (ia - iac) - 1;
      workspace[i1] += AeqTrans[ia - 1] * lambda[ix];
    }

    ix++;
  }

  ix = 3;
  for (iac = 0; iac <= 110; iac += 22) {
    i = iac + nVar;
    for (ia = iac + 1; ia <= i; ia++) {
      i1 = (ia - iac) - 1;
      workspace[i1] += AineqTrans[ia - 1] * lambda[ix];
    }

    ix++;
  }

  i = (unsigned char)mLB;
  for (ix = 0; ix < i; ix++) {
    i1 = finiteLB[ix];
    workspace[i1 - 1] -= lambda[ix + 9];
  }

  for (ix = 0; ix < 9; ix++) {
    i = finiteUB[ix];
    workspace[i - 1] += lambda[((unsigned char)mLB + ix) + 9];
  }
}

/*
 * Arguments    : c_struct_T *TrialState
 *                f_struct_T *MeritFunction
 *                const i_coder_internal_stickyStruct *FcnEvaluator
 *                d_struct_T *FiniteDifferences
 *                struct_T *memspace
 *                e_struct_T *WorkingSet
 *                double Hessian[81]
 *                g_struct_T *QRManager
 *                h_struct_T *CholManager
 *                i_struct_T *QPObjective
 * Return Type  : void
 */
static void b_driver(c_struct_T *TrialState, f_struct_T *MeritFunction, const
                     i_coder_internal_stickyStruct *FcnEvaluator, d_struct_T
                     *FiniteDifferences, struct_T *memspace, e_struct_T
                     *WorkingSet, double Hessian[81], g_struct_T *QRManager,
                     h_struct_T *CholManager, i_struct_T *QPObjective)
{
  static const signed char b_Hessian[81] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const char qpoptions_SolverName[7] = { 'f', 'm', 'i', 'n', 'c', 'o',
    'n' };

  j_struct_T Flags;
  k_struct_T expl_temp;
  int i;
  int i1;
  int iCol;
  int ia;
  int ix;
  int mLB;
  int nVar_tmp;
  int qpoptions_MaxIterations;
  memset(&QPObjective->grad[0], 0, 22U * sizeof(double));
  memset(&QPObjective->Hx[0], 0, 21U * sizeof(double));
  QPObjective->hasLinear = true;
  QPObjective->nvar = 9;
  QPObjective->beta = 0.0;
  QPObjective->rho = 0.0;
  QPObjective->objtype = 3;
  QPObjective->prev_objtype = 3;
  QPObjective->prev_nvar = 0;
  QPObjective->prev_hasLinear = false;
  QPObjective->gammaScalar = 0.0;
  CholManager->ndims = 0;
  CholManager->info = 0;
  CholManager->ConvexCheck = true;
  CholManager->workspace_ = rtInf;
  memset(&CholManager->FMat[0], 0, 1600U * sizeof(double));
  memset(&QRManager->QR[0], 0, 1600U * sizeof(double));
  memset(&QRManager->Q[0], 0, 1600U * sizeof(double));
  QRManager->mrows = 0;
  QRManager->ncols = 0;
  memset(&QRManager->jpvt[0], 0, 40U * sizeof(int));
  memset(&QRManager->tau[0], 0, 40U * sizeof(double));
  QRManager->minRowCol = 0;
  QRManager->usedPivoting = false;
  for (i = 0; i < 81; i++) {
    Hessian[i] = b_Hessian[i];
  }

  nVar_tmp = WorkingSet->nVar;
  mLB = WorkingSet->sizes[3];
  iCol = WorkingSet->nVar;
  ix = WorkingSet->sizes[3] + 15;
  if (iCol >= ix) {
    ix = iCol;
  }

  qpoptions_MaxIterations = 10 * ix;
  TrialState->steplength = 1.0;
  test_exit(MeritFunction, WorkingSet, TrialState, &Flags.fevalOK, &Flags.done,
            &Flags.stepAccepted, &Flags.failedLineSearch, &Flags.stepType);
  iCol = -1;
  i = (unsigned char)nVar_tmp;
  for (ix = 0; ix < 3; ix++) {
    memcpy(&TrialState->JacCeqTrans_old[iCol + 1], &WorkingSet->Aeq[iCol + 1],
           (unsigned int)i * sizeof(double));
    iCol += 22;
  }

  TrialState->sqpFval_old = TrialState->sqpFval;
  for (ix = 0; ix < 9; ix++) {
    TrialState->xstarsqp_old[ix] = TrialState->xstarsqp[ix];
    TrialState->grad_old[ix] = TrialState->grad[ix];
  }

  for (ix = 0; ix < 6; ix++) {
    TrialState->cIneq_old[ix] = TrialState->cIneq[ix];
  }

  TrialState->cEq_old[0] = TrialState->cEq[0];
  TrialState->cEq_old[1] = TrialState->cEq[1];
  TrialState->cEq_old[2] = TrialState->cEq[2];
  if (!Flags.done) {
    TrialState->sqpIterations = 1;
  }

  while (!Flags.done) {
    double d;
    int i2;
    while (!(Flags.stepAccepted || Flags.failedLineSearch)) {
      double constrViolationIneq;
      if (Flags.stepType != 3) {
        updateWorkingSetForNewQP(TrialState->xstarsqp, WorkingSet,
          TrialState->cIneq, TrialState->cEq, mLB);
      }

      expl_temp.MaxIterations = qpoptions_MaxIterations;
      for (i1 = 0; i1 < 7; i1++) {
        expl_temp.SolverName[i1] = qpoptions_SolverName[i1];
      }

      step(&Flags, Hessian, TrialState, MeritFunction, memspace, WorkingSet,
           QRManager, CholManager, QPObjective, &expl_temp);
      if (Flags.stepAccepted) {
        for (iCol = 0; iCol < i; iCol++) {
          TrialState->xstarsqp[iCol] += TrialState->delta_x[iCol];
        }

        TrialState->sqpFval = evalObjAndConstr
          (&FcnEvaluator->next.next.next.next.next.next.next.next.value.workspace,
           TrialState->xstarsqp, TrialState->cEq, &iCol);
        Flags.fevalOK = (iCol == 1);
        TrialState->FunctionEvaluations++;
        computeLinearResiduals(TrialState->xstarsqp, nVar_tmp, TrialState->cIneq,
          WorkingSet->Aineq);
        if (Flags.fevalOK) {
          constrViolationIneq = 0.0;
          for (ix = 0; ix < 6; ix++) {
            d = TrialState->cIneq[ix];
            if (d > 0.0) {
              constrViolationIneq += d;
            }
          }

          MeritFunction->phiFullStep = TrialState->sqpFval +
            MeritFunction->penaltyParam * (((fabs(TrialState->cEq[0]) + fabs
            (TrialState->cEq[1])) + fabs(TrialState->cEq[2])) +
            constrViolationIneq);
        } else {
          MeritFunction->phiFullStep = rtInf;
        }
      }

      if ((Flags.stepType == 1) && Flags.stepAccepted && Flags.fevalOK &&
          (MeritFunction->phi < MeritFunction->phiFullStep) &&
          (TrialState->sqpFval < TrialState->sqpFval_old)) {
        Flags.stepType = 3;
        Flags.stepAccepted = false;
      } else {
        double alpha;
        bool evalWellDefined;
        bool socTaken;
        if ((Flags.stepType == 3) && Flags.stepAccepted) {
          socTaken = true;
        } else {
          socTaken = false;
        }

        evalWellDefined = Flags.fevalOK;
        i1 = WorkingSet->nVar;
        alpha = 1.0;
        iCol = 1;
        constrViolationIneq = MeritFunction->phiFullStep;
        if (i1 - 1 >= 0) {
          memcpy(&TrialState->searchDir[0], &TrialState->delta_x[0], (unsigned
                  int)i1 * sizeof(double));
        }

        int exitg1;
        do {
          exitg1 = 0;
          if (TrialState->FunctionEvaluations < 900) {
            if (evalWellDefined && (constrViolationIneq <= MeritFunction->phi +
                                    alpha * 0.0001 * MeritFunction->phiPrimePlus))
            {
              exitg1 = 1;
            } else {
              bool exitg2;
              bool tooSmallX;
              alpha *= 0.7;
              i2 = (unsigned char)i1;
              for (ix = 0; ix < i2; ix++) {
                TrialState->delta_x[ix] = alpha * TrialState->xstar[ix];
              }

              if (socTaken) {
                constrViolationIneq = alpha * alpha;
                if ((i1 >= 1) && (!(constrViolationIneq == 0.0))) {
                  for (ix = 0; ix < i1; ix++) {
                    TrialState->delta_x[ix] += constrViolationIneq *
                      TrialState->socDirection[ix];
                  }
                }
              }

              tooSmallX = true;
              ix = 0;
              exitg2 = false;
              while ((!exitg2) && (ix <= (unsigned char)i1 - 1)) {
                if (1.0E-6 * fmax(1.0, fabs(TrialState->xstarsqp[ix])) <= fabs
                    (TrialState->delta_x[ix])) {
                  tooSmallX = false;
                  exitg2 = true;
                } else {
                  ix++;
                }
              }

              if (tooSmallX) {
                iCol = -2;
                exitg1 = 1;
              } else {
                for (ix = 0; ix < i2; ix++) {
                  TrialState->xstarsqp[ix] = TrialState->xstarsqp_old[ix] +
                    TrialState->delta_x[ix];
                }

                TrialState->sqpFval = evalObjAndConstr
                  (&FcnEvaluator->next.next.next.next.next.next.next.next.value.workspace,
                   TrialState->xstarsqp, TrialState->cEq, &ix);
                computeLinearResiduals(TrialState->xstarsqp, i1,
                  TrialState->cIneq, WorkingSet->Aineq);
                TrialState->FunctionEvaluations++;
                evalWellDefined = (ix == 1);
                if (evalWellDefined) {
                  constrViolationIneq = 0.0;
                  for (ix = 0; ix < 6; ix++) {
                    d = TrialState->cIneq[ix];
                    if (d > 0.0) {
                      constrViolationIneq += d;
                    }
                  }

                  constrViolationIneq = TrialState->sqpFval +
                    MeritFunction->penaltyParam * (((fabs(TrialState->cEq[0]) +
                    fabs(TrialState->cEq[1])) + fabs(TrialState->cEq[2])) +
                    constrViolationIneq);
                } else {
                  constrViolationIneq = rtInf;
                }
              }
            }
          } else {
            iCol = 0;
            exitg1 = 1;
          }
        } while (exitg1 == 0);

        Flags.fevalOK = evalWellDefined;
        TrialState->steplength = alpha;
        if (iCol > 0) {
          Flags.stepAccepted = true;
        } else {
          Flags.failedLineSearch = true;
        }
      }
    }

    if (Flags.stepAccepted && (!Flags.failedLineSearch)) {
      for (ix = 0; ix < i; ix++) {
        TrialState->xstarsqp[ix] = TrialState->xstarsqp_old[ix] +
          TrialState->delta_x[ix];
      }

      i1 = (unsigned char)(mLB + 18);
      for (ix = 0; ix < i1; ix++) {
        d = TrialState->lambdasqp[ix];
        d += TrialState->steplength * (TrialState->lambda[ix] - d);
        TrialState->lambdasqp[ix] = d;
      }

      TrialState->sqpFval_old = TrialState->sqpFval;
      for (ix = 0; ix < 9; ix++) {
        TrialState->xstarsqp_old[ix] = TrialState->xstarsqp[ix];
        TrialState->grad_old[ix] = TrialState->grad[ix];
      }

      for (ix = 0; ix < 6; ix++) {
        TrialState->cIneq_old[ix] = TrialState->cIneq[ix];
      }

      TrialState->cEq_old[0] = TrialState->cEq[0];
      TrialState->cEq_old[1] = TrialState->cEq[1];
      TrialState->cEq_old[2] = TrialState->cEq[2];
      computeFiniteDifferences(FiniteDifferences, TrialState->sqpFval,
        TrialState->xstarsqp, TrialState->grad);
      TrialState->FunctionEvaluations += FiniteDifferences->numEvals;
      TrialState->sqpFval = evalObjAndConstrAndDerivatives
        (&FcnEvaluator->next.next.next.next.next.next.next.next.value.workspace,
         TrialState->xstarsqp, TrialState->cEq, WorkingSet->Aeq, &iCol);
      Flags.fevalOK = (iCol == 1);
    } else {
      TrialState->sqpFval = TrialState->sqpFval_old;
      memcpy(&TrialState->xstarsqp[0], &TrialState->xstarsqp_old[0], 9U * sizeof
             (double));
      for (ix = 0; ix < 6; ix++) {
        TrialState->cIneq[ix] = TrialState->cIneq_old[ix];
      }

      TrialState->cEq[0] = TrialState->cEq_old[0];
      TrialState->cEq[1] = TrialState->cEq_old[1];
      TrialState->cEq[2] = TrialState->cEq_old[2];
    }

    b_test_exit(&Flags, memspace, MeritFunction, WorkingSet, TrialState,
                QRManager);
    if ((!Flags.done) && Flags.stepAccepted) {
      Flags.stepAccepted = false;
      Flags.stepType = 1;
      Flags.failedLineSearch = false;
      memcpy(&TrialState->delta_gradLag[0], &TrialState->grad[0], (unsigned int)
             i * sizeof(double));
      if (nVar_tmp >= 1) {
        for (ix = 0; ix < nVar_tmp; ix++) {
          TrialState->delta_gradLag[ix] -= TrialState->grad_old[ix];
        }
      }

      ix = 0;
      for (iCol = 0; iCol <= 44; iCol += 22) {
        i1 = iCol + nVar_tmp;
        for (ia = iCol + 1; ia <= i1; ia++) {
          i2 = (ia - iCol) - 1;
          TrialState->delta_gradLag[i2] += WorkingSet->Aeq[ia - 1] *
            TrialState->lambdasqp[ix];
        }

        ix++;
      }

      ix = 0;
      for (iCol = 0; iCol <= 44; iCol += 22) {
        i1 = iCol + nVar_tmp;
        for (ia = iCol + 1; ia <= i1; ia++) {
          i2 = (ia - iCol) - 1;
          TrialState->delta_gradLag[i2] += TrialState->JacCeqTrans_old[ia - 1] *
            -TrialState->lambdasqp[ix];
        }

        ix++;
      }

      iCol = -1;
      for (ix = 0; ix < 3; ix++) {
        memcpy(&TrialState->JacCeqTrans_old[iCol + 1], &WorkingSet->Aeq[iCol + 1],
               (unsigned int)i * sizeof(double));
        iCol += 22;
      }

      BFGSUpdate(nVar_tmp, Hessian, TrialState->delta_x,
                 TrialState->delta_gradLag, memspace->workspace_float);
      TrialState->sqpIterations++;
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                const double x[22]
 * Return Type  : double
 */
static double b_maxConstraintViolation(e_struct_T *obj, const double x[22])
{
  double v;
  int idx;
  int idxUB;
  if (obj->probType == 2) {
    v = c_maxConstraintViolation_AMats_(obj, x);
  } else {
    v = d_maxConstraintViolation_AMats_(obj, x);
  }

  idxUB = (unsigned char)obj->sizes[3];
  for (idx = 0; idx < idxUB; idx++) {
    int idxLB;
    idxLB = obj->indexLB[idx] - 1;
    v = fmax(v, -x[idxLB] - obj->lb[idxLB]);
  }

  for (idx = 0; idx < 9; idx++) {
    idxUB = obj->indexUB[idx] - 1;
    v = fmax(v, x[idxUB] - obj->ub[idxUB]);
  }

  return v;
}

/*
 * Arguments    : j_struct_T *Flags
 *                struct_T *memspace
 *                f_struct_T *MeritFunction
 *                e_struct_T *WorkingSet
 *                c_struct_T *TrialState
 *                g_struct_T *QRManager
 * Return Type  : void
 */
static void b_test_exit(j_struct_T *Flags, struct_T *memspace, f_struct_T
  *MeritFunction, e_struct_T *WorkingSet, c_struct_T *TrialState, g_struct_T
  *QRManager)
{
  double optimRelativeFactor;
  double s;
  double smax;
  int b_i;
  int fullRank_R;
  int i;
  int idx_max;
  int ix;
  int j;
  int mLB;
  int nVar;
  int rankR;
  bool dxTooSmall;
  bool exitg1;
  bool isFeasible;
  nVar = WorkingSet->nVar;
  mLB = WorkingSet->sizes[3];
  i = (unsigned char)(WorkingSet->sizes[3] + 18);
  memcpy(&TrialState->lambdaStopTest[0], &TrialState->lambdasqp[0], (unsigned
          int)i * sizeof(double));
  computeGradLag(TrialState->gradLag, WorkingSet->nVar, TrialState->grad,
                 WorkingSet->Aineq, WorkingSet->Aeq, WorkingSet->indexLB,
                 WorkingSet->sizes[3], WorkingSet->indexUB,
                 TrialState->lambdaStopTest);
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = fabs(TrialState->grad[0]);
      for (ix = 2; ix <= nVar; ix++) {
        s = fabs(TrialState->grad[ix - 1]);
        if (s > smax) {
          idx_max = ix;
          smax = s;
        }
      }
    }
  }

  optimRelativeFactor = fmax(1.0, fabs(TrialState->grad[idx_max - 1]));
  if (rtIsInf(optimRelativeFactor)) {
    optimRelativeFactor = 1.0;
  }

  smax = fmax(fmax(fmax(0.0, fabs(TrialState->cEq[0])), fabs(TrialState->cEq[1])),
              fabs(TrialState->cEq[2]));
  for (rankR = 0; rankR < 6; rankR++) {
    smax = fmax(smax, TrialState->cIneq[rankR]);
  }

  ix = (unsigned char)WorkingSet->sizes[3];
  for (rankR = 0; rankR < ix; rankR++) {
    smax = fmax(smax, -1.0 - TrialState->xstarsqp[WorkingSet->indexLB[rankR] - 1]);
  }

  for (rankR = 0; rankR < 9; rankR++) {
    idx_max = WorkingSet->indexUB[rankR] - 1;
    smax = fmax(smax, TrialState->xstarsqp[idx_max] - (double)iv[idx_max]);
  }

  MeritFunction->nlpPrimalFeasError = smax;
  if (TrialState->sqpIterations == 0) {
    MeritFunction->feasRelativeFactor = fmax(1.0, smax);
  }

  isFeasible = (smax <= 1.0E-6 * MeritFunction->feasRelativeFactor);
  dxTooSmall = true;
  smax = 0.0;
  ix = (unsigned char)WorkingSet->nVar;
  rankR = 0;
  exitg1 = false;
  while ((!exitg1) && (rankR <= ix - 1)) {
    dxTooSmall = ((!rtIsInf(TrialState->gradLag[rankR])) && (!rtIsNaN
      (TrialState->gradLag[rankR])));
    if (!dxTooSmall) {
      exitg1 = true;
    } else {
      smax = fmax(smax, fabs(TrialState->gradLag[rankR]));
      rankR++;
    }
  }

  MeritFunction->nlpDualFeasError = smax;
  if (!dxTooSmall) {
    Flags->done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    double nlpComplErrorTmp;
    MeritFunction->nlpComplError = computeComplError(TrialState->xstarsqp,
      TrialState->cIneq, WorkingSet->indexLB, WorkingSet->sizes[3],
      WorkingSet->indexUB, TrialState->lambdaStopTest);
    MeritFunction->firstOrderOpt = fmax(smax, MeritFunction->nlpComplError);
    if (TrialState->sqpIterations > 1) {
      b_computeGradLag(memspace->workspace_float, WorkingSet->nVar,
                       TrialState->grad, WorkingSet->Aineq, WorkingSet->Aeq,
                       WorkingSet->indexLB, WorkingSet->sizes[3],
                       WorkingSet->indexUB, TrialState->lambdaStopTestPrev);
      s = 0.0;
      rankR = 0;
      while ((rankR <= ix - 1) && ((!rtIsInf(memspace->workspace_float[rankR])) &&
              (!rtIsNaN(memspace->workspace_float[rankR])))) {
        s = fmax(s, fabs(memspace->workspace_float[rankR]));
        rankR++;
      }

      nlpComplErrorTmp = computeComplError(TrialState->xstarsqp,
        TrialState->cIneq, WorkingSet->indexLB, WorkingSet->sizes[3],
        WorkingSet->indexUB, TrialState->lambdaStopTestPrev);
      if ((s < smax) && (nlpComplErrorTmp < MeritFunction->nlpComplError)) {
        MeritFunction->nlpDualFeasError = s;
        MeritFunction->nlpComplError = nlpComplErrorTmp;
        MeritFunction->firstOrderOpt = fmax(s, nlpComplErrorTmp);
        memcpy(&TrialState->lambdaStopTest[0], &TrialState->lambdaStopTestPrev[0],
               (unsigned int)i * sizeof(double));
      } else {
        memcpy(&TrialState->lambdaStopTestPrev[0], &TrialState->lambdaStopTest[0],
               (unsigned int)i * sizeof(double));
      }
    } else {
      memcpy(&TrialState->lambdaStopTestPrev[0], &TrialState->lambdaStopTest[0],
             (unsigned int)i * sizeof(double));
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
          dxTooSmall = true;
          rankR = 0;
          exitg1 = false;
          while ((!exitg1) && (rankR <= ix - 1)) {
            if (1.0E-6 * fmax(1.0, fabs(TrialState->xstarsqp[rankR])) <= fabs
                (TrialState->delta_x[rankR])) {
              dxTooSmall = false;
              exitg1 = true;
            } else {
              rankR++;
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
                updateWorkingSetForNewQP(TrialState->xstarsqp, WorkingSet,
                  TrialState->cIneq, TrialState->cEq, WorkingSet->sizes[3]);
                if (nActiveConstr - 1 >= 0) {
                  memset(&TrialState->lambda[0], 0, (unsigned int)nActiveConstr *
                         sizeof(double));
                }

                factorQRE(QRManager, WorkingSet->ATwset, nVar, nActiveConstr);
                computeQ_(QRManager, QRManager->mrows);
                memset(&memspace->workspace_float[0], 0, (unsigned int)ix *
                       sizeof(double));
                ix = 40 * (nVar - 1) + 1;
                for (rankR = 1; rankR <= ix; rankR += 40) {
                  smax = 0.0;
                  idx_max = rankR + nVar;
                  for (fullRank_R = rankR; fullRank_R < idx_max; fullRank_R++) {
                    smax += QRManager->Q[fullRank_R - 1] * TrialState->
                      grad[fullRank_R - rankR];
                  }

                  idx_max = div_nde_s32_floor(rankR - 1, 40);
                  memspace->workspace_float[idx_max] -= smax;
                }

                if (nVar >= nActiveConstr) {
                  idx_max = nVar;
                } else {
                  idx_max = nActiveConstr;
                }

                smax = fabs(QRManager->QR[0]) * fmin(1.4901161193847656E-8,
                  (double)idx_max * 2.2204460492503131E-16);
                if (nVar <= nActiveConstr) {
                  fullRank_R = nVar;
                } else {
                  fullRank_R = nActiveConstr;
                }

                rankR = 0;
                idx_max = 0;
                while ((rankR < fullRank_R) && (fabs(QRManager->QR[idx_max]) >
                        smax)) {
                  rankR++;
                  idx_max += 41;
                }

                if (rankR != 0) {
                  for (j = rankR; j >= 1; j--) {
                    idx_max = (j + (j - 1) * 40) - 1;
                    memspace->workspace_float[j - 1] /= QRManager->QR[idx_max];
                    for (b_i = 0; b_i <= j - 2; b_i++) {
                      ix = (j - b_i) - 2;
                      memspace->workspace_float[ix] -= memspace->
                        workspace_float[j - 1] * QRManager->QR[(idx_max - b_i) -
                        1];
                    }
                  }
                }

                if (nActiveConstr <= fullRank_R) {
                  fullRank_R = nActiveConstr;
                }

                for (rankR = 0; rankR < fullRank_R; rankR++) {
                  TrialState->lambda[QRManager->jpvt[rankR] - 1] =
                    memspace->workspace_float[rankR];
                }

                sortLambdaQP(TrialState->lambda, WorkingSet->nActiveConstr,
                             WorkingSet->sizes, WorkingSet->isActiveIdx,
                             WorkingSet->Wid, WorkingSet->Wlocalidx,
                             memspace->workspace_float);
                b_computeGradLag(memspace->workspace_float, nVar,
                                 TrialState->grad, WorkingSet->Aineq,
                                 WorkingSet->Aeq, WorkingSet->indexLB, mLB,
                                 WorkingSet->indexUB, TrialState->lambda);
                smax = 0.0;
                rankR = 0;
                while ((rankR <= (unsigned char)nVar - 1) && ((!rtIsInf
                         (memspace->workspace_float[rankR])) && (!rtIsNaN
                         (memspace->workspace_float[rankR])))) {
                  smax = fmax(smax, fabs(memspace->workspace_float[rankR]));
                  rankR++;
                }

                s = computeComplError(TrialState->xstarsqp, TrialState->cIneq,
                                      WorkingSet->indexLB, mLB,
                                      WorkingSet->indexUB, TrialState->lambda);
                nlpComplErrorTmp = fmax(smax, s);
                if (nlpComplErrorTmp <= fmax(MeritFunction->nlpDualFeasError,
                     MeritFunction->nlpComplError)) {
                  MeritFunction->nlpDualFeasError = smax;
                  MeritFunction->nlpComplError = s;
                  MeritFunction->firstOrderOpt = nlpComplErrorTmp;
                  memcpy(&TrialState->lambdaStopTest[0], &TrialState->lambda[0],
                         (unsigned int)i * sizeof(double));
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
          if (TrialState->sqpIterations >= 100) {
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
 * Arguments    : const double x[18]
 *                double y[6]
 * Return Type  : void
 */
static void b_vecnorm(const double x[18], double y[6])
{
  int j;
  for (j = 0; j < 6; j++) {
    double absxk;
    double scale;
    double t;
    double yv;
    scale = 3.3121686421112381E-170;
    absxk = fabs(x[j]);
    if (absxk > 3.3121686421112381E-170) {
      yv = 1.0;
      scale = absxk;
    } else {
      t = absxk / 3.3121686421112381E-170;
      yv = t * t;
    }

    absxk = fabs(x[j + 6]);
    if (absxk > scale) {
      t = scale / absxk;
      yv = yv * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      yv += t * t;
    }

    absxk = fabs(x[j + 12]);
    if (absxk > scale) {
      t = scale / absxk;
      yv = yv * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      yv += t * t;
    }

    y[j] = scale * sqrt(yv);
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
 *                const double A[1600]
 *                int ia0
 *                const double B[880]
 *                double C[1600]
 * Return Type  : void
 */
static void b_xgemm(int m, int n, int k, const double A[1600], int ia0, const
                    double B[880], double C[1600])
{
  int cr;
  int ic;
  int w;
  if ((m != 0) && (n != 0)) {
    int br;
    int i;
    int i1;
    int lastColC;
    lastColC = 40 * (n - 1);
    for (cr = 0; cr <= lastColC; cr += 40) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        memset(&C[i + -1], 0, (unsigned int)((i1 - i) + 1) * sizeof(double));
      }
    }

    br = -1;
    for (cr = 0; cr <= lastColC; cr += 40) {
      int ar;
      ar = ia0;
      i = cr + 1;
      i1 = cr + m;
      for (ic = i; ic <= i1; ic++) {
        double temp;
        temp = 0.0;
        for (w = 0; w < k; w++) {
          temp += A[(w + ar) - 1] * B[(w + br) + 1];
        }

        C[ic - 1] += temp;
        ar += 40;
      }

      br += 40;
    }
  }
}

/*
 * Arguments    : int n
 *                const double x[1600]
 *                int ix0
 * Return Type  : double
 */
static double b_xnrm2(int n, const double x[1600], int ix0)
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      double scale;
      int kend;
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
 * Arguments    : e_struct_T *obj
 *                const double x[22]
 * Return Type  : double
 */
static double c_maxConstraintViolation_AMats_(e_struct_T *obj, const double x[22])
{
  double c;
  double v;
  int ia;
  int iac;
  int k;
  v = 0.0;
  for (k = 0; k < 6; k++) {
    obj->maxConstrWorkspace[k] = -obj->bineq[k];
  }

  for (iac = 0; iac <= 110; iac += 22) {
    c = 0.0;
    k = iac + 9;
    for (ia = iac + 1; ia <= k; ia++) {
      c += obj->Aineq[ia - 1] * x[(ia - iac) - 1];
    }

    k = div_nde_s32_floor(iac, 22);
    obj->maxConstrWorkspace[k] += c;
  }

  for (k = 0; k < 6; k++) {
    c = obj->maxConstrWorkspace[k] - x[k + 9];
    obj->maxConstrWorkspace[k] = c;
    v = fmax(v, c);
  }

  obj->maxConstrWorkspace[0] = -obj->beq[0];
  obj->maxConstrWorkspace[1] = -obj->beq[1];
  obj->maxConstrWorkspace[2] = -obj->beq[2];
  for (iac = 0; iac <= 44; iac += 22) {
    c = 0.0;
    k = iac + 9;
    for (ia = iac + 1; ia <= k; ia++) {
      c += obj->Aeq[ia - 1] * x[(ia - iac) - 1];
    }

    k = div_nde_s32_floor(iac, 22);
    obj->maxConstrWorkspace[k] += c;
  }

  c = (obj->maxConstrWorkspace[0] - x[15]) + x[18];
  obj->maxConstrWorkspace[0] = c;
  v = fmax(v, fabs(c));
  c = (obj->maxConstrWorkspace[1] - x[16]) + x[19];
  obj->maxConstrWorkspace[1] = c;
  v = fmax(v, fabs(c));
  c = (obj->maxConstrWorkspace[2] - x[17]) + x[20];
  obj->maxConstrWorkspace[2] = c;
  return fmax(v, fabs(c));
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
 *                const double x[22]
 * Return Type  : double
 */
static double c_xnrm2(int n, const double x[22])
{
  double y;
  int k;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[0]);
    } else {
      double scale;
      scale = 3.3121686421112381E-170;
      for (k = 0; k < n; k++) {
        double absxk;
        absxk = fabs(x[k]);
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
 * Arguments    : const double vec[3]
 * Return Type  : int
 */
static int checkVectorNonFinite(const double vec[3])
{
  int idx_current;
  int status;
  bool allFinite;
  status = 1;
  allFinite = true;
  idx_current = 0;
  while (allFinite && (idx_current + 1 <= 3)) {
    allFinite = ((!rtIsInf(vec[idx_current])) && (!rtIsNaN(vec[idx_current])));
    idx_current++;
  }

  if (!allFinite) {
    idx_current--;
    if (rtIsNaN(vec[idx_current])) {
      status = -3;
    } else if (vec[idx_current] < 0.0) {
      status = -1;
    } else {
      status = -2;
    }
  }

  return status;
}

/*
 * Arguments    : const double xCurrent[9]
 *                const double cIneq[6]
 *                const int finiteLB[22]
 *                int mLB
 *                const int finiteUB[22]
 *                const double lambda[40]
 * Return Type  : double
 */
static double computeComplError(const double xCurrent[9], const double cIneq[6],
  const int finiteLB[22], int mLB, const int finiteUB[22], const double lambda
  [40])
{
  double nlpComplError;
  double ubDelta;
  double ubLambda;
  int i;
  int idx;
  nlpComplError = 0.0;
  for (idx = 0; idx < 6; idx++) {
    ubLambda = lambda[idx + 3];
    ubDelta = cIneq[idx];
    nlpComplError = fmax(nlpComplError, fmin(fabs(ubDelta * ubLambda), fmin(fabs
      (ubDelta), ubLambda)));
  }

  i = (unsigned char)mLB;
  for (idx = 0; idx < i; idx++) {
    ubLambda = xCurrent[finiteLB[idx] - 1] - -1.0;
    ubDelta = lambda[idx + 9];
    nlpComplError = fmax(nlpComplError, fmin(fabs(ubLambda * ubDelta), fmin(fabs
      (ubLambda), ubDelta)));
  }

  for (idx = 0; idx < 9; idx++) {
    i = finiteUB[idx];
    ubDelta = (double)iv[i - 1] - xCurrent[i - 1];
    ubLambda = lambda[(mLB + idx) + 9];
    nlpComplError = fmax(nlpComplError, fmin(fabs(ubDelta * ubLambda), fmin(fabs
      (ubDelta), ubLambda)));
  }

  return nlpComplError;
}

/*
 * Arguments    : d_struct_T *obj
 *                double fCurrent
 *                double xk[9]
 *                double gradf[22]
 * Return Type  : bool
 */
static bool computeFiniteDifferences(d_struct_T *obj, double fCurrent, double
  xk[9], double gradf[22])
{
  int i;
  int idx;
  int xk_tmp;
  bool evalOK;
  bool exitg1;
  evalOK = true;
  obj->numEvals = 0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx < 9)) {
    double A_wrench[18];
    double b_xk[9];
    double dv1[9];
    double b_obj[6];
    double dv2[6];
    double dv[3];
    double b_deltaX;
    double d;
    double d1;
    double deltaX;
    double temp;
    double ubDiff;
    signed char i1;
    signed char i2;
    signed char i3;
    bool guard1;
    bool modifiedStep;
    deltaX = 1.4901161193847656E-8 * (1.0 - 2.0 * (double)(xk[idx] < 0.0)) *
      fmax(fabs(xk[idx]), 1.0);
    modifiedStep = false;
    if (xk[idx] >= -1.0) {
      i = iv[idx];
      if (xk[idx] <= i) {
        ubDiff = xk[idx] + deltaX;
        if ((ubDiff > i) || (ubDiff < -1.0)) {
          deltaX = -deltaX;
          modifiedStep = true;
          ubDiff = xk[idx] + deltaX;
          if ((ubDiff > i) || (ubDiff < -1.0)) {
            ubDiff = (double)i - xk[idx];
            if (xk[idx] - -1.0 <= ubDiff) {
              deltaX = -(xk[idx] - -1.0);
            } else {
              deltaX = ubDiff;
            }
          }
        }
      }
    }

    b_deltaX = deltaX;
    temp = xk[idx];
    xk[idx] += deltaX;

    /* 'mso_forces:62' @(W) objective(W, center, N, T_max, T_min) */
    /*  Solver */
    /*  目标函数 */
    /* 'mso_forces:16' W = reshape(W, [3, length(W) / 3]); */
    /* 'mso_forces:17' W = W ./ vecnorm(W, 2, 1); */
    /* 'mso_forces:19' A_wrench = [eye(N); -eye(N)] * pinv(W); */
    vecnorm(xk, dv);
    for (i = 0; i < 3; i++) {
      ubDiff = dv[i];
      b_xk[3 * i] = xk[3 * i] / ubDiff;
      xk_tmp = 3 * i + 1;
      b_xk[xk_tmp] = xk[xk_tmp] / ubDiff;
      xk_tmp = 3 * i + 2;
      b_xk[xk_tmp] = xk[xk_tmp] / ubDiff;
    }

    pinv(b_xk, dv1);
    for (i = 0; i < 6; i++) {
      i1 = iv1[i];
      i2 = iv1[i + 6];
      i3 = iv1[i + 12];
      for (xk_tmp = 0; xk_tmp < 3; xk_tmp++) {
        A_wrench[i + 6 * xk_tmp] = ((double)i1 * dv1[3 * xk_tmp] + (double)i2 *
          dv1[3 * xk_tmp + 1]) + (double)i3 * dv1[3 * xk_tmp + 2];
      }
    }

    /* 'mso_forces:20' b_wrench = [T_max; -T_min]; */
    /* 'mso_forces:22' distances = (b_wrench - A_wrench * p) ./ vecnorm(A_wrench, 2, 2); */
    /*  返回最小距离（负号用于最大化） */
    /* 'mso_forces:25' obj_val = -min(distances); */
    b_vecnorm(A_wrench, dv2);
    b_obj[0] = obj->objfun.workspace.T_max[0];
    b_obj[3] = -obj->objfun.workspace.T_min[0];
    b_obj[1] = obj->objfun.workspace.T_max[1];
    b_obj[4] = -obj->objfun.workspace.T_min[1];
    b_obj[2] = obj->objfun.workspace.T_max[2];
    b_obj[5] = -obj->objfun.workspace.T_min[2];
    ubDiff = obj->objfun.workspace.center[0];
    d = obj->objfun.workspace.center[1];
    d1 = obj->objfun.workspace.center[2];
    for (i = 0; i < 6; i++) {
      b_obj[i] = (b_obj[i] - ((A_wrench[i] * ubDiff + A_wrench[i + 6] * d) +
        A_wrench[i + 12] * d1)) / dv2[i];
    }

    ubDiff = -minimum(b_obj);
    evalOK = ((!rtIsInf(ubDiff)) && (!rtIsNaN(ubDiff)));
    if (evalOK) {
      xk[idx] = temp;
    }

    obj->f_1 = ubDiff;
    obj->numEvals++;
    guard1 = false;
    if (!evalOK) {
      if (!modifiedStep) {
        b_deltaX = -deltaX;
        ubDiff = xk[idx] - deltaX;
        if ((ubDiff >= -1.0) && (ubDiff <= iv[idx])) {
          modifiedStep = true;
        } else {
          modifiedStep = false;
        }

        if ((!obj->hasBounds) || modifiedStep) {
          temp = xk[idx];
          xk[idx] = ubDiff;

          /* 'mso_forces:62' @(W) objective(W, center, N, T_max, T_min) */
          /*  Solver */
          /*  目标函数 */
          /* 'mso_forces:16' W = reshape(W, [3, length(W) / 3]); */
          /* 'mso_forces:17' W = W ./ vecnorm(W, 2, 1); */
          /* 'mso_forces:19' A_wrench = [eye(N); -eye(N)] * pinv(W); */
          vecnorm(xk, dv);
          for (i = 0; i < 3; i++) {
            ubDiff = dv[i];
            b_xk[3 * i] = xk[3 * i] / ubDiff;
            xk_tmp = 3 * i + 1;
            b_xk[xk_tmp] = xk[xk_tmp] / ubDiff;
            xk_tmp = 3 * i + 2;
            b_xk[xk_tmp] = xk[xk_tmp] / ubDiff;
          }

          pinv(b_xk, dv1);
          for (i = 0; i < 6; i++) {
            i1 = iv1[i];
            i2 = iv1[i + 6];
            i3 = iv1[i + 12];
            for (xk_tmp = 0; xk_tmp < 3; xk_tmp++) {
              A_wrench[i + 6 * xk_tmp] = ((double)i1 * dv1[3 * xk_tmp] + (double)
                i2 * dv1[3 * xk_tmp + 1]) + (double)i3 * dv1[3 * xk_tmp + 2];
            }
          }

          /* 'mso_forces:20' b_wrench = [T_max; -T_min]; */
          /* 'mso_forces:22' distances = (b_wrench - A_wrench * p) ./ vecnorm(A_wrench, 2, 2); */
          /*  返回最小距离（负号用于最大化） */
          /* 'mso_forces:25' obj_val = -min(distances); */
          b_vecnorm(A_wrench, dv2);
          b_obj[0] = obj->objfun.workspace.T_max[0];
          b_obj[3] = -obj->objfun.workspace.T_min[0];
          b_obj[1] = obj->objfun.workspace.T_max[1];
          b_obj[4] = -obj->objfun.workspace.T_min[1];
          b_obj[2] = obj->objfun.workspace.T_max[2];
          b_obj[5] = -obj->objfun.workspace.T_min[2];
          ubDiff = obj->objfun.workspace.center[0];
          for (i = 0; i < 6; i++) {
            b_obj[i] = (b_obj[i] - ((A_wrench[i] * ubDiff + A_wrench[i + 6] * d)
              + A_wrench[i + 12] * d1)) / dv2[i];
          }

          ubDiff = -minimum(b_obj);
          evalOK = ((!rtIsInf(ubDiff)) && (!rtIsNaN(ubDiff)));
          if (evalOK) {
            xk[idx] = temp;
          }

          obj->f_1 = ubDiff;
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
      gradf[idx] = (obj->f_1 - fCurrent) / b_deltaX;
      idx++;
    }
  }

  return evalOK;
}

/*
 * Arguments    : const i_struct_T *obj
 *                double workspace[880]
 *                const double f[22]
 *                const double x[22]
 * Return Type  : double
 */
static double computeFval_ReuseHx(const i_struct_T *obj, double workspace[880],
  const double f[22], const double x[22])
{
  double val;
  int k;
  val = 0.0;
  switch (obj->objtype) {
   case 5:
    val = obj->gammaScalar * x[obj->nvar - 1];
    break;

   case 3:
    {
      if (obj->hasLinear) {
        int ixlast;
        ixlast = obj->nvar;
        for (k = 0; k < ixlast; k++) {
          workspace[k] = 0.5 * obj->Hx[k] + f[k];
        }

        if (obj->nvar >= 1) {
          for (k = 0; k < ixlast; k++) {
            val += x[k] * workspace[k];
          }
        }
      } else {
        if (obj->nvar >= 1) {
          int ixlast;
          ixlast = obj->nvar;
          for (k = 0; k < ixlast; k++) {
            val += x[k] * obj->Hx[k];
          }
        }

        val *= 0.5;
      }
    }
    break;

   case 4:
    {
      if (obj->hasLinear) {
        int ixlast;
        ixlast = obj->nvar;
        if (ixlast - 1 >= 0) {
          memcpy(&workspace[0], &f[0], (unsigned int)ixlast * sizeof(double));
        }

        ixlast = 21 - obj->nvar;
        for (k = 0; k < ixlast; k++) {
          workspace[obj->nvar + k] = obj->rho;
        }

        for (k = 0; k < 21; k++) {
          double d;
          d = workspace[k] + 0.5 * obj->Hx[k];
          workspace[k] = d;
          val += x[k] * d;
        }
      } else {
        int ixlast;
        for (k = 0; k < 21; k++) {
          val += x[k] * obj->Hx[k];
        }

        val *= 0.5;
        ixlast = obj->nvar + 1;
        for (k = ixlast; k < 22; k++) {
          val += x[k - 1] * obj->rho;
        }
      }
    }
    break;
  }

  return val;
}

/*
 * Arguments    : double workspace[22]
 *                int nVar
 *                const double grad[22]
 *                const double AineqTrans[132]
 *                const double AeqTrans[66]
 *                const int finiteLB[22]
 *                int mLB
 *                const int finiteUB[22]
 *                const double lambda[40]
 * Return Type  : void
 */
static void computeGradLag(double workspace[22], int nVar, const double grad[22],
  const double AineqTrans[132], const double AeqTrans[66], const int finiteLB[22],
  int mLB, const int finiteUB[22], const double lambda[40])
{
  int i;
  int i1;
  int ia;
  int iac;
  int ix;
  i = (unsigned char)nVar;
  memcpy(&workspace[0], &grad[0], (unsigned int)i * sizeof(double));
  ix = 0;
  for (iac = 0; iac <= 44; iac += 22) {
    i = iac + nVar;
    for (ia = iac + 1; ia <= i; ia++) {
      i1 = (ia - iac) - 1;
      workspace[i1] += AeqTrans[ia - 1] * lambda[ix];
    }

    ix++;
  }

  ix = 3;
  for (iac = 0; iac <= 110; iac += 22) {
    i = iac + nVar;
    for (ia = iac + 1; ia <= i; ia++) {
      i1 = (ia - iac) - 1;
      workspace[i1] += AineqTrans[ia - 1] * lambda[ix];
    }

    ix++;
  }

  i = (unsigned char)mLB;
  for (ix = 0; ix < i; ix++) {
    i1 = finiteLB[ix];
    workspace[i1 - 1] -= lambda[ix + 9];
  }

  for (ix = 0; ix < 9; ix++) {
    i = finiteUB[ix];
    workspace[i - 1] += lambda[((unsigned char)mLB + ix) + 9];
  }
}

/*
 * Arguments    : i_struct_T *obj
 *                const double H[81]
 *                const double f[22]
 *                const double x[22]
 * Return Type  : void
 */
static void computeGrad_StoreHx(i_struct_T *obj, const double H[81], const
  double f[22], const double x[22])
{
  int ia;
  int iac;
  int ixlast;
  int lda;
  switch (obj->objtype) {
   case 5:
    {
      int i;
      i = obj->nvar;
      if (i - 2 >= 0) {
        memset(&obj->grad[0], 0, (unsigned int)(i - 1) * sizeof(double));
      }

      obj->grad[obj->nvar - 1] = obj->gammaScalar;
    }
    break;

   case 3:
    {
      int i;
      ixlast = obj->nvar - 1;
      lda = obj->nvar;
      if (obj->nvar != 0) {
        int ix;
        if (ixlast >= 0) {
          memset(&obj->Hx[0], 0, (unsigned int)(ixlast + 1) * sizeof(double));
        }

        ix = 0;
        i = obj->nvar * ixlast + 1;
        for (iac = 1; lda < 0 ? iac >= i : iac <= i; iac += lda) {
          int i1;
          i1 = iac + ixlast;
          for (ia = iac; ia <= i1; ia++) {
            int i2;
            i2 = ia - iac;
            obj->Hx[i2] += H[ia - 1] * x[ix];
          }

          ix++;
        }
      }

      i = obj->nvar;
      if (i - 1 >= 0) {
        memcpy(&obj->grad[0], &obj->Hx[0], (unsigned int)i * sizeof(double));
      }

      if (obj->hasLinear && (obj->nvar >= 1)) {
        ixlast = obj->nvar;
        for (lda = 0; lda < ixlast; lda++) {
          obj->grad[lda] += f[lda];
        }
      }
    }
    break;

   case 4:
    {
      int i;
      int i1;
      ixlast = obj->nvar - 1;
      lda = obj->nvar;
      if (obj->nvar != 0) {
        int ix;
        if (ixlast >= 0) {
          memset(&obj->Hx[0], 0, (unsigned int)(ixlast + 1) * sizeof(double));
        }

        ix = 0;
        i = obj->nvar * (obj->nvar - 1) + 1;
        for (iac = 1; lda < 0 ? iac >= i : iac <= i; iac += lda) {
          i1 = iac + ixlast;
          for (ia = iac; ia <= i1; ia++) {
            int i2;
            i2 = ia - iac;
            obj->Hx[i2] += H[ia - 1] * x[ix];
          }

          ix++;
        }
      }

      i = obj->nvar + 1;
      for (ixlast = i; ixlast < 22; ixlast++) {
        obj->Hx[ixlast - 1] = obj->beta * x[ixlast - 1];
      }

      memcpy(&obj->grad[0], &obj->Hx[0], 21U * sizeof(double));
      if (obj->hasLinear && (obj->nvar >= 1)) {
        ixlast = obj->nvar;
        for (lda = 0; lda < ixlast; lda++) {
          obj->grad[lda] += f[lda];
        }
      }

      if (21 - obj->nvar >= 1) {
        ixlast = obj->nvar;
        i = 20 - obj->nvar;
        for (lda = 0; lda <= i; lda++) {
          i1 = ixlast + lda;
          obj->grad[i1] += obj->rho;
        }
      }
    }
    break;
  }
}

/*
 * Arguments    : const double x[9]
 *                int nVar
 *                double workspaceIneq[6]
 *                const double AineqT[132]
 * Return Type  : void
 */
static void computeLinearResiduals(const double x[9], int nVar, double
  workspaceIneq[6], const double AineqT[132])
{
  int ia;
  int iac;
  int k;
  for (k = 0; k < 6; k++) {
    workspaceIneq[k] = -0.0;
  }

  for (iac = 0; iac <= 110; iac += 22) {
    double c;
    c = 0.0;
    k = iac + nVar;
    for (ia = iac + 1; ia <= k; ia++) {
      c += AineqT[ia - 1] * x[(ia - iac) - 1];
    }

    k = div_nde_s32_floor(iac, 22);
    workspaceIneq[k] += c;
  }
}

/*
 * Arguments    : g_struct_T *obj
 *                int nrows
 * Return Type  : void
 */
static void computeQ_(g_struct_T *obj, int nrows)
{
  double work[40];
  int i;
  int iQR0;
  int ia;
  int idx;
  int k;
  int m;
  int n;
  k = obj->minRowCol;
  for (idx = 0; idx < k; idx++) {
    iQR0 = 40 * idx + idx;
    n = obj->mrows - idx;
    if (n - 2 >= 0) {
      memcpy(&obj->Q[iQR0 + 1], &obj->QR[iQR0 + 1], (unsigned int)(n - 1) *
             sizeof(double));
    }
  }

  m = obj->mrows;
  if (nrows >= 1) {
    int itau;
    for (idx = k; idx < nrows; idx++) {
      ia = idx * 40;
      memset(&obj->Q[ia], 0, (unsigned int)m * sizeof(double));
      obj->Q[ia + idx] = 1.0;
    }

    itau = obj->minRowCol - 1;
    memset(&work[0], 0, 40U * sizeof(double));
    for (i = obj->minRowCol; i >= 1; i--) {
      int b_i;
      int iaii;
      iaii = i + (i - 1) * 40;
      if (i < nrows) {
        int lastc;
        int lastv;
        obj->Q[iaii - 1] = 1.0;
        idx = iaii + 40;
        if (obj->tau[itau] != 0.0) {
          bool exitg2;
          lastv = m - i;
          iQR0 = (iaii + m) - i;
          while ((lastv + 1 > 0) && (obj->Q[iQR0 - 1] == 0.0)) {
            lastv--;
            iQR0--;
          }

          lastc = (nrows - i) - 1;
          exitg2 = false;
          while ((!exitg2) && (lastc + 1 > 0)) {
            int exitg1;
            iQR0 = (iaii + lastc * 40) + 40;
            ia = iQR0;
            do {
              exitg1 = 0;
              if (ia <= iQR0 + lastv) {
                if (obj->Q[ia - 1] != 0.0) {
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
          lastv = -1;
          lastc = -1;
        }

        if (lastv + 1 > 0) {
          double c;
          if (lastc + 1 != 0) {
            if (lastc >= 0) {
              memset(&work[0], 0, (unsigned int)(lastc + 1) * sizeof(double));
            }

            b_i = (iaii + 40 * lastc) + 40;
            for (n = idx; n <= b_i; n += 40) {
              c = 0.0;
              k = n + lastv;
              for (ia = n; ia <= k; ia++) {
                c += obj->Q[ia - 1] * obj->Q[((iaii + ia) - n) - 1];
              }

              iQR0 = div_nde_s32_floor((n - iaii) - 40, 40);
              work[iQR0] += c;
            }
          }

          if (!(-obj->tau[itau] == 0.0)) {
            iQR0 = iaii;
            for (idx = 0; idx <= lastc; idx++) {
              c = work[idx];
              if (c != 0.0) {
                c *= -obj->tau[itau];
                b_i = iQR0 + 40;
                k = lastv + iQR0;
                for (n = b_i; n <= k + 40; n++) {
                  obj->Q[n - 1] += obj->Q[((iaii + n) - iQR0) - 41] * c;
                }
              }

              iQR0 += 40;
            }
          }
        }
      }

      if (i < m) {
        iQR0 = iaii + 1;
        b_i = (iaii + m) - i;
        for (k = iQR0; k <= b_i; k++) {
          obj->Q[k - 1] *= -obj->tau[itau];
        }
      }

      obj->Q[iaii - 1] = 1.0 - obj->tau[itau];
      b_i = (unsigned char)(i - 1);
      for (idx = 0; idx < b_i; idx++) {
        obj->Q[(iaii - idx) - 2] = 0.0;
      }

      itau--;
    }
  }
}

/*
 * Arguments    : const double H[81]
 *                c_struct_T *solution
 *                struct_T *memspace
 *                const g_struct_T *qrmanager
 *                h_struct_T *cholmanager
 *                const i_struct_T *objective
 *                bool alwaysPositiveDef
 * Return Type  : void
 */
static void compute_deltax(const double H[81], c_struct_T *solution, struct_T
  *memspace, const g_struct_T *qrmanager, h_struct_T *cholmanager, const
  i_struct_T *objective, bool alwaysPositiveDef)
{
  int b_i;
  int idx;
  int ix;
  int jA;
  int jjA;
  int mNull_tmp;
  int nVar_tmp;
  nVar_tmp = qrmanager->mrows - 1;
  mNull_tmp = qrmanager->mrows - qrmanager->ncols;
  if (mNull_tmp <= 0) {
    if (nVar_tmp >= 0) {
      memset(&solution->searchDir[0], 0, (unsigned int)(nVar_tmp + 1) * sizeof
             (double));
    }
  } else {
    for (idx = 0; idx <= nVar_tmp; idx++) {
      solution->searchDir[idx] = -objective->grad[idx];
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
              jjA = (nVar_tmp + 1) * idx;
              jA = 40 * idx;
              for (ix = 0; ix <= nVar_tmp; ix++) {
                cholmanager->FMat[jA + ix] = H[jjA + ix];
              }
            }

            cholmanager->info = xpotrf(qrmanager->mrows, cholmanager->FMat);
          } else {
            cholmanager->ndims = qrmanager->mrows;
            for (idx = 0; idx <= nVar_tmp; idx++) {
              jjA = qrmanager->mrows * idx;
              jA = 40 * idx;
              for (ix = 0; ix <= nVar_tmp; ix++) {
                cholmanager->FMat[jA + ix] = H[jjA + ix];
              }
            }

            fullColLDL2_(cholmanager, qrmanager->mrows);
            if (cholmanager->ConvexCheck) {
              idx = 0;
              int exitg1;
              do {
                exitg1 = 0;
                if (idx <= nVar_tmp) {
                  if (cholmanager->FMat[idx + 40 * idx] <= 0.0) {
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
            solve(cholmanager, solution->searchDir);
          } else {
            int i;
            int nVars;
            nVars = cholmanager->ndims - 2;
            if (cholmanager->ndims != 0) {
              for (idx = 0; idx <= nVars + 1; idx++) {
                jjA = idx + idx * 40;
                i = nVars - idx;
                for (b_i = 0; b_i <= i; b_i++) {
                  ix = (idx + b_i) + 1;
                  solution->searchDir[ix] -= solution->searchDir[idx] *
                    cholmanager->FMat[(jjA + b_i) + 1];
                }
              }
            }

            nVars = cholmanager->ndims;
            for (idx = 0; idx < nVars; idx++) {
              solution->searchDir[idx] /= cholmanager->FMat[idx + 40 * idx];
            }

            if (cholmanager->ndims != 0) {
              for (idx = nVars; idx >= 1; idx--) {
                double smax;
                jA = (idx - 1) * 40;
                smax = solution->searchDir[idx - 1];
                i = idx + 1;
                for (b_i = nVars; b_i >= i; b_i--) {
                  smax -= cholmanager->FMat[(jA + b_i) - 1] *
                    solution->searchDir[b_i - 1];
                }

                solution->searchDir[idx - 1] = smax;
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
              jjA = nVars * idx;
              jA = 40 * idx;
              for (ix = 0; ix < nVars; ix++) {
                cholmanager->FMat[jA + ix] = H[jjA + ix];
              }
            }

            cholmanager->info = xpotrf(objective->nvar, cholmanager->FMat);
            if (cholmanager->info != 0) {
              solution->state = -6;
            } else {
              double smax;
              int i;
              solve(cholmanager, solution->searchDir);
              smax = 1.0 / objective->beta;
              jjA = objective->nvar + 1;
              i = qrmanager->mrows;
              for (ix = jjA; ix <= i; ix++) {
                solution->searchDir[ix - 1] *= smax;
              }
            }
          }
        }
        break;
      }
    } else {
      int nullStartIdx_tmp;
      nullStartIdx_tmp = 40 * qrmanager->ncols + 1;
      if (objective->objtype == 5) {
        for (idx = 0; idx < mNull_tmp; idx++) {
          memspace->workspace_float[idx] = -qrmanager->Q[nVar_tmp + 40 *
            (qrmanager->ncols + idx)];
        }

        if (qrmanager->mrows != 0) {
          int i;
          memset(&solution->searchDir[0], 0, (unsigned int)(nVar_tmp + 1) *
                 sizeof(double));
          ix = 0;
          i = nullStartIdx_tmp + 40 * (mNull_tmp - 1);
          for (jjA = nullStartIdx_tmp; jjA <= i; jjA += 40) {
            int nVars;
            nVars = jjA + nVar_tmp;
            for (idx = jjA; idx <= nVars; idx++) {
              jA = idx - jjA;
              solution->searchDir[jA] += qrmanager->Q[idx - 1] *
                memspace->workspace_float[ix];
            }

            ix++;
          }
        }
      } else {
        int i;
        int nVars;
        if (objective->objtype == 3) {
          xgemm(qrmanager->mrows, mNull_tmp, qrmanager->mrows, H,
                qrmanager->mrows, qrmanager->Q, nullStartIdx_tmp,
                memspace->workspace_float);
          b_xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows, qrmanager->Q,
                  nullStartIdx_tmp, memspace->workspace_float, cholmanager->FMat);
        } else if (alwaysPositiveDef) {
          nVars = qrmanager->mrows;
          xgemm(objective->nvar, mNull_tmp, objective->nvar, H, objective->nvar,
                qrmanager->Q, nullStartIdx_tmp, memspace->workspace_float);
          i = objective->nvar + 1;
          for (jA = 0; jA < mNull_tmp; jA++) {
            for (jjA = i; jjA <= nVars; jjA++) {
              memspace->workspace_float[(jjA + 40 * jA) - 1] = objective->beta *
                qrmanager->Q[(jjA + 40 * (jA + qrmanager->ncols)) - 1];
            }
          }

          b_xgemm(mNull_tmp, mNull_tmp, qrmanager->mrows, qrmanager->Q,
                  nullStartIdx_tmp, memspace->workspace_float, cholmanager->FMat);
        }

        if (alwaysPositiveDef) {
          cholmanager->ndims = mNull_tmp;
          cholmanager->info = xpotrf(mNull_tmp, cholmanager->FMat);
        } else {
          cholmanager->ndims = mNull_tmp;
          fullColLDL2_(cholmanager, mNull_tmp);
          if (cholmanager->ConvexCheck) {
            idx = 0;
            int exitg1;
            do {
              exitg1 = 0;
              if (idx <= mNull_tmp - 1) {
                if (cholmanager->FMat[idx + 40 * idx] <= 0.0) {
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
          if (qrmanager->mrows != 0) {
            memset(&memspace->workspace_float[0], 0, (unsigned int)mNull_tmp *
                   sizeof(double));
            i = nullStartIdx_tmp + 40 * (mNull_tmp - 1);
            for (jjA = nullStartIdx_tmp; jjA <= i; jjA += 40) {
              smax = 0.0;
              nVars = jjA + nVar_tmp;
              for (idx = jjA; idx <= nVars; idx++) {
                smax += qrmanager->Q[idx - 1] * objective->grad[idx - jjA];
              }

              nVars = div_nde_s32_floor(jjA - nullStartIdx_tmp, 40);
              memspace->workspace_float[nVars] -= smax;
            }
          }

          if (alwaysPositiveDef) {
            nVars = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (idx = 0; idx < nVars; idx++) {
                jA = idx * 40;
                smax = memspace->workspace_float[idx];
                for (b_i = 0; b_i < idx; b_i++) {
                  smax -= cholmanager->FMat[jA + b_i] *
                    memspace->workspace_float[b_i];
                }

                memspace->workspace_float[idx] = smax / cholmanager->FMat[jA +
                  idx];
              }
            }

            if (cholmanager->ndims != 0) {
              for (idx = nVars; idx >= 1; idx--) {
                jjA = (idx + (idx - 1) * 40) - 1;
                memspace->workspace_float[idx - 1] /= cholmanager->FMat[jjA];
                for (b_i = 0; b_i <= idx - 2; b_i++) {
                  ix = (idx - b_i) - 2;
                  memspace->workspace_float[ix] -= memspace->workspace_float[idx
                    - 1] * cholmanager->FMat[(jjA - b_i) - 1];
                }
              }
            }
          } else {
            nVars = cholmanager->ndims - 2;
            if (cholmanager->ndims != 0) {
              for (idx = 0; idx <= nVars + 1; idx++) {
                jjA = idx + idx * 40;
                i = nVars - idx;
                for (b_i = 0; b_i <= i; b_i++) {
                  ix = (idx + b_i) + 1;
                  memspace->workspace_float[ix] -= memspace->workspace_float[idx]
                    * cholmanager->FMat[(jjA + b_i) + 1];
                }
              }
            }

            nVars = cholmanager->ndims;
            for (idx = 0; idx < nVars; idx++) {
              memspace->workspace_float[idx] /= cholmanager->FMat[idx + 40 * idx];
            }

            if (cholmanager->ndims != 0) {
              for (idx = nVars; idx >= 1; idx--) {
                jA = (idx - 1) * 40;
                smax = memspace->workspace_float[idx - 1];
                i = idx + 1;
                for (b_i = nVars; b_i >= i; b_i--) {
                  smax -= cholmanager->FMat[(jA + b_i) - 1] *
                    memspace->workspace_float[b_i - 1];
                }

                memspace->workspace_float[idx - 1] = smax;
              }
            }
          }

          if (qrmanager->mrows != 0) {
            memset(&solution->searchDir[0], 0, (unsigned int)(nVar_tmp + 1) *
                   sizeof(double));
            ix = 0;
            i = nullStartIdx_tmp + 40 * (mNull_tmp - 1);
            for (jjA = nullStartIdx_tmp; jjA <= i; jjA += 40) {
              nVars = jjA + nVar_tmp;
              for (idx = jjA; idx <= nVars; idx++) {
                jA = idx - jjA;
                solution->searchDir[jA] += qrmanager->Q[idx - 1] *
                  memspace->workspace_float[ix];
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
 * Arguments    : double workspace[880]
 *                c_struct_T *solution
 *                const i_struct_T *objective
 *                const g_struct_T *qrmanager
 * Return Type  : void
 */
static void compute_lambda(double workspace[880], c_struct_T *solution, const
  i_struct_T *objective, const g_struct_T *qrmanager)
{
  int ia;
  int iac;
  int idx;
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
          idxQR = qrmanager->mrows + 40 * (qrmanager->ncols - 1);
          while ((idx > qrmanager->mrows) && (fabs(qrmanager->QR[idxQR - 1]) >=
                  c)) {
            idx--;
            idxQR -= 40;
          }

          nonDegenerate = (idx == qrmanager->mrows);
          if (nonDegenerate) {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }

        if (guard2) {
          idxQR = idx + 40 * (idx - 1);
          while ((idx >= 1) && (fabs(qrmanager->QR[idxQR - 1]) >= c)) {
            idx--;
            idxQR -= 41;
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
      if (qrmanager->mrows != 0) {
        memset(&workspace[0], 0, (unsigned int)nActiveConstr_tmp * sizeof(double));
        idx = 40 * (qrmanager->ncols - 1) + 1;
        for (iac = 1; iac <= idx; iac += 40) {
          c = 0.0;
          idxQR = iac + qrmanager->mrows;
          for (ia = iac; ia < idxQR; ia++) {
            c += qrmanager->Q[ia - 1] * objective->grad[ia - iac];
          }

          idxQR = div_nde_s32_floor(iac - 1, 40);
          workspace[idxQR] += c;
        }
      }

      for (iac = nActiveConstr_tmp; iac >= 1; iac--) {
        idx = (iac + (iac - 1) * 40) - 1;
        workspace[iac - 1] /= qrmanager->QR[idx];
        for (ia = 0; ia <= iac - 2; ia++) {
          idxQR = (iac - ia) - 2;
          workspace[idxQR] -= workspace[iac - 1] * qrmanager->QR[(idx - ia) - 1];
        }
      }

      for (idx = 0; idx < nActiveConstr_tmp; idx++) {
        solution->lambda[idx] = -workspace[idx];
      }
    }
  }
}

/*
 * Arguments    : int x[40]
 *                int xLen
 *                int workspace[40]
 *                int xMin
 *                int xMax
 * Return Type  : void
 */
static void countsort(int x[40], int xLen, int workspace[40], int xMin, int xMax)
{
  int idx;
  int idxFill;
  if ((xLen > 1) && (xMax > xMin)) {
    int i;
    int idxEnd;
    int idxStart;
    i = xMax - xMin;
    if (i >= 0) {
      memset(&workspace[0], 0, (unsigned int)(i + 1) * sizeof(int));
    }

    for (idx = 0; idx < xLen; idx++) {
      idxStart = x[idx] - xMin;
      workspace[idxStart]++;
    }

    for (idx = 2; idx <= i + 1; idx++) {
      workspace[idx - 1] += workspace[idx - 2];
    }

    idxStart = 1;
    idxEnd = workspace[0];
    for (idx = 0; idx < i; idx++) {
      for (idxFill = idxStart; idxFill <= idxEnd; idxFill++) {
        x[idxFill - 1] = idx + xMin;
      }

      idxStart = workspace[idx] + 1;
      idxEnd = workspace[idx + 1];
    }

    for (idx = idxStart; idx <= idxEnd; idx++) {
      x[idx - 1] = xMax;
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                const double x[22]
 * Return Type  : double
 */
static double d_maxConstraintViolation_AMats_(e_struct_T *obj, const double x[22])
{
  double c;
  double v;
  int ia;
  int iac;
  int k;
  v = 0.0;
  for (k = 0; k < 6; k++) {
    obj->maxConstrWorkspace[k] = -obj->bineq[k];
  }

  for (iac = 0; iac <= 110; iac += 22) {
    c = 0.0;
    k = iac + obj->nVar;
    for (ia = iac + 1; ia <= k; ia++) {
      c += obj->Aineq[ia - 1] * x[(ia - iac) - 1];
    }

    k = div_nde_s32_floor(iac, 22);
    obj->maxConstrWorkspace[k] += c;
  }

  for (k = 0; k < 6; k++) {
    v = fmax(v, obj->maxConstrWorkspace[k]);
  }

  obj->maxConstrWorkspace[0] = -obj->beq[0];
  obj->maxConstrWorkspace[1] = -obj->beq[1];
  obj->maxConstrWorkspace[2] = -obj->beq[2];
  for (iac = 0; iac <= 44; iac += 22) {
    c = 0.0;
    k = iac + obj->nVar;
    for (ia = iac + 1; ia <= k; ia++) {
      c += obj->Aeq[ia - 1] * x[(ia - iac) - 1];
    }

    k = div_nde_s32_floor(iac, 22);
    obj->maxConstrWorkspace[k] += c;
  }

  v = fmax(v, fabs(obj->maxConstrWorkspace[0]));
  v = fmax(v, fabs(obj->maxConstrWorkspace[1]));
  return fmax(v, fabs(obj->maxConstrWorkspace[2]));
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
 * Arguments    : g_struct_T *obj
 *                int idx
 * Return Type  : void
 */
static void deleteColMoveEnd(g_struct_T *obj, int idx)
{
  double s;
  double temp_tmp;
  int b_k;
  int i;
  int k;
  if (obj->usedPivoting) {
    i = 1;
    while ((i <= obj->ncols) && (obj->jpvt[i - 1] != idx)) {
      i++;
    }

    idx = i;
  }

  if (idx >= obj->ncols) {
    obj->ncols--;
  } else {
    int b_i;
    int u0;
    b_i = obj->ncols - 1;
    obj->jpvt[idx - 1] = obj->jpvt[b_i];
    i = obj->minRowCol;
    for (k = 0; k < i; k++) {
      obj->QR[k + 40 * (idx - 1)] = obj->QR[k + 40 * b_i];
    }

    obj->ncols = b_i;
    u0 = obj->mrows;
    i = obj->ncols;
    if (u0 <= i) {
      i = u0;
    }

    obj->minRowCol = i;
    if (idx < obj->mrows) {
      double c;
      double temp;
      int QRk0;
      int b_temp_tmp;
      int endIdx;
      int n;
      u0 = obj->mrows - 1;
      endIdx = obj->ncols;
      if (u0 <= endIdx) {
        endIdx = u0;
      }

      k = endIdx;
      i = 40 * (idx - 1);
      while (k >= idx) {
        b_i = k + i;
        temp_tmp = obj->QR[b_i];
        c = xrotg(&obj->QR[b_i - 1], &temp_tmp, &s);
        obj->QR[b_i] = temp_tmp;
        b_i = 40 * (k - 1);
        obj->QR[k + b_i] = 0.0;
        QRk0 = k + 40 * idx;
        n = obj->ncols - idx;
        if (n >= 1) {
          for (b_k = 0; b_k < n; b_k++) {
            b_temp_tmp = QRk0 + b_k * 40;
            temp_tmp = obj->QR[b_temp_tmp - 1];
            temp = c * temp_tmp + s * obj->QR[b_temp_tmp];
            obj->QR[b_temp_tmp] = c * obj->QR[b_temp_tmp] - s * temp_tmp;
            obj->QR[b_temp_tmp - 1] = temp;
          }
        }

        n = obj->mrows;
        for (b_k = 0; b_k < n; b_k++) {
          b_temp_tmp = b_i + b_k;
          temp_tmp = obj->Q[b_temp_tmp + 40];
          temp = c * obj->Q[b_temp_tmp] + s * temp_tmp;
          obj->Q[b_temp_tmp + 40] = c * temp_tmp - s * obj->Q[b_temp_tmp];
          obj->Q[b_temp_tmp] = temp;
        }

        k--;
      }

      b_i = idx + 1;
      for (k = b_i; k <= endIdx; k++) {
        u0 = 40 * (k - 1);
        i = k + u0;
        temp_tmp = obj->QR[i];
        c = xrotg(&obj->QR[i - 1], &temp_tmp, &s);
        obj->QR[i] = temp_tmp;
        QRk0 = k * 41;
        n = obj->ncols - k;
        if (n >= 1) {
          for (b_k = 0; b_k < n; b_k++) {
            b_temp_tmp = QRk0 + b_k * 40;
            temp_tmp = obj->QR[b_temp_tmp - 1];
            temp = c * temp_tmp + s * obj->QR[b_temp_tmp];
            obj->QR[b_temp_tmp] = c * obj->QR[b_temp_tmp] - s * temp_tmp;
            obj->QR[b_temp_tmp - 1] = temp;
          }
        }

        n = obj->mrows;
        for (b_k = 0; b_k < n; b_k++) {
          b_temp_tmp = u0 + b_k;
          temp_tmp = obj->Q[b_temp_tmp + 40];
          temp = c * obj->Q[b_temp_tmp] + s * temp_tmp;
          obj->Q[b_temp_tmp + 40] = c * temp_tmp - s * obj->Q[b_temp_tmp];
          obj->Q[b_temp_tmp] = temp;
        }
      }
    }
  }
}

/*
 * Arguments    : int numerator
 *                int denominator
 * Return Type  : int
 */
static int div_nde_s32_floor(int numerator, int denominator)
{
  int quotient;
  if (((numerator < 0) != (denominator < 0)) && (numerator % denominator != 0))
  {
    quotient = -1;
  } else {
    quotient = 0;
  }

  quotient += numerator / denominator;
  return quotient;
}

/*
 * Arguments    : const double H[81]
 *                const double f[22]
 *                c_struct_T *solution
 *                struct_T *memspace
 *                e_struct_T *workingset
 *                g_struct_T *qrmanager
 *                h_struct_T *cholmanager
 *                i_struct_T *objective
 *                const k_struct_T *options
 *                const k_struct_T *runTimeOptions
 * Return Type  : void
 */
static void driver(const double H[81], const double f[22], c_struct_T *solution,
                   struct_T *memspace, e_struct_T *workingset, g_struct_T
                   *qrmanager, h_struct_T *cholmanager, i_struct_T *objective,
                   const k_struct_T *options, const k_struct_T *runTimeOptions)
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
    i = (unsigned char)workingset->sizes[3];
    for (idx = 0; idx < i; idx++) {
      if (workingset->isActiveConstr[(workingset->isActiveIdx[3] + idx) - 1]) {
        solution->xstar[workingset->indexLB[idx] - 1] = -workingset->
          lb[workingset->indexLB[idx] - 1];
      }
    }

    for (idx = 0; idx < 9; idx++) {
      if (workingset->isActiveConstr[(workingset->isActiveIdx[4] + idx) - 1]) {
        solution->xstar[workingset->indexUB[idx] - 1] = workingset->
          ub[workingset->indexUB[idx] - 1];
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
    solution->maxConstr = b_maxConstraintViolation(workingset, solution->xstar);
    if (solution->maxConstr > 1.0E-6) {
      phaseone(H, f, solution, memspace, workingset, qrmanager, cholmanager,
               objective, options->SolverName, runTimeOptions);
      if (solution->state != 0) {
        solution->maxConstr = b_maxConstraintViolation(workingset,
          solution->xstar);
        if (solution->maxConstr > 1.0E-6) {
          memset(&solution->lambda[0], 0, 40U * sizeof(double));
          d = 0.0;
          switch (objective->objtype) {
           case 5:
            d = objective->gammaScalar * solution->xstar[objective->nvar - 1];
            break;

           case 3:
            linearForm_(objective->hasLinear, objective->nvar,
                        memspace->workspace_float, H, f, solution->xstar);
            if (objective->nvar >= 1) {
              ret = objective->nvar;
              for (idx_local = 0; idx_local < ret; idx_local++) {
                d += solution->xstar[idx_local] * memspace->
                  workspace_float[idx_local];
              }
            }
            break;

           case 4:
            linearForm_(objective->hasLinear, objective->nvar,
                        memspace->workspace_float, H, f, solution->xstar);
            i = objective->nvar + 1;
            for (idx = i; idx < 22; idx++) {
              memspace->workspace_float[idx - 1] = 0.5 * objective->beta *
                solution->xstar[idx - 1] + objective->rho;
            }

            for (idx_local = 0; idx_local < 21; idx_local++) {
              d += solution->xstar[idx_local] * memspace->
                workspace_float[idx_local];
            }
            break;
          }

          solution->fstar = d;
          solution->state = -2;
        } else {
          if (solution->maxConstr > 0.0) {
            if (ret - 1 >= 0) {
              memcpy(&solution->searchDir[0], &solution->xstar[0], (unsigned int)
                     ret * sizeof(double));
            }

            PresolveWorkingSet(solution, memspace, workingset, qrmanager);
            minLambda = b_maxConstraintViolation(workingset, solution->xstar);
            if (minLambda >= solution->maxConstr) {
              solution->maxConstr = minLambda;
              if (ret - 1 >= 0) {
                memcpy(&solution->xstar[0], &solution->searchDir[0], (unsigned
                        int)ret * sizeof(double));
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
    computeGrad_StoreHx(objective, H, f, solution->xstar);
    solution->fstar = computeFval_ReuseHx(objective, memspace->workspace_float,
      f, solution->xstar);
    if (solution->iterations < runTimeOptions->MaxIterations) {
      solution->state = -5;
    } else {
      solution->state = 0;
    }

    memset(&solution->lambda[0], 0, 40U * sizeof(double));
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
            squareQ_appendCol(qrmanager, workingset->ATwset, 22 *
                              (workingset->nActiveConstr - 1) + 1);
            break;

           case -1:
            deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
            break;

           default:
            factorQR(qrmanager, workingset->ATwset, i, workingset->nActiveConstr);
            computeQ_(qrmanager, qrmanager->mrows);
            break;
          }

          ret = memcmp(&options->SolverName[0], &b[0], 7);
          compute_deltax(H, solution, memspace, qrmanager, cholmanager,
                         objective, ret == 0);
          if (solution->state != -5) {
            exitg1 = 1;
          } else if ((c_xnrm2(i, solution->searchDir) < 1.0E-6) ||
                     (workingset->nActiveConstr >= i)) {
            guard4 = true;
          } else {
            minLambda = feasibleratiotest(solution->xstar, solution->searchDir,
              memspace->workspace_float, workingset->nVar, workingset->Aineq,
              workingset->bineq, workingset->lb, workingset->ub,
              workingset->indexLB, workingset->indexUB, workingset->sizes,
              workingset->isActiveIdx, workingset->isActiveConstr,
              workingset->nWConstr, TYPE == 5, &updateFval, &i1, &idx_local);
            if (updateFval) {
              switch (i1) {
               case 3:
                workingset->nWConstr[2]++;
                workingset->isActiveConstr[(workingset->isActiveIdx[2] +
                  idx_local) - 2] = true;
                workingset->nActiveConstr++;
                i1 = workingset->nActiveConstr - 1;
                workingset->Wid[i1] = 3;
                workingset->Wlocalidx[i1] = idx_local;
                ret = 22 * (idx_local - 1);
                activeSetChangeID = 22 * i1;
                i2 = workingset->nVar;
                for (idx = 0; idx < i2; idx++) {
                  workingset->ATwset[activeSetChangeID + idx] =
                    workingset->Aineq[ret + idx];
                }

                workingset->bwset[i1] = workingset->bineq[idx_local - 1];
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
                if (c_xnrm2(objective->nvar, solution->searchDir) > 100.0 *
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
                solution->xstar[idx_local] += minLambda * solution->
                  searchDir[idx_local];
              }
            }

            computeGrad_StoreHx(objective, H, f, solution->xstar);
            guard3 = true;
          }
        } else {
          i1 = (unsigned char)i;
          memset(&solution->searchDir[0], 0, (unsigned int)i1 * sizeof(double));
          guard4 = true;
        }

        if (guard4) {
          compute_lambda(memspace->workspace_float, solution, objective,
                         qrmanager);
          if ((solution->state != -7) || (workingset->nActiveConstr > i)) {
            ret = 0;
            minLambda = 0.0;
            i1 = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
            i2 = workingset->nActiveConstr;
            for (idx = i1; idx <= i2; idx++) {
              d = solution->lambda[idx - 1];
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
                solution->lambda[ret - 1] = solution->lambda
                  [workingset->nActiveConstr];
              }

              solution->lambda[workingset->nActiveConstr] = 0.0;
            }
          } else {
            ret = workingset->nActiveConstr;
            activeSetChangeID = 0;
            globalActiveConstrIdx = workingset->nActiveConstr;
            subProblemChanged = true;
            removeConstr(workingset, workingset->nActiveConstr);
            solution->lambda[ret - 1] = 0.0;
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
              solution->xstar);
            minLambda = solution->maxConstr;
            if (objective->objtype == 5) {
              minLambda = solution->maxConstr - solution->xstar[objective->nvar
                - 1];
            }

            if (minLambda > 1.0E-6) {
              if (ret - 1 >= 0) {
                memcpy(&solution->searchDir[0], &solution->xstar[0], (unsigned
                        int)ret * sizeof(double));
              }

              updateFval = feasibleX0ForWorkingSet(memspace->workspace_float,
                solution->searchDir, workingset, qrmanager);
              if ((!updateFval) && (solution->state != 0)) {
                solution->state = -2;
              }

              activeSetChangeID = 0;
              minLambda = b_maxConstraintViolation(workingset,
                solution->searchDir);
              if (minLambda < solution->maxConstr) {
                if (ret - 1 >= 0) {
                  memcpy(&solution->xstar[0], &solution->searchDir[0], (unsigned
                          int)ret * sizeof(double));
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
            memspace->workspace_float, f, solution->xstar);
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }
}

/*
 * Arguments    : const b_struct_T *c_obj_next_next_next_next_next_
 *                const double x[9]
 *                double Ceq_workspace[3]
 *                int *status
 * Return Type  : double
 */
static double evalObjAndConstr(const b_struct_T *c_obj_next_next_next_next_next_,
  const double x[9], double Ceq_workspace[3], int *status)
{
  double A_wrench[18];
  double b_x[9];
  double b_y[9];
  double d_obj_next_next_next_next_next_[6];
  double dv[6];
  double y[3];
  double d;
  double d1;
  double fval;
  int x_tmp;
  int xi;
  bool b;

  /* 'mso_forces:62' @(W) objective(W, center, N, T_max, T_min) */
  /*  Solver */
  /*  目标函数 */
  /* 'mso_forces:16' W = reshape(W, [3, length(W) / 3]); */
  /* 'mso_forces:17' W = W ./ vecnorm(W, 2, 1); */
  /* 'mso_forces:19' A_wrench = [eye(N); -eye(N)] * pinv(W); */
  vecnorm(x, y);
  for (xi = 0; xi < 3; xi++) {
    d = y[xi];
    b_x[3 * xi] = x[3 * xi] / d;
    x_tmp = 3 * xi + 1;
    b_x[x_tmp] = x[x_tmp] / d;
    x_tmp = 3 * xi + 2;
    b_x[x_tmp] = x[x_tmp] / d;
  }

  pinv(b_x, b_y);
  for (xi = 0; xi < 6; xi++) {
    signed char i;
    signed char i1;
    signed char i2;
    i = iv1[xi];
    i1 = iv1[xi + 6];
    i2 = iv1[xi + 12];
    for (x_tmp = 0; x_tmp < 3; x_tmp++) {
      A_wrench[xi + 6 * x_tmp] = ((double)i * b_y[3 * x_tmp] + (double)i1 * b_y
        [3 * x_tmp + 1]) + (double)i2 * b_y[3 * x_tmp + 2];
    }
  }

  /* 'mso_forces:20' b_wrench = [T_max; -T_min]; */
  /* 'mso_forces:22' distances = (b_wrench - A_wrench * p) ./ vecnorm(A_wrench, 2, 2); */
  /*  返回最小距离（负号用于最大化） */
  /* 'mso_forces:25' obj_val = -min(distances); */
  b_vecnorm(A_wrench, dv);
  d_obj_next_next_next_next_next_[0] = c_obj_next_next_next_next_next_->T_max[0];
  d_obj_next_next_next_next_next_[3] = -c_obj_next_next_next_next_next_->T_min[0];
  d_obj_next_next_next_next_next_[1] = c_obj_next_next_next_next_next_->T_max[1];
  d_obj_next_next_next_next_next_[4] = -c_obj_next_next_next_next_next_->T_min[1];
  d_obj_next_next_next_next_next_[2] = c_obj_next_next_next_next_next_->T_max[2];
  d_obj_next_next_next_next_next_[5] = -c_obj_next_next_next_next_next_->T_min[2];
  d = c_obj_next_next_next_next_next_->center[0];
  fval = c_obj_next_next_next_next_next_->center[1];
  d1 = c_obj_next_next_next_next_next_->center[2];
  for (xi = 0; xi < 6; xi++) {
    d_obj_next_next_next_next_next_[xi] = (d_obj_next_next_next_next_next_[xi] -
      ((A_wrench[xi] * d + A_wrench[xi + 6] * fval) + A_wrench[xi + 12] * d1)) /
      dv[xi];
  }

  fval = -minimum(d_obj_next_next_next_next_next_);
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
    /* 'mso_forces:84' @(W) constraint(W) */
    /*  非线性约束：单位向量约束 */
    /* 'mso_forces:29' c = []; */
    /* 'mso_forces:30' ceq = sum(reshape(W, 3, length(W) / 3).^2, 1) - 1; */
    for (x_tmp = 0; x_tmp < 9; x_tmp++) {
      d = x[x_tmp];
      b_y[x_tmp] = d * d;
    }

    /* 'mso_forces:32' if nargout > 2 */
    for (xi = 0; xi < 3; xi++) {
      x_tmp = xi * 3;
      d = (b_y[x_tmp] + b_y[x_tmp + 1]) + b_y[x_tmp + 2];
      Ceq_workspace[xi] = d - 1.0;
      y[xi] = d - 1.0;
    }

    *status = checkVectorNonFinite(y);
  }

  return fval;
}

/*
 * Arguments    : const b_struct_T *c_obj_next_next_next_next_next_
 *                const double x[9]
 *                double Ceq_workspace[3]
 *                double JacEqTrans_workspace[66]
 *                int *status
 * Return Type  : double
 */
static double evalObjAndConstrAndDerivatives(const b_struct_T
  *c_obj_next_next_next_next_next_, const double x[9], double Ceq_workspace[3],
  double JacEqTrans_workspace[66], int *status)
{
  double A_wrench[18];
  double b_x[9];
  double dv[9];
  double d_obj_next_next_next_next_next_[6];
  double dv1[6];
  double Ceq_tmp[3];
  double d;
  double d1;
  double fval;
  int col;
  int row;
  bool allFinite;

  /* 'mso_forces:62' @(W) objective(W, center, N, T_max, T_min) */
  /*  Solver */
  /*  目标函数 */
  /* 'mso_forces:16' W = reshape(W, [3, length(W) / 3]); */
  /* 'mso_forces:17' W = W ./ vecnorm(W, 2, 1); */
  /* 'mso_forces:19' A_wrench = [eye(N); -eye(N)] * pinv(W); */
  vecnorm(x, Ceq_tmp);
  for (col = 0; col < 3; col++) {
    fval = Ceq_tmp[col];
    b_x[3 * col] = x[3 * col] / fval;
    row = 3 * col + 1;
    b_x[row] = x[row] / fval;
    row = 3 * col + 2;
    b_x[row] = x[row] / fval;
  }

  pinv(b_x, dv);
  for (col = 0; col < 6; col++) {
    signed char i;
    signed char i1;
    signed char i2;
    i = iv1[col];
    i1 = iv1[col + 6];
    i2 = iv1[col + 12];
    for (row = 0; row < 3; row++) {
      A_wrench[col + 6 * row] = ((double)i * dv[3 * row] + (double)i1 * dv[3 *
        row + 1]) + (double)i2 * dv[3 * row + 2];
    }
  }

  /* 'mso_forces:20' b_wrench = [T_max; -T_min]; */
  /* 'mso_forces:22' distances = (b_wrench - A_wrench * p) ./ vecnorm(A_wrench, 2, 2); */
  /*  返回最小距离（负号用于最大化） */
  /* 'mso_forces:25' obj_val = -min(distances); */
  b_vecnorm(A_wrench, dv1);
  d_obj_next_next_next_next_next_[0] = c_obj_next_next_next_next_next_->T_max[0];
  d_obj_next_next_next_next_next_[3] = -c_obj_next_next_next_next_next_->T_min[0];
  d_obj_next_next_next_next_next_[1] = c_obj_next_next_next_next_next_->T_max[1];
  d_obj_next_next_next_next_next_[4] = -c_obj_next_next_next_next_next_->T_min[1];
  d_obj_next_next_next_next_next_[2] = c_obj_next_next_next_next_next_->T_max[2];
  d_obj_next_next_next_next_next_[5] = -c_obj_next_next_next_next_next_->T_min[2];
  fval = c_obj_next_next_next_next_next_->center[0];
  d = c_obj_next_next_next_next_next_->center[1];
  d1 = c_obj_next_next_next_next_next_->center[2];
  for (col = 0; col < 6; col++) {
    d_obj_next_next_next_next_next_[col] = (d_obj_next_next_next_next_next_[col]
      - ((A_wrench[col] * fval + A_wrench[col + 6] * d) + A_wrench[col + 12] *
         d1)) / dv1[col];
  }

  fval = -minimum(d_obj_next_next_next_next_next_);
  *status = 1;
  allFinite = rtIsNaN(fval);
  if (rtIsInf(fval) || allFinite) {
    if (allFinite) {
      *status = -3;
    } else if (fval < 0.0) {
      *status = -1;
    } else {
      *status = -2;
    }
  }

  if (*status == 1) {
    double JacEqTrans_tmp[27];
    mso_forces_anonFcn2(x, Ceq_tmp, JacEqTrans_tmp);
    Ceq_workspace[0] = Ceq_tmp[0];
    Ceq_workspace[1] = Ceq_tmp[1];
    Ceq_workspace[2] = Ceq_tmp[2];
    for (row = 0; row < 9; row++) {
      JacEqTrans_workspace[row] = JacEqTrans_tmp[row];
      JacEqTrans_workspace[row + 22] = JacEqTrans_tmp[row + 9];
      JacEqTrans_workspace[row + 44] = JacEqTrans_tmp[row + 18];
    }

    *status = checkVectorNonFinite(Ceq_workspace);
    if (*status == 1) {
      int idx_mat;
      allFinite = true;
      row = -1;
      col = -1;
      while (allFinite && (col + 2 <= 3)) {
        row = -1;
        while (allFinite && (row + 2 <= 9)) {
          idx_mat = (row + 22 * (col + 1)) + 1;
          allFinite = ((!rtIsInf(JacEqTrans_workspace[idx_mat])) && (!rtIsNaN
            (JacEqTrans_workspace[idx_mat])));
          row++;
        }

        col++;
      }

      if (!allFinite) {
        idx_mat = row + 22 * col;
        if (rtIsNaN(JacEqTrans_workspace[idx_mat])) {
          *status = -3;
        } else if (JacEqTrans_workspace[idx_mat] < 0.0) {
          *status = -1;
        } else {
          *status = -2;
        }
      }
    }
  }

  return fval;
}

/*
 * Arguments    : g_struct_T *obj
 *                const double A[880]
 *                int mrows
 *                int ncols
 * Return Type  : void
 */
static void factorQR(g_struct_T *obj, const double A[880], int mrows, int ncols)
{
  int i;
  int idx;
  int k;
  bool guard1;
  i = mrows * ncols;
  guard1 = false;
  if (i > 0) {
    for (idx = 0; idx < ncols; idx++) {
      int ix0;
      int iy0;
      ix0 = 22 * idx;
      iy0 = 40 * idx;
      i = (unsigned char)mrows;
      for (k = 0; k < i; k++) {
        obj->QR[iy0 + k] = A[ix0 + k];
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
      obj->jpvt[idx] = idx + 1;
    }

    if (mrows <= ncols) {
      i = mrows;
    } else {
      i = ncols;
    }

    obj->minRowCol = i;
    memset(&obj->tau[0], 0, 40U * sizeof(double));
    if (i >= 1) {
      qrf(obj->QR, mrows, ncols, i, obj->tau);
    }
  }
}

/*
 * Arguments    : g_struct_T *obj
 *                const double A[880]
 *                int mrows
 *                int ncols
 * Return Type  : void
 */
static void factorQRE(g_struct_T *obj, const double A[880], int mrows, int ncols)
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
      ix0 = 22 * idx;
      iy0 = 40 * idx;
      y = (unsigned char)mrows;
      for (k = 0; k < y; k++) {
        obj->QR[iy0 + k] = A[ix0 + k];
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
    xgeqp3(obj->QR, mrows, ncols, obj->jpvt, obj->tau);
  }
}

/*
 * Arguments    : double workspace[880]
 *                double xCurrent[22]
 *                e_struct_T *workingset
 *                g_struct_T *qrmanager
 * Return Type  : bool
 */
static bool feasibleX0ForWorkingSet(double workspace[880], double xCurrent[22],
  e_struct_T *workingset, g_struct_T *qrmanager)
{
  double B[880];
  int ar;
  int b_i;
  int ix;
  int j;
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
      int rankQR;
      i = (unsigned char)nVar;
      for (minmn = 0; minmn < i; minmn++) {
        ix = 40 * minmn;
        for (ar = 0; ar < mWConstr; ar++) {
          qrmanager->QR[ar + ix] = workingset->ATwset[minmn + 22 * ar];
        }

        qrmanager->jpvt[minmn] = 0;
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
        xgeqp3(qrmanager->QR, mWConstr, nVar, qrmanager->jpvt, qrmanager->tau);
      }

      computeQ_(qrmanager, qrmanager->mrows);
      rankQR = 0;
      ix = qrmanager->mrows;
      minmn = qrmanager->ncols;
      if (ix <= minmn) {
        minmn = ix;
      }

      if (minmn > 0) {
        ar = qrmanager->ncols;
        if (ix >= ar) {
          ar = ix;
        }

        tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)ar) *
          fabs(qrmanager->QR[0]);
        while ((rankQR < minmn) && (!(fabs(qrmanager->QR[rankQR + 40 * rankQR]) <=
                 tol))) {
          rankQR++;
        }
      }

      for (minmn = 0; minmn < mWConstr; minmn++) {
        d = workingset->bwset[minmn];
        workspace[minmn] = d;
        workspace[minmn + 40] = d;
      }

      i1 = 22 * (mWConstr - 1) + 1;
      for (minmn = 1; minmn <= i1; minmn += 22) {
        tol = 0.0;
        i2 = minmn + nVar;
        for (ix = minmn; ix < i2; ix++) {
          tol += workingset->ATwset[ix - 1] * xCurrent[ix - minmn];
        }

        i2 = div_nde_s32_floor(minmn - 1, 22);
        workspace[i2] -= tol;
      }

      memcpy(&B[0], &workspace[0], 880U * sizeof(double));
      for (j = 0; j <= 40; j += 40) {
        i1 = j + 1;
        i2 = j + nVar;
        if (i1 <= i2) {
          memset(&workspace[i1 + -1], 0, (unsigned int)((i2 - i1) + 1) * sizeof
                 (double));
        }
      }

      minmn = -1;
      for (j = 0; j <= 40; j += 40) {
        ar = -1;
        i1 = j + 1;
        i2 = j + nVar;
        for (b_i = i1; b_i <= i2; b_i++) {
          tol = 0.0;
          for (ix = 0; ix < mWConstr; ix++) {
            tol += qrmanager->Q[(ix + ar) + 1] * B[(ix + minmn) + 1];
          }

          workspace[b_i - 1] += tol;
          ar += 40;
        }

        minmn += 40;
      }

      for (j = 0; j < 2; j++) {
        ix = 40 * j - 1;
        for (k = rankQR; k >= 1; k--) {
          minmn = 40 * (k - 1) - 1;
          i1 = k + ix;
          d = workspace[i1];
          if (d != 0.0) {
            workspace[i1] = d / qrmanager->QR[k + minmn];
            i2 = (unsigned char)(k - 1);
            for (b_i = 0; b_i < i2; b_i++) {
              ar = (b_i + ix) + 1;
              workspace[ar] -= workspace[i1] * qrmanager->QR[(b_i + minmn) + 1];
            }
          }
        }
      }

      i1 = rankQR + 1;
      for (b_i = i1; b_i <= nVar; b_i++) {
        workspace[b_i - 1] = 0.0;
        workspace[b_i + 39] = 0.0;
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace[qrmanager->jpvt[b_i] + 79] = workspace[b_i];
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace[b_i] = workspace[b_i + 80];
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace[qrmanager->jpvt[b_i] + 79] = workspace[b_i + 40];
      }

      for (b_i = 0; b_i < i; b_i++) {
        workspace[b_i + 40] = workspace[b_i + 80];
      }
    } else {
      int i1;
      int rankQR;
      if (mWConstr - 1 >= 0) {
        memset(&qrmanager->jpvt[0], 0, (unsigned int)mWConstr * sizeof(int));
      }

      factorQRE(qrmanager, workingset->ATwset, nVar, mWConstr);
      computeQ_(qrmanager, qrmanager->minRowCol);
      rankQR = 0;
      ix = qrmanager->mrows;
      minmn = qrmanager->ncols;
      if (ix <= minmn) {
        minmn = ix;
      }

      if (minmn > 0) {
        ar = qrmanager->ncols;
        if (ix >= ar) {
          ar = ix;
        }

        tol = fmin(1.4901161193847656E-8, 2.2204460492503131E-15 * (double)ar) *
          fabs(qrmanager->QR[0]);
        while ((rankQR < minmn) && (!(fabs(qrmanager->QR[rankQR + 40 * rankQR]) <=
                 tol))) {
          rankQR++;
        }
      }

      for (minmn = 0; minmn < mWConstr; minmn++) {
        ix = (qrmanager->jpvt[minmn] - 1) * 22;
        tol = 0.0;
        i = (unsigned char)nVar;
        for (k = 0; k < i; k++) {
          tol += workingset->ATwset[ix + k] * xCurrent[k];
        }

        d = workingset->bwset[qrmanager->jpvt[minmn] - 1];
        workspace[minmn] = d - tol;
        workspace[minmn + 40] = d;
      }

      i = (unsigned char)rankQR;
      for (j = 0; j < 2; j++) {
        ix = 40 * j;
        for (b_i = 0; b_i < i; b_i++) {
          minmn = 40 * b_i;
          ar = b_i + ix;
          tol = workspace[ar];
          for (k = 0; k < b_i; k++) {
            tol -= qrmanager->QR[k + minmn] * workspace[k + ix];
          }

          workspace[ar] = tol / qrmanager->QR[b_i + minmn];
        }
      }

      memcpy(&B[0], &workspace[0], 880U * sizeof(double));
      for (j = 0; j <= 40; j += 40) {
        i = j + 1;
        i1 = j + nVar;
        if (i <= i1) {
          memset(&workspace[i + -1], 0, (unsigned int)((i1 - i) + 1) * sizeof
                 (double));
        }
      }

      minmn = 1;
      for (j = 0; j <= 40; j += 40) {
        ar = -1;
        i = minmn + rankQR;
        for (ix = minmn; ix < i; ix++) {
          int i2;
          i1 = j + 1;
          i2 = j + nVar;
          for (b_i = i1; b_i <= i2; b_i++) {
            workspace[b_i - 1] += B[ix - 1] * qrmanager->Q[(ar + b_i) - j];
          }

          ar += 40;
        }

        minmn += 40;
      }
    }

    minmn = 0;
    int exitg1;
    do {
      exitg1 = 0;
      if (minmn <= (unsigned char)nVar - 1) {
        if (rtIsInf(workspace[minmn]) || rtIsNaN(workspace[minmn])) {
          nonDegenerateWset = false;
          exitg1 = 1;
        } else {
          d = workspace[minmn + 40];
          if (rtIsInf(d) || rtIsNaN(d)) {
            nonDegenerateWset = false;
            exitg1 = 1;
          } else {
            minmn++;
          }
        }
      } else {
        double v;
        for (k = 0; k < nVar; k++) {
          workspace[k] += xCurrent[k];
        }

        if (workingset->probType == 2) {
          v = 0.0;
          for (k = 0; k < 6; k++) {
            workingset->maxConstrWorkspace[k] = -workingset->bineq[k];
          }

          for (minmn = 0; minmn <= 110; minmn += 22) {
            tol = 0.0;
            i = minmn + 9;
            for (ix = minmn + 1; ix <= i; ix++) {
              tol += workingset->Aineq[ix - 1] * workspace[(ix - minmn) - 1];
            }

            i = div_nde_s32_floor(minmn, 22);
            workingset->maxConstrWorkspace[i] += tol;
          }

          for (minmn = 0; minmn < 6; minmn++) {
            d = workingset->maxConstrWorkspace[minmn] - workspace[minmn + 9];
            workingset->maxConstrWorkspace[minmn] = d;
            v = fmax(v, d);
          }

          workingset->maxConstrWorkspace[0] = workingset->beq[0];
          workingset->maxConstrWorkspace[1] = workingset->beq[1];
          workingset->maxConstrWorkspace[2] = workingset->beq[2];
          xgemv(9, workingset->Aeq, workspace, workingset->maxConstrWorkspace);
          d = (workingset->maxConstrWorkspace[0] - workspace[15]) + workspace[18];
          workingset->maxConstrWorkspace[0] = d;
          v = fmax(v, fabs(d));
          d = (workingset->maxConstrWorkspace[1] - workspace[16]) + workspace[19];
          workingset->maxConstrWorkspace[1] = d;
          v = fmax(v, fabs(d));
          d = (workingset->maxConstrWorkspace[2] - workspace[17]) + workspace[20];
          workingset->maxConstrWorkspace[2] = d;
          v = fmax(v, fabs(d));
        } else {
          v = 0.0;
          for (k = 0; k < 6; k++) {
            workingset->maxConstrWorkspace[k] = -workingset->bineq[k];
          }

          for (minmn = 0; minmn <= 110; minmn += 22) {
            tol = 0.0;
            i = minmn + workingset->nVar;
            for (ix = minmn + 1; ix <= i; ix++) {
              tol += workingset->Aineq[ix - 1] * workspace[(ix - minmn) - 1];
            }

            i = div_nde_s32_floor(minmn, 22);
            workingset->maxConstrWorkspace[i] += tol;
          }

          for (minmn = 0; minmn < 6; minmn++) {
            v = fmax(v, workingset->maxConstrWorkspace[minmn]);
          }

          workingset->maxConstrWorkspace[0] = workingset->beq[0];
          workingset->maxConstrWorkspace[1] = workingset->beq[1];
          workingset->maxConstrWorkspace[2] = workingset->beq[2];
          xgemv(workingset->nVar, workingset->Aeq, workspace,
                workingset->maxConstrWorkspace);
          v = fmax(v, fabs(workingset->maxConstrWorkspace[0]));
          v = fmax(v, fabs(workingset->maxConstrWorkspace[1]));
          v = fmax(v, fabs(workingset->maxConstrWorkspace[2]));
        }

        i = (unsigned char)workingset->sizes[3];
        for (minmn = 0; minmn < i; minmn++) {
          ix = workingset->indexLB[minmn] - 1;
          v = fmax(v, -workspace[ix] - workingset->lb[ix]);
        }

        for (minmn = 0; minmn < 9; minmn++) {
          ix = workingset->indexUB[minmn] - 1;
          v = fmax(v, workspace[ix] - workingset->ub[ix]);
        }

        tol = maxConstraintViolation(workingset, workspace);
        if ((v <= 2.2204460492503131E-16) || (v < tol)) {
          i = (unsigned char)nVar;
          memcpy(&xCurrent[0], &workspace[0], (unsigned int)i * sizeof(double));
        } else {
          i = (unsigned char)nVar;
          memcpy(&xCurrent[0], &workspace[40], (unsigned int)i * sizeof(double));
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return nonDegenerateWset;
}

/*
 * Arguments    : const double solution_xstar[22]
 *                const double solution_searchDir[22]
 *                double workspace[880]
 *                int workingset_nVar
 *                const double workingset_Aineq[132]
 *                const double workingset_bineq[6]
 *                const double workingset_lb[22]
 *                const double workingset_ub[22]
 *                const int workingset_indexLB[22]
 *                const int workingset_indexUB[22]
 *                const int workingset_sizes[5]
 *                const int workingset_isActiveIdx[6]
 *                const bool workingset_isActiveConstr[40]
 *                const int workingset_nWConstr[5]
 *                bool isPhaseOne
 *                bool *newBlocking
 *                int *constrType
 *                int *constrIdx
 * Return Type  : double
 */
static double feasibleratiotest(const double solution_xstar[22], const double
  solution_searchDir[22], double workspace[880], int workingset_nVar, const
  double workingset_Aineq[132], const double workingset_bineq[6], const double
  workingset_lb[22], const double workingset_ub[22], const int
  workingset_indexLB[22], const int workingset_indexUB[22], const int
  workingset_sizes[5], const int workingset_isActiveIdx[6], const bool
  workingset_isActiveConstr[40], const int workingset_nWConstr[5], bool
  isPhaseOne, bool *newBlocking, int *constrType, int *constrIdx)
{
  double alpha;
  double c;
  double denomTol;
  double phaseOneCorrectionP;
  double phaseOneCorrectionX;
  double ratio;
  int i;
  int ia;
  int k;
  alpha = 1.0E+30;
  *newBlocking = false;
  *constrType = 0;
  *constrIdx = 0;
  denomTol = 2.2204460492503131E-13 * c_xnrm2(workingset_nVar,
    solution_searchDir);
  if (workingset_nWConstr[2] < 6) {
    for (k = 0; k < 6; k++) {
      workspace[k] = -workingset_bineq[k];
    }

    for (k = 0; k <= 110; k += 22) {
      c = 0.0;
      i = k + workingset_nVar;
      for (ia = k + 1; ia <= i; ia++) {
        c += workingset_Aineq[ia - 1] * solution_xstar[(ia - k) - 1];
      }

      i = div_nde_s32_floor(k, 22);
      workspace[i] += c;
    }

    for (k = 0; k < 6; k++) {
      workspace[k + 40] = 0.0;
    }

    for (k = 0; k <= 110; k += 22) {
      c = 0.0;
      i = k + workingset_nVar;
      for (ia = k + 1; ia <= i; ia++) {
        c += workingset_Aineq[ia - 1] * solution_searchDir[(ia - k) - 1];
      }

      i = div_nde_s32_floor(k, 22) + 40;
      workspace[i] += c;
    }

    for (ia = 0; ia < 6; ia++) {
      phaseOneCorrectionX = workspace[ia + 40];
      if ((phaseOneCorrectionX > denomTol) && (!workingset_isActiveConstr
           [(workingset_isActiveIdx[2] + ia) - 1])) {
        c = workspace[ia];
        c = fmin(fabs(c), 1.0E-6 - c) / phaseOneCorrectionX;
        if (c < alpha) {
          alpha = c;
          *constrType = 3;
          *constrIdx = ia + 1;
          *newBlocking = true;
        }
      }
    }
  }

  if (workingset_nWConstr[3] < workingset_sizes[3]) {
    phaseOneCorrectionX = (double)isPhaseOne * solution_xstar[workingset_nVar -
      1];
    phaseOneCorrectionP = (double)isPhaseOne *
      solution_searchDir[workingset_nVar - 1];
    i = (unsigned char)(workingset_sizes[3] - 1);
    for (ia = 0; ia < i; ia++) {
      k = workingset_indexLB[ia];
      c = -solution_searchDir[k - 1] - phaseOneCorrectionP;
      if ((c > denomTol) && (!workingset_isActiveConstr[(workingset_isActiveIdx
            [3] + ia) - 1])) {
        ratio = (-solution_xstar[k - 1] - workingset_lb[k - 1]) -
          phaseOneCorrectionX;
        c = fmin(fabs(ratio), 1.0E-6 - ratio) / c;
        if (c < alpha) {
          alpha = c;
          *constrType = 4;
          *constrIdx = ia + 1;
          *newBlocking = true;
        }
      }
    }

    i = workingset_indexLB[workingset_sizes[3] - 1] - 1;
    phaseOneCorrectionX = -solution_searchDir[i];
    if ((phaseOneCorrectionX > denomTol) && (!workingset_isActiveConstr
         [(workingset_isActiveIdx[3] + workingset_sizes[3]) - 2])) {
      ratio = -solution_xstar[i] - workingset_lb[i];
      c = fmin(fabs(ratio), 1.0E-6 - ratio) / phaseOneCorrectionX;
      if (c < alpha) {
        alpha = c;
        *constrType = 4;
        *constrIdx = workingset_sizes[3];
        *newBlocking = true;
      }
    }
  }

  if (workingset_nWConstr[4] < 9) {
    phaseOneCorrectionX = (double)isPhaseOne * solution_xstar[workingset_nVar -
      1];
    phaseOneCorrectionP = (double)isPhaseOne *
      solution_searchDir[workingset_nVar - 1];
    for (ia = 0; ia < 9; ia++) {
      i = workingset_indexUB[ia];
      c = solution_searchDir[i - 1] - phaseOneCorrectionP;
      if ((c > denomTol) && (!workingset_isActiveConstr[(workingset_isActiveIdx
            [4] + ia) - 1])) {
        ratio = (solution_xstar[i - 1] - workingset_ub[i - 1]) -
          phaseOneCorrectionX;
        c = fmin(fabs(ratio), 1.0E-6 - ratio) / c;
        if (c < alpha) {
          alpha = c;
          *constrType = 5;
          *constrIdx = ia + 1;
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
 * Arguments    : const double fun_workspace_center[3]
 *                const double fun_workspace_T_max[3]
 *                const double fun_workspace_T_min[3]
 *                const double x0[9]
 *                const double Aineq[54]
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
                      fun_workspace_T_max[3], const double fun_workspace_T_min[3],
                      const double x0[9], const double Aineq[54], double x[9],
                      double *exitflag, double *output_iterations, double
                      *output_funcCount, char output_algorithm[3], double
                      *output_constrviolation, double *output_stepsize, double
                      *output_lssteplength, double *output_firstorderopt)
{
  static const signed char b_iv[22] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0 };

  static const signed char b_iv1[6] = { 1, 0, 3, 6, 9, 9 };

  static const signed char iv5[6] = { 1, 0, 3, 6, 10, 9 };

  static const signed char iv6[6] = { 1, 0, 3, 6, 21, 9 };

  static const signed char iv7[6] = { 1, 0, 3, 6, 22, 9 };

  static const signed char iv2[5] = { 0, 3, 6, 10, 9 };

  static const signed char iv3[5] = { 0, 3, 6, 21, 9 };

  static const signed char iv4[5] = { 0, 3, 6, 22, 9 };

  static const signed char obj_tmp[5] = { 0, 3, 6, 9, 9 };

  c_struct_T TrialState;
  d_struct_T FiniteDifferences;
  e_struct_T WorkingSet;
  f_struct_T MeritFunction;
  g_struct_T QRManager;
  h_struct_T CholManager;
  i_coder_internal_stickyStruct FcnEvaluator;
  i_struct_T QPObjective;
  struct_T memspace;
  double c;
  double fval;
  double scale;
  int i;
  int i1;
  int iEq0;
  int idx;
  signed char b_i;
  output_algorithm[0] = 's';
  output_algorithm[1] = 'q';
  output_algorithm[2] = 'p';
  TrialState.sqpFval_old = 0.0;
  TrialState.sqpIterations = 0;
  TrialState.sqpExitFlag = 0;
  memset(&TrialState.lambdasqp[0], 0, 40U * sizeof(double));
  TrialState.steplength = 1.0;
  memset(&TrialState.delta_x[0], 0, 22U * sizeof(double));
  TrialState.fstar = 0.0;
  memset(&TrialState.lambda[0], 0, 40U * sizeof(double));
  TrialState.state = 0;
  TrialState.maxConstr = 0.0;
  TrialState.iterations = 0;
  memcpy(&TrialState.xstarsqp[0], &x0[0], 9U * sizeof(double));
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.center[0]
    = fun_workspace_center[0];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_max[0] =
    fun_workspace_T_max[0];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_min[0] =
    fun_workspace_T_min[0];
  FiniteDifferences.objfun.workspace.center[0] = fun_workspace_center[0];
  FiniteDifferences.objfun.workspace.T_max[0] = fun_workspace_T_max[0];
  FiniteDifferences.objfun.workspace.T_min[0] = fun_workspace_T_min[0];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.center[1]
    = fun_workspace_center[1];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_max[1] =
    fun_workspace_T_max[1];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_min[1] =
    fun_workspace_T_min[1];
  FiniteDifferences.objfun.workspace.center[1] = fun_workspace_center[1];
  FiniteDifferences.objfun.workspace.T_max[1] = fun_workspace_T_max[1];
  FiniteDifferences.objfun.workspace.T_min[1] = fun_workspace_T_min[1];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.center[2]
    = fun_workspace_center[2];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_max[2] =
    fun_workspace_T_max[2];
  FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace.T_min[2] =
    fun_workspace_T_min[2];
  FiniteDifferences.objfun.workspace.center[2] = fun_workspace_center[2];
  FiniteDifferences.objfun.workspace.T_max[2] = fun_workspace_T_max[2];
  FiniteDifferences.objfun.workspace.T_min[2] = fun_workspace_T_min[2];
  FiniteDifferences.f_1 = 0.0;
  FiniteDifferences.numEvals = 0;
  FiniteDifferences.hasBounds = true;
  WorkingSet.nVar = 9;
  memset(&WorkingSet.Aineq[0], 0, 132U * sizeof(double));
  for (i = 0; i < 6; i++) {
    WorkingSet.bineq[i] = 0.0;
  }

  memset(&WorkingSet.Aeq[0], 0, 66U * sizeof(double));
  WorkingSet.beq[0] = 0.0;
  WorkingSet.beq[1] = 0.0;
  WorkingSet.beq[2] = 0.0;
  for (i = 0; i < 22; i++) {
    WorkingSet.lb[i] = 0.0;
    WorkingSet.ub[i] = 0.0;
    b_i = b_iv[i];
    WorkingSet.indexLB[i] = b_i;
    WorkingSet.indexUB[i] = b_i;
  }

  WorkingSet.mEqRemoved = 0;
  WorkingSet.indexEqRemoved[0] = 0;
  WorkingSet.indexEqRemoved[1] = 0;
  WorkingSet.indexEqRemoved[2] = 0;
  memset(&WorkingSet.ATwset[0], 0, 880U * sizeof(double));
  WorkingSet.nActiveConstr = 0;
  memset(&WorkingSet.bwset[0], 0, 40U * sizeof(double));
  memset(&WorkingSet.maxConstrWorkspace[0], 0, 40U * sizeof(double));
  memset(&WorkingSet.Wid[0], 0, 40U * sizeof(int));
  memset(&WorkingSet.Wlocalidx[0], 0, 40U * sizeof(int));
  for (i = 0; i < 40; i++) {
    WorkingSet.isActiveConstr[i] = false;
  }

  WorkingSet.probType = 3;
  WorkingSet.mConstr = 27;
  for (i = 0; i < 5; i++) {
    WorkingSet.nWConstr[i] = 0;
    b_i = obj_tmp[i];
    WorkingSet.sizes[i] = b_i;
    WorkingSet.sizesNormal[i] = b_i;
    WorkingSet.sizesPhaseOne[i] = iv2[i];
    WorkingSet.sizesRegularized[i] = iv3[i];
    WorkingSet.sizesRegPhaseOne[i] = iv4[i];
  }

  for (i = 0; i < 6; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i] = b_iv1[i];
  }

  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (i = 0; i < 6; i++) {
    i1 = WorkingSet.isActiveIdxRegPhaseOne[i];
    WorkingSet.isActiveIdx[i] = i1;
    WorkingSet.isActiveIdxNormal[i] = i1;
    WorkingSet.isActiveIdxRegPhaseOne[i] = iv5[i];
  }

  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (i = 0; i < 6; i++) {
    WorkingSet.isActiveIdxPhaseOne[i] = WorkingSet.isActiveIdxRegPhaseOne[i];
    WorkingSet.isActiveIdxRegPhaseOne[i] = iv6[i];
  }

  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (i = 0; i < 6; i++) {
    WorkingSet.isActiveIdxRegularized[i] = WorkingSet.isActiveIdxRegPhaseOne[i];
    WorkingSet.isActiveIdxRegPhaseOne[i] = iv7[i];
  }

  for (i = 0; i < 5; i++) {
    WorkingSet.isActiveIdxRegPhaseOne[i + 1] +=
      WorkingSet.isActiveIdxRegPhaseOne[i];
  }

  for (i = 0; i < 6; i++) {
    for (iEq0 = 0; iEq0 < 9; iEq0++) {
      WorkingSet.Aineq[iEq0 + 22 * i] = Aineq[i + 6 * iEq0];
    }
  }

  for (idx = 0; idx < 9; idx++) {
    i1 = WorkingSet.indexLB[idx];
    TrialState.xstarsqp[i1 - 1] = fmax(TrialState.xstarsqp[i1 - 1], -1.0);
  }

  for (idx = 0; idx < 9; idx++) {
    i1 = WorkingSet.indexUB[idx];
    TrialState.xstarsqp[i1 - 1] = fmin(TrialState.xstarsqp[i1 - 1], iv[i1 - 1]);
  }

  TrialState.sqpFval = evalObjAndConstrAndDerivatives
    (&FcnEvaluator.next.next.next.next.next.next.next.next.value.workspace,
     TrialState.xstarsqp, TrialState.cEq, WorkingSet.Aeq, &i);
  computeFiniteDifferences(&FiniteDifferences, TrialState.sqpFval,
    TrialState.xstarsqp, TrialState.grad);
  TrialState.FunctionEvaluations = FiniteDifferences.numEvals + 1;
  for (i = 0; i < 6; i++) {
    TrialState.cIneq[i] = -0.0;
  }

  for (i = 0; i <= 110; i += 22) {
    c = 0.0;
    i1 = i + 9;
    for (iEq0 = i + 1; iEq0 <= i1; iEq0++) {
      c += WorkingSet.Aineq[iEq0 - 1] * TrialState.xstarsqp[(iEq0 - i) - 1];
    }

    i1 = div_nde_s32_floor(i, 22);
    TrialState.cIneq[i1] += c;
  }

  i = 0;
  iEq0 = 0;
  for (idx = 0; idx < 3; idx++) {
    c = TrialState.cEq[idx];
    WorkingSet.beq[idx] = -c;
    WorkingSet.bwset[idx] = -c;
    memcpy(&WorkingSet.ATwset[i], &WorkingSet.Aeq[iEq0], 9U * sizeof(double));
    iEq0 = i + 22;
    i += 22;
  }

  for (idx = 0; idx < 6; idx++) {
    WorkingSet.bineq[idx] = -TrialState.cIneq[idx];
  }

  for (idx = 0; idx < 9; idx++) {
    WorkingSet.lb[WorkingSet.indexLB[idx] - 1] = x0[WorkingSet.indexLB[idx] - 1]
      + 1.0;
    WorkingSet.ub[WorkingSet.indexUB[idx] - 1] = (double)
      iv[WorkingSet.indexUB[idx] - 1] - x0[WorkingSet.indexUB[idx] - 1];
  }

  setProblemType(&WorkingSet, 3);
  i = WorkingSet.isActiveIdx[2];
  for (idx = i; idx < 41; idx++) {
    WorkingSet.isActiveConstr[idx - 1] = false;
  }

  WorkingSet.nWConstr[0] = 0;
  WorkingSet.nWConstr[1] = 3;
  WorkingSet.nWConstr[2] = 0;
  WorkingSet.nWConstr[3] = 0;
  WorkingSet.nWConstr[4] = 0;
  WorkingSet.nActiveConstr = 3;
  MeritFunction.initFval = TrialState.sqpFval;
  MeritFunction.penaltyParam = 1.0;
  MeritFunction.threshold = 0.0001;
  MeritFunction.nPenaltyDecreases = 0;
  MeritFunction.linearizedConstrViol = 0.0;
  scale = 0.0;
  i1 = WorkingSet.nVar;
  for (i = 0; i < 3; i++) {
    WorkingSet.Wid[i] = 2;
    WorkingSet.Wlocalidx[i] = i + 1;
    WorkingSet.isActiveConstr[i] = true;
    iEq0 = 22 * i;
    memcpy(&WorkingSet.ATwset[iEq0], &WorkingSet.Aeq[iEq0], (unsigned int)i1 *
           sizeof(double));
    WorkingSet.bwset[i] = WorkingSet.beq[i];
    scale += fabs(TrialState.cEq[i]);
  }

  MeritFunction.initConstrViolationEq = scale;
  scale = 0.0;
  for (idx = 0; idx < 6; idx++) {
    c = TrialState.cIneq[idx];
    if (c > 0.0) {
      scale += c;
    }
  }

  double Hessian[81];
  MeritFunction.initConstrViolationIneq = scale;
  MeritFunction.phi = 0.0;
  MeritFunction.phiPrimePlus = 0.0;
  MeritFunction.phiFullStep = 0.0;
  MeritFunction.feasRelativeFactor = 0.0;
  MeritFunction.nlpPrimalFeasError = 0.0;
  MeritFunction.nlpDualFeasError = 0.0;
  MeritFunction.nlpComplError = 0.0;
  MeritFunction.firstOrderOpt = 0.0;
  b_driver(&TrialState, &MeritFunction, &FcnEvaluator, &FiniteDifferences,
           &memspace, &WorkingSet, Hessian, &QRManager, &CholManager,
           &QPObjective);
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
    absxk = fabs(TrialState.delta_x[i]);
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
 * Arguments    : h_struct_T *obj
 *                int NColsRemain
 * Return Type  : void
 */
static void fullColLDL2_(h_struct_T *obj, int NColsRemain)
{
  int ijA;
  int j;
  int jA;
  int k;
  for (k = 0; k < NColsRemain; k++) {
    double alpha1;
    double y;
    int LD_diagOffset;
    int i;
    int offset1;
    int subMatrixDim;
    LD_diagOffset = 41 * k;
    alpha1 = -1.0 / obj->FMat[LD_diagOffset];
    subMatrixDim = (NColsRemain - k) - 2;
    offset1 = LD_diagOffset + 2;
    y = obj->workspace_;
    for (jA = 0; jA <= subMatrixDim; jA++) {
      y = obj->FMat[(LD_diagOffset + jA) + 1];
    }

    obj->workspace_ = y;
    if (!(alpha1 == 0.0)) {
      jA = LD_diagOffset;
      for (j = 0; j <= subMatrixDim; j++) {
        if (y != 0.0) {
          double temp;
          int i1;
          temp = y * alpha1;
          i = jA + 42;
          i1 = subMatrixDim + jA;
          for (ijA = i; ijA <= i1 + 42; ijA++) {
            obj->FMat[ijA - 1] += y * temp;
          }
        }

        jA += 40;
      }
    }

    alpha1 = 1.0 / obj->FMat[LD_diagOffset];
    i = LD_diagOffset + subMatrixDim;
    for (jA = offset1; jA <= i + 2; jA++) {
      obj->FMat[jA - 1] *= alpha1;
    }
  }
}

/*
 * Arguments    : bool obj_hasLinear
 *                int obj_nvar
 *                double workspace[880]
 *                const double H[81]
 *                const double f[22]
 *                const double x[22]
 * Return Type  : void
 */
static void linearForm_(bool obj_hasLinear, int obj_nvar, double workspace[880],
  const double H[81], const double f[22], const double x[22])
{
  int ia;
  int iac;
  int ix;
  ix = 0;
  if (obj_hasLinear) {
    if (obj_nvar - 1 >= 0) {
      memcpy(&workspace[0], &f[0], (unsigned int)obj_nvar * sizeof(double));
    }

    ix = 1;
  }

  if (obj_nvar != 0) {
    int i;
    if ((ix != 1) && (obj_nvar - 1 >= 0)) {
      memset(&workspace[0], 0, (unsigned int)obj_nvar * sizeof(double));
    }

    ix = 0;
    i = obj_nvar * (obj_nvar - 1) + 1;
    for (iac = 1; obj_nvar < 0 ? iac >= i : iac <= i; iac += obj_nvar) {
      double c;
      int i1;
      c = 0.5 * x[ix];
      i1 = iac + obj_nvar;
      for (ia = iac; ia < i1; ia++) {
        int i2;
        i2 = ia - iac;
        workspace[i2] += H[ia - 1] * c;
      }

      ix++;
    }
  }
}

/*
 * Arguments    : e_struct_T *obj
 *                const double x[880]
 * Return Type  : double
 */
static double maxConstraintViolation(e_struct_T *obj, const double x[880])
{
  double v;
  int i;
  int ia;
  int idxLB;
  if (obj->probType == 2) {
    double c;
    v = 0.0;
    for (idxLB = 0; idxLB < 6; idxLB++) {
      obj->maxConstrWorkspace[idxLB] = -obj->bineq[idxLB];
    }

    for (idxLB = 0; idxLB <= 110; idxLB += 22) {
      c = 0.0;
      i = idxLB + 9;
      for (ia = idxLB + 1; ia <= i; ia++) {
        c += obj->Aineq[ia - 1] * x[(ia - idxLB) + 39];
      }

      i = div_nde_s32_floor(idxLB, 22);
      obj->maxConstrWorkspace[i] += c;
    }

    for (ia = 0; ia < 6; ia++) {
      c = obj->maxConstrWorkspace[ia] - x[ia + 49];
      obj->maxConstrWorkspace[ia] = c;
      v = fmax(v, c);
    }

    obj->maxConstrWorkspace[0] = -obj->beq[0];
    obj->maxConstrWorkspace[1] = -obj->beq[1];
    obj->maxConstrWorkspace[2] = -obj->beq[2];
    for (idxLB = 0; idxLB <= 44; idxLB += 22) {
      c = 0.0;
      i = idxLB + 9;
      for (ia = idxLB + 1; ia <= i; ia++) {
        c += obj->Aeq[ia - 1] * x[(ia - idxLB) + 39];
      }

      i = div_nde_s32_floor(idxLB, 22);
      obj->maxConstrWorkspace[i] += c;
    }

    c = (obj->maxConstrWorkspace[0] - x[55]) + x[58];
    obj->maxConstrWorkspace[0] = c;
    v = fmax(v, fabs(c));
    c = (obj->maxConstrWorkspace[1] - x[56]) + x[59];
    obj->maxConstrWorkspace[1] = c;
    v = fmax(v, fabs(c));
    c = (obj->maxConstrWorkspace[2] - x[57]) + x[60];
    obj->maxConstrWorkspace[2] = c;
    v = fmax(v, fabs(c));
  } else {
    double c;
    v = 0.0;
    for (idxLB = 0; idxLB < 6; idxLB++) {
      obj->maxConstrWorkspace[idxLB] = -obj->bineq[idxLB];
    }

    for (idxLB = 0; idxLB <= 110; idxLB += 22) {
      c = 0.0;
      i = idxLB + obj->nVar;
      for (ia = idxLB + 1; ia <= i; ia++) {
        c += obj->Aineq[ia - 1] * x[(ia - idxLB) + 39];
      }

      i = div_nde_s32_floor(idxLB, 22);
      obj->maxConstrWorkspace[i] += c;
    }

    for (ia = 0; ia < 6; ia++) {
      v = fmax(v, obj->maxConstrWorkspace[ia]);
    }

    obj->maxConstrWorkspace[0] = -obj->beq[0];
    obj->maxConstrWorkspace[1] = -obj->beq[1];
    obj->maxConstrWorkspace[2] = -obj->beq[2];
    for (idxLB = 0; idxLB <= 44; idxLB += 22) {
      c = 0.0;
      i = idxLB + obj->nVar;
      for (ia = idxLB + 1; ia <= i; ia++) {
        c += obj->Aeq[ia - 1] * x[(ia - idxLB) + 39];
      }

      i = div_nde_s32_floor(idxLB, 22);
      obj->maxConstrWorkspace[i] += c;
    }

    v = fmax(v, fabs(obj->maxConstrWorkspace[0]));
    v = fmax(v, fabs(obj->maxConstrWorkspace[1]));
    v = fmax(v, fabs(obj->maxConstrWorkspace[2]));
  }

  i = (unsigned char)obj->sizes[3];
  for (ia = 0; ia < i; ia++) {
    idxLB = obj->indexLB[ia];
    v = fmax(v, -x[idxLB + 39] - obj->lb[idxLB - 1]);
  }

  for (ia = 0; ia < 9; ia++) {
    idxLB = obj->indexUB[ia];
    v = fmax(v, x[idxLB + 39] - obj->ub[idxLB - 1]);
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
 * Arguments    : e_struct_T *obj
 * Return Type  : void
 */
static void modifyOverheadPhaseOne_(e_struct_T *obj)
{
  int i;
  int idx;
  int idxEq;
  int idxEq_tmp;
  idxEq_tmp = obj->nVar - 1;
  obj->Aeq[idxEq_tmp] = 0.0;
  i = 22 * (obj->isActiveIdx[1] - 1);
  obj->ATwset[idxEq_tmp + i] = 0.0;
  idxEq = obj->nVar + 21;
  obj->Aeq[idxEq] = 0.0;
  obj->ATwset[idxEq + i] = 0.0;
  idxEq = obj->nVar + 43;
  obj->Aeq[idxEq] = 0.0;
  obj->ATwset[idxEq + i] = 0.0;
  for (idx = 0; idx < 6; idx++) {
    obj->Aineq[(obj->nVar + 22 * idx) - 1] = -1.0;
  }

  obj->indexLB[obj->sizes[3] - 1] = obj->nVar;
  obj->lb[idxEq_tmp] = 1.0E-5;
  idxEq = obj->isActiveIdx[2];
  i = obj->nActiveConstr;
  for (idx = idxEq; idx <= i; idx++) {
    obj->ATwset[(obj->nVar + 22 * (idx - 1)) - 1] = -1.0;
  }

  idxEq = obj->isActiveIdx[4];
  if (obj->nWConstr[4] > 0) {
    for (idx = 8; idx >= 0; idx--) {
      i = idxEq + idx;
      obj->isActiveConstr[i] = obj->isActiveConstr[i - 1];
    }
  } else {
    obj->isActiveConstr[obj->isActiveIdx[4] + 8] = false;
  }

  obj->isActiveConstr[obj->isActiveIdx[4] - 1] = false;
}

/*
 * @(W)
 *
 * Arguments    : const double W[9]
 *                double varargout_2[3]
 *                double varargout_4[27]
 * Return Type  : void
 */
static void mso_forces_anonFcn2(const double W[9], double varargout_2[3], double
  varargout_4[27])
{
  double y[9];
  int k;

  /* 'mso_forces:84' @(W) constraint(W) */
  /*  非线性约束：单位向量约束 */
  /* 'mso_forces:29' c = []; */
  /* 'mso_forces:30' ceq = sum(reshape(W, 3, length(W) / 3).^2, 1) - 1; */
  for (k = 0; k < 9; k++) {
    double d;
    d = W[k];
    y[k] = d * d;
  }

  sum(y, varargout_2);
  varargout_2[0]--;
  varargout_2[1]--;
  varargout_2[2]--;

  /* 'mso_forces:32' if nargout > 2 */
  /* 'mso_forces:33' gradc = []; */
  /* 'mso_forces:34' gradceq = zeros(length(W), length(W) / 3); */
  memset(&varargout_4[0], 0, 27U * sizeof(double));

  /* 'mso_forces:35' for idx = 1:length(W) / 3 */
  for (k = 0; k < 3; k++) {
    int varargout_4_tmp;
    signed char i;

    /* 'mso_forces:36' gradceq(3*idx-2:3*idx, idx) = 2 * W(3*idx-2:3*idx); */
    i = (signed char)(3 * (k + 1));
    varargout_4_tmp = i + 9 * k;
    varargout_4[varargout_4_tmp - 3] = 2.0 * W[i - 3];
    varargout_4[varargout_4_tmp - 2] = 2.0 * W[i - 2];
    varargout_4[varargout_4_tmp - 1] = 2.0 * W[i - 1];
  }
}

/*
 * Arguments    : const double H[81]
 *                const double f[22]
 *                c_struct_T *solution
 *                struct_T *memspace
 *                e_struct_T *workingset
 *                g_struct_T *qrmanager
 *                h_struct_T *cholmanager
 *                i_struct_T *objective
 *                const char options_SolverName[7]
 *                const k_struct_T *runTimeOptions
 * Return Type  : void
 */
static void phaseone(const double H[81], const double f[22], c_struct_T
                     *solution, struct_T *memspace, e_struct_T *workingset,
                     g_struct_T *qrmanager, h_struct_T *cholmanager, i_struct_T *
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
  solution->xstar[workingset->nVar] = solution->maxConstr + 1.0;
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
    workingset->isActiveConstr[(workingset->isActiveIdx[workingset->
      Wid[idx_global - 1] - 1] + workingset->Wlocalidx[idx_global - 1]) - 2] =
      false;
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
  computeGrad_StoreHx(objective, H, f, solution->xstar);
  solution->fstar = computeFval_ReuseHx(objective, memspace->workspace_float, f,
    solution->xstar);
  solution->state = -5;
  memset(&solution->lambda[0], 0, 40U * sizeof(double));
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
          squareQ_appendCol(qrmanager, workingset->ATwset, 22 *
                            (workingset->nActiveConstr - 1) + 1);
          break;

         case -1:
          deleteColMoveEnd(qrmanager, globalActiveConstrIdx);
          break;

         default:
          factorQR(qrmanager, workingset->ATwset, i, workingset->nActiveConstr);
          computeQ_(qrmanager, qrmanager->mrows);
          break;
        }

        ret = memcmp(&options_SolverName[0], &b[0], 7);
        compute_deltax(H, solution, memspace, qrmanager, cholmanager, objective,
                       ret == 0);
        if (solution->state != -5) {
          exitg1 = 1;
        } else if ((c_xnrm2(i, solution->searchDir) < 1.4901161193847657E-10) ||
                   (workingset->nActiveConstr >= i)) {
          guard2 = true;
        } else {
          minLambda = feasibleratiotest(solution->xstar, solution->searchDir,
            memspace->workspace_float, workingset->nVar, workingset->Aineq,
            workingset->bineq, workingset->lb, workingset->ub,
            workingset->indexLB, workingset->indexUB, workingset->sizes,
            workingset->isActiveIdx, workingset->isActiveConstr,
            workingset->nWConstr, true, &updateFval, &idx_global, &idxMinLambda);
          if (updateFval) {
            switch (idx_global) {
             case 3:
              workingset->nWConstr[2]++;
              workingset->isActiveConstr[(workingset->isActiveIdx[2] +
                idxMinLambda) - 2] = true;
              workingset->nActiveConstr++;
              idx_global = workingset->nActiveConstr - 1;
              workingset->Wid[idx_global] = 3;
              workingset->Wlocalidx[idx_global] = idxMinLambda;
              ret = 22 * (idxMinLambda - 1);
              iAw0 = 22 * idx_global;
              i1 = workingset->nVar;
              for (idx = 0; idx < i1; idx++) {
                workingset->ATwset[iAw0 + idx] = workingset->Aineq[ret + idx];
              }

              workingset->bwset[idx_global] = workingset->bineq[idxMinLambda - 1];
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
              if (c_xnrm2(objective->nvar, solution->searchDir) > 100.0 *
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
              solution->xstar[ret] += minLambda * solution->searchDir[ret];
            }
          }

          computeGrad_StoreHx(objective, H, f, solution->xstar);
          updateFval = true;
          guard1 = true;
        }
      } else {
        idx_global = (unsigned char)i;
        memset(&solution->searchDir[0], 0, (unsigned int)idx_global * sizeof
               (double));
        guard2 = true;
      }

      if (guard2) {
        compute_lambda(memspace->workspace_float, solution, objective, qrmanager);
        if ((solution->state != -7) || (workingset->nActiveConstr > i)) {
          idxMinLambda = 0;
          minLambda = 0.0;
          idx_global = (workingset->nWConstr[0] + workingset->nWConstr[1]) + 1;
          i1 = workingset->nActiveConstr;
          for (idx = idx_global; idx <= i1; idx++) {
            double d;
            d = solution->lambda[idx - 1];
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
              solution->lambda[idxMinLambda - 1] = solution->lambda
                [workingset->nActiveConstr];
            }

            solution->lambda[workingset->nActiveConstr] = 0.0;
          }
        } else {
          idx_global = workingset->nActiveConstr;
          activeSetChangeID = 0;
          globalActiveConstrIdx = workingset->nActiveConstr;
          subProblemChanged = true;
          removeConstr(workingset, workingset->nActiveConstr);
          solution->lambda[idx_global - 1] = 0.0;
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
            solution->xstar);
          minLambda = solution->maxConstr;
          if (objective->objtype == 5) {
            minLambda = solution->maxConstr - solution->xstar[objective->nvar -
              1];
          }

          if (minLambda > 1.0E-6) {
            bool nonDegenerateWset;
            if (iAw0 - 1 >= 0) {
              memcpy(&solution->searchDir[0], &solution->xstar[0], (unsigned int)
                     iAw0 * sizeof(double));
            }

            nonDegenerateWset = feasibleX0ForWorkingSet
              (memspace->workspace_float, solution->searchDir, workingset,
               qrmanager);
            if ((!nonDegenerateWset) && (solution->state != 0)) {
              solution->state = -2;
            }

            activeSetChangeID = 0;
            minLambda = b_maxConstraintViolation(workingset, solution->searchDir);
            if ((minLambda < solution->maxConstr) && (iAw0 - 1 >= 0)) {
              memcpy(&solution->xstar[0], &solution->searchDir[0], (unsigned int)
                     iAw0 * sizeof(double));
            }
          }
        }

        if (updateFval) {
          solution->fstar = computeFval_ReuseHx(objective,
            memspace->workspace_float, f, solution->xstar);
          if ((solution->fstar < 1.0E-6) && ((solution->state != 0) ||
               (objective->objtype != 5))) {
            solution->state = 2;
          }
        }
      }
    } else {
      if (!updateFval) {
        solution->fstar = computeFval_ReuseHx(objective,
          memspace->workspace_float, f, solution->xstar);
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (workingset->isActiveConstr[(workingset->isActiveIdx[3] + workingset->
       sizes[3]) - 2]) {
    bool exitg2;
    idx = 4;
    exitg2 = false;
    while ((!exitg2) && (idx <= workingset->nActiveConstr)) {
      if ((workingset->Wid[idx - 1] == 4) && (workingset->Wlocalidx[idx - 1] ==
           workingset->sizes[3])) {
        removeConstr(workingset, idx);
        exitg2 = true;
      } else {
        idx++;
      }
    }
  }

  ret = workingset->nActiveConstr;
  while ((ret > 3) && (ret > nVar_tmp)) {
    removeConstr(workingset, ret);
    ret--;
  }

  solution->maxConstr = solution->xstar[nVar_tmp];
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
 * Arguments    : double A[1600]
 *                int m
 *                int n
 *                int nfxd
 *                double tau[40]
 * Return Type  : void
 */
static void qrf(double A[1600], int m, int n, int nfxd, double tau[40])
{
  double work[40];
  double atmp;
  int b_i;
  int i;
  memset(&tau[0], 0, 40U * sizeof(double));
  memset(&work[0], 0, 40U * sizeof(double));
  i = (unsigned char)nfxd;
  for (b_i = 0; b_i < i; b_i++) {
    double d;
    int ii;
    int mmi;
    ii = b_i * 40 + b_i;
    mmi = m - b_i;
    if (b_i + 1 < m) {
      atmp = A[ii];
      d = xzlarfg(mmi, &atmp, A, ii + 2);
      tau[b_i] = d;
      A[ii] = atmp;
    } else {
      d = 0.0;
      tau[b_i] = 0.0;
    }

    if (b_i + 1 < n) {
      atmp = A[ii];
      A[ii] = 1.0;
      xzlarf(mmi, (n - b_i) - 1, ii + 1, d, A, ii + 41, work);
      A[ii] = atmp;
    }
  }
}

/*
 * Arguments    : const double Hessian[81]
 *                const double grad[22]
 *                c_struct_T *TrialState
 *                f_struct_T *MeritFunction
 *                struct_T *memspace
 *                e_struct_T *WorkingSet
 *                g_struct_T *QRManager
 *                h_struct_T *CholManager
 *                i_struct_T *QPObjective
 *                k_struct_T *qpoptions
 * Return Type  : void
 */
static void relaxed(const double Hessian[81], const double grad[22], c_struct_T *
                    TrialState, f_struct_T *MeritFunction, struct_T *memspace,
                    e_struct_T *WorkingSet, g_struct_T *QRManager, h_struct_T
                    *CholManager, i_struct_T *QPObjective, k_struct_T *qpoptions)
{
  double beta;
  double d;
  double rho;
  double s;
  double smax;
  int iac;
  int idx;
  int idx_max;
  int mLBOrig;
  int nVarOrig;
  bool b_tf;
  bool tf;
  nVarOrig = WorkingSet->nVar;
  beta = 0.0;
  idx_max = (unsigned char)WorkingSet->nVar;
  for (idx = 0; idx < idx_max; idx++) {
    beta += Hessian[idx + 9 * idx];
  }

  beta /= (double)WorkingSet->nVar;
  if (TrialState->sqpIterations <= 1) {
    mLBOrig = QPObjective->nvar;
    if (QPObjective->nvar < 1) {
      idx_max = 0;
    } else {
      idx_max = 1;
      if (QPObjective->nvar > 1) {
        smax = fabs(grad[0]);
        for (idx = 2; idx <= mLBOrig; idx++) {
          s = fabs(grad[idx - 1]);
          if (s > smax) {
            idx_max = idx;
            smax = s;
          }
        }
      }
    }

    rho = 100.0 * fmax(1.0, fabs(grad[idx_max - 1]));
  } else {
    mLBOrig = WorkingSet->mConstr;
    idx_max = 1;
    smax = fabs(TrialState->lambdasqp[0]);
    for (idx = 2; idx <= mLBOrig; idx++) {
      s = fabs(TrialState->lambdasqp[idx - 1]);
      if (s > smax) {
        idx_max = idx;
        smax = s;
      }
    }

    rho = fabs(TrialState->lambdasqp[idx_max - 1]);
  }

  QPObjective->nvar = WorkingSet->nVar;
  QPObjective->beta = beta;
  QPObjective->rho = rho;
  QPObjective->hasLinear = true;
  QPObjective->objtype = 4;
  setProblemType(WorkingSet, 2);
  mLBOrig = WorkingSet->sizes[3] - 5;
  for (idx = 0; idx < 6; idx++) {
    memspace->workspace_float[idx] = -WorkingSet->bineq[idx];
  }

  for (iac = 0; iac <= 110; iac += 22) {
    smax = 0.0;
    idx_max = iac + nVarOrig;
    for (idx = iac + 1; idx <= idx_max; idx++) {
      smax += WorkingSet->Aineq[idx - 1] * TrialState->xstar[(idx - iac) - 1];
    }

    idx_max = div_nde_s32_floor(iac, 22);
    memspace->workspace_float[idx_max] += smax;
  }

  for (idx = 0; idx < 6; idx++) {
    d = memspace->workspace_float[idx];
    TrialState->xstar[nVarOrig + idx] = (double)(d > 0.0) * d;
  }

  memspace->workspace_float[0] = WorkingSet->beq[0];
  memspace->workspace_float[1] = WorkingSet->beq[1];
  memspace->workspace_float[2] = WorkingSet->beq[2];
  memspace->workspace_float[0] = -memspace->workspace_float[0];
  memspace->workspace_float[1] = -memspace->workspace_float[1];
  memspace->workspace_float[2] = -memspace->workspace_float[2];
  for (iac = 0; iac <= 44; iac += 22) {
    smax = 0.0;
    idx_max = iac + nVarOrig;
    for (idx = iac + 1; idx <= idx_max; idx++) {
      smax += WorkingSet->Aeq[idx - 1] * TrialState->xstar[(idx - iac) - 1];
    }

    idx_max = div_nde_s32_floor(iac, 22);
    memspace->workspace_float[idx_max] += smax;
  }

  if (memspace->workspace_float[0] <= 0.0) {
    TrialState->xstar[nVarOrig + 6] = 0.0;
    TrialState->xstar[nVarOrig + 9] = -memspace->workspace_float[0];
    addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig);
    if (memspace->workspace_float[0] >= -1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 3);
    }
  } else {
    TrialState->xstar[nVarOrig + 6] = memspace->workspace_float[0];
    TrialState->xstar[nVarOrig + 9] = 0.0;
    addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 3);
    if (memspace->workspace_float[0] <= 1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig);
    }
  }

  if (memspace->workspace_float[1] <= 0.0) {
    TrialState->xstar[nVarOrig + 7] = 0.0;
    TrialState->xstar[nVarOrig + 10] = -memspace->workspace_float[1];
    addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 1);
    if (memspace->workspace_float[1] >= -1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 4);
    }
  } else {
    TrialState->xstar[nVarOrig + 7] = memspace->workspace_float[1];
    TrialState->xstar[nVarOrig + 10] = 0.0;
    addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 4);
    if (memspace->workspace_float[1] <= 1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 1);
    }
  }

  if (memspace->workspace_float[2] <= 0.0) {
    TrialState->xstar[nVarOrig + 8] = 0.0;
    TrialState->xstar[nVarOrig + 11] = -memspace->workspace_float[2];
    addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 2);
    if (memspace->workspace_float[2] >= -1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 5);
    }
  } else {
    TrialState->xstar[nVarOrig + 8] = memspace->workspace_float[2];
    TrialState->xstar[nVarOrig + 11] = 0.0;
    addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 5);
    if (memspace->workspace_float[2] <= 1.0E-6) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, mLBOrig + 2);
    }
  }

  mLBOrig = qpoptions->MaxIterations;
  qpoptions->MaxIterations = (qpoptions->MaxIterations + WorkingSet->nVar) -
    nVarOrig;
  driver(Hessian, grad, TrialState, memspace, WorkingSet, QRManager, CholManager,
         QPObjective, qpoptions, qpoptions);
  qpoptions->MaxIterations = mLBOrig;
  mLBOrig = (WorkingSet->isActiveIdx[3] + WorkingSet->sizes[3]) - 13;
  memspace->workspace_int[0] = WorkingSet->isActiveConstr[mLBOrig + 6];
  memspace->workspace_int[3] = WorkingSet->isActiveConstr[mLBOrig + 9];
  tf = WorkingSet->isActiveConstr[mLBOrig + 7];
  b_tf = WorkingSet->isActiveConstr[mLBOrig + 10];
  memspace->workspace_int[1] = tf;
  memspace->workspace_int[4] = b_tf;
  iac = ((WorkingSet->isActiveConstr[mLBOrig + 6] + WorkingSet->
          isActiveConstr[mLBOrig + 9]) + tf) + b_tf;
  tf = WorkingSet->isActiveConstr[mLBOrig + 8];
  b_tf = WorkingSet->isActiveConstr[mLBOrig + 11];
  memspace->workspace_int[2] = tf;
  memspace->workspace_int[5] = b_tf;
  iac = (iac + tf) + b_tf;
  for (idx = 0; idx < 6; idx++) {
    tf = WorkingSet->isActiveConstr[mLBOrig + idx];
    memspace->workspace_int[idx + 6] = tf;
    iac += tf;
  }

  if (TrialState->state != -6) {
    double penaltyParamTrial;
    mLBOrig = nVarOrig + 1;
    smax = 0.0;
    if (21 - nVarOrig >= 1) {
      for (idx = mLBOrig; idx < 22; idx++) {
        smax += fabs(TrialState->xstar[idx - 1]);
      }
    }

    s = 0.0;
    if (21 - nVarOrig >= 1) {
      idx_max = (unsigned char)(21 - nVarOrig);
      for (idx = 0; idx < idx_max; idx++) {
        d = TrialState->xstar[nVarOrig + idx];
        s += d * d;
      }
    }

    beta = (TrialState->fstar - rho * smax) - beta / 2.0 * s;
    penaltyParamTrial = MeritFunction->penaltyParam;
    smax = 0.0;
    for (idx = 0; idx < 6; idx++) {
      d = TrialState->cIneq[idx];
      if (d > 0.0) {
        smax += d;
      }
    }

    rho = ((fabs(TrialState->cEq[0]) + fabs(TrialState->cEq[1])) + fabs
           (TrialState->cEq[2])) + smax;
    smax = MeritFunction->linearizedConstrViol;
    s = 0.0;
    if (21 - nVarOrig >= 1) {
      for (idx = mLBOrig; idx < 22; idx++) {
        s += fabs(TrialState->xstar[idx - 1]);
      }
    }

    MeritFunction->linearizedConstrViol = s;
    smax = (rho + smax) - s;
    if ((smax > 2.2204460492503131E-16) && (beta > 0.0)) {
      if (TrialState->sqpFval == 0.0) {
        d = 1.0;
      } else {
        d = 1.5;
      }

      penaltyParamTrial = d * beta / smax;
    }

    if (penaltyParamTrial < MeritFunction->penaltyParam) {
      MeritFunction->phi = TrialState->sqpFval + penaltyParamTrial * rho;
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
          rho;
      }
    } else {
      MeritFunction->penaltyParam = fmax(penaltyParamTrial, 1.0E-10);
      MeritFunction->phi = TrialState->sqpFval + MeritFunction->penaltyParam *
        rho;
    }

    MeritFunction->phiPrimePlus = fmin(beta - MeritFunction->penaltyParam * rho,
      0.0);
    mLBOrig = WorkingSet->isActiveIdx[1] - 1;
    for (idx = 0; idx < 3; idx++) {
      if ((memspace->workspace_int[idx] != 0) && (memspace->workspace_int[idx +
           3] != 0)) {
        tf = true;
      } else {
        tf = false;
      }

      idx_max = mLBOrig + idx;
      TrialState->lambda[idx_max] *= (double)tf;
    }

    mLBOrig = WorkingSet->isActiveIdx[2];
    idx_max = WorkingSet->nActiveConstr;
    for (idx = mLBOrig; idx <= idx_max; idx++) {
      if (WorkingSet->Wid[idx - 1] == 3) {
        TrialState->lambda[idx - 1] *= (double)memspace->
          workspace_int[WorkingSet->Wlocalidx[idx - 1] + 5];
      }
    }
  }

  mLBOrig = WorkingSet->sizes[3] - 12;
  idx = WorkingSet->nActiveConstr;
  while ((idx > 3) && (iac > 0)) {
    if ((WorkingSet->Wid[idx - 1] == 4) && (WorkingSet->Wlocalidx[idx - 1] >
         mLBOrig)) {
      idx_max = WorkingSet->nActiveConstr - 1;
      smax = TrialState->lambda[idx_max];
      TrialState->lambda[idx_max] = 0.0;
      TrialState->lambda[idx - 1] = smax;
      removeConstr(WorkingSet, idx);
      iac--;
    }

    idx--;
  }

  QPObjective->nvar = nVarOrig;
  QPObjective->hasLinear = true;
  QPObjective->objtype = 3;
  setProblemType(WorkingSet, 3);
  sortLambdaQP(TrialState->lambda, WorkingSet->nActiveConstr, WorkingSet->sizes,
               WorkingSet->isActiveIdx, WorkingSet->Wid, WorkingSet->Wlocalidx,
               memspace->workspace_float);
}

/*
 * Arguments    : e_struct_T *obj
 *                int idx_global
 * Return Type  : void
 */
static void removeConstr(e_struct_T *obj, int idx_global)
{
  int TYPE_tmp;
  int idx;
  TYPE_tmp = obj->Wid[idx_global - 1] - 1;
  obj->isActiveConstr[(obj->isActiveIdx[TYPE_tmp] + obj->Wlocalidx[idx_global -
                       1]) - 2] = false;
  if (idx_global < obj->nActiveConstr) {
    int i;
    int i1;
    i = obj->nActiveConstr - 1;
    obj->Wid[idx_global - 1] = obj->Wid[i];
    obj->Wlocalidx[idx_global - 1] = obj->Wlocalidx[i];
    i1 = (unsigned char)obj->nVar;
    for (idx = 0; idx < i1; idx++) {
      obj->ATwset[idx + 22 * (idx_global - 1)] = obj->ATwset[idx + 22 * i];
    }

    obj->bwset[idx_global - 1] = obj->bwset[i];
  }

  obj->nActiveConstr--;
  obj->nWConstr[TYPE_tmp]--;
}

/*
 * Arguments    : e_struct_T *obj
 *                int PROBLEM_TYPE
 * Return Type  : void
 */
static void setProblemType(e_struct_T *obj, int PROBLEM_TYPE)
{
  int i;
  int idx;
  int idx_col;
  switch (PROBLEM_TYPE) {
   case 3:
    {
      obj->nVar = 9;
      obj->mConstr = 27;
      if (obj->nWConstr[4] > 0) {
        int idxUpperExisting;
        idxUpperExisting = obj->isActiveIdx[4] - 2;
        for (idx = 0; idx < 9; idx++) {
          i = (idxUpperExisting + idx) + 1;
          obj->isActiveConstr[(obj->isActiveIdxNormal[4] + idx) - 1] =
            obj->isActiveConstr[i];
          obj->isActiveConstr[i] = false;
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
    obj->mConstr = 28;
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
      obj->nVar = 21;
      obj->mConstr = 39;
      for (i = 0; i < 5; i++) {
        obj->sizes[i] = obj->sizesRegularized[i];
      }

      if (obj->probType != 4) {
        int colOffsetATw;
        int i1;
        int idxUpperExisting;
        for (idx_col = 0; idx_col < 6; idx_col++) {
          idxUpperExisting = 22 * idx_col - 1;
          i = idx_col + 9;
          if (i >= 10) {
            memset(&obj->Aineq[idxUpperExisting + 10], 0, (unsigned int)(i - 9) *
                   sizeof(double));
          }

          obj->Aineq[(idx_col + idxUpperExisting) + 10] = -1.0;
          i = idx_col + 11;
          memset(&obj->Aineq[i + idxUpperExisting], 0, (unsigned int)
                 (((idxUpperExisting - i) - idxUpperExisting) + 22) * sizeof
                 (double));
        }

        for (idx_col = 0; idx_col < 3; idx_col++) {
          int colOffsetAeq;
          idxUpperExisting = idx_col + 15;
          colOffsetAeq = 22 * idx_col - 1;
          colOffsetATw = colOffsetAeq + 22 * (obj->isActiveIdx[1] - 1);
          for (idx = 0; idx < 6; idx++) {
            obj->Aeq[(idx + colOffsetAeq) + 10] = 0.0;
            obj->ATwset[(idx + colOffsetATw) + 10] = 0.0;
          }

          if (idxUpperExisting >= 16) {
            memset(&obj->Aeq[colOffsetAeq + 16], 0, (unsigned int)
                   (idxUpperExisting - 15) * sizeof(double));
            memset(&obj->ATwset[colOffsetATw + 16], 0, (unsigned int)
                   (idxUpperExisting - 15) * sizeof(double));
          }

          i = idx_col + colOffsetAeq;
          obj->Aeq[i + 16] = -1.0;
          idx = idx_col + colOffsetATw;
          obj->ATwset[idx + 16] = -1.0;
          i1 = idx_col + 17;
          if (i1 <= 18) {
            memset(&obj->Aeq[i1 + colOffsetAeq], 0, (unsigned int)
                   (((colOffsetAeq - i1) - colOffsetAeq) + 19) * sizeof(double));
            memset(&obj->ATwset[i1 + colOffsetATw], 0, (unsigned int)
                   (((colOffsetATw - i1) - colOffsetATw) + 19) * sizeof(double));
          }

          i1 = idx_col + 18;
          if (i1 >= 19) {
            memset(&obj->Aeq[colOffsetAeq + 19], 0, (unsigned int)(i1 - 18) *
                   sizeof(double));
            memset(&obj->ATwset[colOffsetATw + 19], 0, (unsigned int)(i1 - 18) *
                   sizeof(double));
          }

          obj->Aeq[i + 19] = 1.0;
          obj->ATwset[idx + 19] = 1.0;
          i = idx_col + 20;
          if (i <= 21) {
            memset(&obj->Aeq[i + colOffsetAeq], 0, (unsigned int)(((colOffsetAeq
                      - i) - colOffsetAeq) + 22) * sizeof(double));
            memset(&obj->ATwset[i + colOffsetATw], 0, (unsigned int)
                   (((colOffsetATw - i) - colOffsetATw) + 22) * sizeof(double));
          }
        }

        idxUpperExisting = 9;
        for (idx = 0; idx < 12; idx++) {
          idxUpperExisting++;
          obj->indexLB[idx + 9] = idxUpperExisting;
        }

        if (obj->nWConstr[4] > 0) {
          for (idx = 0; idx < 9; idx++) {
            obj->isActiveConstr[obj->isActiveIdxRegularized[4] + idx] =
              obj->isActiveConstr[(obj->isActiveIdx[4] + idx) - 1];
          }
        }

        i = obj->isActiveIdx[4];
        idx = obj->isActiveIdxRegularized[4];
        if (i <= idx - 1) {
          memset(&obj->isActiveConstr[i + -1], 0, (unsigned int)(idx - i) *
                 sizeof(bool));
        }

        memset(&obj->lb[9], 0, 12U * sizeof(double));
        idxUpperExisting = obj->isActiveIdx[2];
        i = obj->nActiveConstr;
        for (idx_col = idxUpperExisting; idx_col <= i; idx_col++) {
          colOffsetATw = 22 * (idx_col - 1) - 1;
          if (obj->Wid[idx_col - 1] == 3) {
            idx = obj->Wlocalidx[idx_col - 1];
            i1 = idx + 8;
            if (i1 >= 10) {
              memset(&obj->ATwset[colOffsetATw + 10], 0, (unsigned int)(i1 - 9) *
                     sizeof(double));
            }

            obj->ATwset[(idx + colOffsetATw) + 9] = -1.0;
            idx += 10;
            if (idx <= 21) {
              memset(&obj->ATwset[idx + colOffsetATw], 0, (unsigned int)
                     (((colOffsetATw - idx) - colOffsetATw) + 22) * sizeof
                     (double));
            }
          } else {
            memset(&obj->ATwset[colOffsetATw + 10], 0, 12U * sizeof(double));
          }
        }
      }

      for (i = 0; i < 6; i++) {
        obj->isActiveIdx[i] = obj->isActiveIdxRegularized[i];
      }
    }
    break;

   default:
    obj->nVar = 22;
    obj->mConstr = 40;
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
 *                const double grad[22]
 *                c_struct_T *TrialState
 *                struct_T *memspace
 *                e_struct_T *WorkingSet
 *                g_struct_T *QRManager
 *                h_struct_T *CholManager
 *                i_struct_T *QPObjective
 *                const k_struct_T *qpoptions
 * Return Type  : bool
 */
static bool soc(const double Hessian[81], const double grad[22], c_struct_T
                *TrialState, struct_T *memspace, e_struct_T *WorkingSet,
                g_struct_T *QRManager, h_struct_T *CholManager, i_struct_T
                *QPObjective, const k_struct_T *qpoptions)
{
  double c;
  int i;
  int i1;
  int idx;
  int idxIneqOffset;
  int idx_Aineq;
  int idx_Partition;
  int idx_lower;
  int idx_upper;
  int nVar_tmp;
  int nWIneq_old;
  int nWLower_old;
  int nWUpper_old;
  bool success;
  nWIneq_old = WorkingSet->nWConstr[2];
  nWLower_old = WorkingSet->nWConstr[3];
  nWUpper_old = WorkingSet->nWConstr[4];
  nVar_tmp = WorkingSet->nVar;
  i = (unsigned char)WorkingSet->nVar;
  memcpy(&TrialState->xstarsqp[0], &TrialState->xstarsqp_old[0], (unsigned int)i
         * sizeof(double));
  memcpy(&TrialState->socDirection[0], &TrialState->xstar[0], (unsigned int)i *
         sizeof(double));
  memcpy(&TrialState->lambdaStopTest[0], &TrialState->lambda[0], 40U * sizeof
         (double));
  idxIneqOffset = WorkingSet->isActiveIdx[2];
  WorkingSet->beq[0] = -TrialState->cEq[0];
  WorkingSet->beq[1] = -TrialState->cEq[1];
  WorkingSet->beq[2] = -TrialState->cEq[2];
  for (idx_Aineq = 0; idx_Aineq <= 44; idx_Aineq += 22) {
    c = 0.0;
    i1 = idx_Aineq + WorkingSet->nVar;
    for (idx_lower = idx_Aineq + 1; idx_lower <= i1; idx_lower++) {
      c += WorkingSet->Aeq[idx_lower - 1] * TrialState->searchDir[(idx_lower -
        idx_Aineq) - 1];
    }

    i1 = div_nde_s32_floor(idx_Aineq, 22);
    WorkingSet->beq[i1] += c;
  }

  WorkingSet->bwset[0] = WorkingSet->beq[0];
  WorkingSet->bwset[1] = WorkingSet->beq[1];
  WorkingSet->bwset[2] = WorkingSet->beq[2];
  for (idx = 0; idx < 6; idx++) {
    WorkingSet->bineq[idx] = -TrialState->cIneq[idx];
  }

  for (idx_Aineq = 0; idx_Aineq <= 110; idx_Aineq += 22) {
    c = 0.0;
    i1 = idx_Aineq + WorkingSet->nVar;
    for (idx_lower = idx_Aineq + 1; idx_lower <= i1; idx_lower++) {
      c += WorkingSet->Aineq[idx_lower - 1] * TrialState->searchDir[(idx_lower -
        idx_Aineq) - 1];
    }

    i1 = div_nde_s32_floor(idx_Aineq, 22);
    WorkingSet->bineq[i1] += c;
  }

  idx_Aineq = 1;
  idx_lower = 7;
  idx_upper = WorkingSet->sizes[3] + 7;
  i1 = WorkingSet->nActiveConstr;
  for (idx = idxIneqOffset; idx <= i1; idx++) {
    switch (WorkingSet->Wid[idx - 1]) {
     case 3:
      idx_Partition = idx_Aineq;
      idx_Aineq++;
      WorkingSet->bwset[idx - 1] = WorkingSet->bineq[WorkingSet->Wlocalidx[idx -
        1] - 1];
      break;

     case 4:
      idx_Partition = idx_lower;
      idx_lower++;
      break;

     default:
      idx_Partition = idx_upper;
      idx_upper++;
      break;
    }

    TrialState->workingset_old[idx_Partition - 1] = WorkingSet->Wlocalidx[idx -
      1];
  }

  memcpy(&TrialState->xstar[0], &TrialState->xstarsqp[0], (unsigned int)i *
         sizeof(double));
  driver(Hessian, grad, TrialState, memspace, WorkingSet, QRManager, CholManager,
         QPObjective, qpoptions, qpoptions);
  while ((WorkingSet->mEqRemoved > 0) && (WorkingSet->indexEqRemoved
          [WorkingSet->mEqRemoved - 1] >= 1)) {
    addAeqConstr(WorkingSet, WorkingSet->indexEqRemoved[WorkingSet->mEqRemoved -
                 1]);
    WorkingSet->mEqRemoved--;
  }

  i = (unsigned char)nVar_tmp;
  for (idx = 0; idx < i; idx++) {
    double oldDirIdx;
    c = TrialState->socDirection[idx];
    oldDirIdx = c;
    c = TrialState->xstar[idx] - c;
    TrialState->socDirection[idx] = c;
    TrialState->xstar[idx] = oldDirIdx;
  }

  success = (c_xnrm2(nVar_tmp, TrialState->socDirection) <= 2.0 * c_xnrm2
             (nVar_tmp, TrialState->xstar));
  idx_Partition = WorkingSet->sizes[3];
  c = -TrialState->cEq[0];
  WorkingSet->beq[0] = c;
  WorkingSet->bwset[0] = c;
  c = -TrialState->cEq[1];
  WorkingSet->beq[1] = c;
  WorkingSet->bwset[1] = c;
  c = -TrialState->cEq[2];
  WorkingSet->beq[2] = c;
  WorkingSet->bwset[2] = c;
  for (idx = 0; idx < 6; idx++) {
    WorkingSet->bineq[idx] = -TrialState->cIneq[idx];
  }

  if (!success) {
    idx_Aineq = WorkingSet->nWConstr[0] + WorkingSet->nWConstr[1];
    idx_lower = idx_Aineq + 1;
    i = WorkingSet->nActiveConstr;
    for (idx_upper = idx_lower; idx_upper <= i; idx_upper++) {
      WorkingSet->isActiveConstr[(WorkingSet->isActiveIdx[WorkingSet->
        Wid[idx_upper - 1] - 1] + WorkingSet->Wlocalidx[idx_upper - 1]) - 2] =
        false;
    }

    WorkingSet->nWConstr[2] = 0;
    WorkingSet->nWConstr[3] = 0;
    WorkingSet->nWConstr[4] = 0;
    WorkingSet->nActiveConstr = idx_Aineq;
    for (idx = 0; idx < nWIneq_old; idx++) {
      idx_Aineq = TrialState->workingset_old[idx];
      WorkingSet->nWConstr[2]++;
      WorkingSet->isActiveConstr[(WorkingSet->isActiveIdx[2] + idx_Aineq) - 2] =
        true;
      WorkingSet->nActiveConstr++;
      i = WorkingSet->nActiveConstr - 1;
      WorkingSet->Wid[i] = 3;
      WorkingSet->Wlocalidx[i] = idx_Aineq;
      idx_lower = 22 * (idx_Aineq - 1);
      idx_upper = 22 * i;
      i1 = WorkingSet->nVar;
      for (nVar_tmp = 0; nVar_tmp < i1; nVar_tmp++) {
        WorkingSet->ATwset[idx_upper + nVar_tmp] = WorkingSet->Aineq[idx_lower +
          nVar_tmp];
      }

      WorkingSet->bwset[i] = WorkingSet->bineq[idx_Aineq - 1];
    }

    for (idx = 0; idx < nWLower_old; idx++) {
      addBoundToActiveSetMatrix_(WorkingSet, 4, TrialState->workingset_old[idx +
        6]);
    }

    for (idx = 0; idx < nWUpper_old; idx++) {
      addBoundToActiveSetMatrix_(WorkingSet, 5, TrialState->workingset_old[(idx
        + idx_Partition) + 6]);
    }

    memcpy(&TrialState->lambda[0], &TrialState->lambdaStopTest[0], 40U * sizeof
           (double));
  } else {
    sortLambdaQP(TrialState->lambda, WorkingSet->nActiveConstr,
                 WorkingSet->sizes, WorkingSet->isActiveIdx, WorkingSet->Wid,
                 WorkingSet->Wlocalidx, memspace->workspace_float);
  }

  return success;
}

/*
 * Arguments    : const h_struct_T *obj
 *                double rhs[22]
 * Return Type  : void
 */
static void solve(const h_struct_T *obj, double rhs[22])
{
  int i;
  int j;
  int jA;
  int n_tmp;
  n_tmp = obj->ndims;
  if (obj->ndims != 0) {
    for (j = 0; j < n_tmp; j++) {
      double temp;
      jA = j * 40;
      temp = rhs[j];
      for (i = 0; i < j; i++) {
        temp -= obj->FMat[jA + i] * rhs[i];
      }

      rhs[j] = temp / obj->FMat[jA + j];
    }
  }

  if (obj->ndims != 0) {
    for (j = n_tmp; j >= 1; j--) {
      jA = (j + (j - 1) * 40) - 1;
      rhs[j - 1] /= obj->FMat[jA];
      for (i = 0; i <= j - 2; i++) {
        int ix;
        ix = (j - i) - 2;
        rhs[ix] -= rhs[j - 1] * obj->FMat[(jA - i) - 1];
      }
    }
  }
}

/*
 * Arguments    : double lambda[40]
 *                int WorkingSet_nActiveConstr
 *                const int WorkingSet_sizes[5]
 *                const int WorkingSet_isActiveIdx[6]
 *                const int WorkingSet_Wid[40]
 *                const int WorkingSet_Wlocalidx[40]
 *                double workspace[880]
 * Return Type  : void
 */
static void sortLambdaQP(double lambda[40], int WorkingSet_nActiveConstr, const
  int WorkingSet_sizes[5], const int WorkingSet_isActiveIdx[6], const int
  WorkingSet_Wid[40], const int WorkingSet_Wlocalidx[40], double workspace[880])
{
  if (WorkingSet_nActiveConstr != 0) {
    int idx;
    int idxOffset;
    int mAll;
    mAll = WorkingSet_sizes[3] + 17;
    idx = (unsigned char)(WorkingSet_sizes[3] + 18);
    memcpy(&workspace[0], &lambda[0], (unsigned int)idx * sizeof(double));
    if (mAll >= 0) {
      memset(&lambda[0], 0, (unsigned int)(mAll + 1) * sizeof(double));
    }

    mAll = 0;
    idx = 0;
    while ((idx + 1 <= WorkingSet_nActiveConstr) && (WorkingSet_Wid[idx] <= 2))
    {
      if (WorkingSet_Wid[idx] == 1) {
        idxOffset = 1;
      } else {
        idxOffset = WorkingSet_isActiveIdx[1];
      }

      lambda[(idxOffset + WorkingSet_Wlocalidx[idx]) - 2] = workspace[mAll];
      mAll++;
      idx++;
    }

    while (idx + 1 <= WorkingSet_nActiveConstr) {
      switch (WorkingSet_Wid[idx]) {
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

      lambda[(idxOffset + WorkingSet_Wlocalidx[idx]) - 2] = workspace[mAll];
      mAll++;
      idx++;
    }
  }
}

/*
 * Arguments    : g_struct_T *obj
 *                const double vec[880]
 *                int iv0
 * Return Type  : void
 */
static void squareQ_appendCol(g_struct_T *obj, const double vec[880], int iv0)
{
  double b_temp_tmp;
  double c;
  double s;
  int Qk0;
  int i;
  int idx;
  int iyend;
  int k;
  int temp_tmp;
  iyend = obj->mrows;
  Qk0 = obj->ncols + 1;
  if (iyend <= Qk0) {
    Qk0 = iyend;
  }

  obj->minRowCol = Qk0;
  Qk0 = 40 * obj->ncols;
  if (obj->mrows != 0) {
    iyend = Qk0 + obj->mrows;
    if (Qk0 + 1 <= iyend) {
      memset(&obj->QR[Qk0], 0, (unsigned int)(iyend - Qk0) * sizeof(double));
    }

    i = 40 * (obj->mrows - 1) + 1;
    for (iyend = 1; iyend <= i; iyend += 40) {
      c = 0.0;
      temp_tmp = iyend + obj->mrows;
      for (idx = iyend; idx < temp_tmp; idx++) {
        c += obj->Q[idx - 1] * vec[((iv0 + idx) - iyend) - 1];
      }

      temp_tmp = Qk0 + div_nde_s32_floor(iyend - 1, 40);
      obj->QR[temp_tmp] += c;
    }
  }

  obj->ncols++;
  i = obj->ncols - 1;
  obj->jpvt[i] = obj->ncols;
  for (idx = obj->mrows - 2; idx + 2 > obj->ncols; idx--) {
    temp_tmp = idx + 40 * i;
    b_temp_tmp = obj->QR[temp_tmp + 1];
    c = xrotg(&obj->QR[temp_tmp], &b_temp_tmp, &s);
    obj->QR[temp_tmp + 1] = b_temp_tmp;
    Qk0 = 40 * idx;
    iyend = obj->mrows;
    if (obj->mrows >= 1) {
      for (k = 0; k < iyend; k++) {
        double temp;
        temp_tmp = Qk0 + k;
        b_temp_tmp = obj->Q[temp_tmp + 40];
        temp = c * obj->Q[temp_tmp] + s * b_temp_tmp;
        obj->Q[temp_tmp + 40] = c * b_temp_tmp - s * obj->Q[temp_tmp];
        obj->Q[temp_tmp] = temp;
      }
    }
  }
}

/*
 * Arguments    : j_struct_T *stepFlags
 *                double Hessian[81]
 *                c_struct_T *TrialState
 *                f_struct_T *MeritFunction
 *                struct_T *memspace
 *                e_struct_T *WorkingSet
 *                g_struct_T *QRManager
 *                h_struct_T *CholManager
 *                i_struct_T *QPObjective
 *                k_struct_T *qpoptions
 * Return Type  : void
 */
static void step(j_struct_T *stepFlags, double Hessian[81], c_struct_T
                 *TrialState, f_struct_T *MeritFunction, struct_T *memspace,
                 e_struct_T *WorkingSet, g_struct_T *QRManager, h_struct_T
                 *CholManager, i_struct_T *QPObjective, k_struct_T *qpoptions)
{
  e_struct_T obj;
  double dv[22];
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
    idxStartIneq = (unsigned char)WorkingSet->nVar;
    memcpy(&TrialState->xstar[0], &TrialState->xstarsqp[0], (unsigned int)
           idxStartIneq * sizeof(double));
  } else if (nVar >= 0) {
    memcpy(&TrialState->searchDir[0], &TrialState->xstar[0], (unsigned int)(nVar
            + 1) * sizeof(double));
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
        driver(Hessian, TrialState->grad, TrialState, memspace, WorkingSet,
               QRManager, CholManager, QPObjective, qpoptions, qpoptions);
        if (WorkingSet->probType == 2) {
          obj = *WorkingSet;
          constrViolationIneq = c_maxConstraintViolation_AMats_(&obj,
            TrialState->xstar);
        } else {
          obj = *WorkingSet;
          constrViolationIneq = d_maxConstraintViolation_AMats_(&obj,
            TrialState->xstar);
        }

        idxStartIneq = (unsigned char)WorkingSet->sizes[3];
        for (iH0 = 0; iH0 < idxStartIneq; iH0++) {
          constrViolationIneq = fmax(constrViolationIneq, -TrialState->
            xstar[obj.indexLB[iH0] - 1] - obj.lb[obj.indexLB[iH0] - 1]);
        }

        for (iH0 = 0; iH0 < 9; iH0++) {
          constrViolationIneq = fmax(constrViolationIneq, TrialState->
            xstar[obj.indexUB[iH0] - 1] - obj.ub[obj.indexUB[iH0] - 1]);
        }

        if ((TrialState->state > 0) || ((TrialState->state == 0) &&
             (constrViolationIneq <= 1.0E-6))) {
          double constrViolation;
          double penaltyParamTrial;
          penaltyParamTrial = MeritFunction->penaltyParam;
          constrViolationIneq = 0.0;
          for (iH0 = 0; iH0 < 6; iH0++) {
            linearizedConstrViolPrev = TrialState->cIneq[iH0];
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

        sortLambdaQP(TrialState->lambda, WorkingSet->nActiveConstr,
                     WorkingSet->sizes, WorkingSet->isActiveIdx, WorkingSet->Wid,
                     WorkingSet->Wlocalidx, memspace->workspace_float);
        nonlinEqRemoved = (WorkingSet->mEqRemoved > 0);
        while ((WorkingSet->mEqRemoved > 0) && (WorkingSet->
                indexEqRemoved[WorkingSet->mEqRemoved - 1] >= 1)) {
          addAeqConstr(WorkingSet, WorkingSet->indexEqRemoved
                       [WorkingSet->mEqRemoved - 1]);
          WorkingSet->mEqRemoved--;
        }

        if (nonlinEqRemoved) {
          WorkingSet->Wlocalidx[0] = 1;
          WorkingSet->Wlocalidx[1] = 2;
          WorkingSet->Wlocalidx[2] = 3;
        }

        if (stepFlags->stepType != 2) {
          if (nVar >= 0) {
            memcpy(&TrialState->delta_x[0], &TrialState->xstar[0], (unsigned int)
                   (nVar + 1) * sizeof(double));
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
        WorkingSet->isActiveConstr[(WorkingSet->isActiveIdx[WorkingSet->
          Wid[idx_global - 1] - 1] + WorkingSet->Wlocalidx[idx_global - 1]) - 2]
          = false;
      }

      WorkingSet->nWConstr[2] = 0;
      WorkingSet->nWConstr[3] = 0;
      WorkingSet->nWConstr[4] = 0;
      WorkingSet->nActiveConstr = iH0;
      memcpy(&dv[0], &TrialState->xstar[0], 22U * sizeof(double));
      idxStartIneq = (unsigned char)WorkingSet->sizes[3];
      for (iH0 = 0; iH0 < idxStartIneq; iH0++) {
        linearizedConstrViolPrev = WorkingSet->lb[WorkingSet->indexLB[iH0] - 1];
        if (-dv[WorkingSet->indexLB[iH0] - 1] > linearizedConstrViolPrev) {
          dv[WorkingSet->indexLB[iH0] - 1] = (WorkingSet->ub[WorkingSet->
            indexLB[iH0] - 1] - linearizedConstrViolPrev) / 2.0;
        }
      }

      for (iH0 = 0; iH0 < 9; iH0++) {
        linearizedConstrViolPrev = WorkingSet->ub[WorkingSet->indexUB[iH0] - 1];
        if (dv[WorkingSet->indexUB[iH0] - 1] > linearizedConstrViolPrev) {
          dv[WorkingSet->indexUB[iH0] - 1] = (linearizedConstrViolPrev -
            WorkingSet->lb[WorkingSet->indexUB[iH0] - 1]) / 2.0;
        }
      }

      memcpy(&TrialState->xstar[0], &dv[0], 22U * sizeof(double));
      relaxed(Hessian, TrialState->grad, TrialState, MeritFunction, memspace,
              WorkingSet, QRManager, CholManager, QPObjective, qpoptions);
      if (nVar >= 0) {
        memcpy(&TrialState->delta_x[0], &TrialState->xstar[0], (unsigned int)
               (nVar + 1) * sizeof(double));
      }

      guard1 = true;
      break;

     default:
      memcpy(&dv[0], &TrialState->grad[0], 22U * sizeof(double));
      checkBoundViolation = soc(Hessian, dv, TrialState, memspace, WorkingSet,
        QRManager, CholManager, QPObjective, qpoptions);
      stepFlags->stepAccepted = checkBoundViolation;
      if (stepFlags->stepAccepted && (TrialState->state != -6)) {
        idxStartIneq = (unsigned char)(nVar + 1);
        for (iH0 = 0; iH0 < idxStartIneq; iH0++) {
          TrialState->delta_x[iH0] = TrialState->xstar[iH0] +
            TrialState->socDirection[iH0];
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
          constrViolationIneq = fmax(constrViolationIneq, fabs(TrialState->
            grad[iH0]));
          linearizedConstrViolPrev = fmax(linearizedConstrViolPrev, fabs
            (TrialState->xstar[iH0]));
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
          idxStartIneq = 7 - idxEndIneq;
          if (idxStartIneq >= 0) {
            memset(&Hessian[iH0 + 1], 0, (unsigned int)(idxStartIneq + 1) *
                   sizeof(double));
          }
        }
      }
    }
  } while (exitg1 == 0);

  if (checkBoundViolation) {
    memcpy(&dv[0], &TrialState->delta_x[0], 22U * sizeof(double));
    idxStartIneq = (unsigned char)WorkingSet->sizes[3];
    for (iH0 = 0; iH0 < idxStartIneq; iH0++) {
      constrViolationIneq = dv[WorkingSet->indexLB[iH0] - 1];
      linearizedConstrViolPrev = (TrialState->xstarsqp[WorkingSet->indexLB[iH0]
        - 1] + constrViolationIneq) - -1.0;
      if (linearizedConstrViolPrev < 0.0) {
        dv[WorkingSet->indexLB[iH0] - 1] = constrViolationIneq -
          linearizedConstrViolPrev;
        TrialState->xstar[WorkingSet->indexLB[iH0] - 1] -=
          linearizedConstrViolPrev;
      }
    }

    for (iH0 = 0; iH0 < 9; iH0++) {
      constrViolationIneq = dv[WorkingSet->indexUB[iH0] - 1];
      linearizedConstrViolPrev = ((double)iv[WorkingSet->indexUB[iH0] - 1] -
        TrialState->xstarsqp[WorkingSet->indexUB[iH0] - 1]) -
        constrViolationIneq;
      if (linearizedConstrViolPrev < 0.0) {
        dv[WorkingSet->indexUB[iH0] - 1] = constrViolationIneq +
          linearizedConstrViolPrev;
        TrialState->xstar[WorkingSet->indexUB[iH0] - 1] +=
          linearizedConstrViolPrev;
      }
    }

    memcpy(&TrialState->delta_x[0], &dv[0], 22U * sizeof(double));
  }
}

/*
 * Arguments    : const double x[9]
 *                double y[3]
 * Return Type  : void
 */
static void sum(const double x[9], double y[3])
{
  int xi;
  for (xi = 0; xi < 3; xi++) {
    int xpageoffset;
    xpageoffset = xi * 3;
    y[xi] = (x[xpageoffset] + x[xpageoffset + 1]) + x[xpageoffset + 2];
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
 * Arguments    : f_struct_T *MeritFunction
 *                const e_struct_T *WorkingSet
 *                c_struct_T *TrialState
 *                bool *Flags_fevalOK
 *                bool *Flags_done
 *                bool *Flags_stepAccepted
 *                bool *Flags_failedLineSearch
 *                int *Flags_stepType
 * Return Type  : bool
 */
static bool test_exit(f_struct_T *MeritFunction, const e_struct_T *WorkingSet,
                      c_struct_T *TrialState, bool *Flags_fevalOK, bool
                      *Flags_done, bool *Flags_stepAccepted, bool
                      *Flags_failedLineSearch, int *Flags_stepType)
{
  double s;
  double smax;
  int i;
  int idx_max;
  int k;
  int nVar;
  bool Flags_gradOK;
  bool exitg1;
  bool isFeasible;
  *Flags_fevalOK = true;
  *Flags_done = false;
  *Flags_stepAccepted = false;
  *Flags_failedLineSearch = false;
  *Flags_stepType = 1;
  nVar = WorkingSet->nVar;
  i = (unsigned char)(WorkingSet->sizes[3] + 18);
  memcpy(&TrialState->lambdaStopTest[0], &TrialState->lambdasqp[0], (unsigned
          int)i * sizeof(double));
  computeGradLag(TrialState->gradLag, WorkingSet->nVar, TrialState->grad,
                 WorkingSet->Aineq, WorkingSet->Aeq, WorkingSet->indexLB,
                 WorkingSet->sizes[3], WorkingSet->indexUB,
                 TrialState->lambdaStopTest);
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = fabs(TrialState->grad[0]);
      for (k = 2; k <= nVar; k++) {
        s = fabs(TrialState->grad[k - 1]);
        if (s > smax) {
          idx_max = k;
          smax = s;
        }
      }
    }
  }

  s = fmax(1.0, fabs(TrialState->grad[idx_max - 1]));
  if (rtIsInf(s)) {
    s = 1.0;
  }

  smax = fmax(fmax(fmax(0.0, fabs(TrialState->cEq[0])), fabs(TrialState->cEq[1])),
              fabs(TrialState->cEq[2]));
  for (idx_max = 0; idx_max < 6; idx_max++) {
    smax = fmax(smax, TrialState->cIneq[idx_max]);
  }

  nVar = (unsigned char)WorkingSet->sizes[3];
  for (idx_max = 0; idx_max < nVar; idx_max++) {
    smax = fmax(smax, -1.0 - TrialState->xstarsqp[WorkingSet->indexLB[idx_max] -
                1]);
  }

  for (idx_max = 0; idx_max < 9; idx_max++) {
    nVar = WorkingSet->indexUB[idx_max] - 1;
    smax = fmax(smax, TrialState->xstarsqp[nVar] - (double)iv[nVar]);
  }

  MeritFunction->nlpPrimalFeasError = smax;
  MeritFunction->feasRelativeFactor = fmax(1.0, smax);
  isFeasible = (smax <= 1.0E-6 * MeritFunction->feasRelativeFactor);
  Flags_gradOK = true;
  smax = 0.0;
  nVar = (unsigned char)WorkingSet->nVar;
  idx_max = 0;
  exitg1 = false;
  while ((!exitg1) && (idx_max <= nVar - 1)) {
    Flags_gradOK = ((!rtIsInf(TrialState->gradLag[idx_max])) && (!rtIsNaN
      (TrialState->gradLag[idx_max])));
    if (!Flags_gradOK) {
      exitg1 = true;
    } else {
      smax = fmax(smax, fabs(TrialState->gradLag[idx_max]));
      idx_max++;
    }
  }

  MeritFunction->nlpDualFeasError = smax;
  if (!Flags_gradOK) {
    *Flags_done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    MeritFunction->nlpComplError = 0.0;
    MeritFunction->firstOrderOpt = smax;
    memcpy(&TrialState->lambdaStopTestPrev[0], &TrialState->lambdaStopTest[0],
           (unsigned int)i * sizeof(double));
    if (isFeasible && (smax <= 1.0E-6 * s)) {
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
 *                e_struct_T *WorkingSet
 *                const double cIneq[6]
 *                const double cEq[3]
 *                int mLB
 * Return Type  : void
 */
static void updateWorkingSetForNewQP(const double xk[9], e_struct_T *WorkingSet,
  const double cIneq[6], const double cEq[3], int mLB)
{
  int b_i;
  int i;
  int iEq0;
  int idx;
  int iw0;
  iw0 = 0;
  iEq0 = 0;
  for (idx = 0; idx < 3; idx++) {
    double d;
    d = cEq[idx];
    WorkingSet->beq[idx] = -d;
    WorkingSet->bwset[idx] = -d;
    i = WorkingSet->nVar;
    for (b_i = 0; b_i < i; b_i++) {
      WorkingSet->ATwset[iw0 + b_i] = WorkingSet->Aeq[iEq0 + b_i];
    }

    iEq0 = iw0 + 22;
    iw0 += 22;
  }

  for (idx = 0; idx < 6; idx++) {
    WorkingSet->bineq[idx] = -cIneq[idx];
  }

  i = (unsigned char)mLB;
  for (idx = 0; idx < i; idx++) {
    WorkingSet->lb[WorkingSet->indexLB[idx] - 1] = xk[WorkingSet->indexLB[idx] -
      1] + 1.0;
  }

  for (idx = 0; idx < 9; idx++) {
    WorkingSet->ub[WorkingSet->indexUB[idx] - 1] = (double)iv
      [WorkingSet->indexUB[idx] - 1] - xk[WorkingSet->indexUB[idx] - 1];
  }

  if (WorkingSet->nActiveConstr > 3) {
    i = WorkingSet->nActiveConstr;
    for (idx = 4; idx <= i; idx++) {
      switch (WorkingSet->Wid[idx - 1]) {
       case 4:
        WorkingSet->bwset[idx - 1] = WorkingSet->lb[WorkingSet->
          indexLB[WorkingSet->Wlocalidx[idx - 1] - 1] - 1];
        break;

       case 5:
        WorkingSet->bwset[idx - 1] = WorkingSet->ub[WorkingSet->
          indexUB[WorkingSet->Wlocalidx[idx - 1] - 1] - 1];
        break;

       default:
        WorkingSet->bwset[idx - 1] = WorkingSet->bineq[WorkingSet->Wlocalidx[idx
          - 1] - 1];
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
 *                const double B[1600]
 *                int ib0
 *                double C[880]
 * Return Type  : void
 */
static void xgemm(int m, int n, int k, const double A[81], int lda, const double
                  B[1600], int ib0, double C[880])
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
    lastColC = 40 * (n - 1);
    for (cr = 0; cr <= lastColC; cr += 40) {
      i = cr + 1;
      i1 = cr + m;
      if (i <= i1) {
        memset(&C[i + -1], 0, (unsigned int)((i1 - i) + 1) * sizeof(double));
      }
    }

    for (cr = 0; cr <= lastColC; cr += 40) {
      int ar;
      ar = -1;
      i = br + k;
      for (ib = br; ib < i; ib++) {
        int i2;
        i1 = cr + 1;
        i2 = cr + m;
        for (ic = i1; ic <= i2; ic++) {
          C[ic - 1] += B[ib - 1] * A[(ar + ic) - cr];
        }

        ar += lda;
      }

      br += 40;
    }
  }
}

/*
 * Arguments    : int m
 *                const double A[66]
 *                const double x[880]
 *                double y[40]
 * Return Type  : void
 */
static void xgemv(int m, const double A[66], const double x[880], double y[40])
{
  int ia;
  int iac;
  y[0] = -y[0];
  y[1] = -y[1];
  y[2] = -y[2];
  for (iac = 0; iac <= 44; iac += 22) {
    double c;
    int i;
    c = 0.0;
    i = iac + m;
    for (ia = iac + 1; ia <= i; ia++) {
      c += A[ia - 1] * x[(ia - iac) - 1];
    }

    i = div_nde_s32_floor(iac, 22);
    y[i] += c;
  }
}

/*
 * Arguments    : double A[1600]
 *                int m
 *                int n
 *                int jpvt[40]
 *                double tau[40]
 * Return Type  : void
 */
static void xgeqp3(double A[1600], int m, int n, int jpvt[40], double tau[40])
{
  double vn1[40];
  double vn2[40];
  double work[40];
  double temp;
  int b_i;
  int k;
  int minmn_tmp;
  int pvt;
  if (m <= n) {
    minmn_tmp = m;
  } else {
    minmn_tmp = n;
  }

  memset(&tau[0], 0, 40U * sizeof(double));
  if (minmn_tmp < 1) {
    for (pvt = 0; pvt < n; pvt++) {
      jpvt[pvt] = pvt + 1;
    }
  } else {
    int i;
    int ix;
    int iy;
    int nfxd;
    int temp_tmp;
    nfxd = 0;
    for (pvt = 0; pvt < n; pvt++) {
      if (jpvt[pvt] != 0) {
        nfxd++;
        if (pvt + 1 != nfxd) {
          ix = pvt * 40;
          iy = (nfxd - 1) * 40;
          for (k = 0; k < m; k++) {
            temp_tmp = ix + k;
            temp = A[temp_tmp];
            i = iy + k;
            A[temp_tmp] = A[i];
            A[i] = temp;
          }

          jpvt[pvt] = jpvt[nfxd - 1];
          jpvt[nfxd - 1] = pvt + 1;
        } else {
          jpvt[pvt] = pvt + 1;
        }
      } else {
        jpvt[pvt] = pvt + 1;
      }
    }

    if (nfxd > minmn_tmp) {
      nfxd = minmn_tmp;
    }

    qrf(A, m, n, nfxd, tau);
    if (nfxd < minmn_tmp) {
      double d;
      memset(&work[0], 0, 40U * sizeof(double));
      memset(&vn1[0], 0, 40U * sizeof(double));
      memset(&vn2[0], 0, 40U * sizeof(double));
      i = nfxd + 1;
      for (pvt = i; pvt <= n; pvt++) {
        d = b_xnrm2(m - nfxd, A, (nfxd + (pvt - 1) * 40) + 1);
        vn1[pvt - 1] = d;
        vn2[pvt - 1] = d;
      }

      for (b_i = i; b_i <= minmn_tmp; b_i++) {
        double s;
        int ii;
        int ip1;
        int mmi;
        int nmi;
        ip1 = b_i + 1;
        iy = (b_i - 1) * 40;
        ii = (iy + b_i) - 1;
        nmi = (n - b_i) + 1;
        mmi = m - b_i;
        if (nmi < 1) {
          nfxd = -2;
        } else {
          nfxd = -1;
          if (nmi > 1) {
            temp = fabs(vn1[b_i - 1]);
            for (k = 2; k <= nmi; k++) {
              s = fabs(vn1[(b_i + k) - 2]);
              if (s > temp) {
                nfxd = k - 2;
                temp = s;
              }
            }
          }
        }

        pvt = b_i + nfxd;
        if (pvt + 1 != b_i) {
          ix = pvt * 40;
          for (k = 0; k < m; k++) {
            temp_tmp = ix + k;
            temp = A[temp_tmp];
            nfxd = iy + k;
            A[temp_tmp] = A[nfxd];
            A[nfxd] = temp;
          }

          nfxd = jpvt[pvt];
          jpvt[pvt] = jpvt[b_i - 1];
          jpvt[b_i - 1] = nfxd;
          vn1[pvt] = vn1[b_i - 1];
          vn2[pvt] = vn2[b_i - 1];
        }

        if (b_i < m) {
          temp = A[ii];
          d = xzlarfg(mmi + 1, &temp, A, ii + 2);
          tau[b_i - 1] = d;
          A[ii] = temp;
        } else {
          d = 0.0;
          tau[b_i - 1] = 0.0;
        }

        if (b_i < n) {
          temp = A[ii];
          A[ii] = 1.0;
          xzlarf(mmi + 1, nmi - 1, ii + 1, d, A, ii + 41, work);
          A[ii] = temp;
        }

        for (pvt = ip1; pvt <= n; pvt++) {
          nfxd = b_i + (pvt - 1) * 40;
          d = vn1[pvt - 1];
          if (d != 0.0) {
            temp = fabs(A[nfxd - 1]) / d;
            temp = 1.0 - temp * temp;
            if (temp < 0.0) {
              temp = 0.0;
            }

            s = d / vn2[pvt - 1];
            s = temp * (s * s);
            if (s <= 1.4901161193847656E-8) {
              if (b_i < m) {
                d = b_xnrm2(mmi, A, nfxd + 1);
                vn1[pvt - 1] = d;
                vn2[pvt - 1] = d;
              } else {
                vn1[pvt - 1] = 0.0;
                vn2[pvt - 1] = 0.0;
              }
            } else {
              vn1[pvt - 1] = d * sqrt(temp);
            }
          }
        }
      }
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
 *                double A[1600]
 * Return Type  : int
 */
static int xpotrf(int n, double A[1600])
{
  int ia;
  int iac;
  int info;
  int j;
  int nmj;
  bool exitg1;
  info = 0;
  j = 0;
  exitg1 = false;
  while ((!exitg1) && (j <= n - 1)) {
    double c;
    double ssq;
    int idxA1j;
    int idxAjj;
    idxA1j = j * 40;
    idxAjj = idxA1j + j;
    ssq = 0.0;
    if (j >= 1) {
      for (nmj = 0; nmj < j; nmj++) {
        c = A[idxA1j + nmj];
        ssq += c * c;
      }
    }

    ssq = A[idxAjj] - ssq;
    if (ssq > 0.0) {
      ssq = sqrt(ssq);
      A[idxAjj] = ssq;
      if (j + 1 < n) {
        int i;
        int ia0;
        int idxAjjp1;
        nmj = (n - j) - 2;
        ia0 = idxA1j + 41;
        idxAjjp1 = idxAjj + 41;
        if ((j != 0) && (nmj + 1 != 0)) {
          i = (idxA1j + 40 * nmj) + 41;
          for (iac = ia0; iac <= i; iac += 40) {
            int i1;
            c = 0.0;
            i1 = iac + j;
            for (ia = iac; ia < i1; ia++) {
              c += A[ia - 1] * A[(idxA1j + ia) - iac];
            }

            i1 = (idxAjj + div_nde_s32_floor((iac - idxA1j) - 41, 40) * 40) + 40;
            A[i1] -= c;
          }
        }

        ssq = 1.0 / ssq;
        i = (idxAjj + 40 * nmj) + 41;
        for (nmj = idxAjjp1; nmj <= i; nmj += 40) {
          A[nmj - 1] *= ssq;
        }
      }

      j++;
    } else {
      A[idxAjj] = ssq;
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
 *                double C[1600]
 *                int ic0
 *                double work[40]
 * Return Type  : void
 */
static void xzlarf(int m, int n, int iv0, double tau, double C[1600], int ic0,
                   double work[40])
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
    while ((lastv > 0) && (C[i - 2] == 0.0)) {
      lastv--;
      i--;
    }

    lastc = n - 1;
    exitg2 = false;
    while ((!exitg2) && (lastc + 1 > 0)) {
      int exitg1;
      i = ic0 + lastc * 40;
      ia = i;
      do {
        exitg1 = 0;
        if (ia <= (i + lastv) - 1) {
          if (C[ia - 1] != 0.0) {
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
    lastc = -1;
  }

  if (lastv > 0) {
    double c;
    int b_i;
    if (lastc + 1 != 0) {
      if (lastc >= 0) {
        memset(&work[0], 0, (unsigned int)(lastc + 1) * sizeof(double));
      }

      b_i = ic0 + 40 * lastc;
      for (iac = ic0; iac <= b_i; iac += 40) {
        c = 0.0;
        i = iac + lastv;
        for (ia = iac; ia < i; ia++) {
          c += C[ia - 1] * C[((iv0 + ia) - iac) - 1];
        }

        i = div_nde_s32_floor(iac - ic0, 40);
        work[i] += c;
      }
    }

    if (!(-tau == 0.0)) {
      i = ic0;
      for (iac = 0; iac <= lastc; iac++) {
        c = work[iac];
        if (c != 0.0) {
          c *= -tau;
          b_i = lastv + i;
          for (ia = i; ia < b_i; ia++) {
            C[ia - 1] += C[((iv0 + ia) - i) - 1] * c;
          }
        }

        i += 40;
      }
    }
  }
}

/*
 * Arguments    : int n
 *                double *alpha1
 *                double x[1600]
 *                int ix0
 * Return Type  : double
 */
static double xzlarfg(int n, double *alpha1, double x[1600], int ix0)
{
  double tau;
  int k;
  tau = 0.0;
  if (n > 0) {
    double xnorm;
    xnorm = b_xnrm2(n - 1, x, ix0);
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
            x[k - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while ((fabs(xnorm) < 1.0020841800044864E-292) && (knt < 20));

        a_tmp = fabs(*alpha1);
        xnorm = fabs(b_xnrm2(n - 1, x, ix0));
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
          x[k - 1] *= a_tmp;
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
          x[k - 1] *= a_tmp;
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
 * function [best_result, best_radius, best_exitflag] = mso_forces(center, T_min, T_max)
 *
 * MultiStart Optimize Forces
 *  INPUT
 *  center (3x1): 球心
 *  T_min (Nx1): 最小张力
 *  T_max (Nx1): 最大张力
 *
 * Arguments    : const double center[3]
 *                const double T_min[3]
 *                const double T_max[3]
 *                double best_result[9]
 *                double *best_radius
 *                double *best_exitflag
 * Return Type  : void
 */
void mso_forces(const double center[3], const double T_min[3], const double
                T_max[3], double best_result[9], double *best_radius, double
                *best_exitflag)
{
  double A[54];
  double qxy_lb[9];
  double qxy_ub[9];
  double b_expl_temp;
  double c_expl_temp;
  double d_expl_temp;
  double e_expl_temp;
  double f_expl_temp;
  double g_expl_temp;
  double q_tmp;
  int A_tmp;
  int i;

  /* 'mso_forces:8' if ~isvector(T_max) || ~isvector(T_min) */
  /* 'mso_forces:12' N = length(T_max); */
  /*  不等式约束 */
  /* 'mso_forces:42' psi = deg2rad(360) / N * (0:N-1)'; */
  /* 'mso_forces:43' range = deg2rad(180 / N); */
  /* 'mso_forces:44' qxy_ub = [sin(psi + range) -cos(psi + range) zeros(N, 1)]; */
  /* 'mso_forces:45' qxy_lb = [sin(psi - range) -cos(psi - range) zeros(N, 1)]; */
  qxy_ub[0] = 0.8660254037844386;
  qxy_ub[3] = -0.50000000000000011;
  qxy_ub[6] = 0.0;
  qxy_lb[0] = -0.8660254037844386;
  qxy_lb[3] = -0.50000000000000011;
  qxy_lb[6] = 0.0;
  qxy_ub[1] = 1.2246467991473532E-16;
  qxy_ub[4] = 1.0;
  qxy_ub[7] = 0.0;
  qxy_lb[1] = 0.8660254037844386;
  qxy_lb[4] = -0.50000000000000011;
  qxy_lb[7] = 0.0;
  qxy_ub[2] = -0.866025403784439;
  qxy_ub[5] = -0.49999999999999933;
  qxy_ub[8] = 0.0;
  qxy_lb[2] = 1.2246467991473532E-16;
  qxy_lb[5] = 1.0;
  qxy_lb[8] = 0.0;

  /* 'mso_forces:46' A = zeros([N * 2, N * 3]); */
  memset(&A[0], 0, 54U * sizeof(double));

  /* 'mso_forces:47' for i = 1:N */
  for (i = 0; i < 3; i++) {
    signed char b_i;

    /* 'mso_forces:48' A(i, 3*i-2:3*i) = -qxy_ub(i, :); */
    b_i = (signed char)(3 * (i + 1));

    /* 'mso_forces:49' A(i + N, 3*i-2:3*i) = qxy_lb(i, :); */
    A_tmp = i + 6 * (b_i - 3);
    A[A_tmp] = -qxy_ub[i];
    A[A_tmp + 3] = qxy_lb[i];
    A_tmp = i + 6 * (b_i - 2);
    A[A_tmp] = -qxy_ub[i + 3];
    A[A_tmp + 3] = qxy_lb[i + 3];
    A_tmp = i + 6 * (b_i - 1);
    A[A_tmp] = -0.0;
    A[A_tmp + 3] = 0.0;
  }

  /* 'mso_forces:51' b = zeros(N * 2, 1); */
  /*  等式约束     */
  /* 'mso_forces:54' Aeq = []; */
  /* 'mso_forces:55' beq = []; */
  /*  上下限 */
  /* 'mso_forces:58' ub = repmat([1 1 0], 1, N); */
  /* 'mso_forces:59' lb = -ones(1, 3 * N); */
  /*  目标函数 */
  /* 'mso_forces:62' f_objective = @(W) objective(W, center, N, T_max, T_min); */
  /*  优化选项 */
  /* 'mso_forces:65' options = optimoptions('fmincon', 'Display', 'off', ... */
  /* 'mso_forces:66'                         'Algorithm', 'sqp', 'SpecifyConstraintGradient', true, ... */
  /* 'mso_forces:67'                         'MaxIterations', 100); */
  /*  MultiStart */
  /* 'mso_forces:70' theta_list = deg2rad([10, 20, 30, 40, 50, 60, 70]); */
  /* 'mso_forces:71' q = zeros([3 * N, 1]); */
  /* 'mso_forces:73' best_radius = -nan; */
  *best_radius = rtNaN;

  /* 'mso_forces:74' best_result = zeros(size(q)); */
  memset(&best_result[0], 0, 9U * sizeof(double));

  /* 'mso_forces:75' best_exitflag = 0; */
  *best_exitflag = 0.0;

  /* 'mso_forces:77' for i = 1:length(theta_list) */
  qxy_ub[1] = -0.0;
  for (i = 0; i < 7; i++) {
    double b_qxy_ub[9];
    double q_tmp_tmp;

    /* 'mso_forces:78' for j = 1:N */
    q_tmp_tmp = 0.17453292519943295 * (double)i + 0.17453292519943295;
    q_tmp = cos(q_tmp_tmp);
    q_tmp_tmp = sin(q_tmp_tmp);

    /* 'mso_forces:79' q(3*j-2:3*j) = [-cos(psi(j)) * cos(theta_list(i)) -sin(psi(j)) * cos(theta_list(i)) sin(theta_list(i))]; */
    qxy_ub[0] = -q_tmp;
    qxy_ub[2] = q_tmp_tmp;

    /* 'mso_forces:79' q(3*j-2:3*j) = [-cos(psi(j)) * cos(theta_list(i)) -sin(psi(j)) * cos(theta_list(i)) sin(theta_list(i))]; */
    qxy_ub[3] = 0.49999999999999978 * q_tmp;
    qxy_ub[4] = -0.86602540378443871 * q_tmp;
    qxy_ub[5] = q_tmp_tmp;

    /* 'mso_forces:79' q(3*j-2:3*j) = [-cos(psi(j)) * cos(theta_list(i)) -sin(psi(j)) * cos(theta_list(i)) sin(theta_list(i))]; */
    qxy_ub[6] = 0.50000000000000044 * q_tmp;
    qxy_ub[7] = 0.86602540378443837 * q_tmp;
    qxy_ub[8] = q_tmp_tmp;

    /*  优化求解 */
    /* 'mso_forces:83' [result, maxfval, exitflag, output] = fmincon(f_objective, -q, A, b, Aeq, beq, ... */
    /* 'mso_forces:84'                         lb, ub, @(W) constraint(W), options); */
    for (A_tmp = 0; A_tmp < 9; A_tmp++) {
      b_qxy_ub[A_tmp] = -qxy_ub[A_tmp];
    }

    char expl_temp[3];
    q_tmp_tmp = fmincon(center, T_max, T_min, b_qxy_ub, A, qxy_lb, &q_tmp,
                        &b_expl_temp, &c_expl_temp, expl_temp, &d_expl_temp,
                        &e_expl_temp, &f_expl_temp, &g_expl_temp);

    /*  输出结果 */
    /*  disp_result(i, result, maxfval, exitflag, output); */
    /*  保存最优结果 */
    /* 'mso_forces:90' if ~isfinite(best_radius) || best_radius < -maxfval */
    if (rtIsInf(*best_radius) || rtIsNaN(*best_radius) || (*best_radius <
         -q_tmp_tmp)) {
      /* 'mso_forces:91' best_radius = -maxfval; */
      *best_radius = -q_tmp_tmp;

      /* 'mso_forces:92' best_result = -result; */
      for (A_tmp = 0; A_tmp < 9; A_tmp++) {
        best_result[A_tmp] = -qxy_lb[A_tmp];
      }

      /* 'mso_forces:93' best_exitflag = exitflag; */
      *best_exitflag = q_tmp;
    }
  }
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void mso_forces_initialize(void)
{
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void mso_forces_terminate(void)
{
}

/*
 * File trailer for mso_forces.c
 *
 * [EOF]
 */
