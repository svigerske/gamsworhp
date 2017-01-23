/*******************************************************************************
* GAMS / WORHP Interface
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "gmomcc.h"
#include "gevmcc.h"

#include "worhp.h"

/* Objective function */
void UserF(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int objMinMaxFac);
/* Function of constraints */
void UserG(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt);
/* Gradient structure */
int StructDF(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt);
/* Gradient of objective function */
void UserDF(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int objMinMaxFac);
/* Jacobian structure */
int StructDG(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** DGrowStartInit, int** DGcolInit, double** DGvalInit, double** DGvalDenseInit);
/* Jacobian of constraints */
void UserDG(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** DGrowStartInit, int** DGcolInit, double** DGvalInit, double** DGvalDenseInit);
/* Hessian structure */
int StructHM(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** HMrowInit, int** HMcolInit, double** HMvalInit, int** HMpermInit, int* HMdimMiss);
/* Hessian of Lagrangian */
void UserHM(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** HMrowInit, int** HMcolInit, double** HMvalInit, int** HMpermInit, int* HMdimMiss, int objMinMaxFac);

int main(int argc, char** argv)
{
   gmoHandle_t gmo = NULL;
   gevHandle_t gev = NULL;
   int rc = EXIT_FAILURE;
   char buffer[1024];
   int status;
   int* HMrowInit;
   int* HMcolInit;
   int* HMpermInit;
   double* HMvalInit;
   int HMdimMiss;
   int* DGcolInit;
   int* DGrowStartInit;
   double* DGvalInit;
   double* DGvalDenseInit;
   double clockStart;
   double objMinMaxFac;

   /* WORHP data structures */
   OptVar    opt;
   Workspace wsp;
   Params    par;
   Control   cnt;

   if (argc < 2)
   {
      printf("usage: %s <cntrlfile>\n", argv[0]);
      return 1;
   }

   /*
   * GAMS initialize GMO and GEV libraries
   */
   if (!gmoCreate(&gmo, buffer, sizeof(buffer)) || !gevCreate(&gev, buffer, sizeof(buffer)))
   {
      fprintf(stderr, "%s\n", buffer);
      goto TERMINATE;
   }

   /*
   * GAMS load control file
   */
   if (gevInitEnvironmentLegacy(gev, argv[1]))
   {
      fprintf(stderr, "Could not load control file %s\n", argv[1]);
      goto TERMINATE;
   }

   /*
   * GAMS let gmo know about gev
   */
   if (gmoRegisterEnvironment(gmo, gev, buffer))
   {
      fprintf(stderr, "Error registering GAMS Environment: %s\n", buffer);
      goto TERMINATE;
   }

   /*
   * GAMS load instance data
   */
   if (gmoLoadDataLegacy(gmo, buffer))
   {
      fprintf(stderr, "Could not load model data.\n");
      goto TERMINATE;
   }

   sprintf(buffer, "This is WORHP %d.%d." WORHP_PATCH ".", WORHP_MAJOR, WORHP_MINOR);
   gevLogStat(gev, buffer);

   /*
   * GAMS general
   */
   gmoObjStyleSet(gmo, gmoObjType_Fun);
   gmoObjReformSet(gmo, 1);
   gmoIndexBaseSet(gmo, 0);
   objMinMaxFac = (gmoSense(gmo) == gmoObj_Max) ? -1.0 : 1.0;

   /*
   * WORHP check Version
   */
   CHECK_WORHP_VERSION

   /*
   * WORHP pre-initialization
   */
   WorhpPreInit(&opt, &wsp, &par, &cnt);

   /*
   * WORHP parameter initialization
   */
   InitParams(&status, &par);
   WorhpSetIntParam(&par, "MaxIter", gevGetIntOpt(gev, gevIterLim));
   WorhpSetDoubleParam(&par, "Timeout", gevGetDblOpt(gev, gevResLim));
   if (gmoNLM(gmo) == 0) {
      par.ScaledKKT = false;
      par.UserHM = false;
      par.BFGSmethod = 2;
      par.BFGSmaxblockSize = 1;
      par.BFGSminblockSize = 1;
   }

   /*
   * problem size
   */
   opt.n = gmoN(gmo);
   opt.m = gmoM(gmo);
   wsp.DF.nnz = gmoObjNZ(gmo);
   wsp.DG.nnz = WorhpMatrix_Dont_Allocate;
   wsp.HM.nnz = WorhpMatrix_Dont_Allocate;

   /*
   * WORHP data structure initialization
   */
   WorhpInit(&opt, &wsp, &par, &cnt);
   if (cnt.status != FirstCall)
   {
      gevLogStat(gev, "Error: WORHP Initialisation failed.");
      return EXIT_FAILURE;
   }

   /*
   * initial point
   */
   gmoGetVarL(gmo, opt.X);
   gmoGetEquM(gmo, opt.Mu);
   gmoGetVarM(gmo, opt.Lambda);

   for (int i = 0; i < opt.m; ++i)
      opt.Mu[i] *= -1.0;

   /*
   * bounds
   */
   gmoGetVarLower(gmo, opt.XL);
   gmoGetVarUpper(gmo, opt.XU);
   gmoGetRhs(gmo, opt.GU);

   for (int i = 0; i < opt.m; ++i)
   {
      switch (gmoGetEquTypeOne(gmo, i))
      {
         case gmoequ_E:
            opt.GL[i] = opt.GU[i];
            break;
         case gmoequ_G:
            opt.GL[i] = opt.GU[i];
            opt.GU[i] = par.Infty;
            break;
         case gmoequ_L:
            opt.GL[i] = -par.Infty;
            break;
         case gmoequ_N:
            opt.GL[i] = -par.Infty;
            opt.GU[i] = par.Infty;
            break;
         case gmoequ_C:
            gevLogStat(gev, "Error: Conic constraints not supported.");
            return EXIT_FAILURE;
            break;
         case gmoequ_B:
            gevLogStat(gev, "Error: Logic constraints not supported.");
            return EXIT_FAILURE;
            break;
         default:
            gevLogStat(gev, "Error: Unsupported equation type.");
            return EXIT_FAILURE;
            break;
      }
   }

   /*
   * WORHP set structure of gradient
   */
   if (wsp.DF.NeedStructure)
      StructDF(&gmo, &gev, &opt, &wsp, &par, &cnt);

   /*
   * WORHP set structure of jacobian
   */
   StructDG(&gmo, &gev, &opt, &wsp, &par, &cnt, &DGrowStartInit, &DGcolInit, &DGvalInit, &DGvalDenseInit);

   /*
   * WORHP set structure of hessian
   */
   if (par.UserHM || par.FidifHM || par.BFGSmethod > 1)
      StructHM(&gmo, &gev, &opt, &wsp, &par, &cnt, &HMrowInit, &HMcolInit, &HMvalInit, &HMpermInit, &HMdimMiss);

   /*
   * WORHP Reverse Communication loop.
   */
   clockStart = gevTimeDiffStart(gev);
   while (cnt.status < TerminateSuccess && cnt.status > TerminateError)
   {
      /*
      * WORHP's main routine.
      */
      if (GetUserAction(&cnt, callWorhp))
      {
         Worhp(&opt, &wsp, &par, &cnt);
         /* No DoneUserAction! */
      }

      /*
      * Show iteration output.
      */
      if (GetUserAction(&cnt, iterOutput))
      {
         IterationOutput(&opt, &wsp, &par, &cnt);
         DoneUserAction(&cnt, iterOutput);
      }

      /*
      * Evaluate the objective function.
      */
      if (GetUserAction(&cnt, evalF))
      {
         UserF(&gmo, &gev, &opt, &wsp, &par, &cnt, objMinMaxFac);
         DoneUserAction(&cnt, evalF);
      }

      /*
      * Evaluate the constraints.
      */
      if (GetUserAction(&cnt, evalG))
      {
         UserG(&gmo, &gev, &opt, &wsp, &par, &cnt);
         DoneUserAction(&cnt, evalG);
      }

      /*
      * Evaluate the gradient of the objective function.
      */
      if (GetUserAction(&cnt, evalDF))
      {
         UserDF(&gmo, &gev, &opt, &wsp, &par, &cnt, objMinMaxFac);
         DoneUserAction(&cnt, evalDF);
      }

      /*
      * Evaluate the Jacobian of the constraints.
      */
      if (GetUserAction(&cnt, evalDG))
      {
         UserDG(&gmo, &gev, &opt, &wsp, &par, &cnt, &DGrowStartInit, &DGcolInit, &DGvalInit, &DGvalDenseInit);
         DoneUserAction(&cnt, evalDG);
      }

      /*
      * Evaluate the Hessian matrix of the Lagrange function (L = f + mu*g)
      */
      if (GetUserAction(&cnt, evalHM))
      {
         UserHM(&gmo, &gev, &opt, &wsp, &par, &cnt, &HMrowInit, &HMcolInit, &HMvalInit, &HMpermInit, &HMdimMiss, objMinMaxFac);
         DoneUserAction(&cnt, evalHM);
      }

      /*
      * Use finite differences with RC to determine derivatives
      */
      if (GetUserAction(&cnt, fidif))
      {
         WorhpFidif(&opt, &wsp, &par, &cnt);
         /* No DoneUserAction! */
      }
   }

   /*
   * WORHP translate status flag into a meaningful message.
   */
   StatusMsg(&opt, &wsp, &par, &cnt);

   /*
   * GAMS set solution
   */
   gmoSetHeadnTail(gmo, gmoHiterused, wsp.MajorIter);
   gmoSetHeadnTail(gmo, gmoHresused, gevTimeDiffStart(gev) - clockStart);
   /* swap sign of constraints marginals if minimization */
   if (gmoSense(gmo) == gmoObj_Min)
      for (int i = 0; i < gmoM(gmo); ++i)
         opt.Mu[i] = -opt.Mu[i];
   gmoSetSolution(gmo, opt.X, opt.Lambda, opt.Mu, opt.G);
   switch (cnt.status)
   {
      /* successful terminations */
      case OptimalSolution:
      case SearchDirectionZero:
      case SearchDirectionSmall:
      case StationaryPointFound:
      case AcceptableSolution:
      case AcceptablePrevious:
      case FritzJohn:
      case NotDiffable:
      case LowPassFilterOptimal:
      case LowPassFilterAcceptable:
      case OptimalSolutionBoxEqual:
      case AcceptableSolutionConstantF:
      case AcceptablePreviousConstantF:
      case OptimalSolutionConstantF:
      case AcceptableSolutionSKKT:
#if WORHP_MAJOR >= 2
      case AcceptableSolutionScaled:
      case AcceptablePreviousScaled:
#endif
         gmoModelStatSet(gmo, (gmoObjNLNZ(gmo) || gmoNLNZ(gmo)) ? gmoModelStat_OptimalLocal : gmoModelStat_OptimalGlobal);
         gmoSolveStatSet(gmo, gmoSolveStat_Normal);
         break;
      case FeasibleSolution:
         gmoModelStatSet(gmo, gmoModelStat_Feasible);
         gmoSolveStatSet(gmo, gmoSolveStat_Normal);
         break;

      case Unbounded:
         gmoModelStatSet(gmo, gmoModelStat_Unbounded);
         gmoSolveStatSet(gmo, gmoSolveStat_Normal);
         break;

      /* unsuccessful terminations */
      case InitError:
      case DataError:
      case QPerror:
      case FDError:
      case LinearSolverFailed:
         gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
         gmoSolveStatSet(gmo, gmoSolveStat_InternalErr);
         break;

      case MaxCalls:
      case MaxIter:
#if WORHP_MAJOR >= 2
      case MaxIterUnscaled:
#endif
         gmoModelStatSet(gmo, wsp.Feasible ? gmoModelStat_Feasible : gmoModelStat_InfeasibleIntermed);
         gmoSolveStatSet(gmo, gmoSolveStat_Iteration);
         break;
      case Timeout:
         gmoModelStatSet(gmo, wsp.Feasible ? gmoModelStat_Feasible : gmoModelStat_InfeasibleIntermed);
         gmoSolveStatSet(gmo, gmoSolveStat_Resource);
         break;

      case TerminatedByUser:
      case TerminatedByCheckFD:
         gmoModelStatSet(gmo, wsp.Feasible ? gmoModelStat_Feasible : gmoModelStat_InfeasibleIntermed);
         gmoSolveStatSet(gmo, gmoSolveStat_User);
         break;

      case MinimumStepsize:
      case RegularizationFailed:
         gmoModelStatSet(gmo, wsp.Feasible ? gmoModelStat_Feasible : gmoModelStat_InfeasibleIntermed);
         gmoSolveStatSet(gmo, gmoSolveStat_Solver);
         break;

      case ProblemInfeasible:
      case LocalInfeas:
#if WORHP_MAJOR >= 2
      case DivergingPenaltyObj:
      case DivergingPenaltyFeas:
      case LocalInfeasOptimal:
#endif
         gmoModelStatSet(gmo, (gmoObjNLNZ(gmo) || gmoNLNZ(gmo)) ? gmoModelStat_InfeasibleLocal : gmoModelStat_InfeasibleGlobal);
         gmoSolveStatSet(gmo, gmoSolveStat_Normal);
         break;
      case evalsNaN:
      case TooBig:
      case FunctionErrorF:
      case FunctionErrorG:
      case FunctionErrorDF:
      case FunctionErrorDG:
      case FunctionErrorHM:
         gmoModelStatSet(gmo, gmoSolveStat_EvalError);
         gmoSolveStatSet(gmo, gmoModelStat_ErrorNoSolution);
         break;
#if WORHP_MAJOR >= 2
      case DivergingPrimal:
      case DivergingDual:
#endif
         gmoModelStatSet(gmo, gmoModelStat_Unbounded);
         gmoSolveStatSet(gmo, gmoSolveStat_Normal);
         break;
      case LicenseError:
      case LicenseWarnExpiryDays:
         gmoModelStatSet(gmo, gmoModelStat_LicenseError);
         gmoSolveStatSet(gmo, gmoSolveStat_License);
         break;
   }
   gmoUnloadSolutionLegacy(gmo);
   gmoHessUnload(gmo);

   /*
   * WORHP deallocation
   */
   WorhpFree(&opt, &wsp, &par, &cnt);
   free(HMrowInit);
   free(HMcolInit);
   free(HMvalInit);
   free(HMpermInit);
   free(DGrowStartInit);
   free(DGcolInit);
   free(DGvalInit);
   free(DGvalDenseInit);

   rc = EXIT_SUCCESS;

TERMINATE:
   if(gmo != NULL)
      gmoFree(&gmo);
   if(gev != NULL)
      gevFree(&gev);

   return rc;
}

/* Objective function */
void UserF(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int objMinMaxFac)
{
   int numerr;
   int rc;

   rc = gmoEvalFuncObj(*gmo, opt->X, &opt->F, &numerr);
   if (rc != 0)
   {
      char buffer[255];
      sprintf(buffer, "Objective could not be evaluated. Error code: %d\n", rc);
      gevLogStatPChar(*gev, buffer);
   }

   opt->F *= wsp->ScaleObj;
   opt->F *= objMinMaxFac;
}

/* Function of constraints */
void UserG(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
   int numerr;
   int rc;

   for (int j = 0; j < opt->m; ++j)
   {
      rc = gmoEvalFunc(*gmo, j, opt->X, &opt->G[j], &numerr);
      if (rc != 0)
      {
         char buffer[255];
         sprintf(buffer, "Constraints could not be evaluated. Error code: %d\n", rc);
         gevLogStatPChar(*gev, buffer);
      }
   }
}

/* Gradient structure */
int StructDF(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt)
{
   UserDF(gmo, gev, opt, wsp, par, cnt, 1.0);
   return 0;
}

/* Gradient of objective function */
void UserDF(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int objMinMaxFac)
{
   int nz;
   int nlnz;
   int rc;

   rc = gmoGetObjSparse(*gmo, wsp->DF.row, wsp->DF.val, NULL, &nz, &nlnz);
   if (rc != 0)
   {
      char buffer[255];
      sprintf(buffer, "Gradient of objective could not be evaluated. Error code: %d\n", rc);
      gevLogStatPChar(*gev, buffer);
   }

   /* adapt coordinate storage format to WORHP */
   for (int i = 0; i < wsp->DF.nnz; ++i)
   {
      wsp->DF.row[i] += 1;
      wsp->DF.val[i] *= wsp->ScaleObj;
      wsp->DF.val[i] *= objMinMaxFac;
   }
}

/* Jacobian structure */
int StructDG(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** DGrowStartInit, int** DGcolInit, double** DGvalInit, double** DGvalDenseInit)
{
   int status;

   wsp->DG.nnz = gmoNZ(*gmo);
   *DGcolInit = (int*) malloc(wsp->DG.nnz * sizeof(int));
   *DGrowStartInit = (int*) malloc((opt->m+1) * sizeof(int));
   *DGvalInit = (double*) malloc(wsp->DG.nnz * sizeof(double));
   *DGvalDenseInit = (double*) malloc(opt->n * sizeof(double));

   gmoGetMatrixRow(*gmo, *DGrowStartInit, *DGcolInit, *DGvalInit, NULL);

   /* Tell Worhp to initialise the permutation vector for jacobian */
   wsp->DG.dim_perm = wsp->DG.nnz;

   /* Initialising the gradient */
   status = InitWorhpMatrix(&wsp->DG, "DG", 0, par->MatrixCC, par->MatrixCC);
   if (status != OK)
   {
      gevLogStat(*gev, "Error: Could not allocate DG structure");
      return -1;
   }

   /* init row and col structure and adapt to WORHP format*/
   for (int i = 0; i < opt->m; ++i)
      for (int j = (*DGrowStartInit)[i]; j < (*DGrowStartInit)[i+1]; ++j)
         wsp->DG.row[j] = i+1;
   for (int i = 0; i < wsp->DG.nnz; ++i)
      wsp->DG.col[i] = (*DGcolInit)[i] + 1;

   /* Tell Worhp to sort the jacobian */
   SortWorhpMatrix(&wsp->DG);

   return 0;
}

/* Jacobian of constraints */
void UserDG(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** DGrowStartInit, int** DGcolInit, double** DGvalInit, double** DGvalDenseInit)
{
   int numerr;
   int rc;
   double g;
   double DGx;

   int k = 0;
   for (int i = 0; i < opt->m; ++i)
   {
      rc = gmoEvalGrad(*gmo, i, opt->X, &g, *DGvalDenseInit, &DGx, &numerr);
      if (rc != 0)
      {
         char buffer[255];
         sprintf(buffer, "Jacobian of constraints could not be evaluated. Error code: %d\n", rc);
         gevLogStatPChar(*gev, buffer);
      }
      for ( ; k < (*DGrowStartInit)[i+1]; ++k)
      {
         (*DGvalInit)[k] = (*DGvalDenseInit)[(*DGcolInit)[k]];
      }
   }

   for (int i = 0; i < wsp->DG.nnz; ++i)
   {
      wsp->DG.val[i] = (*DGvalInit)[wsp->DG.perm[i]-1];
   }
}

/* Hessian structure */
int StructHM(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** HMrowInit, int** HMcolInit, double** HMvalInit, int** HMpermInit, int* HMdimMiss)
{
   int do2dir, dohess, status;

   /* enable hessian evaluation */
   do2dir = 0;
   dohess = 1;
   gmoHessLoad(*gmo, 0, &do2dir, &dohess);
   if (!dohess)
   {
      gevLogStat(*gev, "Error: Failed to initialize Hessian structure.");
      return -1;
   }

   /* load hessian structure */
   wsp->HM.nnz = gmoHessLagNz(*gmo);
   *HMrowInit = (int*) malloc(wsp->HM.nnz * sizeof(int));
   *HMcolInit = (int*) malloc(wsp->HM.nnz * sizeof(int));
   *HMvalInit = (double*) malloc(wsp->HM.nnz * sizeof(double));
   gmoHessLagStruct(*gmo, *HMrowInit, *HMcolInit);

   /* Tell Worhp to initialise the permutation vector for hessian */
   wsp->HM.dim_perm = wsp->HM.nnz;

   /* Initialising the hessian, while taking care of our relaxation variables.
    * Extend must be specified in this case. */
   status = InitWorhpMatrix(&wsp->HM, "HM", wsp->RelaxNvar, par->MatrixCC, par->MatrixCC);
   if (status != OK)
   {
      gevLogStat(*gev, "Error: Could not allocate HM structure");
      return -1;
   }

   /* init row and col structure and adapt to WORHP format*/
   for (int i = 0; i < wsp->HM.nnz; ++i)
   {
      wsp->HM.row[i] = (*HMcolInit)[i] + 1;
      wsp->HM.col[i] = (*HMrowInit)[i] + 1;
   }

   /* Tell Worhp to sort the hessian */
   SortWorhpMatrix(&wsp->HM);

   /* determine number of missing elements */
   *HMdimMiss = 0;
   if (wsp->HM.nnz > opt->n)
   {
      if (wsp->HM.row[wsp->HM.nnz - opt->n] != 1 || wsp->HM.col[wsp->HM.nnz - opt->n] != 1)
      {
         int i;
         for (i = wsp->HM.nnz - opt->n + 1; i < wsp->HM.nnz; ++i)
            if (wsp->HM.row[i] == wsp->HM.col[i])
               break;
         *HMdimMiss = opt->n - (wsp->HM.nnz - (i+1) +1);
      }
   }
   else
   {
      if (wsp->HM.nnz == 0)
      {
         *HMdimMiss = opt->n;
      }
      else
      {
         int i;
         for (i = 0; i < wsp->HM.nnz; ++i)
            if (!(i >= wsp->HM.nnz))
               if (wsp->HM.row[i] == wsp->HM.col[i])
                  break;
         *HMdimMiss = opt->n - (wsp->HM.nnz - i);
      }
   }
   if (*HMdimMiss)
   {
      /* increase memory for hessian */
      *HMpermInit = (int*) malloc((opt->n - *HMdimMiss) * sizeof(int));
      wsp->HM.nnz += *HMdimMiss;
      wsp->HM.val = (double*) wRealloc(wsp->HM.val, (wsp->HM.nnz + wsp->RelaxNvar) * sizeof(double));
      wsp->HM.row = (mat_int*) wRealloc(wsp->HM.row, (wsp->HM.nnz + wsp->RelaxNvar) * sizeof(mat_int));
      wsp->HM.col = (mat_int*) wRealloc(wsp->HM.col, (wsp->HM.nnz + wsp->RelaxNvar) * sizeof(mat_int));
      wsp->HM.dim_row = wsp->HM.nnz + wsp->RelaxNvar;
      wsp->HM.dim_val = wsp->HM.nnz + wsp->RelaxNvar;
      wsp->HM.dim_col = wsp->HM.nnz + wsp->RelaxNvar;

      /* Worhp will resize hessian during optimisation, we need to set default
       * again (happened in WorhpInit) */
      wsp->HM.nnzDefault = wsp->HM.nnz;
      wsp->HM.nRowDefault = wsp->HM.nRow;
      wsp->HM.nColDefault = wsp->HM.nCol;

      /* adapt coordinate storage format */
      int indexCount = 1;
      int changeIndex = 1;
      for (int i = wsp->HM.nnz - opt->n; i < wsp->HM.nnz; ++i)
      {
         if (wsp->HM.row[i] != indexCount || wsp->HM.col[i] != indexCount)
         {
            for (int j = wsp->HM.nnz-1; j > i; --j)
            {
               wsp->HM.row[j] = wsp->HM.row[j-1];
               wsp->HM.col[j] = wsp->HM.col[j-1];
            }
            wsp->HM.row[i] = indexCount;
            wsp->HM.col[i] = indexCount;
            wsp->HM.val[i] = 0.0;
         }
         else
         {
            (*HMpermInit)[changeIndex++ - 1] = i;
         }
         ++indexCount;
      }
   }

   return 0;
}

/* Hessian of Lagrangian */
void UserHM(gmoHandle_t *gmo, gevHandle_t *gev, OptVar *opt, Workspace *wsp, Params *par, Control *cnt, int** HMrowInit, int** HMcolInit, double** HMvalInit, int** HMpermInit, int* HMdimMiss, int objMinMaxFac)
{
   int numerr;
   int nnz_init;
   int rc;

   nnz_init = gmoHessLagNz(*gmo);
   rc = gmoHessLagValue(*gmo, opt->X, opt->Mu, *HMvalInit, wsp->ScaleObj*objMinMaxFac, -1.0, &numerr);
   if (rc != 0)
   {
      char buffer[255];
      sprintf(buffer, "Hessian could not be evaluated. Error code: %d\n", rc);
      gevLogStatPChar(*gev, buffer);
   }

   if (nnz_init < wsp->HM.nnz)
   {
      for (int k = 0; k < (nnz_init - (opt->n - *HMdimMiss)); ++k)
      {
         wsp->HM.val[k] = (*HMvalInit)[wsp->HM.perm[k]-1];
      }
      int j = 0;
      for (int k = (nnz_init - (opt->n - *HMdimMiss)); k < nnz_init; ++k)
      {
         wsp->HM.val[(*HMpermInit)[j]] = (*HMvalInit)[wsp->HM.perm[k]-1];
      }
   }
   else
   {
      for (int k = 0; k < nnz_init; ++k)
      {
         wsp->HM.val[k] = (*HMvalInit)[wsp->HM.perm[k]-1];
      }
   }
}
