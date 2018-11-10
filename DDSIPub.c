/*  Authors:           Andreas M"arkert, Ralf Gollmer
    Copyright to:      University of Duisburg-Essen
Language:          C
License:
This file is part of DDSIP.
Description:
This file contains procedures to evaluate upper bounds (for suggestions
for the first-stage variables produced via heuristics).

DDSIP is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

DDSIP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with DDSIP; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include <DDSIP.h>
#include <DDSIPconst.h>

#define CHECKINTEGERCUT
//#define CHECKIDENTICAL
//#define DEBUG

static int DDSIP_PrintModFileUb (int);
static int DDSIP_WarmUb (void);
static int DDSIP_GetCpxSolution (int, double *, double *, double *);

//==========================================================================
// Print lp file
int
DDSIP_PrintModFileUb (int scen)
{
    // Debugging information (core files) are written
    if ((DDSIP_param->files > 2 && DDSIP_bb->nonode == 1) || (DDSIP_param->files > 3 && DDSIP_param->files < 5))
    {
        int status;
        char fname[DDSIP_ln_fname];
        if (DDSIP_bb->DDSIP_step == eev)
            sprintf (fname, "%s/eev_sc%d%s", DDSIP_outdir, scen + 1, DDSIP_param->coretype);
        else
            sprintf (fname, "%s/ub_sc%d_n%d_gc%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, DDSIP_bb->uboutcnt++, DDSIP_param->coretype);
        status = CPXwriteprob (DDSIP_env, DDSIP_lp, fname, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write problem\n");
            return status;
        }
        else if (DDSIP_param->outlev > 21)
        {
#ifdef DEBUG
            printf ("  LP file %s written.\n",fname);
#endif
            fprintf (DDSIP_bb->moreoutfile, "  LP file %s written.\n",fname);
        }
    }

    return 0;
}

//==========================================================================
// Warm starts
int
DDSIP_WarmUb ()
{
    int status;
    if ((status = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 2)
    {
        status  = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 1, status-2);
    }
    if (DDSIP_param->hot)
    {
        // ADVIND must be 2 when using supplied start values or the last incumbent
        status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 2);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
            return status;
        }
    }				// end if (DDSIP_param->hot>2)
    return 0;
}

//==========================================================================
// Objective function contribution
int
DDSIP_Contrib (double *mipx, int scen)
{
    int j, status;

    double *obj = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar),
                                          "obj (UB)");
    status = CPXgetobj (DDSIP_env, DDSIP_lp, obj, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    for (j = 0; j < DDSIP_bb->novar; j++)
    {
        if (!scen)
            DDSIP_bb->objcontrib[j] = DDSIP_data->prob[scen] * obj[j] * mipx[j];
        else
            DDSIP_bb->objcontrib[j] += DDSIP_data->prob[scen] * obj[j] * mipx[j];
    }

    DDSIP_Free ((void **) &(obj));

    return 0;
}

//==========================================================================
// Objective function contribution in Lower Bound
int
DDSIP_Contrib_LB (double *mipx, int scen)
{
    int j, status;
    double h;

    double *obj = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar),
                                          "obj (UB)");
    status = CPXgetobj (DDSIP_env, DDSIP_lp, obj, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }
    //reset cost coeffs for first stage vars to their original value
    if (DDSIP_param->cb)
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            obj[DDSIP_bb->firstindex[j]] = DDSIP_bb->cost[j];
        }

    (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[scen] = 0.;
    for (j = 0; j < DDSIP_bb->novar; j++)
    {
        h = obj[j] * mipx[j];
        if (!scen)
            DDSIP_bb->objcontrib[j] = DDSIP_data->prob[scen] * h;
        else
            DDSIP_bb->objcontrib[j] += DDSIP_data->prob[scen] * h;
        (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[scen] += h;
    }
#ifdef DEBUG
    if (DDSIP_param->outlev > 30)
    {
        if (DDSIP_param->cb)
        {
            printf ("  -- Scenario %d, original obj=%g,\t first stage:", scen + 1, (DDSIP_node[DDSIP_bb->curnode]->ref_scenobj)[scen]);
        }
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            printf (" %.14g", mipx[DDSIP_bb->firstindex[j]]);
        }
        printf ("--\n");
    }
#endif

    DDSIP_Free ((void **) &(obj));

    return 0;
}

//==========================================================================
int
DDSIP_GetCpxSolution (int mipstatus, double *objval, double *bobjval, double *mipx)
{
    int j, status;

    status = CPXgetobjval (DDSIP_env, DDSIP_lp, objval);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective value \n");
        return status;
    }

    if (mipstatus == CPXMIP_OPTIMAL)
    {
        *bobjval = *objval;
    }
    else
    {
        status = CPXgetbestobjval (DDSIP_env, DDSIP_lp, bobjval);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to get value of best remaining node\n");
            return status;
        }
    }

    status = CPXgetx (DDSIP_env, DDSIP_lp, mipx, 0, DDSIP_bb->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get solution (UBVal) \n");
        return status;
    }
    // Returns sometimes rubbish, don't know why..
    for (j = 0; j < DDSIP_bb->novar; j++)
        if (DDSIP_Equal (mipx[j], 0.0))
            mipx[j] = 0.0;

    return status;
}

//==========================================================================
// We check if the suggested solution has already occured (return 1)
int
DDSIP_SolChk (double* cutViolation)
{
    int i, cnt = 0;
    double lhs;
    cutpool_t *currentCut;
    sug_t *tmp;

    *cutViolation = 0.;
    if (DDSIP_param->outlev >= 50)
        // Print output for debugging purposes even for multiple suggestions
    {
        int ih;
        fprintf (DDSIP_bb->moreoutfile, "\nSuggested first stage solution:\n");

        for (i = 0; i < DDSIP_bb->firstvar; i++)
        {
            if (DDSIP_bb->firsttype[i] == 'B' || DDSIP_bb->firsttype[i] == 'I' || DDSIP_bb->firsttype[i] == 'N')
            {
                ih = (int) floor ((DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] + 0.5);
                fprintf (DDSIP_bb->moreoutfile, " %20d,", ih);
            }
            else
                fprintf (DDSIP_bb->moreoutfile, " %20.15g,", (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i]);

            if (!((i + 1) % 5))
                fprintf (DDSIP_bb->moreoutfile, "\n");
        }
        fprintf (DDSIP_bb->moreoutfile, "\n");
    }
    // Check for multiple suggestion
    if (DDSIP_bb->sug[DDSIP_bb->curnode])
    {
        tmp = (DDSIP_bb->sug[DDSIP_bb->curnode]);
        do
        {
            cnt = 0;
            while (cnt < DDSIP_bb->firstvar && DDSIP_Equal ((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval[cnt], tmp->firstval[cnt]))
            {
                cnt++;
            }
            // solution occurred previously ?
            if (cnt == DDSIP_bb->firstvar)
            {
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "... multiple suggestion (same node).\n");
                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 10)
                    printf ("... multiple suggestion (same node).\n");
                DDSIP_bb->skip = 4;
                return 0;
            }
            tmp = tmp->next;
        }
        while (tmp);
    }
    // in case of asd the target (expected value) most probably has changed from the former nodes. An identical proposal has to be checked, too
    if (DDSIP_param->riskmod != 3)
    {
        for (i = 0; i < DDSIP_bb->nonode; i++)
        {
            if (DDSIP_bb->curnode != i && DDSIP_param->nodelim + 2 != i)
            {
                if (DDSIP_bb->sug[i])
                {
                    tmp = DDSIP_bb->sug[i];
                    do
                    {
                        cnt = 0;
                        while (cnt < DDSIP_bb->firstvar && DDSIP_Equal ((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval[cnt], tmp->firstval[cnt]))
                        {
                            cnt++;
                        }
                        // solution occurred previously ?
                        if (cnt == DDSIP_bb->firstvar)
                        {
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "... multiple suggestion (node %d).\n", i);
                            if (DDSIP_param->cpxscr || DDSIP_param->outlev > 10)
                                printf ("... multiple suggestion (node %d).\n", i);
                            DDSIP_bb->skip = 4;
                            return 0;
                        }
                        tmp = tmp->next;
                    }
                    while (tmp);
                }
            }
        }
    }

    // if cuts have been inserted test for violation
    if (DDSIP_bb->cutpool)
    {
        currentCut = DDSIP_bb->cutpool;
        while (currentCut)
        {
            lhs = 0.;
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                lhs += (DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval[i] * currentCut->matval[i];
            }
            if (lhs < currentCut->rhs - 1.e-9*fabs(currentCut->rhs))
            {
                *cutViolation = currentCut->rhs - lhs;
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "... violates cut %d, violation %g.\n", currentCut->number, *cutViolation);
                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 10)
                    printf ("... violates cut %d.\n", currentCut->number);
                if (!(DDSIP_param->testOtherScens) || DDSIP_bb->curnode > 6 /* || DDSIP_param->heuristic < 3 */)
                {
                    DDSIP_bb->skip = 4;
                    return 0;
                }
                break;
            }
            currentCut = currentCut->prev;
        }
    }

    // a new solution
    return 1;
}

//==========================================================================
// Functions checks feasibility of the solution gained by Heuristics and
// calculates a new objective value

int
DDSIP_UpperBound (int nrScenarios, int feasCheckOnly)
{
    int status, scen, iscen, mipstatus = 0, nr;
    int i, j, k, Bi, Bs, fs, prematureStop = 0;
    int *index;
    int wall_hrs, wall_mins,cpu_hrs, cpu_mins;

    double tmpbestvalue = 0., tmpfeasbound = 0., rest_bound, tmprisk = 0., tmprisk4 = -DDSIP_infty, tmpprob = 0.;
    double security_factor = 1.0-2.e-14, bobjval, objval, time_start, time_end, time_lap, wall_secs, cpu_secs, gap, meanGap;

    double *mipx, *values;
    double *subsol;

    double we, wr, d, mipgap, oldviol = DDSIP_infty, viol;
    double rhs;
    //  double         *Tx;

    double * sort_array;


    // if the user has supplied a start point, there is no useful information for additional variables
    if (DDSIP_bb->DDSIP_step == adv || DDSIP_bb->DDSIP_step == eev)
    {
        DDSIP_bb->curnode = 0;
        if (DDSIP_bb->DDSIP_step == adv && abs(DDSIP_param->riskmod) == 4)
        {
            (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[DDSIP_data->firstvar] = 1.e10;
        }
    }

    // Solution occured before ?
    if (DDSIP_bb->skip != -1)
        if (!DDSIP_SolChk (&oldviol))
            return 0;

    DDSIP_bb->skip = 0;

    sort_array = (double *) DDSIP_Alloc(sizeof(double), DDSIP_param->scenarios, "sort_array(UpperBound)");
    if (DDSIP_param->outlev >= DDSIP_suggested_first_stage_outlev && DDSIP_param->outlev < 50)
        // Print output
    {
        int ih;
        fprintf (DDSIP_bb->moreoutfile, "\nSuggested first stage solution:\n");

        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
            {
                ih = (int) floor ((DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j] + 0.5);
                fprintf (DDSIP_bb->moreoutfile, " %20d", ih);
            }
            else
                fprintf (DDSIP_bb->moreoutfile, " %20.15g", (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[j]);

            if (!((j + 1) % 5))
                fprintf (DDSIP_bb->moreoutfile, "\n");
        }
        fprintf (DDSIP_bb->moreoutfile, "\n");
    }

    mipx = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar), "mipx(UpperBound)");
    subsol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "subsol(UpperBound)");
    //Tx = (double *) DDSIP_Alloc(sizeof(double), DDSIP_bb->firstcon + DDSIP_bb->seccon,"(UpperBound)");

    if (!feasCheckOnly)
        DDSIP_bb->UBIters++;

    if (DDSIP_param->outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        if (feasCheckOnly)
            fprintf (DDSIP_bb->moreoutfile, "Checking first-stage solution for possible cuts:\n");
        else
            fprintf (DDSIP_bb->moreoutfile, "Evaluating first-stage solution (Upper bounds):\n");
    }

    status = DDSIP_SetCpxPara (DDSIP_param->cpxnoub, DDSIP_param->cpxubisdbl, DDSIP_param->cpxubwhich, DDSIP_param->cpxubwhat);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to set CPLEX parameters (UB) \n");
        goto TERMINATE;
    }
    if (!DDSIP_param->scalarization)
    {
        // Delete risk objective and retain objective of expected value problem
        if (DDSIP_param->riskvar)
            DDSIP_DeleteRiskObj ();

        // Some of the dual risk measures (WC, TVaR) require additional first-stage variables which need not to be
        // fixed. Therefore ...
        fs = DDSIP_bb->firstvar;
        if ((abs (DDSIP_param->riskmod) == 4 || abs (DDSIP_param->riskmod) == 5) && !DDSIP_param->riskalg)
            fs = DDSIP_bb->firstvar - 1;
    }
    else
        fs = DDSIP_bb->firstvar;

    // Change objective to the original one
    index = (int *) DDSIP_Alloc (sizeof (int), DDSIP_data->novar, "index,Newobj");
    for (j = 0; j < DDSIP_data->novar; j++)
        index[j] = j;
    status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_data->novar, index, DDSIP_data->cost + DDSIP_param->scenarios * DDSIP_param->stoccost);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change objective \n");
        goto TERMINATE;
    }
    DDSIP_Free ((void **) &(index));
    /* */
    // Heuristic returns suggested solution
    // Fix first stage variables to these values
    // Lower bounds
    values = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "values(UpperBound)");
    for(j=0; j<DDSIP_bb->firstvar; j++)
    {
        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
            values[j] = (int) floor((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval[j] + 0.5);
        else
        {
            values[j] = DDSIP_Dmax (DDSIP_bb->lborg[j],(DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval[j]);
            values[j] = DDSIP_Dmin (DDSIP_bb->uborg[j],values[j]);
        }
        for (k = 0; k < DDSIP_bb->curbdcnt; k++)
        {
            if (DDSIP_bb->curind[k] == j)
            {
                values[j] = DDSIP_Dmax (values[j], DDSIP_bb->curlb[k]);
                values[j] = DDSIP_Dmin (values[j], DDSIP_bb->curub[k]);
            }
        }
    }
    status = CPXchgbds (DDSIP_env, DDSIP_lp, fs, DDSIP_bb->firstindex, DDSIP_bb->lbident, values);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change bounds \n");
        goto TERMINATE;
    }
    // Upper bounds - first-stage variables are fixed, so values are the same
    status = CPXchgbds (DDSIP_env, DDSIP_lp, fs, DDSIP_bb->firstindex, DDSIP_bb->ubident, values);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change bounds \n");
        goto TERMINATE;
    }

    // prepare for decision about stopping: rest_bound is the expectation of all the lower bounds in this node
    if (DDSIP_bb->DDSIP_step == adv || DDSIP_bb->DDSIP_step == eev || feasCheckOnly)
        rest_bound = -DDSIP_infty;
    else if(DDSIP_param->cb)
    {
        rest_bound = DDSIP_node[DDSIP_bb->curnode]->BoundNoLag;
    }
    else
    {
        rest_bound = DDSIP_node[DDSIP_bb->curnode]->bound;
    }

    // UpperBound single-scenario problems
    meanGap = tmpprob = tmpbestvalue = 0.;
    nr = feasCheckOnly ? 1 : DDSIP_param->scenarios;
    for (iscen = 0; iscen < nr; iscen++)
    {
        scen = DDSIP_bb->ub_scen_order[iscen];

        tmpprob += DDSIP_data->prob[scen];

        status = DDSIP_ChgProb (scen, 0);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change problem \n");
            goto TERMINATE;
        }
        if (!feasCheckOnly)
        {
            // Calculate Tx+h
            //      if (DDSIP_param->outlev)  status = RecFun(param,bb,Tx);
    
            if (DDSIP_param->cpxubscr)
            {
                printf ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
                printf ("Calculating objective value for scenario problem %d in node %d (heuristic %d).....\n", scen + 1,DDSIP_bb->curnode,DDSIP_param->heuristic);
            }
            // Debugging information
            DDSIP_PrintModFileUb (scen);
    
            // Warm start
            if (DDSIP_param->hot)
            {
                status = DDSIP_WarmUb ();
                if (status)
                    goto TERMINATE;
            }
    
#ifdef DEBUG
            // query time limit amd mip rel. gap parameters
            if (DDSIP_param->cpxubscr || DDSIP_param->outlev > 21)
            {
                status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
                printf ("   -- 1st optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
            }
#endif
            //
            time_start = DDSIP_GetCpuTime ();
            // Optimize
            DDSIP_bb->scenUBIters++;
            status = CPXmipopt (DDSIP_env, DDSIP_lp);
            if (DDSIP_Error (status))
            {
                fprintf (stderr, "ERROR: Failed to optimize (UB)\n");
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to optimize (UB), status= %d\n", status);
                goto TERMINATE;
            }
            mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
            if (!DDSIP_Infeasible (mipstatus))
            {
                if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                {
                    fprintf (stderr, "ERROR: Query of CPLEX mip gap from 1st optimization failed (UpperBound) \n");
                    fprintf (stderr, "       CPXgetstat returned: %d\n",mipstatus);
                    mipgap = 1.e+30;
                }
                time_lap = DDSIP_GetCpuTime ();
#ifdef DEBUG
                if (DDSIP_param->cpxubscr ||  DDSIP_param->outlev > 11)
                {
                    j = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
                    printf ("      UB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,j,time_lap-time_start);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile,"      UB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,j,time_lap-time_start);
                }
#endif
                if (DDSIP_param->watchkappa)
                {
                    double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                    status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                              CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                    if (status)
                    {
                        fprintf (stderr, "Failed to obtain Kappa information. \n");
                    }
                    else
                    {
                        printf("maximal Kappa value        = %g\n",maxkappaval);
                        if (stablekappaval < 1.)
                        {
                            if (stablekappaval)
                                printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                            if (suspiciouskappaval)
                                printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                            if (unstablekappaval)
                                printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                            if (illposedkappaval)
                                printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                        }
                        if (DDSIP_param->outlev)
                        {
                            fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                            if (stablekappaval < 1.)
                            {
                                if (stablekappaval)
                                    fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                if (suspiciouskappaval)
                                    fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                if (unstablekappaval)
                                    fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                if (illposedkappaval)
                                    fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                            }
                        }
                    }
                }
                // If even with the lower bounds we would reach a greater value we may stop here and save some time
                // reduce the sum of the bounds rest_bound by the term for the current scenario
                if (DDSIP_param->cb)
                    rest_bound -= (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag)[scen] * DDSIP_data->prob[scen];
                else
                    rest_bound -= (DDSIP_node[DDSIP_bb->curnode]->subbound)[scen] * DDSIP_data->prob[scen];
                status = CPXgetbestobjval(DDSIP_env, DDSIP_lp, &bobjval);
    
                if (DDSIP_param->prematureStop && DDSIP_bb->bestvalue < DDSIP_infty && !DDSIP_param->riskalg && iscen < DDSIP_param->scenarios - 1 &&
                    !DDSIP_param->scalarization && DDSIP_bb->DDSIP_step != dual)
                {
                    // for risk models: premature stop requires also partial risk value
                    if (DDSIP_param->riskmod > 0)
                    {
                        // Expected excess
                        if (DDSIP_param->riskmod == 1)
                        {
                            tmprisk  += DDSIP_data->prob[scen] * DDSIP_Dmax (bobjval - DDSIP_param->risktarget, 0);
                        }
                        // Excess probability
                        if (DDSIP_param->riskmod == 2)
                        {
                            if (bobjval > DDSIP_param->risktarget)
                                tmprisk += DDSIP_data->prob[scen];
                        }
                        // Worst case cost
                        if (DDSIP_param->riskmod == 4)
                        {
                            if (bobjval > tmprisk4)
                                tmprisk4 = bobjval;
                            tmprisk = tmpprob * tmprisk4;
                        }
                        if (DDSIP_param->outlev > 50 && ((DDSIP_param->riskmod == 1)||(DDSIP_param->riskmod == 2)||(DDSIP_param->riskmod == 4)))
                            fprintf(DDSIP_bb->moreoutfile, " for premature stop: riskweight*tmprisk = %g\n", DDSIP_param->riskweight*tmprisk);
                    }
                    if (tmpbestvalue + DDSIP_param->riskweight*tmprisk + bobjval * DDSIP_data->prob[scen] + rest_bound > DDSIP_bb->bestvalue + DDSIP_param->accuracy)
                    {
                        if (!(DDSIP_bb->heurval < DDSIP_infty))
                            DDSIP_bb->skip = 100 + iscen;
                        time (&DDSIP_bb->cur_time);
                        time_end = DDSIP_GetCpuTime ();
                        time_start = time_end-time_start;
                        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                        DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
                        status = DDSIP_GetCpxSolution (mipstatus, &objval, &bobjval, mipx);
                        if (status)
                            objval = DDSIP_infty;
                        gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
                        meanGap += DDSIP_data->prob[scen] * gap;
                        if (DDSIP_param->outlev)
                        {
                            fprintf (DDSIP_bb->moreoutfile,
                                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g\t(%9.4g%%)\tStatus=%3.0d\t%3dh %02d:%02.0f cpu %3dh %02d:%05.2f (%6.2f)\n",
                                     iscen + 1, scen + 1, objval, bobjval, gap, mipstatus,
                                     wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start);
                            fprintf (DDSIP_bb->moreoutfile,
                                     "Lower bound for suggested solution yields expected value already greater than the best known\n (after %d scenarios reached %.16g, plus bound for the remaining scenarios: %.16g)\n", iscen + 1,tmpbestvalue + bobjval * DDSIP_data->prob[scen] + DDSIP_param->riskweight*tmprisk,tmpbestvalue + bobjval * DDSIP_data->prob[scen] + DDSIP_param->riskweight*tmprisk +rest_bound);
                        }
                        prematureStop = 1;
                        if (iscen > DDSIP_bb->shifts + 4)
                        {
                            cpu_secs = sort_array[DDSIP_bb->shifts];
                            for (wall_hrs = DDSIP_bb->shifts+1; wall_hrs <= iscen; wall_hrs++)
                            {
                                cpu_secs = DDSIP_Dmax (cpu_secs, sort_array[wall_hrs]);
                            }
                            wall_secs = sort_array[0];
                            for (wall_hrs = 1; wall_hrs <= iscen; wall_hrs++)
                            {
                                wall_secs += sort_array[wall_hrs];
                            }
                            wall_secs /= (1. + iscen);
                            viol = 13.5 - 5.5*(iscen - DDSIP_bb->shifts)/(DDSIP_param->scenarios - DDSIP_bb->shifts + 1.);
#ifdef DEBUG
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "### max time %g, mean %g, max>%g*mean: %d ###\n", cpu_secs, wall_secs, viol, cpu_secs > viol*wall_secs);
#endif
                            if (cpu_secs > viol*wall_secs)
                            {
                                cpu_hrs = 0;
                                for (wall_hrs = DDSIP_bb->shifts; wall_hrs <= iscen-cpu_hrs; wall_hrs++)
                                {
                                    if (sort_array[wall_hrs] > viol*wall_secs)
                                    {
                                         cpu_mins = DDSIP_bb->ub_scen_order[wall_hrs];
                                         cpu_secs = sort_array[wall_hrs];
#ifdef DEBUG
                                         if (DDSIP_param->outlev)
                                             fprintf (DDSIP_bb->moreoutfile, "### shifting scenario %d with time %g to the end of ub_scen_order ###\n", cpu_mins+1, sort_array[wall_hrs]);
#endif
                                         for (wall_mins = wall_hrs+1; wall_mins < DDSIP_param->scenarios; wall_mins++)
                                         {
                                             DDSIP_bb->ub_scen_order[wall_mins-1] = DDSIP_bb->ub_scen_order[wall_mins];
                                             sort_array[wall_mins-1] = sort_array[wall_mins];
                                         }
                                         DDSIP_bb->ub_scen_order[DDSIP_param->scenarios-1] = cpu_mins; 
                                         sort_array[DDSIP_param->scenarios-1] = cpu_secs; 
                                         cpu_hrs++;
                                         wall_hrs--;
                                    }
                                }
                            }
                        }
                        // in the first nodes check all the remaining scenarios whether they give rise to cut
                        if (DDSIP_param->alwaysBendersCuts && (DDSIP_param->testOtherScens || DDSIP_bb->curnode < 3) && DDSIP_param->heuristic > 3 )
                        {
                            time_start = time_lap = DDSIP_GetCpuTime ();
                            CPXLPptr     DDSIP_dual_lp  = NULL;
                            cutpool_t * newCut;
                            for (Bi = iscen + 1; Bi < DDSIP_param->scenarios; Bi++)
                            {
                                time_lap = DDSIP_GetCpuTime ();
                                Bs = DDSIP_bb->ub_scen_order[Bi];
                                status = DDSIP_ChgProb (Bs, 0);
                                if (status)
                                {
                                    fprintf (stderr, "ERROR: Failed to change problem \n");
                                    goto TERMINATE;
                                }
                                DDSIP_dual_lp = CPXcloneprob (DDSIP_env, DDSIP_lp, &status);
                                if (status)
                                {
                                    fprintf (stderr, "ERROR: Failed to clone problem for Benders.\n");
                                    fprintf (DDSIP_outfile, "ERROR: Failed to clone problem for Benders.\n");
                                    DDSIP_param->addBendersCuts = 0;
                                }
                                else
                                {
                                    int numRows = CPXgetnumrows (DDSIP_env, DDSIP_dual_lp);
                                    status = CPXchgprobtype (DDSIP_env, DDSIP_dual_lp, CPXPROB_LP);
                                    if (status)
                                    {
                                        fprintf (stderr, "ERROR: Failed to change type of problem for Benders to LP.\n");
                                        fprintf (DDSIP_outfile, "ERROR: Failed to change type of problem for Benders to LP.\n");
                                        return status;
                                    }
                                    status = CPXsetintparam (DDSIP_env, CPX_PARAM_PREIND, CPX_OFF);
                                    if (status)
                                    {
                                        fprintf (stderr, "ERROR: Failed to switch off preprocessing for Benders.\n");
                                        fprintf (DDSIP_outfile, "ERROR: Failed to switch off preprocessing for Benders.\n");
                                    }
                                    status = CPXdualopt (DDSIP_env, DDSIP_dual_lp);
                                    if (status)
                                    {
                                        fprintf (stderr, "ERROR: dualopt of LP for Benders failed.\n");
                                        fprintf (DDSIP_outfile, "ERROR: dualopt of LP for Benders failed.\n");
                                        return status;
                                    }
                                    status = CPXgetstat (DDSIP_env, DDSIP_dual_lp);
                                    time_end = DDSIP_GetCpuTime ();
                                    if (status == CPX_STAT_INFEASIBLE)
                                    {
                                        double *ray = (double *) DDSIP_Alloc (sizeof (double), numRows, "ray (UpperBound)");
                                        status = CPXdualfarkas (DDSIP_env, DDSIP_dual_lp, ray, &viol);
                                        if (status)
                                        {
                                            if (DDSIP_param->outlev > 21)
                                            {
                                                if (status == CPXERR_NOT_DUAL_UNBOUNDED)
                                                    fprintf (DDSIP_bb->moreoutfile," ------------ LP infeasible for scen. %3d, dualfarkas returned CPXERR_NOT_DUAL_UNBOUNDED (%6.2f) ---------------\n", Bs+1, time_end - time_lap);
                                                else
                                                    fprintf (DDSIP_bb->moreoutfile," ------------ LP infeasible for scen. %3d, dualfarkas returned %d (%6.2f) ---------------\n", Bs+1, status, time_end - time_lap);
                                            }
                                            time_lap = DDSIP_GetCpuTime ();
                                        }
                                        else
                                        {
                                            if (DDSIP_param->outlev > 21)
                                                fprintf (DDSIP_bb->moreoutfile," ------------ LP infeasible for scen. %3d, violation %g, oldviol %g (%6.2f) ---------------\n", Bs+1, viol, oldviol, time_end - time_lap);
                                        
                                        }
                                        time_lap = DDSIP_GetCpuTime ();
                                        if (!status && (viol > 1e-3) && (fabs(viol - oldviol)/viol > 1.1e-8))
                                        {
                                            int rmatbeg;
                                            int *rmatind = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->firstvar + 1), "rmatind(UpperBound)");
                                            double *rmatval = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_data->firstvar + 1), "rmatval (UpperBound)");
                                            char sense;
                                            char **rowname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(UpperBound)");
                                            char *rowstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(UpperBound)");
                                            int *cmatbeg = (int *) DDSIP_Alloc (sizeof (int*), 1, "cmatbeg(UpperBound)");
                                            int *cmatind = (int *) DDSIP_Alloc (sizeof (int), numRows, "cmatind(UpperBound)");
                                            double lhs;
                                            double *cmatval = (double *) DDSIP_Alloc (sizeof (double), numRows, "cmatval (UpperBound)");
                                            int col_nonzeros, surplus, ii;

                                            oldviol = viol;

                                            for (k = 0; k < DDSIP_bb->seccon; k++)
                                            {
                                                if (fabs (ray[DDSIP_bb->secondrowind[k]]) > 1.e-16)
                                                    break;
                                            }
                                            if (k < DDSIP_bb->seccon)
                                            {
                                                rhs = 1.;
                                                sense = 'G';
                                                rmatbeg = 0;
                                                for (k = 0; k < DDSIP_data->firstvar; k++)
                                                {
                                                    status = CPXgetcols (DDSIP_env, DDSIP_lp, &col_nonzeros, cmatbeg, cmatind, cmatval, DDSIP_bb->nocon+DDSIP_bb->cutCntr+1, &surplus,
                                                                         DDSIP_bb->firstindex[k], DDSIP_bb->firstindex[k]);
                                                    if (status)
                                                    {
                                                        fprintf(stdout, "ERROR getting column %d, surplus= %d\n", DDSIP_bb->firstindex[k], surplus);
                                                        fprintf(DDSIP_outfile, "ERROR getting column %d, surplus= %d\n", DDSIP_bb->firstindex[k], surplus);
                                                        exit(111);
                                                    }
                                                    rmatind[k] = DDSIP_bb->firstindex[k];
                                                    for (i = 0; i < numRows; i++)
                                                    {
                                                        if (fabs (ray[i]) > 1.e-16)
                                                        {
                                                            for (ii = 0; ii < col_nonzeros; ii++)
                                                            {
                                                                if (cmatind[ii] == i)
                                                                {
                                                                    rmatval[k] += ray[i]*cmatval[ii];
                                                                    break;
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                lhs = 0.;
                                                i = 0;
                                                for (k = 0; k < DDSIP_data->firstvar; k++)
                                                {
                                                    if (rmatval[k])
                                                    {
                                                        lhs += rmatval[k]*values[k];
                                                        i++;
                                                    }
                                                }
                                                if (i > 1)
                                                    rhs = lhs + viol*security_factor;
                                                else
                                                    rhs = lhs + viol;
#ifdef CHECKIDENTICAL
                                                // Check whether cut with identical matval was already added.
                                                // current one may be dominated!
                                                // This is possible with continued check of the same solution for other scens
                                                i = 1;
                                                newCut = DDSIP_bb->cutpool;
                                                while (newCut)
                                                {
                                                    for (k = 0; k < DDSIP_data->firstvar; k++)
                                                    {
                                                        if (!DDSIP_Equal (rmatval[k], newCut->matval[k]))
                                                        {
#ifdef DEBUG
////////////////////
if (DDSIP_param->outlev /*&& DDSIP_bb->curnode > 25*/)
{
  fprintf(DDSIP_bb->moreoutfile, "### 1 ### check for identical cut: Cut no. %d is different, index %d: %22.16g - %22.16g = %g \n", newCut->number, k, rmatval[k], newCut->matval[k],rmatval[k] - newCut->matval[k]);
}
////////////////////
#endif
                                                            break;
                                                        }
                                                    }
                                                    if (k < DDSIP_data->firstvar)
                                                        newCut = newCut->prev;
                                                    else
                                                    {
////////////////////
if (DDSIP_param->outlev > 21)
{
  fprintf(DDSIP_bb->moreoutfile, "### 1 ### check for identical cut: Cut no. %d is identical, rhs was: %g, now: %g ###\n", newCut->number, newCut->rhs, rhs);
}
////////////////////
                                                        if (rhs > newCut->rhs)
                                                            newCut->rhs = rhs;
                                                        i = 0;
                                                        break;
                                                    }
                                                }
#endif
                                                if (i)
                                                {
                                                    DDSIP_bb->cutNumber++;
                                                    DDSIP_bb->cutCntr++;
                                                    rowname[0] = rowstore;
                                                    sprintf (rowstore, "DDSIPBendersCut%.04d",DDSIP_bb->cutNumber);
                                                    if (DDSIP_param->outlev)
                                                    {
                                                        fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s  (infeas. scen %2d), violation %g ############\n", rowstore, Bs+1, viol);
                                                        if (DDSIP_param->outlev > 8)
                                                           printf (" ############ adding cut %s  (infeas. scen %2d) ############\n", rowstore, Bs+1);
                                                    }
                                                    if ((status = CPXaddrows(DDSIP_env, DDSIP_lp, 0, 1, DDSIP_data->firstvar, &rhs, &sense, &rmatbeg, rmatind, rmatval, NULL, rowname)))
                                                    {
                                                        fprintf (stdout, "ERROR: Failed to add cut inequality (UpperBound), status=%d \n",status);
                                                        if (DDSIP_param->outlev)
                                                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to add cut inequality (UpperBound) \n");
                                                    }
                                                    else
                                                    {
                                                        DDSIP_bb->cutAdded++;
                                                    }
                                                    // store cut in bb->cutpool
                                                    newCut = (cutpool_t*) DDSIP_Alloc(sizeof(cutpool_t), 1, "cutpool (DDSIP_UpperBound)");
                                                    newCut->prev = DDSIP_bb->cutpool;
                                                    newCut->matval = rmatval;
                                                    newCut->rhs    = rhs;
                                                    newCut->number = DDSIP_bb->cutNumber;
                                                    newCut->Benders = 1;
                                                    DDSIP_bb->cutpool = newCut;
                                                    rmatval = NULL;
                                                    //shift this infeasible scenario to first place, such that next time it is checked first
                                                    if (Bi)
                                                    {
/////////////////////
if(DDSIP_param->outlev)
  fprintf(DDSIP_bb->moreoutfile," ######## shift infeasible scen to top, Bi= %d, shifts= %d\n", Bi, DDSIP_bb->shifts);
/////////////////////
                                                        if (Bi >= DDSIP_bb->shifts)
                                                        {
                                                            DDSIP_bb->shifts++;
/////////////////////
if(DDSIP_param->outlev)
  fprintf(DDSIP_bb->moreoutfile," ##### increased shifts= %d\n", DDSIP_bb->shifts);
/////////////////////
                                                        }
                                                        k = DDSIP_bb->ub_scen_order[Bi];
                                                        for(j = Bi; j>0; j--)
                                                            DDSIP_bb->ub_scen_order[j] = DDSIP_bb->ub_scen_order[j-1];
                                                        DDSIP_bb->ub_scen_order[0] = k;
                                                        iscen = 0;
                                                        scen = k;
                                                    }
                                                }
                                                else
                                                {
                                                    if (DDSIP_param->outlev > 20)
                                                    {
                                                        fprintf (DDSIP_bb->moreoutfile," ############ not adding identical cut ############\n");
                                                    }
                                                }
                                            }
/////////////////////////////////
                                            else if (DDSIP_param->outlev > 20)
                                            {
                                                fprintf (DDSIP_bb->moreoutfile," ####----------ray doesn't contain a coefficient for a second-stage constraint\n");
                                            }
/////////////////////////////////

                                            DDSIP_Free ((void *) &rmatind);
                                            DDSIP_Free ((void *) &rmatval);
                                            DDSIP_Free ((void *) &cmatbeg);
                                            DDSIP_Free ((void *) &cmatind);
                                            DDSIP_Free ((void *) &cmatval);
                                            DDSIP_Free ((void *) &rowname);
                                            DDSIP_Free ((void *) &rowstore);
                                        }
                                        DDSIP_Free ((void *) &ray);
                                        if (!(DDSIP_param->testOtherScens))
                                            break;
                                    }
                                    else
                                    {
                                        if (DDSIP_param->outlev > 21)
                                            fprintf (DDSIP_bb->moreoutfile," ------------ LP   feasible for scen. %3d  (%6.2f) ---------------\n", Bs+1, time_end - time_lap);
                                    }
                                }
                                if (DDSIP_dual_lp != NULL)
                                {
                                    status = CPXfreeprob (DDSIP_env, &DDSIP_dual_lp);
                                    if (status)
                                        fprintf (stderr, "ERROR: CPXfreeprob failed, error code %d\n", status);
                                }
                            }
                            if (DDSIP_param->outlev > 20)
                            {
                                time_end = DDSIP_GetCpuTime ();
                                fprintf (DDSIP_bb->moreoutfile," ------------ total time for checking for cuts  %6.2fs ---------------\n", time_end-time_start);
                            }
                            status = CPXsetintparam (DDSIP_env, CPX_PARAM_PREIND, CPX_ON);
                            if (status)
                            {
                                fprintf (stderr, "ERROR: Failed to switch preprocessing back to on.\n");
                                fprintf (DDSIP_outfile, "ERROR: Failed to switch off preprocessing back to on.\n");
                            }
                        }
                        goto TERMINATE;
                    }
                }
                if (DDSIP_param->cpxnoub2 && status != CPXMIP_OPTIMAL && !feasCheckOnly)
                {
                    // more iterations with different settings
                    status = DDSIP_SetCpxPara (DDSIP_param->cpxnoub2, DDSIP_param->cpxubisdbl2, DDSIP_param->cpxubwhich2, DDSIP_param->cpxubwhat2);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to set CPLEX parameters (UpperBound) \n");
                        goto TERMINATE;
                    }
                    // query time limit amd mip rel. gap parameters
                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
                    status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
                    // continue if desired gap is not reached yet
                    if (mipgap > wr)
                    {
#ifdef DEBUG
                        if (DDSIP_param->cpxubscr || DDSIP_param->outlev > 21)
                        {
                            printf ("   -- 2nd optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                        }
#endif
                        status = CPXmipopt (DDSIP_env, DDSIP_lp);
                        mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
#ifdef DEBUG
                        if (DDSIP_param->cpxubscr || DDSIP_param->outlev > 11)
                        {
                            if (CPXgetmiprelgap(DDSIP_env, DDSIP_lp, &mipgap))
                            {
                                fprintf (stderr, "ERROR: Query of CPLEX mip gap from 2nd optimization failed (UpperBound) \n");
                                fprintf (stderr, "       CPXgetstat returned: %d\n",mipstatus);
                                mipgap = 1.e+30;
                            }
                            else
                            {
                                j = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
                                time_end = DDSIP_GetCpuTime ();
                                printf ("      UB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,j,time_end-time_lap);
                                if (DDSIP_param->outlev)
                                    fprintf (DDSIP_bb->moreoutfile,"      UB: after 2nd optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,j,time_end-time_lap);
                            }
                        }
#endif
                        if (DDSIP_param->watchkappa)
                        {
                            double maxkappaval, stablekappaval, suspiciouskappaval, unstablekappaval, illposedkappaval;
                            status = (CPXgetdblquality(DDSIP_env,DDSIP_lp,&maxkappaval,CPX_KAPPA_MAX) ||
                                      CPXgetdblquality(DDSIP_env,DDSIP_lp,&illposedkappaval,CPX_KAPPA_ILLPOSED) ||
                                      CPXgetdblquality(DDSIP_env,DDSIP_lp,&stablekappaval,CPX_KAPPA_STABLE) ||
                                      CPXgetdblquality(DDSIP_env,DDSIP_lp,&suspiciouskappaval,CPX_KAPPA_SUSPICIOUS) ||
                                      CPXgetdblquality(DDSIP_env,DDSIP_lp,&unstablekappaval,CPX_KAPPA_UNSTABLE));
                            if (status)
                            {
                                fprintf (stderr, "Failed to obtain Kappa information. \n\n");
                            }
                            else
                            {
                                printf("maximal Kappa value        = %g\n",maxkappaval);
                                if (stablekappaval < 1.)
                                {
                                    if (stablekappaval)
                                        printf("numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                    if (suspiciouskappaval)
                                        printf("numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                    if (unstablekappaval)
                                        printf("numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                    if (illposedkappaval)
                                        printf("numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                }
                                if (DDSIP_param->outlev)
                                {
                                    fprintf(DDSIP_bb->moreoutfile, "             maximal Kappa value        = %g\n",maxkappaval);
                                    if (stablekappaval < 1.)
                                    {
                                        if (stablekappaval)
                                            fprintf(DDSIP_bb->moreoutfile, "             numerically stable bases     %7.3f%%\n",stablekappaval*1e2);
                                        if (suspiciouskappaval)
                                            fprintf(DDSIP_bb->moreoutfile, "             numerically suspicious bases %7.3f%%\n",suspiciouskappaval*1e2);
                                        if (unstablekappaval)
                                            fprintf(DDSIP_bb->moreoutfile, "             numerically unstable bases   %7.3f%%\n",unstablekappaval*1e2);
                                        if (illposedkappaval)
                                            fprintf(DDSIP_bb->moreoutfile, "             numerically ill-posed bases  %7.3f%%\n",illposedkappaval*1e2);
                                    }
                                }
                            }
                        }
                    }
                    // reset CPLEX parameters for 1st optimization
                    status = DDSIP_SetCpxPara (DDSIP_param->cpxnoub, DDSIP_param->cpxubisdbl, DDSIP_param->cpxubwhich, DDSIP_param->cpxubwhat);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to set CPLEX parameters (UpperBound) \n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to set CPLEX parameters (UpperBound) \n");
                        goto TERMINATE;
                    }
                }
            }
            if ((k = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 3)
            {
                status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 2, k-1);
            }
            // Infeasible, unbounded .. ?
            if (!DDSIP_Infeasible (mipstatus))
                mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
            else
                DDSIP_bb->skip = -2;
        }
        else // feasCheckOnly
        {
            mipstatus = CPXMIP_INFEASIBLE;
        }
        if (DDSIP_NoSolution (mipstatus))
        {
            time_start = DDSIP_GetCpuTime ();
            if (!feasCheckOnly)
            {
                if (DDSIP_param->outlev)
                {
                    fprintf (DDSIP_bb->moreoutfile,
                        "Suggested solution infeasible for scenario problem %d (Status=%d)\n", scen + 1, mipstatus);
                }
                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 10)
                {
                    printf ("Suggested solution infeasible for scenario problem %d (Status=%d)\n", scen + 1, mipstatus);
                }
            }
#ifdef ADDBENDERSCUTS
            if (DDSIP_bb->DDSIP_step != adv && DDSIP_bb->DDSIP_step != eev &&
                (DDSIP_param->addBendersCuts))
            {
                CPXLPptr     DDSIP_dual_lp  = NULL;
                cutpool_t * newCut;
                int beg, end;
                if (feasCheckOnly)
                {
                    beg = 0;
                    //end = (DDSIP_param->addBendersCuts > 1 || (DDSIP_bb->curnode < 7 && DDSIP_param->testOtherScens) || !nrScenarios) ? DDSIP_param->scenarios : nrScenarios;
                    end = (DDSIP_param->addBendersCuts > 1 || (DDSIP_bb->curnode < 11 && DDSIP_param->testOtherScens) || !nrScenarios) ? DDSIP_param->scenarios : nrScenarios;
                }
                else
                {
                    beg = iscen;
                    end = ((DDSIP_param->addBendersCuts > 1 || (DDSIP_bb->curnode < 8 && DDSIP_param->testOtherScens)) /* && (DDSIP_param->heuristic > 3) */) ? DDSIP_param->scenarios : DDSIP_Imax(iscen+1, DDSIP_bb->shifts);
                }
                for (Bi = beg; Bi < end; Bi++)
                {
                    time_lap = DDSIP_GetCpuTime ();
                    Bs = DDSIP_bb->ub_scen_order[Bi];
                    status = DDSIP_ChgProb (Bs, 0);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to change problem \n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to change problem \n");
                        goto TERMINATE;
                    }
                    DDSIP_dual_lp = CPXcloneprob (DDSIP_env, DDSIP_lp, &status);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to clone problem for Benders.\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to clone problem for Benders.\n");
                        DDSIP_param->addBendersCuts = 0;
                    }
                    else
                    {
                        int numRows = CPXgetnumrows (DDSIP_env, DDSIP_dual_lp);
                        status = CPXchgprobtype (DDSIP_env, DDSIP_dual_lp, CPXPROB_LP);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to change type of problem for Benders to LP.\n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to change type of problem for Benders to LP.\n");
                            return status;
                        }
                        status = CPXsetintparam (DDSIP_env, CPX_PARAM_PREIND, CPX_OFF);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: Failed to switch off preprocessing for Benders.\n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to switch off preprocessing for Benders.\n");
                        }
                        status = CPXdualopt (DDSIP_env, DDSIP_dual_lp);
                        if (status)
                        {
                            fprintf (stderr, "ERROR: dualopt of LP for Benders failed.\n");
                            fprintf (DDSIP_outfile, "ERROR: dualopt of LP for Benders failed.\n");
                            return status;
                        }
                        status = CPXgetstat (DDSIP_env, DDSIP_dual_lp);
                        time_end = DDSIP_GetCpuTime ();
                        if (status == CPX_STAT_INFEASIBLE)
                        {
                            double *ray = (double *) DDSIP_Alloc (sizeof (double), numRows, "ray (UpperBound)");
                            viol = 0.;
                            status = CPXdualfarkas (DDSIP_env, DDSIP_dual_lp, ray, &viol);
                            if (status)
                            {
                                if (DDSIP_param->outlev > 21)
                                {
                                    if (status == CPXERR_NOT_DUAL_UNBOUNDED)
                                        fprintf (DDSIP_bb->moreoutfile," ------------ LP infeasible for scen. %3d, dualfarkas returned CPXERR_NOT_DUAL_UNBOUNDED (%6.2f) ---------------\n", Bs+1, time_end - time_lap);
                                    else
                                        fprintf (DDSIP_bb->moreoutfile," ------------ LP infeasible for scen. %3d, dualfarkas returned %d (%6.2f) ---------------\n", Bs+1, status, time_end - time_lap);
                                }
                                time_lap = DDSIP_GetCpuTime ();
                            }
                            else
                            {
                                if (DDSIP_param->outlev > 21)
                                    fprintf (DDSIP_bb->moreoutfile," ------------ LP infeasible for scen. %3d, violation %g, oldviol %g (%6.2f) ---------------\n", Bs+1, viol, oldviol, time_end - time_lap);
                            }
                            
                            time_lap = DDSIP_GetCpuTime ();
                            if (!status && (viol > 1e-1) && (fabs(viol - oldviol)/viol > 1.1e-8))
                            {
                                int rmatbeg;
                                int *rmatind = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->firstvar + 1), "rmatind(UpperBound)");
                                double *rmatval = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_data->firstvar + 1), "rmatval (UpperBound)");
                                char sense;
                                char **rowname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(UpperBound)");
                                char *rowstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(UpperBound)");
                                int *cmatbeg = (int *) DDSIP_Alloc (sizeof (int*), 1, "cmatbeg(UpperBound)");
                                int *cmatind = (int *) DDSIP_Alloc (sizeof (int), numRows, "cmatind(UpperBound)");
                                double lhs;
                                double *cmatval = (double *) DDSIP_Alloc (sizeof (double), numRows, "cmatval (UpperBound)");
                                int col_nonzeros, surplus, ii;

                                oldviol = viol;

                                for (k = 0; k < DDSIP_bb->seccon; k++)
                                {
                                    if (fabs (ray[DDSIP_bb->secondrowind[k]]) > 1.e-16)
                                        break;
                                }
                                if (k < DDSIP_bb->seccon)
                                {
                                    rhs = 1.;
                                    sense = 'G';
                                    rmatbeg = 0;
                                    for (k = 0; k < DDSIP_data->firstvar; k++)
                                    {
                                        status = CPXgetcols (DDSIP_env, DDSIP_lp, &col_nonzeros, cmatbeg, cmatind, cmatval, DDSIP_bb->nocon+DDSIP_bb->cutCntr+1, &surplus,
                                                             DDSIP_bb->firstindex[k], DDSIP_bb->firstindex[k]);
                                        if (status)
                                        {
                                            fprintf(stdout, "ERROR getting column %d, surplus= %d\n", DDSIP_bb->firstindex[k], surplus);
                                            fprintf(DDSIP_outfile, "ERROR getting column %d, surplus= %d\n", DDSIP_bb->firstindex[k], surplus);
                                            exit(111);
                                        }
                                        rmatind[k] = DDSIP_bb->firstindex[k];
                                        for (i = 0; i < numRows; i++)
                                        {
                                            if (fabs (ray[i]) > 1.e-16)
                                            {
                                                for (ii = 0; ii < col_nonzeros; ii++)
                                                {
                                                    if (cmatind[ii] == i)
                                                    {
                                                        rmatval[k] += ray[i]*cmatval[ii];
                                                        break;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    lhs = 0.;
                                    i = 0;
                                    for (k = 0; k < DDSIP_data->firstvar; k++)
                                    {
                                        if (rmatval[k])
                                        {
                                            lhs += rmatval[k]*values[k];
                                            i++;
                                        }
                                    }
                                    if (i > 1)
                                        rhs = lhs + viol*security_factor;
                                    else
                                        rhs = lhs + viol;
#ifdef CHECKIDENTICAL
                                    // Check whether cut with identical matval was already added.
                                    // current one may be dominated!
                                    // This is possible with continued check of the same solution for other scens
                                    i = 1;
                                    newCut = DDSIP_bb->cutpool;
                                    while (newCut)
                                    {
                                        for (k = 0; k < DDSIP_data->firstvar; k++)
                                        {
                                            if (!DDSIP_Equal (rmatval[k], newCut->matval[k]))
                                            {
#ifdef DEBUG
////////////////////
if (DDSIP_param->outlev /*&& DDSIP_bb->curnode > 25*/)
{
  fprintf(DDSIP_bb->moreoutfile, "### 2 ### check for identical cut: Cut no. %d is different, index %d: %22.16g - %22.16g = %g \n", newCut->number, k, rmatval[k], newCut->matval[k],rmatval[k] - newCut->matval[k]);
}
////////////////////
#endif
                                                break;
                                            }
                                        }
                                        if (k < DDSIP_data->firstvar)
                                            newCut = newCut->prev;
                                        else
                                        {
////////////////////
if (DDSIP_param->outlev > 21)
{
  fprintf(DDSIP_bb->moreoutfile, "### 2 ### check for identical cut: Cut no. %d is identical, rhs was: %g, now: %g ###\n", newCut->number, newCut->rhs, rhs);
}
////////////////////
                                            if (rhs > newCut->rhs)
                                                newCut->rhs = rhs;
                                            i = 0;
                                            break;
                                        }
                                    }
#endif
                                    if (i)
                                    {
                                        DDSIP_bb->cutNumber++;
                                        DDSIP_bb->cutCntr++;
                                        rowname[0] = rowstore;
                                        sprintf (rowstore, "DDSIPBendersCut%.04d",DDSIP_bb->cutNumber);
                                        if (DDSIP_param->outlev)
                                        {
                                            fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s  (infeas. scen %2d), violation %g ############\n", rowstore, Bs+1, viol);
                                            if (DDSIP_param->outlev > 8)
                                                printf (" ############ adding cut %s  (infeas. scen %2d) ############\n", rowstore, Bs+1);
                                        }
                                        if ((status = CPXaddrows(DDSIP_env, DDSIP_lp, 0, 1, DDSIP_data->firstvar, &rhs, &sense, &rmatbeg, rmatind, rmatval, NULL, rowname)))
                                        {
                                            fprintf (stdout, "ERROR: Failed to add cut inequality (UpperBound), status=%d \n",status);
                                            if (DDSIP_param->outlev)
                                                fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to add cut inequality (UpperBound) \n");
                                        }
                                        else
                                        {
                                            DDSIP_bb->cutAdded++;
                                        }
                                        // store cut in bb->cutpool
                                        newCut = (cutpool_t*) DDSIP_Alloc(sizeof(cutpool_t), 1, "cutpool (DDSIP_UpperBound)");
                                        newCut->prev = DDSIP_bb->cutpool;
                                        newCut->matval = rmatval;
                                        newCut->rhs    = rhs;
                                        newCut->number = DDSIP_bb->cutNumber;
                                        newCut->Benders = 1;
                                        DDSIP_bb->cutpool = newCut;
                                        rmatval = NULL;
                                    }
                                    //shift this infeasible scenario to first place, such that next time it is checked first
                                    if (Bi != iscen)
                                    {
                                        if (Bi >= DDSIP_bb->shifts)
                                        {
                                            DDSIP_bb->shifts++;
                                        }
                                        k = DDSIP_bb->ub_scen_order[Bi];
                                        for(j = Bi; j>0; j--)
                                            DDSIP_bb->ub_scen_order[j] = DDSIP_bb->ub_scen_order[j-1];
                                        DDSIP_bb->ub_scen_order[0] = k;
                                        iscen = 0;
                                        scen = k;
                                    }
                                }
/////////////////////////////////
                                else if (DDSIP_param->outlev > 20)
                                {
                                    fprintf (DDSIP_bb->moreoutfile," ####----------ray doesn't contain a coefficient for a second-stage constraint\n");
                                }
/////////////////////////////////

                                DDSIP_Free ((void *) &rmatind);
                                DDSIP_Free ((void *) &rmatval);
                                DDSIP_Free ((void *) &cmatbeg);
                                DDSIP_Free ((void *) &cmatind);
                                DDSIP_Free ((void *) &cmatval);
                                DDSIP_Free ((void *) &rowname);
                                DDSIP_Free ((void *) &rowstore);
                            }
                            DDSIP_Free ((void *) &ray);
                        }
                        else
                        {
                            if (DDSIP_param->outlev > 21)
                                fprintf (DDSIP_bb->moreoutfile," ------------ LP   feasible for scen. %3d  (%6.2f) ---------------\n", Bs+1, time_end - time_lap);
                        }
                    }
                    if (DDSIP_dual_lp != NULL)
                    {
                        status = CPXfreeprob (DDSIP_env, &DDSIP_dual_lp);
                        if (status)
                            fprintf (stderr, "ERROR: CPXfreeprob failed, error code %d\n", status);
                    }
                    status = CPXsetintparam (DDSIP_env, CPX_PARAM_PREIND, CPX_ON);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to switch preprocessing back to on.\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to switch off preprocessing back to on.\n");
                    }
                }
                if (DDSIP_param->outlev > 20)
                {
                    time_end = DDSIP_GetCpuTime ();
                    fprintf (DDSIP_bb->moreoutfile," ------------ total time for checking for Benders cuts  %6.2fs ---------------\n", time_end-time_start);
                }
            }
#endif
#ifdef ADDINTEGERCUTS
            //if all first-stage variables are binary ones, we can add an inequality, cutting off this point
            if (DDSIP_bb->DDSIP_step != adv && DDSIP_bb->DDSIP_step != eev &&
                (DDSIP_param->addIntegerCuts && DDSIP_param->heuristic > 3 /* && !feasCheckOnly */))
            {
                int rmatbeg;
                int *rmatind = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->novar), "rmatind(UpperBound)");
                double *rmatval = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar), "rmatval (UpperBound)");
                char sense;
                char **rowname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(UpperBound)");
                char *rowstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(UpperBound)");
#ifdef CHECKINTEGERCUT
                double objv = DDSIP_infty;
                CPXLPptr     DDSIP_dual_lp  = NULL;
#endif
                rowname[0] = rowstore;
                rhs = 1.;
                sense = 'G';
                rmatbeg = 0;
                for (k = 0; k < DDSIP_bb->firstvar; k++)
                {
                    rmatind[k] = DDSIP_bb->firstindex[k];
                    if (floor((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval[k] + 0.5))
                    {
                        rmatval[k] = -1.;
                        rhs -= 1.;
                    }
                    else
                    {
                        rmatval[k] = 1.;
                    }
                }
                for (; k < DDSIP_bb->novar; k++)
                {
                  rmatind[k] = DDSIP_bb->secondindex[k - DDSIP_data->firstvar];
                  rmatval[k] = 0.;
                }
#ifdef CHECKINTEGERCUT
                // check whether this cut is dominated
                if (DDSIP_param->addBendersCuts  /*&& DDSIP_bb->from_scenario > -1*/)
                {
                    status = DDSIP_RestoreBoundAndType ();
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to change bounds for check, error code %d\n", status);
                        goto TERMINATE;
                    }
                    status = DDSIP_ChgProb (DDSIP_bb->from_scenario, 0);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to change problem for check to scenario %d, error code %d\n", scen + 1, status);
                        goto TERMINATE;
                    }
                    DDSIP_dual_lp = CPXcloneprob (DDSIP_env, DDSIP_lp, &status);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to clone problem for check.\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to clone problem for check.\n");
                        return status;
                    }
                    status = CPXchgprobtype (DDSIP_env, DDSIP_dual_lp, CPXPROB_LP);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to change type of problem for check to LP.\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to change type of problem for check to LP.\n");
                        return status;
                    }
                    status = CPXchgobj (DDSIP_env, DDSIP_dual_lp, DDSIP_bb->novar, rmatind, rmatval);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to update objective coefficients\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to update objective coefficients\n");
                        return status;
                    }
                    status = CPXdualopt (DDSIP_env, DDSIP_dual_lp);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: dualopt of LP for check failed.\n");
                        fprintf (DDSIP_outfile, "ERROR: dualopt of LP for check failed.\n");
                        return status;
                    }
#ifdef DEBUG
/////////////////////////////////
                    fprintf (stderr," return value of dualopt: %d\n", status);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile," return value of dualopt: %d\n", status);
/////////////////////////////////
#endif
                    status = CPXgetstat (DDSIP_env, DDSIP_dual_lp);
                    if (!status)
                    {
                        fprintf (stderr, "ERROR: dualopt of LP for check failed, getstat= %d.\n", status);
                        fprintf (DDSIP_outfile, "ERROR: dualopt of LP for check failed, getstat= %d.\n", status);
                        return status;
                    }
#ifdef DEBUG
/////////////////////////////////
                    fprintf (stderr," getstat for dualopt: %d\n", status);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile," getstat for dualopt: %d\n", status);
/////////////////////////////////
#endif
                    status = CPXgetobjval (DDSIP_env, DDSIP_dual_lp, &objv);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: dualopt: getobjval failed.\n");
                        fprintf (DDSIP_outfile, "ERROR: dualopt: getobjval failed.\n");
                        return status;
                    }
#ifdef DEBUG
/////////////////////////////////
                    fprintf (stderr," return value of getobjval: %d, value=%g\n", status, objv);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile," return value of getobjval: %d, value=%g\n", status, objv);
/////////////////////////////////
#endif
                }
                else
                {
                    objv = -1.e20;
                }
                if (objv < rhs - 3e-2)
                {
                    cutpool_t * newCut;
                    DDSIP_bb->cutNumber++;
                    DDSIP_bb->cutCntr++;
                    sprintf (rowstore, "DDSIPIntegerCut%.04d",DDSIP_bb->cutNumber);
                    if (DDSIP_param->outlev)
                    {
                        if (objv > -1.e20)
                            fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s  (objval = %g < %g) ############\n", rowstore, objv, rhs);
                        else
                            fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s ############\n", rowstore);
                        if (DDSIP_param->outlev > 8)
                            printf (" ############ adding cut %s  (objval = %g < %g) ############\n", rowstore, objv, rhs);
                    }
#else
                    cutpool_t * newCut;
                    DDSIP_bb->cutNumber++;
                    DDSIP_bb->cutCntr++;
                    sprintf (rowstore, "DDSIPIntegerCut%.04d",DDSIP_bb->cutNumber);
                    if (DDSIP_param->outlev)
                    {
                        fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s ############\n", rowstore);
                        if (DDSIP_param->outlev > 8)
                            printf (" ############ adding cut %s ############\n", rowstore);
                    }
#endif
                    if ((status = CPXaddrows(DDSIP_env, DDSIP_lp, 0, 1, DDSIP_data->firstvar, &rhs, &sense, &rmatbeg, rmatind, rmatval, NULL, rowname)))
                    {
                        fprintf (stderr, "ERROR: Failed to add cut inequality (UpperBound), status=%d \n",status);
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to add cut inequality (UpperBound) \n");
                    }
                    else
                    {
                        DDSIP_bb->cutAdded++;
                    }
                    // store cut in bb->cutpool
                    newCut = (cutpool_t*) DDSIP_Alloc(sizeof(cutpool_t), 1, "cutpool (DDSIP_UpperBound)");
                    newCut->prev = DDSIP_bb->cutpool;
                    newCut->matval = rmatval;
                    newCut->rhs    = rhs;
                    newCut->number = DDSIP_bb->cutNumber;
                    newCut->Benders = 0;
                    DDSIP_bb->cutpool = newCut;
                    rmatval = NULL;
#ifdef CHECKINTEGERCUT
                }
                else
                {
                    if (DDSIP_param->outlev > 20)
                    {
                        fprintf (DDSIP_bb->moreoutfile," ############ integer cut is not effective: objval = %g, rhs = %g ############\n", objv, rhs);
                        printf (" ############ integer cut is not effective: objval = %g, rhs = %g ############\n", objv, rhs);
                    }
                }
                if (DDSIP_dual_lp != NULL)
                {
                    status = CPXfreeprob (DDSIP_env, &DDSIP_dual_lp);
                    if (status)
                        fprintf (stderr, "ERROR: CPXfreeprob failed, error code %d\n", status);
                }
#endif
                DDSIP_Free ((void *) &rmatind);
                DDSIP_Free ((void *) &rmatval);
                DDSIP_Free ((void *) &rowname);
                DDSIP_Free ((void *) &rowstore);
            }
#endif
            //shift this infeasible scenario to first place, such that next time it is checked first
            if (iscen || !DDSIP_bb->shifts)
            {
                if (iscen >= DDSIP_bb->shifts)
                    DDSIP_bb->shifts++;
                k = DDSIP_bb->ub_scen_order[iscen];
                cpu_secs = sort_array[iscen];
                for(j = iscen; j>0; j--)
                {
                    DDSIP_bb->ub_scen_order[j] = DDSIP_bb->ub_scen_order[j-1];
                    sort_array[j] = sort_array[j-1];
                }
                DDSIP_bb->ub_scen_order[0] = k;
                sort_array[0] = cpu_secs;
            }
            if (iscen > DDSIP_bb->shifts + 4)
            {
                cpu_secs = sort_array[DDSIP_bb->shifts];
                for (wall_hrs = DDSIP_bb->shifts+1; wall_hrs < iscen; wall_hrs++)
                {
                    cpu_secs = DDSIP_Dmax (cpu_secs, sort_array[wall_hrs]);
                }
                wall_secs = sort_array[0];
                for (wall_hrs = 1; wall_hrs < iscen; wall_hrs++)
                {
                    wall_secs += sort_array[wall_hrs];
                }
                wall_secs /= iscen;
                viol = 14.5 - 5.5*(iscen - DDSIP_bb->shifts)/(DDSIP_param->scenarios - DDSIP_bb->shifts + 1.);
#ifdef DEBUG
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "### max time %g, mean %g, max>%g*mean: %d ###\n", cpu_secs, wall_secs, viol, cpu_secs > viol*wall_secs);
#endif
                if (cpu_secs > viol*wall_secs)
                {
                    cpu_hrs = 0;
                    for (wall_hrs = DDSIP_bb->shifts; wall_hrs < iscen-cpu_hrs; wall_hrs++)
                    {
                        if (sort_array[wall_hrs] > viol*wall_secs)
                        {
                             cpu_mins = DDSIP_bb->ub_scen_order[wall_hrs];
                             cpu_secs = sort_array[wall_hrs];
#ifdef DEBUG
                             if (DDSIP_param->outlev)
                                 fprintf (DDSIP_bb->moreoutfile, "### shifting scenario %d with time %g to the end of ub_scen_order ###\n", cpu_mins+1, sort_array[wall_hrs]);
#endif
                             for (wall_mins = wall_hrs+1; wall_mins < DDSIP_param->scenarios; wall_mins++)
                             {
                                 DDSIP_bb->ub_scen_order[wall_mins-1] = DDSIP_bb->ub_scen_order[wall_mins];
                                 sort_array[wall_mins-1] = sort_array[wall_mins];
                             }
                             DDSIP_bb->ub_scen_order[DDSIP_param->scenarios-1] = cpu_mins; 
                             sort_array[DDSIP_param->scenarios-1] = cpu_secs; 
                             cpu_hrs++;
                             wall_hrs--;
                        }
                    }
                }
            }
            if (feasCheckOnly)
                continue;
            goto TERMINATE;
        }
        if (!feasCheckOnly)
        {
            // Get solution
            status = DDSIP_GetCpxSolution (mipstatus, &objval, &bobjval, mipx);
            if (status)
                goto TERMINATE;

            // Remember objective function contribution
            if (DDSIP_param->outlev > 10)
                DDSIP_Contrib (mipx, scen);

            // Store second stage solutions
            for (j = 0; j < DDSIP_bb->secvar; j++)
                DDSIP_bb->cur_secstage[scen][j] = mipx[DDSIP_bb->secondindex[j]];

            if ( mipstatus == 101 && bobjval > DDSIP_infty )
                bobjval=objval;

            subsol[scen] = objval;

            gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
            meanGap += DDSIP_data->prob[scen] * gap;
            time_end = DDSIP_GetCpuTime ();
            time_start = time_end-time_start;
            // in order to sort the scenarios which take much longer to the end (hopefully often not to be evaluated due to premature stop)
            sort_array[iscen] = time_start;
            // Debugging information
            if (DDSIP_param->outlev)
            {
                time (&DDSIP_bb->cur_time);
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
                fprintf (DDSIP_bb->moreoutfile,
                         "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g\t(%9.4g%%)\tStatus=%3.0d\t%3dh %02d:%02.0f cpu %3dh %02d:%05.2f (%6.2f)\n",
                         iscen + 1, scen + 1, objval, bobjval, gap, mipstatus,
                         wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start);
                if (DDSIP_param->outlev>8)
                {
                    printf ("%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g\t(%9.4g%%)\tStatus=%3.0d\t%3dh %02d:%02.0f cpu %3dh %02d:%05.2f (%6.2f)\n",
                            iscen + 1, scen + 1, objval, bobjval, gap, mipstatus,
                            wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start);
                    if (DDSIP_param->cpxscr)
                        printf ("\n");
                }
            }

            tmpfeasbound += DDSIP_Dmin (objval, bobjval) * DDSIP_data->prob[scen];
            tmpbestvalue += objval * DDSIP_data->prob[scen];

            // If even with the lower bounds we would reach a greater value we may stop here and save some time
            // reduce the sum of the bounds rest_bound by the term for the current scenario
            if (DDSIP_param->prematureStop && DDSIP_bb->DDSIP_step != dual && DDSIP_param->heuristic > 99 && DDSIP_bb->bestvalue < DDSIP_infty && !DDSIP_param->riskalg && iscen < DDSIP_param->scenarios - 1 &&
                    tmpbestvalue + rest_bound > DDSIP_bb->bestvalue + DDSIP_param->accuracy && !DDSIP_param->scalarization)
            {
                if (!(DDSIP_bb->heurval < DDSIP_infty))
                    DDSIP_bb->skip = 100 + iscen;
                if (DDSIP_param->outlev)
                {
                    printf (
                        "Suggested solution yields expected value already greater than the best known\n (after %d scenarios reached %.16g, plus bound for the remaining scenarios: %.16g)\n", iscen + 1,tmpbestvalue,tmpbestvalue+rest_bound);
                    fprintf (DDSIP_bb->moreoutfile,
                             "Suggested solution yields expected value already greater than the best known\n (after %d scenarios reached %.16g, plus bound for the remaining scenarios: %.16g)\n", iscen + 1,tmpbestvalue,tmpbestvalue+rest_bound);
                    // evaluate solution times in order to see whether there are very difficult ones
                }
                if (iscen > DDSIP_bb->shifts + 4)
                {
                    wall_secs = cpu_secs = sort_array[0];
                    for (wall_hrs = 1; wall_hrs <= iscen; wall_hrs++)
                    {
                        wall_secs += sort_array[wall_hrs];
                        cpu_secs = DDSIP_Dmax (cpu_secs, sort_array[wall_hrs]);
                    }
                    viol = 14.5 - 5.5*(iscen - DDSIP_bb->shifts)/(DDSIP_param->scenarios - DDSIP_bb->shifts + 1.);
                    wall_secs /= (0.01 + iscen);
#ifdef DEBUG
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "### max time %g, mean %g, max>%g*mean: %d ###\n", cpu_secs, wall_secs, viol, cpu_secs > 5.2*wall_secs);
#endif
                    if (cpu_secs > viol*wall_secs)
                    {
                        cpu_hrs = 0;
                        for (wall_hrs = DDSIP_bb->shifts; wall_hrs <= iscen-cpu_hrs; wall_hrs++)
                        {
                            if (sort_array[wall_hrs] > viol*wall_secs)
                            {
                                 cpu_mins = DDSIP_bb->ub_scen_order[wall_hrs];
                                 cpu_secs = sort_array[wall_hrs];
#ifdef DEBUG
                                 if (DDSIP_param->outlev)
                                     fprintf (DDSIP_bb->moreoutfile, "### shifting scenario %d with time %g to the end of ub_scen_order ###\n", cpu_mins+1, sort_array[wall_hrs]);
#endif
                                 for (wall_mins = wall_hrs+1; wall_mins < DDSIP_param->scenarios; wall_mins++)
                                 {
                                     DDSIP_bb->ub_scen_order[wall_mins-1] = DDSIP_bb->ub_scen_order[wall_mins];
                                     sort_array[wall_mins-1] = sort_array[wall_mins];
                                 }
                                 DDSIP_bb->ub_scen_order[DDSIP_param->scenarios-1] = cpu_mins; 
                                 sort_array[DDSIP_param->scenarios-1] = cpu_secs; 
                                 cpu_hrs++;
                                 wall_hrs--;
                            }
                        }
                    }
                }
                goto TERMINATE;
            }

            if ((j = CPXgetnummipstarts(DDSIP_env, DDSIP_lp)) > 3)
            {
                status    = CPXdelmipstarts (DDSIP_env, DDSIP_lp, 2, j-1);
            }

            if (DDSIP_param->hot)
            {
                status = CPXsetintparam (DDSIP_env, CPX_PARAM_ADVIND, 2);
                if (status)
                {
                    fprintf (DDSIP_outfile, "ERROR: Failed to set cplex parameter CPX_PARAM_ADVIND.\n");
                    return status;
                }
            }
        }
    }				// end for iscen=..
    if (feasCheckOnly)
        goto TERMINATE;

    if (DDSIP_param->riskmod)
    {
        DDSIP_bb->curexp = 0.0;
        for (j = 0; j < DDSIP_param->scenarios; j++)
            DDSIP_bb->curexp += DDSIP_data->prob[j] * subsol[j];

        // Calculate risk objective
        DDSIP_RiskObjective (subsol);

        if (DDSIP_param->scalarization)
        {
            we = (DDSIP_bb->curexp - DDSIP_param->ref_point[0]) / DDSIP_param->ref_scale[0];
            wr = (DDSIP_bb->currisk - DDSIP_param->ref_point[1]) / DDSIP_param->ref_scale[1];
            d = DDSIP_Dmax (we, wr);
            tmpbestvalue = d + DDSIP_param->ref_eps * (we + wr);
        }
        else
        {
            if (DDSIP_param->riskmod > 0)
                tmpbestvalue = DDSIP_bb->curexp + DDSIP_param->riskweight * DDSIP_bb->currisk;
            else
                tmpbestvalue = DDSIP_bb->currisk;
        }
    }


    DDSIP_bb->heurval = DDSIP_Dmin (tmpbestvalue, DDSIP_bb->heurval);

    DDSIP_bb->meanGapUB = DDSIP_Dmax(DDSIP_bb->meanGapUB, meanGap);

    // Debugging information
    if (DDSIP_param->outlev)
    {
        if (DDSIP_param->riskmod)
            fprintf (DDSIP_bb->moreoutfile, "\tQ_E = %.10g \t Q_R = %.10g\n", DDSIP_bb->curexp, DDSIP_bb->currisk);
        fprintf (DDSIP_bb->moreoutfile, "\tNew suggested upper bound = %.16g      \t(mean MIP gap: %g%%)\n", tmpbestvalue, meanGap);

        if (DDSIP_param->outlev > 30)
        {
            char **colname;
            char *colstore;

            colname = (char **) DDSIP_Alloc (sizeof (char *), (DDSIP_bb->novar), "colname(UpperBound)");
            colstore = (char *) DDSIP_Alloc (sizeof (char), (DDSIP_bb->novar) * DDSIP_ln_varname, "colstore(UpperBound)");

            status =
                CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore,
                               (DDSIP_bb->novar) * DDSIP_ln_varname, &j, 0, DDSIP_bb->novar - 1);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to get column names (UB)\n");
                goto TERMINATE;
            }

            fprintf (DDSIP_bb->moreoutfile, "\nObjective function composition:\n");
            for (j = 0; j < DDSIP_bb->novar; j++)
                if (fabs (DDSIP_bb->objcontrib[j]) > DDSIP_param->accuracy)
                    fprintf (DDSIP_bb->moreoutfile, "\t%s:  %.12g\n", colname[j], DDSIP_bb->objcontrib[j]);

            DDSIP_Free ((void **) &(colstore));
            DDSIP_Free ((void **) &(colname));
        }
    }
    // Update if better solution was found
    if (tmpbestvalue < DDSIP_bb->bestvalue)
    {
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\t(Current best bound, improvement %16.12g, %g%%)\n",DDSIP_bb->bestvalue-tmpbestvalue,1e+2*(DDSIP_bb->bestvalue-tmpbestvalue)/(fabs(DDSIP_bb->bestvalue) + 1.e-16));
        DDSIP_bb->bestsol_in_curnode = 1;
        DDSIP_bb->bestvalue = tmpbestvalue;
        DDSIP_bb->feasbound = tmpfeasbound;
        DDSIP_bb->bestrisk = DDSIP_bb->currisk;
        DDSIP_bb->bestexp = DDSIP_bb->curexp;

        // We store all risk measures
        if (DDSIP_param->riskmod)
            for (j = 0; j < DDSIP_maxrisk; j++)
                DDSIP_bb->bestriskval[j] = DDSIP_bb->curriskval[j];

        for (j = 0; j < DDSIP_bb->firstvar; j++)
            DDSIP_bb->bestsol[j] = ((DDSIP_bb->sug[DDSIP_param->nodelim + 2])->firstval)[j];

        // compare with bound_optimal_node
        if (DDSIP_bb->bestvalue < DDSIP_bb->bound_optimal_node)
        {
            DDSIP_bb->bound_optimal_node = DDSIP_bb->bestvalue;
            DDSIP_bb->found_optimal_node = 0;
        }

        rhs = DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestvalue;
        if (rhs > 0.)
        {
            // something's wrong with the accuracy of solutions - the heuristic gave a better value than the lower bound in the same node!
            printf ("*\t*** WARNING: heuristic value is better than lower bound of the node by %g ***\n", rhs);
            if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_outfile, "*\t*** WARNING: heuristic value is better than lower bound of the current node by %g", rhs);
                fprintf (DDSIP_bb->moreoutfile, "*\t*** WARNING: heuristic value is better than lower bound of the current node by %g", rhs);
                if (rhs/(fabs(tmpbestvalue)+1.e-11) > 1e-9)
                {
                    fprintf (DDSIP_outfile, ", MIP gaps or tolerances probably not small enough ***\n");
                    fprintf (DDSIP_bb->moreoutfile, ", MIP gaps or tolerances probably not small enough ***\n");
                }
                else
                {
                    fprintf (DDSIP_outfile, " ***\n");
                    fprintf (DDSIP_bb->moreoutfile, " ***\n");
                }
            }
            DDSIP_bb->correct_bounding = DDSIP_Dmax (DDSIP_bb->correct_bounding,  rhs);
        }
        // copy second stage solutions of current best
        for (k = 0; k < DDSIP_param->scenarios; k++)
        {
            memcpy(DDSIP_bb->secstage[k], DDSIP_bb->cur_secstage[k], DDSIP_bb->secvar*sizeof(double));
            DDSIP_bb->subsol[k] = subsol[k];
        }
        // Print a line of output
        DDSIP_PrintStateUB (DDSIP_param->heuristic);
        // heuristic was successful
        DDSIP_bb->heurSuccess++;
    }

    if (DDSIP_bb->shifts + 2 < DDSIP_param->scenarios)
    {
        // evaluate solution times in order to shift very difficult ones to the end
        wall_secs = cpu_secs = 0.;
        for (wall_hrs = DDSIP_bb->shifts; wall_hrs < DDSIP_param->scenarios; wall_hrs++)
        {
            wall_secs += sort_array[wall_hrs];
            cpu_secs = DDSIP_Dmax (cpu_secs, sort_array[wall_hrs]);
        }
        wall_secs /= (0.01 + DDSIP_param->scenarios - DDSIP_bb->shifts);
        viol = 13.5 - 4.5*(iscen - DDSIP_bb->shifts)/(DDSIP_param->scenarios - DDSIP_bb->shifts + 1.);
#ifdef DEBUG
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "### max time %g, mean %g, max>%g*mean: %d ###\n", cpu_secs, wall_secs, viol, cpu_secs > viol*wall_secs);
#endif
        if (cpu_secs > viol*wall_secs)
        {
            cpu_hrs = 0;
            for (wall_hrs = DDSIP_bb->shifts; wall_hrs < DDSIP_param->scenarios-cpu_hrs; wall_hrs++)
            {
                if (sort_array[wall_hrs] > viol*wall_secs)
                {
                     cpu_mins = DDSIP_bb->ub_scen_order[wall_hrs];
                     cpu_secs = sort_array[wall_hrs];
#ifdef DEBUG
                     if (DDSIP_param->outlev)
                         fprintf (DDSIP_bb->moreoutfile, "### shifting scenario %d with time %g to the end of ub_scen_order ###\n", cpu_mins+1, sort_array[wall_hrs]);
#endif
                     for (wall_mins = wall_hrs+1; wall_mins < DDSIP_param->scenarios; wall_mins++)
                     {
                         DDSIP_bb->ub_scen_order[wall_mins-1] = DDSIP_bb->ub_scen_order[wall_mins];
                         sort_array[wall_mins-1] = sort_array[wall_mins];
                     }
                     DDSIP_bb->ub_scen_order[DDSIP_param->scenarios-1] = cpu_mins; 
                     sort_array[DDSIP_param->scenarios-1] = cpu_secs; 
                     cpu_hrs++;
                     wall_hrs--;
                }
            }
        }
    }

TERMINATE:

    DDSIP_Free ((void **) &(sort_array));
    // if not only feasibility was tested, add the suggested first-stage to the list of suggested solutions
    if (!feasCheckOnly)
    {
        DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = DDSIP_bb->sug[DDSIP_bb->curnode];
        DDSIP_bb->sug[DDSIP_bb->curnode] = DDSIP_bb->sug[DDSIP_param->nodelim + 2];
        DDSIP_bb->sug[DDSIP_param->nodelim + 2] = NULL;
    }
    // Store optimal value of EEV problem
    if (DDSIP_bb->DDSIP_step == eev)
    {
        if (DDSIP_NoSolution (mipstatus))
            DDSIP_bb->expbest = DDSIP_infty;
        else
            DDSIP_bb->expbest = tmpbestvalue;
    }

    if (DDSIP_bb->DDSIP_step == adv)
    {
        if (DDSIP_NoSolution (mipstatus))
        {
            fprintf (DDSIP_outfile, " SOLUTION infeasible.\n");
            if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_bb->moreoutfile, " SOLUTION infeasible.\n");
                printf ("\t\t Solution infeasible.\n");
            }
        }
        else
        {
            fprintf (DDSIP_outfile, " OBJECTIVE OF START VALUE:     %13.7f\n", tmpbestvalue);
            if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_bb->moreoutfile, " OBJECTIVE OF START VALUE:     %13.7f\n", tmpbestvalue);
                printf ("\t\t Best: %13.7f\n", tmpbestvalue);
            }
        }
    }

    DDSIP_Free ((void**) &values);
    DDSIP_Free ((void **) &(subsol));
    DDSIP_Free ((void **) &(mipx));
    
    status = DDSIP_RestoreBoundAndType ();
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to restore bounds\n");
    }
    if (!prematureStop && (DDSIP_bb->nofront == 1) && (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - fabs(DDSIP_bb->bestvalue)*DDSIP_param->relgap))
    {
        // specified gap is reached, further heuristics not necessary
        status = 100000;
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "### specified gap %g was reached\n", DDSIP_param->relgap);
    }
    return status;
}

void DDSIP_EvaluateScenarioSolutions (int* comb)
{
    int i, i_scen, status = 0;
    double tmpbestheur = DDSIP_infty;
    enum DDSIP_step_t currentStep;

    if ((DDSIP_bb->skip == 3) && (DDSIP_node[DDSIP_bb->curnode]->bound >= DDSIP_bb->bestvalue))
        return;
    currentStep = DDSIP_bb->DDSIP_step;
    DDSIP_bb->heurSuccess = 0;
    DDSIP_bb->DDSIP_step = neobj;
    // Initialize, heurval contains the current heuristic solution
    DDSIP_bb->heurval = DDSIP_infty;
    // sort scenarios according to lower bounds for the initial solution
    if (!(DDSIP_bb->ub_sorted))
    {
        // in order to allow for premature cutoff: sort scenarios according to lower bound in root node in descending order
        double * sort_array;
        sort_array = (double *) DDSIP_Alloc(sizeof(double), DDSIP_param->scenarios, "sort_array(LowerBound)");
        for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
        {
            sort_array[DDSIP_bb->ub_scen_order[i_scen]] = DDSIP_data->prob[DDSIP_bb->ub_scen_order[i_scen]] * (DDSIP_node[DDSIP_bb->curnode]->subbound)[DDSIP_bb->ub_scen_order[i_scen]];
        }
        DDSIP_qsort_ins_D (sort_array, DDSIP_bb->ub_scen_order, DDSIP_bb->shifts, DDSIP_param->scenarios-1);
        // insert the least-cost scenario after the first bb->shifts scenarios - to be sure to check that for feasibility, too
        if (DDSIP_bb->shifts < DDSIP_param->scenarios - 2)
        {
            int j;
            i_scen = DDSIP_bb->ub_scen_order[DDSIP_param->scenarios - 1];
            for (j = DDSIP_param->scenarios - 1; j>DDSIP_bb->shifts; j--)
                DDSIP_bb->ub_scen_order[j] = DDSIP_bb->ub_scen_order[j-1];
            DDSIP_bb->ub_scen_order[DDSIP_bb->shifts] = i_scen;
        }
        DDSIP_bb->ub_sorted = 1;

        if (DDSIP_param->outlev > 20)
        {
            // debug output
            fprintf (DDSIP_bb->moreoutfile,"order of scenarios after sorting ub order (%d shifts)\n", DDSIP_bb->shifts);
            for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
            {
                fprintf(DDSIP_bb->moreoutfile," %3d: Scen %3d  %20.15g =  %10.07g * %20.15g\n",
                        i_scen+1, DDSIP_bb->ub_scen_order[i_scen]+1, sort_array[DDSIP_bb->ub_scen_order[i_scen]],
                        DDSIP_data->prob[DDSIP_bb->ub_scen_order[i_scen]],
                        (DDSIP_node[DDSIP_bb->curnode]->subbound)[DDSIP_bb->ub_scen_order[i_scen]]);
            }
            // debug output
        }
        DDSIP_Free ((void**) &sort_array);
    }
    if (DDSIP_param->heuristic == 100)  	// use subsequently different heuristics in the same node
    {
	    for (i = 1; i < DDSIP_param->heuristic_num; i++)
	    {
		    DDSIP_param->heuristic = (int) floor (DDSIP_param->heuristic_vector[i] + 0.1);
		    if (!DDSIP_Heuristics (comb, DDSIP_param->scenarios, 0))
		    {
			    // Evaluate the proposed first-stage solution (if DDSIP_bb->skip was not set)
			    if (DDSIP_bb->skip != -4 && DDSIP_bb->sug[DDSIP_param->nodelim + 2])
			    {
				    if (!(status = DDSIP_UpperBound (DDSIP_param->scenarios, 0)) || status == 100000)
				    {
					    if (DDSIP_bb->heurval < tmpbestheur)
						    tmpbestheur = DDSIP_bb->heurval;
				    }
				    if (status == 100000)
				    {
					    DDSIP_bb->skip = -5;
					    if (DDSIP_param->interrupt_heur > 0)
						    break;
				    }
			    }
		    }
		    if (DDSIP_killsignal)
		    {
			    DDSIP_bb->skip = -5;
			    break;
		    }
	    }
	    // use (expensive) heuristic 12 using all single-scenario solutions as suggestions
	    if (!(DDSIP_param->interrupt_heur && DDSIP_bb->skip == -5) && (DDSIP_bb->heurSuccess || DDSIP_bb->curnode < 13 || DDSIP_bb->noiter%250 > 247))
	    {
		    DDSIP_param->heuristic = 12;
		    if (!DDSIP_Heuristics (comb, DDSIP_param->scenarios, 0))
		    {
			    // Evaluate the proposed first-stage solution (if DDSIP_bb->skip was not set)
			    if (DDSIP_bb->sug[DDSIP_param->nodelim + 2])
			    {
				    if (!DDSIP_UpperBound (DDSIP_param->scenarios, 0))
				    {
					    if (DDSIP_bb->heurval < tmpbestheur)
						    tmpbestheur = DDSIP_bb->heurval;
				    }
			    }
		    }
		    DDSIP_bb->heurval = tmpbestheur;
	    }
	    DDSIP_param->heuristic = 100;
    }
    else if (DDSIP_param->heuristic == 99)  	// use subsequently different heuristics in the same node
    {
	    for (i = 1; i < DDSIP_param->heuristic_num; i++)
	    {
		    DDSIP_param->heuristic = (int) floor (DDSIP_param->heuristic_vector[i] + 0.1);
		    if (!DDSIP_Heuristics (comb, DDSIP_param->scenarios, 0))
		    {
			    // Evaluate the proposed first-stage solution
			    if (!(status = DDSIP_UpperBound (DDSIP_param->scenarios, 0)) || status == 100000)
			    {
				    if (DDSIP_bb->heurval < tmpbestheur)
					    tmpbestheur = DDSIP_bb->heurval;
			    }
			    if (status == 100000 || DDSIP_killsignal)
			    {
				    DDSIP_bb->skip = -5;
				    if (DDSIP_param->interrupt_heur > 0)
					    break;
			    }
		    }
	    }
	    DDSIP_param->heuristic = 99;
	    DDSIP_bb->heurval = tmpbestheur;
    }
    else
    {
	    if (!DDSIP_Heuristics (comb, DDSIP_param->scenarios, 0))
		    // Evaluate the proposed first-stage solution
		    if (!(status = DDSIP_UpperBound (DDSIP_param->scenarios, 0)) || status == 100000)
		    {
			    if (DDSIP_bb->heurval < tmpbestheur)
				    tmpbestheur = DDSIP_bb->heurval;
		    }
	    if (status == 100000 || DDSIP_killsignal)
	    {
		    DDSIP_bb->skip = -5;
	    }
    }
    DDSIP_bb->DDSIP_step = currentStep;
    DDSIP_bb->heurSuccess = -1;
    return;
}
