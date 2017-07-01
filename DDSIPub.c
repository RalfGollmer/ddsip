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

//int DDSIP_RecFun (double *);
int DDSIP_PrintModFileUb (int);
int DDSIP_WarmUb (void);
int DDSIP_GetCpxSolution (int, double *, double *, double *);

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
        else if (DDSIP_param->outlev > 29)
        {
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
                ih = floor ((DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i] + 0.5);
                fprintf (DDSIP_bb->moreoutfile, " %18d,", ih);
            }
            else
                fprintf (DDSIP_bb->moreoutfile, " %18.16g,", (DDSIP_bb->sug[DDSIP_param->nodelim + 2]->firstval)[i]);

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
            if (lhs < currentCut->rhs - 1.e-7)
            {
                *cutViolation = currentCut->rhs - lhs;
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "... violates cut %d, violation %g.\n", currentCut->number, *cutViolation);
                if (DDSIP_param->cpxscr || DDSIP_param->outlev > 10)
                    printf ("... violates cut %d.\n", currentCut->number);
                //if (DDSIP_param->heuristic < 3)
                if (!(DDSIP_param->testOtherScens) || DDSIP_bb->curnode > 6 || DDSIP_param->heuristic < 3)
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
    DDSIP_bb->sug[DDSIP_param->nodelim + 2]->next = DDSIP_bb->sug[DDSIP_bb->curnode];
    DDSIP_bb->sug[DDSIP_bb->curnode] = DDSIP_bb->sug[DDSIP_param->nodelim + 2];
    DDSIP_bb->sug[DDSIP_param->nodelim + 2] = NULL;
    return 1;
}

//==========================================================================
// Functions checks feasibility of the solution gained by Heuristics and
// calculates a new objective value

int
DDSIP_UpperBound (void)
{
    int status, scen, iscen, mipstatus = 0;
    int i, j, k, Bi, Bs, prematureStop = 0;
    int *index;
    int wall_hrs, wall_mins,cpu_hrs, cpu_mins;

    double tmpbestvalue = 0., tmpfeasbound = 0., rest_bound, tmprisk = 0., tmprisk4 = -DDSIP_infty, tmpprob = 0.;
    double bobjval, objval, fs, time_start, time_end, time_lap, wall_secs, cpu_secs, gap, meanGap;

    double *mipx, *values;
    double **tmpsecsol;
    double *subsol;

    double we, wr, d, mipgap, oldviol = DDSIP_infty, viol;
    //  double         *Tx;

    DDSIP_bb->skip = 0;

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
    if (!DDSIP_SolChk (&oldviol))
        return 0;

    if (DDSIP_param->outlev >= DDSIP_suggested_first_stage_outlev && DDSIP_param->outlev < 50)
        // Print output
    {
        int ih;
        fprintf (DDSIP_bb->moreoutfile, "\nSuggested first stage solution:\n");

        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
            {
                ih = floor ((DDSIP_bb->sug[DDSIP_bb->curnode]->firstval)[j] + 0.5);
                fprintf (DDSIP_bb->moreoutfile, " %18d,", ih);
            }
            else
                fprintf (DDSIP_bb->moreoutfile, " %18.16g,", (DDSIP_bb->sug[DDSIP_bb->curnode]->firstval)[j]);

            if (!((j + 1) % 5))
                fprintf (DDSIP_bb->moreoutfile, "\n");
        }
        fprintf (DDSIP_bb->moreoutfile, "\n");
    }

    mipx = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar), "mipx(UpperBound)");
    tmpsecsol = (double **) DDSIP_Alloc (sizeof (double *), DDSIP_bb->secvar, "tmpsecsol(UpperBound)");
    subsol = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "subsol(UpperBound)");
    //Tx = (double *) DDSIP_Alloc(sizeof(double), DDSIP_bb->firstcon + DDSIP_bb->seccon,"(UpperBound)");

    for (j = 0; j < DDSIP_bb->secvar; j++)
        tmpsecsol[j] = (double *) DDSIP_Alloc (sizeof (double), DDSIP_param->scenarios, "tmpsecsol[j](UpperBound)");

    DDSIP_bb->neobjcnt++;

    if (DDSIP_param->outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
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
    // Heuristic returns suggested solution, SolChck puts it in the first element of list sug[DDSIP_bb->curnode]
    // Fix first stage variables to these values
    // Lower bounds
    values = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "values(UpperBound)");
    for(j=0; j<DDSIP_bb->firstvar; j++)
    {
        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
            values[j] = floor((DDSIP_bb->sug[DDSIP_bb->curnode])->firstval[j] + 0.5);
        else
            values[j] = DDSIP_Dmax (DDSIP_bb->lborg[j],(DDSIP_bb->sug[DDSIP_bb->curnode])->firstval[j]);
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
    // Upper bounds
    for(j=0; j<DDSIP_data->firstvar; j++)
    {
        if (DDSIP_bb->firsttype[j] == 'B' || DDSIP_bb->firsttype[j] == 'I' || DDSIP_bb->firsttype[j] == 'N')
            values[j] = floor((DDSIP_bb->sug[DDSIP_bb->curnode])->firstval[j] + 0.5);
        else
            values[j] = DDSIP_Dmin(DDSIP_bb->uborg[j],(DDSIP_bb->sug[DDSIP_bb->curnode])->firstval[j]);
        for (k = 0; k < DDSIP_bb->curbdcnt; k++)
        {
            if (DDSIP_bb->curind[k] == j)
                values[j] = DDSIP_Dmin (values[j], DDSIP_bb->curub[k]);
        }
    }
    status = CPXchgbds (DDSIP_env, DDSIP_lp, fs, DDSIP_bb->firstindex, DDSIP_bb->ubident, values);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change bounds \n");
        goto TERMINATE;
    }

    // prepare for decision about stopping: rest_bound is the expectation of all the lower bounds in this node
    if (DDSIP_bb->DDSIP_step == adv || DDSIP_bb->DDSIP_step == eev)
        rest_bound = -DDSIP_infty;
    else if(DDSIP_param->cb)
    {
        // CAUTION: subboundNoLag is not really a secure lower bound for the problem without Lagrangean!
        cpu_secs = 0.0;
        for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
        {
            cpu_secs += (DDSIP_node[DDSIP_bb->curnode]->subboundNoLag)[iscen] * DDSIP_data->prob[iscen];
        }
        rest_bound = DDSIP_Dmin(cpu_secs, DDSIP_node[DDSIP_bb->curnode]->bound);
    }
    else
    {
        rest_bound = DDSIP_node[DDSIP_bb->curnode]->bound;
    }

    // UpperBound single-scenario problems
    meanGap = tmpprob = tmpbestvalue = 0.;
    for (iscen = 0; iscen < DDSIP_param->scenarios; iscen++)
    {
        scen = DDSIP_bb->ub_scen_order[iscen];
        if (DDSIP_killsignal)
        {
            status = 0;
            goto TERMINATE;
        }

        tmpprob += DDSIP_data->prob[scen];

        if (DDSIP_param->cpxubscr)
        {
            printf ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n");
            printf ("Calculating objective value for scenario problem %d in node %d (heuristic %d).....\n", scen + 1,DDSIP_bb->curnode,DDSIP_param->heuristic);
        }

        status = DDSIP_ChgProb (scen);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change problem \n");
            goto TERMINATE;
        }
        // Calculate Tx+h
        //      if (DDSIP_param->outlev)  status = RecFun(param,bb,Tx);

        // Debugging information
        DDSIP_PrintModFileUb (scen);

        // Warm start
        if (DDSIP_param->hot)
        {
            status = DDSIP_WarmUb ();
            if (status)
                goto TERMINATE;
        }

        // query time limit amd mip rel. gap parameters
        if (DDSIP_param->cpxubscr || DDSIP_param->outlev > 21)
        {
            status = CPXgetdblparam (DDSIP_env,CPX_PARAM_TILIM,&we);
            status = CPXgetdblparam (DDSIP_env,CPX_PARAM_EPGAP,&wr);
            printf ("   -- 1st optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
        }
        //
        time_start = DDSIP_GetCpuTime ();
        // Optimize
        status = CPXmipopt (DDSIP_env, DDSIP_lp);
        if (DDSIP_Error (status))
        {
            fprintf (stderr, "ERROR: Failed to optimze (UB)\n");
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
            if (DDSIP_param->cpxubscr ||  DDSIP_param->outlev > 11)
            {
                j = CPXgetnodecnt (DDSIP_env,DDSIP_lp);
                printf ("      UB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,j,time_lap-time_start);
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile,"      UB: after 1st optimization: mipgap %% %-12lg %7d nodes  (%6.2fs)\n",mipgap*100.0,j,time_lap-time_start);
            }
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
                rest_bound -= (DDSIP_node[DDSIP_bb->curnode]->subboundNoLag)[scen] * DDSIP_data->prob[scen];
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
                    goto TERMINATE;
                }
            }
            if (DDSIP_param->cpxnoub2 && status != CPXMIP_OPTIMAL)
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
                    if (DDSIP_param->cpxubscr || DDSIP_param->outlev > 21)
                    {
                        printf ("   -- 2nd optimization time limit: %gs, rel. gap: %g%% --\n",we,wr*100.0);
                    }
                    status = CPXmipopt (DDSIP_env, DDSIP_lp);
                    mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);
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

        if (DDSIP_NoSolution (mipstatus))
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
#ifdef ADDBENDERSCUTS
            if (DDSIP_param->addBendersCuts && (DDSIP_param->heuristic > 3 || iscen))
            {
                CPXLPptr     DDSIP_dual_lp  = NULL;
                cutpool_t * newCut;
                int end;
                end = (DDSIP_param->addBendersCuts > 1 || (DDSIP_bb->curnode < 7 && DDSIP_param->testOtherScens)) ? DDSIP_param->scenarios : iscen+1;
                for (Bi = iscen; Bi < end; Bi++)
                {
                    Bs = DDSIP_bb->ub_scen_order[Bi];
                    status = DDSIP_ChgProb (Bs);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to change problem \n");
                        goto TERMINATE;
                    }
                    DDSIP_dual_lp = CPXcloneprob (DDSIP_env, DDSIP_lp, &status);
                    if (status)
                    {
                        printf ("ERROR: Failed to clone problem for Benders.\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to clone problem for Benders.\n");
                        DDSIP_param->addBendersCuts = 0;
                    }
                    else
                    {
                        int numRows = CPXgetnumrows (DDSIP_env, DDSIP_dual_lp);
                        status = CPXchgprobtype (DDSIP_env, DDSIP_dual_lp, CPXPROB_LP);
                        if (status)
                        {
                            printf ("ERROR: Failed to change type of problem for Benders to LP.\n");
                            fprintf (DDSIP_outfile, "ERROR: Failed to change type of problem for Benders to LP.\n");
                            return status;
                        }
                        status = CPXsetintparam (DDSIP_env, CPX_PARAM_PREIND, CPX_OFF);
                        status = CPXdualopt (DDSIP_env, DDSIP_dual_lp);
                        if (status)
                        {
                            printf ("ERROR: dualopt of LP for Benders failed.\n");
                            fprintf (DDSIP_outfile, "ERROR: dualopt of LP for Benders failed.\n");
                            return status;
                        }
                        status = CPXgetstat (DDSIP_env, DDSIP_dual_lp);
                        if (status == CPX_STAT_INFEASIBLE)
                        {
                            double *ray = (double *) DDSIP_Alloc (sizeof (double), numRows, "ray (UpperBound)");
                            status = CPXdualfarkas (DDSIP_env, DDSIP_dual_lp, ray, &viol);
                            if (!status && (viol > 1e-1) && (fabs(viol - oldviol)/viol > 1.1e-8))
                            {
                                int rmatbeg;
                                int *rmatind = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->firstvar + 1), "rmatind(UpperBound)");
                                double *rmatval = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_data->firstvar + 1), "rmatval (UpperBound)");
                                double rhs;
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
                                    if (ray[DDSIP_bb->secondrowind[k]] != 0.)
                                        break;
                                }
                                if (k < DDSIP_bb->seccon)
                                {
                                    DDSIP_bb->cutCntr++;
                                    rowname[0] = rowstore;
                                    sprintf (rowstore, "DDSIPBendersCut%d",DDSIP_bb->cutCntr);
                                    if (DDSIP_param->outlev)
                                        fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s  (from scen %2d), violation %g ############\n", rowstore, Bs+1, viol);
                                    printf (" ############ adding cut %s  (from scen %2d) ############\n", rowstore, Bs+1);
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
                                            if (ray[i] != 0.)
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
                                    for (k = 0; k < DDSIP_data->firstvar; k++)
                                    {
                                        lhs += rmatval[k]*values[k];
                                    }
                                    rhs = lhs + viol*0.999999925;
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
                                    newCut->number = DDSIP_bb->cutCntr;
                                    DDSIP_bb->cutpool = newCut;
                                    rmatval = NULL;
                                }

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
                            if (DDSIP_param->outlev > 23)
                                fprintf (DDSIP_bb->moreoutfile," ------------ LP feasible ---------------\n");
                        }
                    }
                    if (DDSIP_dual_lp != NULL)
                    {
                        status = CPXfreeprob (DDSIP_env, &DDSIP_dual_lp);
                        if (status)
                            printf ("ERROR: CPXfreeprob failed, error code %d\n", status);
                    }
                }
            }
#endif
#ifdef ADDINTEGERCUTS
            //if all first-stage variables are binary ones, we can add an inequality, cutting off this point
            if (DDSIP_param->addIntegerCuts)
            {
                int rmatbeg;
                int *rmatind = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_bb->novar), "rmatind(UpperBound)");
                double *rmatval = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_bb->novar), "rmatval (UpperBound)");
                double rhs, objv = DDSIP_infty;
                char sense;
                char **rowname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(UpperBound)");
                char *rowstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(UpperBound)");
#ifdef CHECKINTEGERCUT
                CPXLPptr     DDSIP_dual_lp  = NULL;
#endif
                rowname[0] = rowstore;
                rhs = 1.;
                sense = 'G';
                rmatbeg = 0;
                for (k = 0; k < DDSIP_bb->firstvar; k++)
                {
                    rmatind[k] = DDSIP_bb->firstindex[k];
                    if (floor((DDSIP_bb->sug[DDSIP_bb->curnode])->firstval[k] + 0.5))
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
                // check whether this cut is dominated
                status = DDSIP_RestoreBoundAndType ();
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to change bounds for check, error code %d\n", status);
                    goto TERMINATE;
                }
                status = DDSIP_ChgProb (scen);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to change problem for check to scenario %d, error code %d\n", scen + 1, status);
                    goto TERMINATE;
                }
#ifdef CHECKINTEGERCUT
                DDSIP_dual_lp = CPXcloneprob (DDSIP_env, DDSIP_lp, &status);
                if (status)
                {
                    printf ("ERROR: Failed to clone problem for check.\n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to clone problem for check.\n");
                }
                else
                {
                    status = CPXchgprobtype (DDSIP_env, DDSIP_dual_lp, CPXPROB_LP);
                    if (status)
                    {
                        printf ("ERROR: Failed to change type of problem for check to LP.\n");
                        fprintf (DDSIP_outfile, "ERROR: Failed to change type of problem for check to LP.\n");
                        return status;
                    }
                    status = CPXchgobj (DDSIP_env, DDSIP_dual_lp, DDSIP_bb->novar, rmatind, rmatval);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to update objective coefficients\n");
                        return status;
                    }
                    status = CPXdualopt (DDSIP_env, DDSIP_dual_lp);
                    if (status)
                    {
                        printf ("ERROR: dualopt of LP for check failed.\n");
                        fprintf (DDSIP_outfile, "ERROR: dualopt of LP for check failed.\n");
                        return status;
                    }
                    status = CPXgetobjval (DDSIP_env, DDSIP_dual_lp, &objv);
                    if (objv < rhs - 3e-2)
                    {
                        cutpool_t * newCut;
                        DDSIP_bb->cutCntr++;
                        sprintf (rowstore, "DDSIPIntegerCut%d",DDSIP_bb->cutCntr);
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s  (objval = %g < %g) ############\n", rowstore, objv, rhs);
                        printf (" ############ adding cut %s  (objval = %g < %g) ############\n", rowstore, objv, rhs);
#else
                        cutpool_t * newCut;
                        DDSIP_bb->cutCntr++;
                        sprintf (rowstore, "DDSIPIntegerCut%d",DDSIP_bb->cutCntr);
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile," ############ adding cut %s ############\n", rowstore);
                        printf (" ############ adding cut %s  (objval = %g < %g) ############\n", rowstore, objv, rhs);
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
                        newCut->number = DDSIP_bb->cutCntr;
                        DDSIP_bb->cutpool = newCut;
                        rmatval = NULL;
#ifdef CHECKINTEGERCUT
                    }
                    else
                    {
                        if (DDSIP_param->outlev > 20)
                            fprintf (DDSIP_bb->moreoutfile," ############ integer cut is not effective: objval = %g, rhs = %g ############\n", objv, rhs);
                        printf (" ############ integer cut is not effective: objval = %g, rhs = %g ############\n", objv, rhs);
                    }
                }
                if (DDSIP_dual_lp != NULL)
                {
                    status = CPXfreeprob (DDSIP_env, &DDSIP_dual_lp);
                    if (status)
                        printf ("ERROR: CPXfreeprob failed, error code %d\n", status);
                }
#endif
                DDSIP_Free ((void *) &rmatind);
                DDSIP_Free ((void *) &rmatval);
                DDSIP_Free ((void *) &rowname);
                DDSIP_Free ((void *) &rowstore);
            }
#endif
            //shift this infeasible scenario to first place, such that next time it is checked first
            if (iscen >= DDSIP_bb->shifts)
                DDSIP_bb->shifts++;
            k = DDSIP_bb->ub_scen_order[iscen];
            for(j = iscen; j>0; j--)
                DDSIP_bb->ub_scen_order[j] = DDSIP_bb->ub_scen_order[j-1];
            DDSIP_bb->ub_scen_order[0] = k;
            goto TERMINATE;
        }
        // Get solution
        status = DDSIP_GetCpxSolution (mipstatus, &objval, &bobjval, mipx);
        if (status)
            goto TERMINATE;

        // Remember objective function contribution
        if (DDSIP_param->outlev > 10)
            DDSIP_Contrib (mipx, scen);

        // Store second stage solutions
        for (j = 0; j < DDSIP_bb->secvar; j++)
            tmpsecsol[j][scen] = mipx[DDSIP_bb->secondindex[j]];

        if ( mipstatus == 101 && bobjval > DDSIP_infty )
            bobjval=objval;

        subsol[scen] = objval;

        gap = 100.0*(objval-bobjval)/(fabs(objval)+1e-4);
        meanGap += DDSIP_data->prob[scen] * gap;
        // Debugging information
        if (DDSIP_param->outlev)
        {
            time (&DDSIP_bb->cur_time);
            time_end = DDSIP_GetCpuTime ();
            time_start = time_end-time_start;
            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
            DDSIP_translate_time (time_end,&cpu_hrs,&cpu_mins,&cpu_secs);
            fprintf (DDSIP_bb->moreoutfile,
                     "%4d Scenario %4.0d:  Best=%-20.14g\tBound=%-20.14g\t(%9.4g%%)\tStatus=%3.0d\t%3dh %02d:%02.0f cpu %3dh %02d:%05.2f (%6.2f)\n",
                     iscen + 1, scen + 1, objval, bobjval, gap, mipstatus,
                     wall_hrs,wall_mins,wall_secs,cpu_hrs,cpu_mins,cpu_secs, time_start);
            if (DDSIP_param->outlev>7)
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
    }				// end for iscen=..

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
        DDSIP_bb->bestvalue = tmpbestvalue;
        DDSIP_bb->feasbound = tmpfeasbound;
        DDSIP_bb->bestrisk = DDSIP_bb->currisk;
        DDSIP_bb->bestexp = DDSIP_bb->curexp;

        // We store all risk measures
        if (DDSIP_param->riskmod)
            for (j = 0; j < DDSIP_maxrisk; j++)
                DDSIP_bb->bestriskval[j] = DDSIP_bb->curriskval[j];

        for (j = 0; j < DDSIP_bb->firstvar; j++)
            DDSIP_bb->bestsol[j] = ((DDSIP_bb->sug[DDSIP_bb->curnode])->firstval)[j];

        if (tmpbestvalue < DDSIP_node[DDSIP_bb->curnode]->bound)
        {
            // something's wrong with the accuracy of solutions - the heuristic gave a better value than the lower bound in the same node!
          printf ("*\t*** WARNING: heuristic value is better than lower bound of the node by %g, MIP gaps or tolerances probably not small enough ***\n", DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestvalue);
            if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_outfile, "*\t*** WARNING: heuristic value is better than lower bound of the current node by %g, MIP gaps or tolerances probably not small enough ***\n", DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestvalue);
                fprintf (DDSIP_bb->moreoutfile, "*\t*** WARNING: heuristic value is better than lower bound of the current node by %g, MIP gaps or tolerances probably not small enough ***\n", DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestvalue);
            }
            DDSIP_bb->correct_bounding = DDSIP_Dmax (DDSIP_bb->correct_bounding,  DDSIP_node[DDSIP_bb->curnode]->bound - tmpbestvalue);
        }

        for (j = 0; j < DDSIP_bb->secvar; j++)
            for (k = 0; k < DDSIP_param->scenarios; k++)
                DDSIP_bb->secstage[j][k] = tmpsecsol[j][k];
        for (j = 0; j < DDSIP_param->scenarios; j++)
            DDSIP_bb->subsol[j] = subsol[j];
        // Print a line of output
        DDSIP_PrintStateUB (DDSIP_param->heuristic);
        // heuristic was successful
        DDSIP_bb->heurSuccess++;
    }


TERMINATE:

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
            printf ("\t\t Solution infeasible.\n");
            fprintf (DDSIP_outfile, " SOLUTION infeasible.\n");
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, " SOLUTION infeasible.\n");
        }
        else
        {
            printf ("\t\t Best: %13.7f\n", tmpbestvalue);
            fprintf (DDSIP_outfile, " OBJECTIVE OF START VALUE:     %13.7f\n", tmpbestvalue);
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, " OBJECTIVE OF START VALUE:     %13.7f\n", tmpbestvalue);
        }
    }

    DDSIP_Free ((void**) &values);
    DDSIP_Free ((void **) &(subsol));
    for (j = 0; j < DDSIP_bb->secvar; j++)
        DDSIP_Free ((void **) &(tmpsecsol[j]));
    DDSIP_Free ((void **) &(tmpsecsol));
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
        fprintf (stderr, "specified gap was reached.\n");
    }
    return status;
}

void DDSIP_EvaluateScenarioSolutions (void)
{
    int comb, i, i_scen, status = 0;
    double tmpbestheur = DDSIP_infty;

    // Initialize, heurval contains the current heuristic solution
    DDSIP_bb->heurval = DDSIP_infty;
    // sort scenarios according to lower bounds for the initial solution
    if (!DDSIP_bb->curnode)
    {
        // in order to allow for premature cutoff: sort scenarios according to lower bound in root node in descending order
        double * sort_array;
        sort_array = (double *) DDSIP_Alloc(sizeof(double), DDSIP_param->scenarios, "sort_array(LowerBound)");
        for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
        {
            sort_array[i_scen] = DDSIP_data->prob[DDSIP_bb->ub_scen_order[i_scen]] * (DDSIP_node[DDSIP_bb->curnode]->subbound)[DDSIP_bb->ub_scen_order[i_scen]];
        }
        DDSIP_qsort_ins_D (sort_array, DDSIP_bb->ub_scen_order, DDSIP_bb->shifts, DDSIP_param->scenarios-1);

        if (DDSIP_param->outlev > 30)
        {
            // debug output
            fprintf (DDSIP_bb->moreoutfile,"order of scenarios after sorting ub order (%d shifts)\n", DDSIP_bb->shifts);
            for (i_scen = 0; i_scen < DDSIP_param->scenarios; i_scen++)
            {
                fprintf(DDSIP_bb->moreoutfile," %3d: Scen%3d  %g =  %g * %g\n",
                        i_scen+1, DDSIP_bb->ub_scen_order[i_scen]+1, sort_array[DDSIP_bb->ub_scen_order[i_scen]],
                        DDSIP_data->prob[DDSIP_bb->ub_scen_order[i_scen]],
                        (DDSIP_node[DDSIP_bb->curnode]->subbound)[DDSIP_bb->ub_scen_order[i_scen]]);
            }
            // debug output
        }
        DDSIP_Free ((void**) &sort_array);
    }
    DDSIP_bb->cutAdded = 0;
    if (DDSIP_param->heuristic == 100)  	// use subsequently different heuristics in the same node
    {
        for (i = 1; i < DDSIP_param->heuristic_num; i++)
        {
            DDSIP_param->heuristic = floor (DDSIP_param->heuristic_vector[i] + 0.1);
            if (!DDSIP_Heuristics (&comb))
            {
                // Evaluate the proposed first-stage solution (if DDSIP_bb->skip was not set)
                if (DDSIP_bb->skip != -4 && DDSIP_bb->sug[DDSIP_param->nodelim + 2])
                {
                    if (!(status = DDSIP_UpperBound ()) || status == 100000)
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
        }
        DDSIP_param->heuristic = 12;
        if (!(DDSIP_param->interrupt_heur && (DDSIP_bb->skip == -5)))
        {
            if (!DDSIP_Heuristics (&comb))
            {
                // Evaluate the proposed first-stage solution (if DDSIP_bb->skip was not set)
                if (DDSIP_bb->sug[DDSIP_param->nodelim + 2])
                {
                    if (!DDSIP_UpperBound ())
                    {
                        if (DDSIP_bb->heurval < tmpbestheur)
                            tmpbestheur = DDSIP_bb->heurval;
                    }
                }
            }
        }
        DDSIP_param->heuristic = 100;
        DDSIP_bb->heurval = tmpbestheur;
    }
    else if (DDSIP_param->heuristic == 99)  	// use subsequently different heuristics in the same node
    {
        for (i = 1; i < DDSIP_param->heuristic_num; i++)
        {
            DDSIP_param->heuristic = floor (DDSIP_param->heuristic_vector[i] + 0.1);
            if (!DDSIP_Heuristics (&comb))
            {
                // Evaluate the proposed first-stage solution
                if (!(status = DDSIP_UpperBound ()) || status == 100000)
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
        DDSIP_param->heuristic = 99;
        DDSIP_bb->heurval = tmpbestheur;
    }
    else
    {
        if (!DDSIP_Heuristics (&comb))
            // Evaluate the proposed first-stage solution
            if (!(status = DDSIP_UpperBound ()) || status == 100000)
            {
                if (DDSIP_bb->heurval < tmpbestheur)
                    tmpbestheur = DDSIP_bb->heurval;
            }
            if (status == 100000)
            {
                DDSIP_bb->skip = -5;
            }
    }
    return;
}
