/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
	Description:
	This file contains the implementation of the dual method, i.e. all
	functions required by ConicBundle of C. Helmberg.

	License:
	This file is part of DDSIP.

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

#ifdef CONIC_BUNDLE
// compile only in case ConicBundle ist available

//#define MIX 0
//#define MIX 1
#define MIX 2
#define MIXFACTOR 0.5
#define CANCELLATION
//#define ONLYHEUR12
//#define DEBUG
//#define MDEBUG

#include <DDSIP.h>

#include <DDSIPconst.h>

//==========================================================================
int
DDSIP_NonAnt (void)
{
    int status = 0, scen, i, k, j;

    DDSIP_data->nabeg = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar * DDSIP_param->scenarios, "nabeg(NonAnt)");
    DDSIP_data->nacnt = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar * DDSIP_param->scenarios, "nacnt(NonAnt)");

    if (DDSIP_param->nonant == 3)
    {
        DDSIP_data->naind =
            (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->firstvar * DDSIP_param->scenarios * (DDSIP_param->scenarios - 1), "naind(NonAnt)");
        DDSIP_data->naval =
            (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar * DDSIP_param->scenarios * (DDSIP_param->scenarios - 1), "naval(NonAnt)");
    }
    else
    {
        DDSIP_data->naind = (int *) DDSIP_Alloc (sizeof (int), 2 * DDSIP_bb->firstvar * DDSIP_param->scenarios, "naind(NonAnt)");
        DDSIP_data->naval = (double *) DDSIP_Alloc (sizeof (double), 2 * DDSIP_bb->firstvar * DDSIP_param->scenarios, "naval(NonAnt)");
    }

    // x1=x2, x2=x3 ... x(S-1)=xS
    if (DDSIP_param->nonant == 2)
    {
        k = 0;
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] = k;
                if (scen && scen < DDSIP_param->scenarios - 1)
                    DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i] = 2;
                else
                    DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i] = 1;

                if (scen)
                {
                    DDSIP_data->naind[k] = (scen - 1) * DDSIP_bb->firstvar + i;
                    DDSIP_data->naval[k++] = -1.0;
                }
                if (scen < DDSIP_param->scenarios - 1)
                {
                    DDSIP_data->naind[k] = scen * DDSIP_bb->firstvar + i;
                    DDSIP_data->naval[k++] = 1.0;
                }
            }
    }
    // x_j = sum p_i x_i
    else if (DDSIP_param->nonant == 3)
    {
        k = 0;
        for (scen = 0; scen < DDSIP_param->scenarios; scen++)
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] = k;
                DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i] = DDSIP_param->scenarios - 1;

                for (j = 0; j < DDSIP_param->scenarios-1; j++)
                {
                    DDSIP_data->naind[k] = j * DDSIP_bb->firstvar + i;
                    if (j == scen)
                        DDSIP_data->naval[k++] = 1 - DDSIP_data->prob[scen];
                    else
                        DDSIP_data->naval[k++] = -DDSIP_data->prob[scen];
                }
            }
    }
    // x1=x2, x1=x3 ... x1=xS
    else
    {
        k = 0;
        // First scenario
        for (i = 0; i < DDSIP_bb->firstvar; i++)
        {
            DDSIP_data->nabeg[i] = k;
            DDSIP_data->nacnt[i] = DDSIP_param->scenarios - 1;

            for (scen = 0; scen < DDSIP_param->scenarios - 1; scen++)
            {
                DDSIP_data->naind[k] = scen * DDSIP_bb->firstvar + i;
                DDSIP_data->naval[k++] = 1.0;
            }
        }
        // Other scenarios
        for (i = 0; i < DDSIP_bb->firstvar; i++)
            for (scen = 1; scen < DDSIP_param->scenarios; scen++)
            {
                DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] = k;
                DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i] = 1;
                DDSIP_data->naind[k] = (scen - 1) * DDSIP_bb->firstvar + i;
                DDSIP_data->naval[k++] = -1.0;
            }
    }

#ifdef DEBUG
    if (DDSIP_param->outlev > 97)
    {
        int rows = 0, r, surplus;
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        for (i = 0; i < DDSIP_bb->firstvar * DDSIP_param->scenarios; i++)
            for (k = DDSIP_data->nabeg[i]; k < DDSIP_data->nabeg[i] + DDSIP_data->nacnt[i]; k++)
            {
                if (DDSIP_data->naind[k] > rows)
                    rows = DDSIP_data->naind[k];
            }
        fprintf (DDSIP_bb->moreoutfile, "\nNonanticipativity constraints:\n");
        fprintf (DDSIP_bb->moreoutfile, "( %d rows )\n",rows+1);
        for (r=0; r<=rows; r++)
        {
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                if (CPXgetcolname(DDSIP_env, DDSIP_lp, DDSIP_bb->name_buffer, DDSIP_bb->n_buffer, DDSIP_bb->n_buffer_len, &surplus, DDSIP_bb->firstindex[i], DDSIP_bb->firstindex[i]))
                {
                    fprintf (stderr,"Error in querying columns name for index %d\n",DDSIP_bb->firstindex[i]);
                    exit (1);
                }
                for (scen = 0; scen < DDSIP_param->scenarios; scen++)
                {
                    for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                            j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                    {
                        if (DDSIP_data->naind[j] == r)
                        {
                            if (DDSIP_data->naval[j] >= 0.)
                                fprintf (DDSIP_bb->moreoutfile, " +%g * %s_SC%.2d",DDSIP_data->naval[j],DDSIP_bb->n_buffer,scen);
                            else
                                fprintf (DDSIP_bb->moreoutfile, " %g * %s_SC%.2d",DDSIP_data->naval[j],DDSIP_bb->n_buffer,scen);
                        }
                    }
                }
            }
            fprintf (DDSIP_bb->moreoutfile, " = 0\n");
        }
    }
#endif
    return status;
}

//==========================================================================
// This function is called by Helmbergs conic bundle method to update the subgradient
// and the objective value
// Note: We maximize with respect to the Lagrangian multipliers, Helmberg minimizes
int
DDSIP_DualUpdate (void *function_key, double *dual, double relprec,
                  int max_new_subg, double *objective_value, int *new_subg,
                  double *subgval, double *subgradient, double *primal)
{
    int status, scen, i, j;
    int k;
    double h;
#ifdef CANCELLATION
    double * maxval = NULL, h_min;
#endif

    if (DDSIP_killsignal)
    {
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\nTermination signal received.\n");
        DDSIP_bb->dualObjVal = -DDSIP_infty;
        *new_subg = 0;
        return -1;
    }

    // update node->dual
    memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, dual, sizeof (double) * DDSIP_bb->dimdual);

    if (DDSIP_param->outlev > 25)
    {
        fprintf (DDSIP_bb->moreoutfile, "\nDualUpdate: Current lambda for node %d: (%p)\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
        for (i = 0; i < DDSIP_bb->dimdual; i++)
        {
            fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
            if (!((i+1)%10))
                fprintf (DDSIP_bb->moreoutfile, "\n");
        }
        fprintf(DDSIP_bb->moreoutfile, "\n");
    }
    i = 0;
    // LowerBound scenario problems
ONCE_AGAIN:
    if ((status = DDSIP_CBLowerBound (objective_value, relprec)) > 1 || status == -111)
    {
        if (DDSIP_param->outlev)
        {
            printf ("Failed to solve scenario problems: status = %d .\n", status);
            fprintf (DDSIP_bb->moreoutfile, "############ Failed to solve scenario problems: status = %d .\n", status);
        }
        *new_subg = 0;
        if (status != -111)
            return status;
        else
        {
            if (i)
            {
                *new_subg = 0;
                return -1;
            }
            if (!i && DDSIP_bb->dualitcnt<2 && DDSIP_bb->curnode > 0)
            {
                if (DDSIP_param->outlev)
                {
                    printf ("Possibly due to unbounded variables. Trying zero Lagrange multipliers.\n");
                    fprintf (DDSIP_bb->moreoutfile, "Possibly due to unbounded variables. Trying zero Lagrange multipliers.\n");
                }
                memset (DDSIP_node[DDSIP_bb->curnode]->dual, '\0', sizeof (double) * DDSIP_bb->dimdual);
                i = 1;
                goto ONCE_AGAIN;
            }
            else
            {
                if (DDSIP_param->outlev)
                {
                    printf ("Possibly due to unbounded variables. Trying smaller stepsize = higher weight.\n");
                    fprintf (DDSIP_bb->moreoutfile, "Possibly due to unbounded variables. Trying smaller stepsize = higher weight.\n");
                }
                *new_subg = 0;
                return -111;
            }
        }
    }

    DDSIP_bb->dualitcnt++;

    // Calculate new subgradient H * x
    // Initialization
    memset (subgradient, '\0', sizeof (double) * DDSIP_bb->dimdual); 

    if (DDSIP_param->scalarization && DDSIP_bb->skip != 2)
    {
        if (DDSIP_param->outlev > 29)
        {
            printf ("DDSIP_bb->ref_max= %d, DDSIP_bb->risk:", DDSIP_bb->ref_max);
            for (scen = 0; scen < DDSIP_param->scenarios; scen++)
                printf (" %g|", DDSIP_bb->ref_risk[scen]);
        }
        for (i = 0; i < DDSIP_bb->firstvar; i++)
        {
            for (scen = 0; scen < DDSIP_param->scenarios; scen++)
            {
                for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                        j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                {
                    subgradient[DDSIP_data->naind[j]] -= DDSIP_data->naval[j] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][i] / DDSIP_param->ref_scale[0];
                    if (DDSIP_param->outlev > 29)
                        printf (" Subgradient[%d]: Scenario %d, Variable %d (Werte: %g): -= %g\n", DDSIP_data->naind[j], scen, i,
                                (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][i],
                                DDSIP_data->naval[j] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][i] / DDSIP_param->ref_scale[0]);
                }
            }
        }
        if (DDSIP_param->outlev > 29)
        {
            printf ("\n");
            printf ("   curexp=%g, currisk=%g\n", DDSIP_bb->curexp, DDSIP_bb->currisk);
        }
    }
    else if (DDSIP_bb->skip != 2)
    {
#ifdef CANCELLATION
        maxval = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->dimdual, "maxval(DualUpdate)");
#endif
        if (DDSIP_bb->violations)
        {
            for (i = 0; i < DDSIP_bb->firstvar; i++)
                for (scen = 0; scen < DDSIP_param->scenarios; scen++)
                    for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                            j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                    {
			k = DDSIP_data->naind[j];
                        h = DDSIP_data->naval[j] * (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][i];
                        subgradient[k] -= h;
#ifdef CANCELLATION
                        maxval[k] = DDSIP_Dmax (maxval[k], fabs (h));
#endif
                    }
        }
#ifdef CANCELLATION
        {
            h = 0.;
            h_min = DDSIP_infty;
            i = 0;
            for (k = 0; k < DDSIP_bb->dimdual; k++)
            {
                if (subgradient[k] && fabs (subgradient[k]) < 4.e-15 * maxval[k])
                {
                    i++;
                    h = DDSIP_Dmax (h, fabs(subgradient[k])/(maxval[k]));
                    h_min = DDSIP_Dmin (h_min, fabs(subgradient[k])/(maxval[k]));
                    subgradient[k] = 0.;
                }
            }
            DDSIP_Free ((void **) &maxval);
            if (i && DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, "   ## it %3d: %3d  cancellations, max.ratio= %g, min.ratio= %g\n", DDSIP_bb->dualitcnt, i, h, h_min);
        }
#endif
    }
    //
    subgval[0] = *objective_value;
    if (!(DDSIP_bb->dualdescitcnt))
        DDSIP_bb->dualObjVal =  -(*objective_value);
    if (DDSIP_param->outlev > DDSIP_first_stage_outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n After scenarios solving: (Total it. %d in node %d) ", DDSIP_bb->dualitcnt, DDSIP_bb->curnode);
        fprintf (DDSIP_bb->moreoutfile, "\n SC  VAR     FIRSTSOL         LAMBDA       returned subgradient (=H*FIRSTSOL)\n");
        for (scen = 0; scen < DDSIP_param->scenarios - 1; scen++)
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                fprintf (DDSIP_bb->moreoutfile,
                         "%3d  %3d  %15.10f   %15.10f     %22.14g\n", scen + 1,
                         i + 1, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[scen][i],
                         dual[scen * DDSIP_bb->firstvar + i], subgradient[scen * DDSIP_bb->firstvar + i]);
            }
        if (DDSIP_param->nonant == 3)
        {
            scen = DDSIP_param->scenarios - 1;
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                fprintf (DDSIP_bb->moreoutfile, "%3d  %3d  %15.10f   %15.10f     %22.14g\n",
                         DDSIP_param->scenarios, i + 1, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[DDSIP_param->scenarios - 1][i],
                         dual[scen * DDSIP_bb->firstvar + i], subgradient[scen * DDSIP_bb->firstvar + i]);
            }
        }
        else
        {
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                fprintf (DDSIP_bb->moreoutfile, "%3d  %3d  %15.10f\n",
                         DDSIP_param->scenarios, i + 1, (DDSIP_node[DDSIP_bb->curnode]->first_sol)[DDSIP_param->scenarios - 1][i]);
            }
        }
        fprintf (DDSIP_bb->moreoutfile," Dual objective    = %18.12f\n", *objective_value);
        fprintf (DDSIP_bb->moreoutfile,"-DDSIP_node->bound = %18.12f\n", -DDSIP_node[DDSIP_bb->curnode]->bound);
        printf (" Dual objective    = %18.12f\n", *objective_value);
        printf ("-DDSIP_node->bound = %18.12f\n", -DDSIP_node[DDSIP_bb->curnode]->bound);
        fprintf (DDSIP_bb->moreoutfile, "\n");
        printf ("\n");
    }
//#ifdef DEBUG
    if (DDSIP_param->outlev > 20)
    {
        double h = 0., h1, hmin = DDSIP_infty, hmax = -DDSIP_infty;
        int nonzeros = 0;
        for (scen = 0; scen < DDSIP_param->scenarios - 1; scen++)
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                h1 = fabs(subgradient[scen * DDSIP_bb->firstvar + i]);
                h += h1*h1;
                if (h1)
                {
                    nonzeros++;
                    hmin = DDSIP_Dmin (hmin, h1);
                    hmax = DDSIP_Dmax (hmax, h1);
                }
            }
        h = sqrt (h);
        fprintf (DDSIP_bb->moreoutfile, "   ## 2-norm of subgradient: %.8g,   #nonzeros: %d, abs. values in [%g, %g]\n", h, nonzeros, hmin, hmax);
    }
//#endif
    *new_subg = 1;
    if (DDSIP_killsignal)
    {
        fprintf (DDSIP_outfile, "\nTermination signal received.\n");
        if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\nTermination signal received.\n");
        DDSIP_bb->dualObjVal = -DDSIP_infty;
        DDSIP_bb->skip = 1;
    }
    return 0;
}

//==========================================================================
// Main procedure of Helmbergs conic bundle implementation
int
DDSIP_DualOpt (void)
{
    cb_problemp p;
    int cb_status, status, i_scen, j, cnt, noIncreaseCounter = 0, comb = 3, init_iters = 1;
    double *minfirst;
    double *maxfirst;
    double *center_point;
    double obj, old_obj, start_weight, diff;
    int    wall_hrs, wall_mins, cpu_hrs, cpu_mins, limits_reset, cur_iters, current_maxsteps, repeated_increase = 2, weight_decreases = 0, many_iters = 0, cycleCnt  = 0;
    double wall_secs, cpu_secs, inherited_bound, inhMult_bound = -1e+20, rgap;
    double old_cpxrelgap = 1.e-16, old_cpxtimelim = 1000000., old_cpxrelgap2 = 1.e-16, old_cpxtimelim2 = 1000000., last_weight, next_weight;
    double reduction_factor = 0.50;
    double nfactor;
    minfirst = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
               "minfirst(DualOpt)");
    maxfirst = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar,
               "maxfirst(DualOpt)");
    center_point = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->dimdual,
                           "center_point(DualOpt)");

    DDSIP_bb->DDSIP_step = dual;
    if (DDSIP_bb->curnode >= DDSIP_param->cbBreakIters)
        DDSIP_bb->current_itlim = DDSIP_param->cbitlim;
    else if (!DDSIP_bb->curnode)
        DDSIP_bb->current_itlim = DDSIP_param->cbrootitlim;
    else
        DDSIP_bb->current_itlim = (DDSIP_node[DDSIP_bb->curnode]->depth <= DDSIP_param->cb_depth)?
                                   DDSIP_param->cb_depth_iters:(DDSIP_param->cbitlim+1)/2;
    DDSIP_bb->last_dualitcnt = 0;
    diff = -1.;
    inherited_bound = DDSIP_node[DDSIP_bb->curnode]->bound;
    if (DDSIP_param->outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Invoking ConicBundle...\n");
    }
    memcpy (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual + 3));
    // New cplex parameters
    if (DDSIP_param->cpxnodual)
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxnodual, DDSIP_param->cpxdualisdbl, DDSIP_param->cpxdualwhich, DDSIP_param->cpxdualwhat);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to set CPLEX parameters (Dual) \n");
            DDSIP_Free ((void **) &(minfirst));
            DDSIP_Free ((void **) &(maxfirst));
            DDSIP_Free ((void **) &(center_point));
            return status;
        }
    }

    p = cb_construct_problem (0);
    if (p == 0)
    {
        fprintf (stderr, "ERROR: construct_problem failed\n");
        DDSIP_Free ((void **) &(minfirst));
        DDSIP_Free ((void **) &(maxfirst));
        DDSIP_Free ((void **) &(center_point));
        return 1;
    }
    DDSIP_bb->dualProblem = p;
    cb_clear (p);
    cb_set_defaults (p);
    if (cb_init_problem (p, DDSIP_bb->dimdual, NULL, NULL))
    {
        fprintf (stderr, "ERROR: init_problem failed\n");
        cb_destruct_problem (&p);
        DDSIP_Free ((void **) &(minfirst));
        DDSIP_Free ((void **) &(maxfirst));
        DDSIP_Free ((void **) &(center_point));
        return 1;
    }
    if (cb_add_function (p, (void *) DDSIP_DualUpdate, DDSIP_DualUpdate, 0, 0))
    {
        fprintf (stderr, "ERROR: add DUAL_UPDATE failed\n");
        cb_destruct_problem (&p);
        DDSIP_Free ((void **) &(minfirst));
        DDSIP_Free ((void **) &(maxfirst));
        DDSIP_Free ((void **) &(center_point));
        return 1;
    }
    if (cb_reinit_function_model(p, (void *) DDSIP_DualUpdate))
    {
        fprintf (stderr, "ERROR: reinit_function_model failed\n");
        cb_destruct_problem (&p);
        DDSIP_Free ((void **) &(minfirst));
        DDSIP_Free ((void **) &(maxfirst));
        DDSIP_Free ((void **) &(center_point));
        return 1;
    }
    DDSIP_bb->dualObjVal = -DDSIP_infty;
    old_obj = DDSIP_node[DDSIP_bb->curnode]->bound;

    cb_set_print_level (p, DDSIP_param->cbprint);
    cb_set_term_relprec (p, DDSIP_param->cbrelgap);
    cb_set_eval_limit (p, DDSIP_param->cbtotalitlim);
    cb_set_max_bundlesize (p, (void *) DDSIP_DualUpdate, DDSIP_param->cbbundlesz);
    cb_set_max_new_subgradients (p, (void *) DDSIP_DualUpdate, DDSIP_param->cbmaxsubg);
    cb_set_inner_update_limit 	(p, DDSIP_param->cb_maxsteps + 2);
    
    // The choice of the starting weight influences the convergence and results of Conic Bundle.
    // Is there a generally good choice?
    // More testing has to be done!
    if (DDSIP_bb->curnode)
    {
        start_weight = DDSIP_param->cbfactor*DDSIP_param->cbweight + (1. - DDSIP_param->cbfactor)*DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual];
        if (DDSIP_param->cb_increaseWeight || DDSIP_param->cb_reduceWeight)
        {
            // idea: limit too big start weights - speed up or slow down?
            //       cpu_secs is used for the limit
            cpu_secs = 5.0 * fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 10.;
            if (start_weight > cpu_secs && (!DDSIP_param->riskmod || start_weight > 1.e+3*DDSIP_param->riskweight))
            {
                start_weight = (start_weight-cpu_secs)*0.5  + cpu_secs;
//#############
if(DDSIP_param->outlev)
fprintf(DDSIP_bb->moreoutfile, "### new start_weight= %.8g\n", start_weight);
//#############
            }
            if (DDSIP_param->riskmod == 4)
            {
                if (start_weight < DDSIP_param->riskweight*0.2)
                    start_weight = DDSIP_param->riskweight*0.2;
            }
        }
    }
    else
        start_weight = DDSIP_param->cbweight;
    last_weight = next_weight = start_weight;
    if (start_weight >= 0.)
    {
        cb_set_next_weight (p, start_weight);
    }
    else
    {
        last_weight = next_weight = 12.3456789;
        cb_set_next_weight (p, next_weight);
    }
    DDSIP_bb->dualitcnt     = 0;
    DDSIP_bb->dualdescitcnt = 0;
    DDSIP_bb->cutAdded = 0;
    if (DDSIP_param->outlev)
        fprintf (DDSIP_bb->moreoutfile, "\nInitial dual evaluation\n");

    // paranoid dimension check
    if (cb_get_dim(p) != DDSIP_bb->dimdual)
    {
        fprintf (stderr, "XXX Error in dimension of the conic bundle problem: is %d, should be %d.\n",cb_get_dim(p),DDSIP_bb->dimdual);
        exit (1);
    }
    // save the inherited multipliers
    if (DDSIP_bb->curnode)
    {
        memcpy (DDSIP_bb->startinfo_multipliers, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual + 3));
#ifdef DEBUG
        if (DDSIP_param->outlev > 21 && DDSIP_param->outlev < DDSIP_current_lambda_outlev)
        {
            int i;
            fprintf (DDSIP_bb->moreoutfile, "\nCurrent lambda (inherited from node %.0g) for node %d: (%p)\n", DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual+1], DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
            for (i = 0; i < DDSIP_bb->dimdual; i++)
            {
                fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                if (!((i+1)%10))
                    fprintf (DDSIP_bb->moreoutfile, "\n");
            }
            fprintf(DDSIP_bb->moreoutfile, "\n");
        }
#endif
    }
    // Initialize multipliers (in root node: start multipliers if given, multipliers of father node otherwise)
    if (DDSIP_bb->curnode)
    {
        if (DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_bb->bestvalue - 0.5*(fabs(DDSIP_bb->bestvalue) + 1.e-10)*DDSIP_param->relgap)
        {
            for (cpu_hrs=0; cpu_hrs < DDSIP_param->scenarios; cpu_hrs++)
            {
                if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs])
                {
                    if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs])[DDSIP_bb->firstvar] - 0.9)))
                        for (j = cpu_hrs + 1; cnt && j < DDSIP_param->scenarios; j++)
                        {
                            if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]))
                            {
                                ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] = NULL;
                                cnt--;
                            }
                        }
                    DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs]));
                }
            }
            memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->startinfo_multipliers, sizeof (double) * (DDSIP_bb->dimdual+3));
            if ((status = cb_set_new_center_point (p, DDSIP_bb->startinfo_multipliers)))
            {
                fprintf (stderr, "set_new_center_point returned %d\n", status);
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                cb_destruct_problem (&p);
                DDSIP_Free ((void **) &(minfirst));
                DDSIP_Free ((void **) &(maxfirst));
                DDSIP_Free ((void **) &(center_point));
                return status;
            }
            obj = DDSIP_bb->dualObjVal;
            inhMult_bound = obj;
            if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_outfile, "\n   -- Dual:  Descent    Total  Objective        Weight ----------------------- Bound ------------------------------------------ Wall Time    CPU Time\n");
                printf ("\n   -- Dual:  Descent    Total  Objective        Weight ----------------------- Bound ------------------------------------------ Wall Time    CPU Time\n");
            }
            if (obj > old_obj)
            {
                memcpy (DDSIP_bb->local_bestdual, DDSIP_bb->startinfo_multipliers, sizeof (double) * (DDSIP_bb->dimdual + 3));
                old_obj = obj;
            }
        }
    }
    // in case initial_multiplier is -1, it was only partially set (e.g. the startinfo corresponds to a smaller number of scenarios).
    else // if (!DDSIP_bb->curnode)
    {
        DDSIP_bb->multipliers = 1;
        memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->startinfo_multipliers, sizeof (double) * (DDSIP_bb->dimdual));
        for (cpu_hrs=0; cpu_hrs < DDSIP_param->scenarios; cpu_hrs++)
        {
            if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs])
            {
                if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs])[DDSIP_bb->firstvar] - 0.9)))
                    for (j = cpu_hrs + 1; cnt && j < DDSIP_param->scenarios; j++)
                    {
                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]))
                        {
                            ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j] = NULL;
                            cnt--;
                        }
                    }
                DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[cpu_hrs]));
            }
        }
        if ((status = cb_set_new_center_point (p, DDSIP_bb->startinfo_multipliers)))
        {
            fprintf (stderr, "set_new_center_point returned %d\n", status);
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
            cb_destruct_problem (&p);
            DDSIP_Free ((void **) &(minfirst));
            DDSIP_Free ((void **) &(maxfirst));
            DDSIP_Free ((void **) &(center_point));
            return status;
        }
        obj = DDSIP_bb->dualObjVal;
        inhMult_bound = obj;
        if (obj > old_obj)
        {
            memcpy (DDSIP_bb->local_bestdual, DDSIP_bb->startinfo_multipliers, sizeof (double) * (DDSIP_bb->dimdual));
            memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->startinfo_multipliers, sizeof (double) * (DDSIP_bb->dimdual));
            old_obj = obj;
        }
        if (DDSIP_param->outlev && obj <= DDSIP_bb->bestvalue)
        {
            fprintf (DDSIP_outfile, "\n   -- Dual:  Descent    Total  Objective        Weight ----------------------- Bound ------------------------------------------ Wall Time    CPU Time\n");
            printf ("\n   -- Dual:  Descent    Total  Objective        Weight ----------------------- Bound ------------------------------------------ Wall Time    CPU Time\n");
        }
        if (!DDSIP_bb->initial_multiplier)
        {
            memcpy (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->subbound, DDSIP_param->scenarios*sizeof(double));
            DDSIP_node[DDSIP_bb->curnode]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[0] * DDSIP_data->prob[0];
            for (j = 1; j < DDSIP_param->scenarios; j++)
                DDSIP_node[DDSIP_bb->curnode]->BoundNoLag += DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j] * DDSIP_data->prob[j];
        //////////////////////////////
            if (DDSIP_param->outlev > 21)
            {
//                for (j = 0; j < DDSIP_param->scenarios; j++)
//                {
//                    fprintf (DDSIP_bb->moreoutfile, " ## bound zero mult. scen %3d:  %24.15g\n", j+1, DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j]);
//                }
                fprintf (DDSIP_bb->moreoutfile, " ## -----------------------\n ## node %d  BoundNoLag:         %24.15g\n ## -----------------------\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->BoundNoLag);
            }
        //////////////////////////////
        }
        if (DDSIP_param->outlev)
        {
            diff = obj - inherited_bound;
            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
            time (&DDSIP_bb->cur_time);
            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
            if (diff > 0.)
            {
                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                    time (&DDSIP_bb->cur_time);
                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                    if (DDSIP_bb->initial_multiplier)
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  startinfo mult.   %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight,
                             DDSIP_node[DDSIP_bb->curnode]->bound, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    else
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g       zero mult.   %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight,
                             DDSIP_node[DDSIP_bb->curnode]->bound, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
            }
            else
            {
                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                    time (&DDSIP_bb->cur_time);
                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                    if (DDSIP_bb->initial_multiplier)
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  startinfo mult.   %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight,
                             DDSIP_node[DDSIP_bb->curnode]->bound, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    else
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g       zero mult.   %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight,
                             DDSIP_node[DDSIP_bb->curnode]->bound, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
            }
        }
    }
    DDSIP_bb->multipliers = 1;
    obj = DDSIP_bb->dualObjVal;
    if (obj > old_obj)
        old_obj = obj;
    limits_reset = 0;
    if (obj < DDSIP_node[DDSIP_bb->curnode]->bound - 1.e-5)
    {
        if (DDSIP_param->cb_increaseWeight)
        {
            // start with bigger weight - hopefully this helps
            next_weight *= 10.0;
            cb_set_next_weight (p, next_weight);
        }
        if (DDSIP_param->cb_changetol)
        {
            // the branched problems produced a worse bound than in the father node - reset temporarily cplex tolerance and time limit
            limits_reset = 1;
            if (DDSIP_param->outlev)
            {
                printf("Initial evaluation at branched node gave obj=%g, worse than bound %g from father node - reset CPLEX relgap and/or time limit.\n", obj, DDSIP_node[DDSIP_bb->curnode]->bound);
                fprintf(DDSIP_bb->moreoutfile, "Initial evaluation at branched node gave obj=%g, worse than bound %g from father node - reset CPLEX relgap and/or time limit.\n", obj, DDSIP_node[DDSIP_bb->curnode]->bound);
            }
            for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual; cpu_hrs++)
            {
                if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_TILIM)
                {
                    old_cpxtimelim = DDSIP_param->cpxdualwhat[cpu_hrs];
                    if (DDSIP_param->outlev)
                        printf(" reset: old_timelim= %g, DDSIP_param->cpxdualwhich[%d]= %d, what= %g\n", old_cpxtimelim, cpu_hrs, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                    DDSIP_param->cpxdualwhat[cpu_hrs] = 1.5*old_cpxtimelim;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                        }
                    }
                }
                if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_EPGAP)
                {
                    old_cpxrelgap = DDSIP_param->cpxdualwhat[cpu_hrs];
                    if (DDSIP_param->outlev)
                        printf(" reset: old_relgap = %g, DDSIP_param->cpxdualwhich[%d]= %d, what= %g\n", old_cpxtimelim, cpu_hrs, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                    DDSIP_param->cpxdualwhat[cpu_hrs] = 0.10*old_cpxrelgap;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to reset CPLEX relgap to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                        }
                    }
                }
            }
            for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual2; cpu_hrs++)
            {
                if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_TILIM)
                {
                    old_cpxtimelim2 = DDSIP_param->cpxdualwhat2[cpu_hrs];
                    DDSIP_param->cpxdualwhat2[cpu_hrs] = 1.5*old_cpxtimelim2;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to reset CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("temporarily reset CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                        }
                    }
                }
                if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_EPGAP)
                {
                    old_cpxrelgap2 = DDSIP_param->cpxdualwhat2[cpu_hrs];
                    DDSIP_param->cpxdualwhat2[cpu_hrs] = 0.10*old_cpxrelgap2;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to reset CPLEX relgap 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("temporarily reset CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                        }
                    }
                }
            }
            // reevaluate initial point with changed parameters
            if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
            {
                fprintf (stderr, "set_new_center_point returned %d\n", status);
                cb_destruct_problem (&p);
                DDSIP_Free ((void **) &(minfirst));
                DDSIP_Free ((void **) &(maxfirst));
                DDSIP_Free ((void **) &(center_point));
                return status;
            }
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "set_new_center_point to bestdual, reevaluate with new tolerances and limits\n");
            obj = DDSIP_bb->dualObjVal;
        }
    }

    memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->local_bestdual, sizeof (double) * (DDSIP_bb->dimdual + 3));
    noIncreaseCounter = 0;
    if (DDSIP_bb->curnode <= DDSIP_param->cb_cutnodes)
    {
        if (DDSIP_param->outlev)
        {
            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
            time (&DDSIP_bb->cur_time);
            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
            if (DDSIP_bb->bestvalue < DDSIP_infty)
            {
                if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
                    rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / fabs (DDSIP_bb->bestvalue);
                else
                    rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / (fabs (DDSIP_bb->bestvalue) + DDSIP_param->accuracy);
                rgap = DDSIP_Dmin (rgap, 100.0);
                if (!DDSIP_bb->cutAdded)
                {
                    printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
                else
                {
                    printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
            }
            else
            {
                if (!DDSIP_bb->cutAdded)
                {
                    printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
                else
                {
                    printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %16dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %16dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
            }
        }
        DDSIP_bb->cutAdded = 0;
        if (!DDSIP_killsignal &&  DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_bb->bestvalue)
        {
            DDSIP_EvaluateScenarioSolutions (&comb);
            DDSIP_bb->keepSols = 1;
        }
        if (DDSIP_bb->bestvalue < 0.)
        {
            nfactor = 1.+DDSIP_Dmin(0.9*DDSIP_param->relgap, 1.e-8);
        }
        else
        {
            nfactor = 1.-DDSIP_Dmin(0.9*DDSIP_param->relgap, 1.e-8);
        }
        if (DDSIP_bb->cutAdded)
        {
            // reinit model
            if (cb_reinit_function_model(p, (void *) DDSIP_DualUpdate))
            {
                fprintf (stderr, "ERROR: reinit_function_model failed\n");
                cb_destruct_problem (&p);
                DDSIP_Free ((void **) &(minfirst));
                DDSIP_Free ((void **) &(maxfirst));
                DDSIP_Free ((void **) &(center_point));
                return 1;
            }
            else if (DDSIP_param->outlev > 20)
            {
                fprintf (DDSIP_bb->moreoutfile, "######### cb_reinit_function_model successful #########\n");
            }
            if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
            {
                fprintf (stderr, "set_new_center_point returned %d\n", status);
                cb_destruct_problem (&p);
                DDSIP_Free ((void **) &(minfirst));
                DDSIP_Free ((void **) &(maxfirst));
                DDSIP_Free ((void **) &(center_point));
                return status;
            }
            // use the result of set_center as initial dualObjVal
            if(DDSIP_param->outlev > 10)
                fprintf(DDSIP_bb->moreoutfile," after 1. reinit: currentDualObjVal = %20.14g, old_obj= %20.14g, incr.= %g\n", DDSIP_bb->currentDualObjVal, old_obj, DDSIP_bb->dualObjVal-old_obj);
	    if (inherited_bound <= -DDSIP_infty)
                inherited_bound = obj;

            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                time (&DDSIP_bb->cur_time);
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                if (!DDSIP_bb->curnode && DDSIP_bb->bestvalue > DDSIP_infty)
                {
                    if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
                        rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / fabs (DDSIP_bb->bestvalue);
                    else
                        rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / (fabs (DDSIP_bb->bestvalue) + DDSIP_param->accuracy);
                    rgap = DDSIP_Dmin (rgap, 100.0);
                    if (!DDSIP_bb->cutAdded)
                    {
                        printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                 0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    }
                    else
                    {
                        printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                 0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    }
                }
                else
                {
                    if (!DDSIP_bb->cutAdded)
                    {
                        printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                 0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    }
                    else
                    {
                        printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %16dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %16dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                 0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    }
                }
            }
            DDSIP_bb->cutAdded = 0;
            if (!DDSIP_killsignal)
            {
//#ifdef ONLYHEUR12
//                cpu_hrs = DDSIP_param->heuristic;
//                DDSIP_param->heuristic = 12;
//                DDSIP_EvaluateScenarioSolutions (&comb);
//                DDSIP_param->heuristic = cpu_hrs;
//                DDSIP_bb->keepSols = 1;
//#else
                DDSIP_EvaluateScenarioSolutions (&comb);
                DDSIP_bb->keepSols = 1;
//#endif
                if (DDSIP_param->outlev && DDSIP_bb->cutAdded)
                {
                    fprintf (DDSIP_outfile, "  |%17d%88d cuts\n", DDSIP_bb->dualdescitcnt, DDSIP_bb->cutAdded);
                }
                DDSIP_bb->cutAdded = 0;
            }
            if (DDSIP_bb->cutCntr)
            {
                cnt = 1;
                do
                {
                    DDSIP_bb->cutAdded = 0;
                    // reinit model
                    if (cb_reinit_function_model(p, (void *) DDSIP_DualUpdate))
                    {
                        fprintf (stderr, "ERROR: reinit_function_model failed\n");
                        cb_destruct_problem (&p);
                        DDSIP_Free ((void **) &(minfirst));
                        DDSIP_Free ((void **) &(maxfirst));
                        DDSIP_Free ((void **) &(center_point));
                        return 1;
                    }
                    else if (DDSIP_param->outlev > 20)
                    {
                        fprintf (DDSIP_bb->moreoutfile, "######### cb_reinit_function_model %d successful #########\n", cnt);
                    }
                    if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
                    {
                        fprintf (stderr, "set_new_center_point returned %d\n", status);
                        cb_destruct_problem (&p);
                        DDSIP_Free ((void **) &(minfirst));
                        DDSIP_Free ((void **) &(maxfirst));
                        DDSIP_Free ((void **) &(center_point));
                        return status;
                    }
                    obj = DDSIP_bb->currentDualObjVal;
                    rgap = (obj - old_obj)/(fabs(old_obj) + 2e-15);
                    if (rgap > 1.e-9)
                    {
                        noIncreaseCounter = 0;
//#ifdef DEBUG
                        if(DDSIP_param->outlev > 20)
                            fprintf(DDSIP_bb->moreoutfile," ## BETTER BOUND   currentDualObjVal= %20.14g, %20.14g=old_obj, rel. change: %g -> noIncreaseCounter= %d\n", DDSIP_bb->currentDualObjVal, old_obj, rgap, noIncreaseCounter);
//#endif
                    }
                    else if (rgap < -1.e-14)
                    {
                        noIncreaseCounter += 2;
//#ifdef DEBUG
                        if(DDSIP_param->outlev > 20)
                            fprintf(DDSIP_bb->moreoutfile," ## WORSE BOUND    currentDualObjVal= %20.14g, %20.14g=old_obj, rel. change: %g -> noIncreaseCounter= %d\n", DDSIP_bb->currentDualObjVal, old_obj, rgap, noIncreaseCounter);
//#endif
                    }
//#endif
                    else
                    {
                        noIncreaseCounter++;
//#ifdef DEBUG
                        if(DDSIP_param->outlev > 20)
                            fprintf(DDSIP_bb->moreoutfile," ## ELSE           currentDualObjVal= %20.14g, %20.14g=old_obj, rel. change: %g    noIncreaseCounter= %d\n", DDSIP_bb->currentDualObjVal, old_obj, rgap, noIncreaseCounter);
                    }
//#endif
                    cnt++;
                    if(DDSIP_param->outlev > 10)
                        fprintf(DDSIP_bb->moreoutfile," after %d. reinit: currentDualObjVal = %20.14g, old_obj= %20.14g, incr.= %g,  noIncreaseCounter= %d\n", cnt, DDSIP_bb->currentDualObjVal, old_obj, DDSIP_bb->currentDualObjVal-old_obj, noIncreaseCounter);
                    if (!DDSIP_bb->curnode && !DDSIP_bb->initial_multiplier && rgap > 0.)
                    {
                        memcpy (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->subbound, DDSIP_param->scenarios*sizeof(double));
                        DDSIP_node[DDSIP_bb->curnode]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[0] * DDSIP_data->prob[0];
                        for (j = 1; j < DDSIP_param->scenarios; j++)
                            DDSIP_node[DDSIP_bb->curnode]->BoundNoLag += DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j] * DDSIP_data->prob[j];
                    //////////////////////////////
                        if (DDSIP_param->outlev > 21)
                        {
                            fprintf (DDSIP_bb->moreoutfile, " ## -----------------------\n ## node %d  BoundNoLag:         %24.15g\n ## -----------------------\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->BoundNoLag);
                        }
                    //////////////////////////////
                    }
                    if (cnt > DDSIP_param->numberReinits)
                        break;
                    rgap = 100.;
                    if (DDSIP_param->outlev)
                    {
                        DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                        time (&DDSIP_bb->cur_time);
                        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                        if (!DDSIP_bb->curnode && DDSIP_bb->bestvalue < DDSIP_infty)
                        {
                            if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
                                rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / fabs (DDSIP_bb->bestvalue);
                            else
                                rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / (fabs (DDSIP_bb->bestvalue) + DDSIP_param->accuracy);
                            rgap = DDSIP_Dmin (rgap, 100.0);
                            if (!DDSIP_bb->cutAdded)
                            {
                                printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                        0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                         0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                            }
                            else
                            {
                                printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                        0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                         0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                            }
                        }
                        else
                        {
                            if (!DDSIP_bb->cutAdded)
                            {
                                printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                        0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                         0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                            }
                            else
                            {
                                printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %16dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                        0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %8d cuts %16dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                         0, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->cutAdded, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                            }
                        }
                    }
                    if (DDSIP_bb->currentDualObjVal > old_obj)
                    {
#ifdef DEBUG
                        if(DDSIP_param->outlev > 10)
                            fprintf(DDSIP_bb->moreoutfile," ##  update old_obj: currentDualObjVal= %20.14g, dualObjVal= %20.14g, old_obj= %20.14g\n", DDSIP_bb->currentDualObjVal, DDSIP_bb->dualObjVal, old_obj);
#endif
                        old_obj = DDSIP_bb->currentDualObjVal;
                        obj = DDSIP_bb->dualObjVal;
                    }
                    cpu_mins = DDSIP_bb->cutAdded;
      if (DDSIP_param->outlev > 21)
      {
          fprintf (DDSIP_bb->moreoutfile, "### vor  EvaluateScenarioSolutions %d: cutAdded= %d, cutCntr= %d\n", cnt, DDSIP_bb->cutAdded, DDSIP_bb->cutCntr);
      }
                    if (!DDSIP_killsignal)
                    {
#ifdef ONLYHEUR12
                        cpu_hrs = DDSIP_param->heuristic;
                        DDSIP_param->heuristic = 12;
                        DDSIP_EvaluateScenarioSolutions (&comb);
                        DDSIP_param->heuristic = cpu_hrs;
                        DDSIP_bb->keepSols = 1;
#else
                        DDSIP_EvaluateScenarioSolutions (&comb);
                        DDSIP_bb->keepSols = 1;
#endif
                    }
                    else
                        break;
      if (DDSIP_param->outlev > 21)
      {
          fprintf (DDSIP_bb->moreoutfile, "### nach EvaluateScenarioSolutions %d: cutAdded= %d, cutCntr= %d\n", cnt, DDSIP_bb->cutAdded, DDSIP_bb->cutCntr);
      }
                    if (DDSIP_param->outlev && (cpu_hrs = DDSIP_bb->cutAdded - cpu_mins))
                    {
                        fprintf (DDSIP_outfile, "  |%17d%88d cuts\n", DDSIP_bb->dualdescitcnt, cpu_hrs);
                    }
                    cpu_mins = DDSIP_bb->cutAdded;
                    DDSIP_bb->cutAdded = 0;
                } while (cpu_mins && (((obj - old_obj)/(fabs(obj)+1e-16) > 4.e-12) || (noIncreaseCounter < 6))
                         && cnt < DDSIP_param->numberReinits && rgap > 99.*DDSIP_param->relgap);
                if (DDSIP_param->outlev && DDSIP_bb->cutAdded)
                {
                    fprintf (DDSIP_outfile, "  |%17d%88d cuts\n", DDSIP_bb->dualdescitcnt, DDSIP_bb->cutAdded);
                }
                DDSIP_bb->cutAdded = 0;
                if (!DDSIP_bb->curnode && DDSIP_bb->cutCntr > 1)
                {
                    if (DDSIP_param->outlev > 10)
                    {
                       fprintf (DDSIP_bb->moreoutfile, " ########  total number of cuts added: %3d  #########################\n", DDSIP_bb->cutCntr);
                    }
                    fprintf (DDSIP_outfile, "  | %16d  --- %82d cuts total ---\n", 0, DDSIP_bb->cutCntr);
                    if (DDSIP_param->cb_increaseWeight && obj > inherited_bound + 1.e-3 && DDSIP_param->cbrootitlim > 5)
                    {
                        next_weight *= 5.0;
                        cb_set_next_weight (p, next_weight);
                    }
                }
                if (DDSIP_bb->currentDualObjVal > DDSIP_bb->dualObjVal)
                {
                    DDSIP_bb->dualObjVal = DDSIP_Dmax (old_obj, DDSIP_bb->currentDualObjVal);
                }
                old_obj = obj = DDSIP_bb->dualObjVal;
            }
            last_weight = next_weight;
        }
        init_iters = DDSIP_bb->dualitcnt;
    }
    if (DDSIP_bb->bestvalue < 0.)
    {
        nfactor = 1.+DDSIP_Dmin(0.9*DDSIP_param->relgap, 1.e-8);
    }
    else
    {
        nfactor = 1.-DDSIP_Dmin(0.9*DDSIP_param->relgap, 1.e-8);
    }
    noIncreaseCounter = 0;
    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 2] = 0;
    if (obj > inhMult_bound)
    {
        if (DDSIP_param->outlev > 20)
        {
           fprintf (DDSIP_bb->moreoutfile, " ######## update inhMult_bound  (=%.14g) to obj (=%.14g) #########################\n", inhMult_bound, obj);
        }
        inhMult_bound = obj;
    }
    if (DDSIP_bb->curnode)
    {
        diff = obj - inherited_bound;
        if (obj > inherited_bound)
        {
            if (DDSIP_param->outlev > 20)
            {
               fprintf (DDSIP_bb->moreoutfile, " ######## update inherited_bound (=%.14g) to obj (=%.14g) #########################\n", inherited_bound, obj);
            }
            inherited_bound = obj;
        }
        DDSIP_bb->keepSols = 0;
        if (DDSIP_param->cb_checkBestdual && DDSIP_bb->bestdual_cnt &&
            (DDSIP_bb->bestvalue == DDSIP_infty || diff < 0.97*(DDSIP_bb->bestvalue - DDSIP_bb->bestbound)) &&
            (DDSIP_bb->dualObjVal < DDSIP_bb->bestvalue - DDSIP_Dmax(1.e-11,0.5*DDSIP_param->relgap)*(fabs(DDSIP_bb->bestvalue)+1.e-10)))
        {
            bbest_t * tmp_bestdual, * tmp_previous, * tmp_maxbound = NULL, * tmp_minbound = NULL, * tmp_minprevious = NULL;
            double max_weight, max_bound, min_bound;
#if MIX > 0
            double * tmp_b1 = NULL, * tmp_b2 = NULL;
            double b_b1 = -DDSIP_infty, b_b2 = -DDSIP_infty;
#endif
            tmp_maxbound = tmp_bestdual = tmp_previous = DDSIP_bb->bestdual;
            max_bound = -DDSIP_infty;
            min_bound =  DDSIP_infty;
            max_weight = start_weight;
            cnt = 0;
            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                time (&DDSIP_bb->cur_time);
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                if (diff > 0.)
                {
                        DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                        time (&DDSIP_bb->cur_time);
                        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  mult. node %5g  %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                 DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight,
                                 DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1], DDSIP_bb->dualObjVal, diff,
                                 wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
                else
                {
                        DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                        time (&DDSIP_bb->cur_time);
                        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  mult. node %5g  %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                 DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight,
                                 DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1], DDSIP_bb->dualObjVal, diff,
                                 wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
            }
            old_obj = obj = DDSIP_bb->dualObjVal;
            if (obj > inhMult_bound)
            {
                if (DDSIP_param->outlev > 20)
                {
                   fprintf (DDSIP_bb->moreoutfile, " ######## update inhMult_bound  (=%.14g) to obj (=%.14g) #########################\n", inhMult_bound, obj);
                }
                inhMult_bound = obj;
            }
            // test Lagrange multipliers from DDSIP_bb->bestdual
            while (tmp_bestdual && obj < DDSIP_bb->bestvalue - DDSIP_Dmin (0.3*DDSIP_param->relgap, 2.5e-9)*(fabs(DDSIP_bb->bestvalue)+1.e-12))
            {
                if (tmp_bestdual->node_nr != (int) DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual + 1]
                    && (!DDSIP_MultEqual (DDSIP_bb->local_bestdual, tmp_bestdual->dual))
                    && (!DDSIP_MultEqual (DDSIP_bb->startinfo_multipliers, tmp_bestdual->dual))
                   )
                {
                    memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, tmp_bestdual->dual, sizeof (double) * DDSIP_bb->dimdual);
                    DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual+1] = tmp_bestdual->node_nr;
                    if (DDSIP_param->cb_increaseWeight || DDSIP_param->cb_reduceWeight)
                    {
                        next_weight = last_weight =
                            //DDSIP_Dmin(
                            //10.*fabs(DDSIP_bb->bestbound) + 1.,
                            DDSIP_param->cbfactor*DDSIP_param->cbweight + (1. - DDSIP_param->cbfactor)*
                            last_weight > tmp_bestdual->weight?
                            (0.7*last_weight + 0.3*tmp_bestdual->weight):
                            (0.3*last_weight + 0.7*tmp_bestdual->weight)
                            //  )
                           ;
                        max_weight = DDSIP_Dmax (max_weight, next_weight);
                        cb_set_next_weight (p, next_weight);
                    }
                    else
                    {
                        max_weight = next_weight = last_weight = start_weight;
                        cb_set_next_weight (p, next_weight);
                    }
//#ifdef DEBUG
//                    if (DDSIP_param->outlev > 21 && DDSIP_param->outlev < DDSIP_current_lambda_outlev)
//                    {
//                        int i;
//                        fprintf (DDSIP_bb->moreoutfile, "\nCurrent lambda (mult. node %d) for node %d: (%p)\n", tmp_bestdual->node_nr, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
//                        for (i = 0; i < DDSIP_bb->dimdual; i++)
//                        {
//                            fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
//                            if (!((i+1)%10))
//                                fprintf (DDSIP_bb->moreoutfile, "\n");
//                        }
//                        fprintf(DDSIP_bb->moreoutfile, "\n");
//                    }
//#endif

                    if ((status = cb_set_new_center_point (p, DDSIP_node[DDSIP_bb->curnode]->dual)))
                    {
                        fprintf (stderr, "set_new_center_point returned %d\n", status);
                        cb_destruct_problem (&p);
                        DDSIP_Free ((void **) &(minfirst));
                        DDSIP_Free ((void **) &(maxfirst));
                        DDSIP_Free ((void **) &(center_point));
                        return status;
                    }
                    cnt = 1;
                    obj = DDSIP_bb->currentDualObjVal;
                    if (obj > inhMult_bound && obj > max_bound)
                    {
                        memcpy (DDSIP_bb->local_bestdual, tmp_bestdual->dual, sizeof (double) * (DDSIP_bb->dimdual));
                        DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                        DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = tmp_maxbound->node_nr;
                    }
                    diff = obj - inhMult_bound;
                    tmp_bestdual->mean_diff = (tmp_bestdual->nr_uses * tmp_bestdual->mean_diff + diff)/(tmp_bestdual->nr_uses + 1.);
                    tmp_bestdual->nr_uses += 1.;
//*********************
if (DDSIP_param->outlev > 20)
{
    fprintf (DDSIP_bb->moreoutfile, "  Lagrangian from node %d used %g times, diff= %16.12g, mean_diff= %16.12g\n", tmp_bestdual->node_nr, tmp_bestdual->nr_uses, diff, tmp_bestdual->mean_diff);
}
//*********************
                    if (obj < min_bound)
                    {
                        min_bound = obj;
                        tmp_minbound = tmp_bestdual;
                    }
                    if (diff > 0.)
                    {
#if MIX > 0
                        if (!tmp_b1)
                        {
                            tmp_b1 = tmp_bestdual->dual;
                            b_b1   = obj;
                        }
                        else
                        {
                            if (obj < b_b1)
                            {
                                if (obj > b_b2)
                                {
                                    tmp_b2 = tmp_bestdual->dual;
                                    b_b2   = obj;
                                }
                            }
                            else
                            {
                                tmp_b2 = tmp_b1;
                                b_b2   = b_b1;
                                tmp_b1 = tmp_bestdual->dual;
                                b_b1   = obj;
                            }
                        }
#endif
                        if (DDSIP_param->outlev)
                        {
                            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                            time (&DDSIP_bb->cur_time);
                            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                            fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  mult. node %5d  %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight, tmp_bestdual->node_nr,
                                     obj, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                        }
                        if (obj > max_bound)
                        {
                            max_bound = obj;
                            tmp_maxbound = tmp_bestdual;
                            DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                            DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = tmp_maxbound->node_nr;
                            // put the max bound at the top of the list
                            if (tmp_bestdual != DDSIP_bb->bestdual)
                            {
                                tmp_previous->next = tmp_bestdual->next;
                                tmp_bestdual->next = DDSIP_bb->bestdual;
                                DDSIP_bb->bestdual = tmp_bestdual;
//######################
#ifdef DEBUG
if (DDSIP_param->outlev > 10)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
fprintf (DDSIP_bb->moreoutfile, " ## shifted mult. from node %d to top of bestdual list, #entries: %d\n", tmp1_bestdual->node_nr, DDSIP_bb->bestdual_cnt);
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
#endif
//######################
                                tmp_bestdual = tmp_previous;
                            }
                        }
                        if (max_bound > DDSIP_bb->bestvalue - fabs (DDSIP_bb->bestvalue) * DDSIP_param->accuracy ||
                            //diff > DDSIP_Dmax(1.e-1, 1.e+2*DDSIP_param->relgap) * (fabs(obj) + 1.e-12) ||
                            (DDSIP_bb->bestvalue < DDSIP_infty && (obj - DDSIP_bb->bestbound) > 0.975*(DDSIP_bb->bestvalue - DDSIP_bb->bestbound)))
                        {
                            // cutoff or sufficient increase reached
                            old_obj = obj;
                            if (obj >  DDSIP_bb->bestvalue - 0.5*(fabs(DDSIP_bb->bestvalue) + 1e-9)*DDSIP_param->relgap - 1e-9)
                            {
                                cnt = 0;
                            }
                            else if (DDSIP_Dmax (obj, max_bound) >=  DDSIP_bb->bestvalue - 10.*(fabs(DDSIP_bb->bestvalue) + 1e-10)*DDSIP_param->relgap)
                            {
                                cnt = 2;
                                cb_set_next_weight (p, DDSIP_Dmax (5.*next_weight, 10.));
                            }
                            else
                                cnt = 1;
                            break;
                        }
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                            time (&DDSIP_bb->cur_time);
                            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                            fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  mult. node %5d  %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight, tmp_bestdual->node_nr,
                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                        }
                        cnt = 1;
                        if (obj > max_bound)
                        {
                            max_bound = obj;
                            tmp_maxbound = tmp_bestdual;
                            // put this one at the top of the list
                            if (tmp_bestdual != DDSIP_bb->bestdual)
                            {
                                tmp_previous->next = tmp_bestdual->next;
                                tmp_bestdual->next = DDSIP_bb->bestdual;
                                DDSIP_bb->bestdual = tmp_bestdual;
//######################
#ifdef DEBUG
if (DDSIP_param->outlev > 10)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
fprintf (DDSIP_bb->moreoutfile, " ## shifted mult. from node %d to top of bestdual list, #entries: %d\n", tmp1_bestdual->node_nr, DDSIP_bb->bestdual_cnt);
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
#endif
//######################
                                tmp_bestdual =  tmp_previous;
                            }
                        }
                    }
                }
                else if (tmp_bestdual != DDSIP_bb->bestdual && DDSIP_MultEqual (DDSIP_bb->startinfo_multipliers, tmp_bestdual->dual) && tmp_minbound)
                {
                    if (tmp_bestdual->bound > tmp_minbound->bound)
                    {
                        // put this one at the top of the list
                        if (tmp_bestdual != DDSIP_bb->bestdual)
                        {
                            tmp_previous->next = tmp_bestdual->next;
                            tmp_bestdual->next = DDSIP_bb->bestdual;
                            DDSIP_bb->bestdual = tmp_bestdual;
//######################
//#ifdef DEBUG
if (DDSIP_param->outlev > 10)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
fprintf (DDSIP_bb->moreoutfile, " ## shifted mult. from node %d to top of bestdual list, #entries: %d\n", tmp1_bestdual->node_nr, DDSIP_bb->bestdual_cnt);
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
//#endif
//######################
                            tmp_bestdual =  tmp_previous;
                        }
                    }
                }
                tmp_previous = tmp_bestdual;
                tmp_bestdual = tmp_bestdual->next;
//#ifdef DEBUG
if (DDSIP_param->outlev > 20)
{
fprintf (DDSIP_bb->moreoutfile, " ### tmp_bestdual now: %p, tmp_minbound: %p", tmp_bestdual, tmp_minbound);
if (tmp_minbound)
    fprintf (DDSIP_bb->moreoutfile, ", ->next: %p\n", tmp_minbound->next);
else
    fprintf (DDSIP_bb->moreoutfile, "\n");
}
//#endif
            }
            // if all were tested put the multiplier which gave the minimal diff and then the minimal mean_diff to the end of the list
            if (!tmp_bestdual && tmp_minbound && tmp_minbound->next)
            {
                bbest_t * tmp_meanprevious, * tmp_minmean, * tmp_last;
                double h = DDSIP_infty;
                tmp_previous = tmp_meanprevious = NULL;
                tmp_minmean = tmp_bestdual = DDSIP_bb->bestdual;
                while (tmp_bestdual)
                {
                   if (tmp_bestdual == tmp_minbound)
                       tmp_minprevious = tmp_previous;
                   if (tmp_bestdual->mean_diff < h)
                   {
                       h = tmp_bestdual->mean_diff;
                       tmp_minmean = tmp_bestdual;
                       tmp_meanprevious = tmp_previous;
                   }
                   tmp_previous = tmp_bestdual;
                   tmp_bestdual = tmp_bestdual->next;
                }
                tmp_last = tmp_previous;
                if (tmp_minprevious)
                {
                    if (tmp_minbound->next)
                    {
                        tmp_minprevious->next = tmp_minbound->next;
                        tmp_previous->next = tmp_minbound;
                        tmp_minbound->next = NULL;
                        tmp_last = tmp_minbound;
                    }
                }
                else
                {
                    if (tmp_minbound->next)
                    {
                        DDSIP_bb->bestdual = tmp_minbound->next;
                        tmp_previous->next = tmp_minbound;
                        tmp_minbound->next = NULL;
                        tmp_last = tmp_minbound;
                    }
                }
//######################
//#ifdef DEBUG
if (DDSIP_param->outlev > 20)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
fprintf (DDSIP_bb->moreoutfile, " ## diff of last iter: shifted mult. from node %d to end of bestdual list, #entries: %d\n", tmp_minbound->node_nr, DDSIP_bb->bestdual_cnt);
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, mean_diff %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->mean_diff, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
if (DDSIP_param->outlev)
   fprintf (DDSIP_bb->moreoutfile, " ##   tmp_minmean= %p != tmp_minbound= %p : %d)\n",tmp_minmean, tmp_minbound, tmp_minmean != tmp_minbound);
//#endif
//######################
                if (tmp_minmean != tmp_minbound)
                {
                    h = DDSIP_infty;
                    tmp_previous = tmp_meanprevious = NULL;
                    tmp_minmean = tmp_bestdual = DDSIP_bb->bestdual;
                    while (tmp_bestdual)
                    {
                       if (tmp_bestdual->mean_diff < h)
                       {
                           h = tmp_bestdual->mean_diff;
                           tmp_minmean = tmp_bestdual;
                           tmp_meanprevious = tmp_previous;
                       }
                       tmp_previous = tmp_bestdual;
                       tmp_bestdual = tmp_bestdual->next;
                    }
                    if (tmp_meanprevious)
                    {
                        if (tmp_minmean->next)
                        {
                            tmp_meanprevious->next = tmp_minmean->next;
                            tmp_last->next = tmp_minmean;
                            tmp_minmean->next = NULL;
                        }
                    }
                    else
                    {
                        if (tmp_minmean->next)
                        {
                            DDSIP_bb->bestdual = tmp_minmean->next;
                            tmp_last->next = tmp_minmean;
                            tmp_minmean->next = NULL;
                        }
                    }
//######################
//#ifdef DEBUG
if (DDSIP_param->outlev > 20)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
fprintf (DDSIP_bb->moreoutfile, " ## mean diff: shifted mult. from node %d to end of bestdual list, #entries: %d\n", tmp_minmean->node_nr, DDSIP_bb->bestdual_cnt);
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, mean_diff %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->mean_diff, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
//#endif
//######################
                }
            }
            if (cnt && DDSIP_bb->skip != 2)
            {
                // go back to the multipliers which gave the best bound - either from bestdual or the inherited one
                if (DDSIP_param->cb_increaseWeight || DDSIP_param->cb_reduceWeight)
                {
                    last_weight = next_weight = 0.2*max_weight + 0.8*start_weight;
                    cb_set_next_weight (p, next_weight);
                }
                if (max_bound >= inhMult_bound)
                {
//////////////////////////
if(DDSIP_param->outlev > 20)
    fprintf(DDSIP_bb->moreoutfile, "### max_bound (%g) >= inhMult_bound (%g), tmp_maxbound= %p\n", max_bound, inhMult_bound, tmp_maxbound);
//////////////////////////
                    if (tmp_maxbound)
                    {
                        if (DDSIP_param->cb_test_line && max_bound < DDSIP_bb->bestvalue - DDSIP_Dmin (0.49*DDSIP_param->relgap, 2.5e-9)*(fabs(DDSIP_bb->bestvalue) + 1e-12))
                        {
                            // but first try a point on the line between inherited and best multipliers
                            if (cnt == 2)
                            {
                                for (status=0; status < DDSIP_bb->dimdual; status++)
                                    DDSIP_node[DDSIP_bb->curnode]->dual[status] = 0.00002 * DDSIP_bb->startinfo_multipliers[status]
                                                                                + 0.99998 * tmp_maxbound->dual[status];
                            }
                            else
                            {
#if MIX == 0
                                for (status=0; status < DDSIP_bb->dimdual; status++)
                                        DDSIP_node[DDSIP_bb->curnode]->dual[status] = 0.008 * DDSIP_bb->startinfo_multipliers[status]
                                                                                    + 0.992 * tmp_maxbound->dual[status];
#else
                                if (!tmp_b2)
                                {
                                    for (status=0; status < DDSIP_bb->dimdual; status++)
                                        DDSIP_node[DDSIP_bb->curnode]->dual[status] = 0.010 * DDSIP_bb->startinfo_multipliers[status]
                                                                                    + 0.990 * tmp_maxbound->dual[status];
                                }
                                else
                                {
#if MIX == 1
                                    double h1, h2;
#ifdef DEBUG
if (DDSIP_param->outlev)
    fprintf(DDSIP_bb->moreoutfile, "# MIX 1, factor %g\n", MIXFACTOR);
#endif
                                    h2 = MIXFACTOR*(b_b2 - inhMult_bound)/(b_b1 + b_b2 - 2.*inhMult_bound);
                                    h1 = 1. - h2;
                                    for (status=0; status < DDSIP_bb->dimdual; status++)
                                        DDSIP_node[DDSIP_bb->curnode]->dual[status] = 0.008 * DDSIP_bb->startinfo_multipliers[status]
                                                                                    + 0.992 * (h1*tmp_b1[status]+h2*tmp_b2[status]);
#elif MIX == 2
#ifdef DEBUG
if (DDSIP_param->outlev)
    fprintf(DDSIP_bb->moreoutfile, "# MIX 2\n");
#endif
                                    for (status=0; status < DDSIP_bb->dimdual; status++)
                                        DDSIP_node[DDSIP_bb->curnode]->dual[status] = 0.00002 * DDSIP_bb->startinfo_multipliers[status] +
                                                                                      0.99998 * (0.002 * tmp_b2[status]
                                                                                              + 0.998 * tmp_b1[status]);
#endif
                                }
#endif
                            }
#ifdef DEBUG
                            if (DDSIP_param->outlev > 21 && DDSIP_param->outlev < DDSIP_current_lambda_outlev)
                            {
                                int i;
                                fprintf (DDSIP_bb->moreoutfile, "\nCurrent lambda (line betw. %d) for node %d: (%p)\n", tmp_maxbound->node_nr, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
                                for (i = 0; i < DDSIP_bb->dimdual; i++)
                                {
                                    fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                                    if (!((i+1)%10))
                                        fprintf (DDSIP_bb->moreoutfile, "\n");
                                }
                                fprintf(DDSIP_bb->moreoutfile, "\n");
                            }
#endif
                            if ((status = cb_set_new_center_point (p, DDSIP_node[DDSIP_bb->curnode]->dual)))
                            {
                                fprintf (stderr, "set_new_center_point returned %d\n", status);
                                if (DDSIP_param->outlev)
                                    fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                                cb_destruct_problem (&p);
                                DDSIP_Free ((void **) &(minfirst));
                                DDSIP_Free ((void **) &(maxfirst));
                                DDSIP_Free ((void **) &(center_point));
                                return status;
                                }
                            if (DDSIP_param->outlev)
                            {
                                DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                                time (&DDSIP_bb->cur_time);
                                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                diff = DDSIP_bb->currentDualObjVal - max_bound;
                                if (diff > 0.)
                                {
                                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  + line between    %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                 DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                 DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                }
                                else
                                {
                                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  + line between    %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                 DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                 DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                }
                            }
                            if (DDSIP_bb->currentDualObjVal >= max_bound)
                            {
                                memcpy (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual));
                                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = DDSIP_bb->curnode;
                                cnt = 0;
                            }
                            else
                            {
                                // try a extrapolated point on the line between inherited and best multipliers
                                if (cnt == 2)
                                {
                                    for (status=0; status < DDSIP_bb->dimdual; status++)
                                        DDSIP_node[DDSIP_bb->curnode]->dual[status] = -0.001 * DDSIP_bb->startinfo_multipliers[status]
                                                                                     + 1.001 * tmp_maxbound->dual[status];
                                }
                                else
                                {
#if MIX == 0
                                    for (status=0; status < DDSIP_bb->dimdual; status++)
                                        DDSIP_node[DDSIP_bb->curnode]->dual[status] = -0.012 * DDSIP_bb->startinfo_multipliers[status]
                                                                                     + 1.012 * tmp_maxbound->dual[status];
#else
                                    if (!tmp_b2)
                                    {
                                        for (status=0; status < DDSIP_bb->dimdual; status++)
                                            DDSIP_node[DDSIP_bb->curnode]->dual[status] = -0.01 * DDSIP_bb->startinfo_multipliers[status]
                                                                                         + 1.01 * tmp_maxbound->dual[status];
                                    }
                                    else
                                    {
#if MIX == 1
                                        double h1, h2;
                                        h2 = MIXFACTOR*(b_b2 - inhMult_bound)/(b_b1 + b_b2 - 2.*inhMult_bound);
                                        h1 = 1. - h2;
                                        for (status=0; status < DDSIP_bb->dimdual; status++)
                                            DDSIP_node[DDSIP_bb->curnode]->dual[status] = -0.01 * DDSIP_bb->startinfo_multipliers[status]
                                                                                         + 1.01 * (h1*tmp_b1[status]+h2*tmp_b2[status]);
#elif MIX == 2
                                        for (status=0; status < DDSIP_bb->dimdual; status++)
                                           DDSIP_node[DDSIP_bb->curnode]->dual[status] = -0.00002 * DDSIP_bb->startinfo_multipliers[status]
                                                                                        + 1.00002 * (-0.002 * tmp_b2[status]
                                                                                                  + 1.002 * tmp_b1[status]);
#endif
                                    }
#endif
                                }
#ifdef DEBUG
                                if (DDSIP_param->outlev > 21 && DDSIP_param->outlev < DDSIP_current_lambda_outlev)
                                {
                                    int i;
                                    fprintf (DDSIP_bb->moreoutfile, "\nCurrent lambda (line extr. %d) for node %d: (%p)\n", tmp_maxbound->node_nr, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
                                    for (i = 0; i < DDSIP_bb->dimdual; i++)
                                    {
                                        fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                                        if (!((i+1)%10))
                                            fprintf (DDSIP_bb->moreoutfile, "\n");
                                    }
                                    fprintf(DDSIP_bb->moreoutfile, "\n");
                                }
#endif
                                if ((status = cb_set_new_center_point (p, DDSIP_node[DDSIP_bb->curnode]->dual)))
                                {
                                    fprintf (stderr, "set_new_center_point returned %d\n", status);
                                    if (DDSIP_param->outlev)
                                        fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                                    cb_destruct_problem (&p);
                                    DDSIP_Free ((void **) &(minfirst));
                                    DDSIP_Free ((void **) &(maxfirst));
                                    DDSIP_Free ((void **) &(center_point));
                                    return status;
                                }
                                if (DDSIP_param->outlev)
                                {
                                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                                    time (&DDSIP_bb->cur_time);
                                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                    diff = DDSIP_bb->currentDualObjVal - max_bound;
                                    if (diff > 0.)
                                    {
                                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  + line extrap.    %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                    }
                                    else
                                    {
                                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  + line extrap.    %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                    }
                                }
                                if (DDSIP_bb->currentDualObjVal >= max_bound)
                                {
                                    memcpy (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual));
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = DDSIP_bb->curnode;
                                    cnt = 0;
                                }
                                else
                                {
                                    memcpy (DDSIP_bb->local_bestdual, tmp_maxbound->dual, sizeof (double) * (DDSIP_bb->dimdual));
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = tmp_maxbound->node_nr;
                                    cnt = 1;
                                }
                            }
                        }
                        else 
                        {
                            if (memcmp(DDSIP_bb->local_bestdual, tmp_maxbound->dual, sizeof (double) * (DDSIP_bb->dimdual)))
                            {
                                memcpy (DDSIP_bb->local_bestdual, tmp_maxbound->dual, sizeof (double) * (DDSIP_bb->dimdual));
                                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                                DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = tmp_maxbound->node_nr;
                                cnt = 1;
                            }
                            else
                                cnt = 0;
                        }
                    }
                }
                else
                {
                    cnt = 1;
//////////////////////////
if(DDSIP_param->outlev > 20)
    fprintf(DDSIP_bb->moreoutfile, "### max_bound (%g) < inherited_bound (%g), tmp_maxbound= %p\n", max_bound, inherited_bound, tmp_maxbound);
//////////////////////////
                    if (DDSIP_param->cb_test_line)
                    {
                        if (tmp_maxbound)
                        {
                            {
                                // try a extrapolated point on the line between inherited and best multipliers
                                for (status=0; status < DDSIP_bb->dimdual; status++)
                                    DDSIP_node[DDSIP_bb->curnode]->dual[status] =  1.006 * DDSIP_bb->startinfo_multipliers[status]
                                                                                 - 0.006 * tmp_maxbound->dual[status];
                                    //DDSIP_node[DDSIP_bb->curnode]->dual[status] =  1.002 * DDSIP_bb->startinfo_multipliers[status]
                                    //                                             - 0.002 * tmp_maxbound->dual[status];
#ifdef DEBUG
                                    if (DDSIP_param->outlev > 21 && DDSIP_param->outlev < DDSIP_current_lambda_outlev)
                                    {
                                        int i;
                                        fprintf (DDSIP_bb->moreoutfile, "\nCurrent lambda (line extr. %d) for node %d: (%p)\n", tmp_maxbound->node_nr, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
                                        for (i = 0; i < DDSIP_bb->dimdual; i++)
                                        {
                                            fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                                            if (!((i+1)%10))
                                                fprintf (DDSIP_bb->moreoutfile, "\n");
                                        }
                                        fprintf(DDSIP_bb->moreoutfile, "\n");
                                    }
#endif
                                if ((status = cb_set_new_center_point (p, DDSIP_node[DDSIP_bb->curnode]->dual)))
                                {
                                    fprintf (stderr, "set_new_center_point returned %d\n", status);
                                    if (DDSIP_param->outlev)
                                        fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                                    cb_destruct_problem (&p);
                                    DDSIP_Free ((void **) &(minfirst));
                                    DDSIP_Free ((void **) &(maxfirst));
                                    DDSIP_Free ((void **) &(center_point));
                                    return status;
                                }
                                if (DDSIP_param->outlev)
                                {
                                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                                    time (&DDSIP_bb->cur_time);
                                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                    diff = DDSIP_bb->currentDualObjVal - inhMult_bound;
                                    if (diff > 0.)
                                    {
                                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  - line extrap.    %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                    }
                                    else
                                    {
                                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  - line extrap.    %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                    }
                                }
                                if (DDSIP_bb->currentDualObjVal >= inhMult_bound)
                                {
                                    memcpy (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual));
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = DDSIP_bb->curnode;
                                    cnt = 0;
                                }
                                else
                                {
                                    cnt = 1;
                                }
                            }
                            if (cnt)
                            {
                                // try a point on the line between inherited and best multipliers
                                for (status=0; status < DDSIP_bb->dimdual; status++)
                                    DDSIP_node[DDSIP_bb->curnode]->dual[status] = 0.992 * DDSIP_bb->startinfo_multipliers[status]
                                                                                + 0.008 * tmp_maxbound->dual[status];
#ifdef DEBUG
                                    if (DDSIP_param->outlev > 21 && DDSIP_param->outlev < DDSIP_current_lambda_outlev)
                                    {
                                        int i;
                                        fprintf (DDSIP_bb->moreoutfile, "\nCurrent lambda (line betw. %d) for node %d: (%p)\n", tmp_maxbound->node_nr, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual);
                                        for (i = 0; i < DDSIP_bb->dimdual; i++)
                                        {
                                            fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                                            if (!((i+1)%10))
                                                fprintf (DDSIP_bb->moreoutfile, "\n");
                                        }
                                        fprintf(DDSIP_bb->moreoutfile, "\n");
                                    }
#endif
                                if ((status = cb_set_new_center_point (p, DDSIP_node[DDSIP_bb->curnode]->dual)))
                                {
                                    fprintf (stderr, "set_new_center_point returned %d\n", status);
                                    if (DDSIP_param->outlev)
                                        fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                                    cb_destruct_problem (&p);
                                    DDSIP_Free ((void **) &(minfirst));
                                    DDSIP_Free ((void **) &(maxfirst));
                                    DDSIP_Free ((void **) &(center_point));
                                    return status;
                                }
                                if (DDSIP_param->outlev)
                                {
                                    DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                                    time (&DDSIP_bb->cur_time);
                                    DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                                    diff = DDSIP_bb->currentDualObjVal - inhMult_bound;
                                    if (diff > 0.)
                                    {
                                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  - line between    %-20.14g incr.  %-12.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                    }
                                    else
                                    {
                                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  - line between    %-20.14g diff. %-13.7g %10dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                                     DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                                     DDSIP_bb->dualObjVal, diff, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                                    }
                                }
                                if (DDSIP_bb->currentDualObjVal >= inhMult_bound)
                                {
                                    memcpy (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual));
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = last_weight;
                                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1] = DDSIP_bb->curnode;
                                    cnt = 0;
                                }
                                else
                                {
                                    cnt = 1;
                                }
                            }
                        }
                    }
                }
                if (cnt)
                {
                    // increase weight
                    if (DDSIP_param->cb_increaseWeight && last_weight < 1.e3 * DDSIP_node[DDSIP_bb->curnode]->bound)
                    {
                        last_weight = next_weight *= 1.2;
                        cb_set_next_weight (p, next_weight);
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile,"########### points on line worse -> increased next weight to %g\n", next_weight);
                    }
                    if (memcmp (DDSIP_bb->local_bestdual, DDSIP_node[DDSIP_bb->curnode]->dual, sizeof (double) * (DDSIP_bb->dimdual)))
                    {
                        if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
                        {
                            fprintf (stderr, "set_new_center_point returned %d\n", status);
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                            cb_destruct_problem (&p);
                            DDSIP_Free ((void **) &(minfirst));
                            DDSIP_Free ((void **) &(maxfirst));
                            DDSIP_Free ((void **) &(center_point));
                            return status;
                        }
                    }
                    old_obj = obj = DDSIP_bb->dualObjVal;
                    if (obj > inhMult_bound)
                    {
                        if (DDSIP_param->outlev > 20)
                        {
                           fprintf (DDSIP_bb->moreoutfile, " ######## update inhMult_bound  (=%.14g) to obj (=%.14g) #########################\n", inhMult_bound, obj);
                        }
                        inhMult_bound = obj;
                    }
                    if (DDSIP_param->outlev)
                    {
                        DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                        time (&DDSIP_bb->cur_time);
                        DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                        fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  mult. node %5g  %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal, last_weight,
                                 DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 1], obj,
                                 wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    }
                }
            }
        }
        
        old_obj = obj = DDSIP_bb->dualObjVal;
        DDSIP_bb->last_dualitcnt = DDSIP_bb->dualitcnt;
        if (DDSIP_param->cb_increaseWeight &&
            ((DDSIP_node[DDSIP_bb->curnode]->depth <= DDSIP_param->cb_depth ||
             ((DDSIP_bb->bestvalue - obj) < DDSIP_Dmin (1.e4*DDSIP_param->relgap, 1.e-5)*(fabs(DDSIP_bb->bestvalue)+1.e-10))) &&
              next_weight < 10.))
        {
            last_weight = next_weight *= 2.5;
            cb_set_next_weight (p, next_weight);
        }
        if (DDSIP_param->outlev)
        {
            if (DDSIP_bb->cutAdded)
            {
                printf ("  | %16d  %86d cuts\n", 0, DDSIP_bb->cutAdded);
                fprintf (DDSIP_outfile, "  | %16d  %86d cuts\n", 0, DDSIP_bb->cutAdded);
            }
            DDSIP_bb->cutAdded = 0;
            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
            time (&DDSIP_bb->cur_time);
            DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
            if (!DDSIP_bb->curnode && DDSIP_bb->bestvalue < DDSIP_infty)
            {
                if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
                    rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / fabs (DDSIP_bb->bestvalue);
                else
                    rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / (fabs (DDSIP_bb->bestvalue) + DDSIP_param->accuracy);
                rgap = DDSIP_Dmin (rgap, 100.0);
                printf ("  | %16d  %7d  %-16.12g %-11.6g  ----------------  %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                        DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  ----------------  %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                         DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
            }
            else
            {
                printf ("  | %16d  %7d  %-16.12g %-11.6g  ----------------  %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                        DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g  ----------------  %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                         DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
            }
        }
    }
    DDSIP_bb->last_dualitcnt = DDSIP_bb->dualitcnt;

    if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - DDSIP_Dmin (0.5*DDSIP_param->relgap, 8e-9)*(fabs(DDSIP_bb->bestvalue) + 1e-14))
    {
        memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->local_bestdual, sizeof (double) * (DDSIP_bb->dimdual + 3));
        if (DDSIP_node[DDSIP_bb->curnode]->bound >= DDSIP_bb->bestvalue)
        {
            DDSIP_bb->cutoff++;
            if (DDSIP_bb->dualitcnt == 1)
               DDSIP_bb->dualdescitcnt = -1;
        }
        //
        if ((DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue)/(fabs(DDSIP_bb->bestvalue)+1.e-10) > DDSIP_param->accuracy)
        {
            if (DDSIP_param->outlev)
                DDSIP_Print2 ("   --------- termination status: cutoff. --------------------------------------------------------------", "\n", 0, 0);
            DDSIP_bb->skip = 3;
        }
        else if (!DDSIP_bb->violations &&
                 DDSIP_Equal (DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue))
        {
            if (DDSIP_param->outlev)
                DDSIP_Print2 ("   --------- termination status: optimal. --------------------------------------------------------------", "\n", 0, 0);
            DDSIP_bb->skip = 3;
            DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
        }
        else
        {
            if (DDSIP_param->outlev)
                DDSIP_Print2 ("   --------- termination status: within relative gap. Number of violations of nonanticipativity: ", " ------------------------------------------\n", 1.*DDSIP_node[DDSIP_bb->curnode]->violations, 1);
            DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
            if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - 1.e-11*fabs(DDSIP_bb->bestvalue))
                DDSIP_bb->skip = -2;
            if ((DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue*nfactor) &&
                (!(DDSIP_bb->found_optimal_node) || (DDSIP_bb->found_optimal_node && DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bound_optimal_node)))
            {
    //////////////////////////////////////////////////////////////////
                   if (DDSIP_param->outlev > 20)
                       fprintf (DDSIP_bb->moreoutfile, "########## setting found_optimal_node: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bestvalue (%20.15g) = %.8g, - bestvalue*nfactor (%20.15g) = %.8g\n",
                                DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,
                                DDSIP_bb->bestvalue*nfactor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*nfactor);
    //////////////////////////////////////////////////////////////////
                DDSIP_bb->found_optimal_node = DDSIP_bb->curnode;
                DDSIP_bb->bound_optimal_node = DDSIP_node[DDSIP_bb->curnode]->bound;
            }
    //////////////////////////////////////////////////////////////////
            else
            {
               if (DDSIP_param->outlev > 20)
                   fprintf (DDSIP_bb->moreoutfile, "########## else: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bound_optimal_node (%20.15g) = %.8g, - bestvalue*nfactor (%20.15g) = %.8g\n",
                            DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bound_optimal_node, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bound_optimal_node,
                            DDSIP_bb->bestvalue*nfactor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*nfactor);
            }
    //////////////////////////////////////////////////////////////////
        }
    }
    else if ((DDSIP_bb->no_reduced_front == 1) && (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - fabs(DDSIP_bb->bestvalue)*DDSIP_param->relgap))
    {
        if (DDSIP_param->outlev)
        {
            DDSIP_Print2 ("   --------- termination status: gap reached. --------------------------------------------------------------", "\n", 0, 0);
        }
        DDSIP_bb->dualitcnt    = 0;
    }
    else
    {
        if (!(DDSIP_bb->curnode) && DDSIP_param->deleteRedundantCuts)
            DDSIP_CheckRedundancy(1);
        noIncreaseCounter = 0;
        many_iters = 0;
////////////////////////////////////////////////////////////////////////////////////
////////// check the cuts - only for that reason the same center is used
//        {
//            int outfiles;
//            outfiles = DDSIP_param->files;
//            DDSIP_bb->keepSols = 0;
//            DDSIP_param->files = 4;
//            if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
//            {
//                fprintf (stderr, "set_new_center_point returned %d\n", status);
//                cb_destruct_problem (&p);
//                DDSIP_Free ((void **) &(minfirst));
//                DDSIP_Free ((void **) &(maxfirst));
//                DDSIP_Free ((void **) &(center_point));
//                return status;
//            }
//            DDSIP_param->files = outfiles;
//        }
////////////////////////////////////////////////////////////////////////////////////
	DDSIP_bb->last_weight = last_weight = next_weight;
        DDSIP_bb->last_dualitcnt = DDSIP_bb->dualitcnt;
        while ((!cb_termination_code (p)) && DDSIP_bb->violations && DDSIP_bb->skip != 2
                && (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time) < DDSIP_param->timelim)
                && DDSIP_bb->dualdescitcnt < DDSIP_bb->current_itlim
                && DDSIP_bb->dualitcnt < DDSIP_param->cbtotalitlim && !(obj > DDSIP_bb->bestvalue - DDSIP_param->accuracy)
                && (DDSIP_node[DDSIP_bb->curnode]->bound < DDSIP_bb->bestvalue - (fabs(DDSIP_bb->bestvalue) + 1.e-12) * DDSIP_Dmax (DDSIP_Dmin (0.3*DDSIP_param->relgap, 4e-9), 2.e-12))
                && cycleCnt < 2)
        {
            if ((DDSIP_bb->no_reduced_front == 1) &&
                    (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - 4.e-1*fabs(DDSIP_bb->bestvalue) * DDSIP_param->relgap))
            {
                break;
            }

            if (DDSIP_killsignal)
            {
                break;
            }
            if (DDSIP_param->outlev && DDSIP_bb->cutAdded)
            {
                fprintf (DDSIP_outfile, "  |%17d%88d cuts\n", DDSIP_bb->dualdescitcnt, DDSIP_bb->cutAdded);
            }
            DDSIP_bb->dualdescitcnt++;
            if (DDSIP_bb->dualdescitcnt <= 1)
                DDSIP_bb->dualObjVal = -DDSIP_infty;
            DDSIP_bb->cutAdded = 0;
            if (DDSIP_param->outlev > 8)
            {
                printf ("\nDescent step %d    next weight %g,  last weight %g\n", DDSIP_bb->dualdescitcnt, next_weight, last_weight);
                fprintf (DDSIP_bb->moreoutfile, "\nDescent step %d    next weight %g\n", DDSIP_bb->dualdescitcnt, next_weight);
            }
            DDSIP_bb->keepSols = 0;

#ifdef DEBUG
            if (DDSIP_param->outlev > 99)
            {
                // check the changes of the center point
                if ((status = cb_get_center (p, center_point)))
                {
                    fprintf (stderr, "get_center returned %d\n", status);
                }
                else
                {
                    fprintf (DDSIP_bb->moreoutfile, " --------------- current center point ---------------\n");
                    for (status=0; status < DDSIP_bb->dimdual; status++)
                    {
                        fprintf (DDSIP_bb->moreoutfile, " %14.8g,", center_point[status]);
                        if (!((status+1)%10))
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    fprintf (DDSIP_bb->moreoutfile, "\n ------------------------------\n");
                }
            }
#endif
            DDSIP_bb->dualObjVal = -DDSIP_infty;

            /* Make a descent step */
            /* status = cb_do_descent_step (p); */
            // /* DLW Dec 2014 - limit the number of null steps in order to avoid loops of reevaluations */
NEXT_TRY:
            current_maxsteps = DDSIP_param->cb_maxsteps + (DDSIP_bb->dualdescitcnt==1?5:0) + ((DDSIP_bb->curnode < 5)?4:0);
            cb_set_inner_update_limit (p, current_maxsteps + 2);
            cb_status = cb_do_maxsteps(p, current_maxsteps); /* DLW Dec 2014 */
            // update dual solution
            cb_get_center (p,DDSIP_node[DDSIP_bb->curnode]->dual);
            if (DDSIP_bb->dualdescitcnt == 1)
                last_weight = DDSIP_bb->last_weight;
            DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual] = cb_get_last_weight(p);
            if (DDSIP_killsignal )
            {
                // Killsignal
                break;
            }
            else if (cb_status)
            {
                if (!DDSIP_bb->curnode)
                    fprintf (DDSIP_outfile, "############### cb_do_maxsteps returned %d, newTry= %d\n", cb_status, DDSIP_bb->newTry);
                if (!DDSIP_bb->newTry || DDSIP_bb->newTry > 5)
                {
                    fprintf (stderr, "############### cb_do_maxsteps returned %d\n", cb_status);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "############### cb_do_maxsteps returned %d\n", cb_status);
                    cb_status = cb_termination_code(p);
                    fprintf (stderr, "############### cb_termination_code is %d\n", cb_status);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "############### cb_termination_code is %d\n", cb_status);
                    if (DDSIP_bb->curnode || !DDSIP_bb->dualdescitcnt)
                    {
                        cb_destruct_problem (&p);
                        // Reset first stage solutions to the ones that gave the best bound
                        for (j = 0; j < DDSIP_param->scenarios; j++)
                        {
                            if (DDSIP_bb->bestfirst[j].first_sol)
                            {
                                if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])
                                {
                                    if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])[DDSIP_bb->firstvar] - 0.9)))
                                        for (i_scen = j + 1; cnt && i_scen < DDSIP_param->scenarios; i_scen++)
                                        {
                                            if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i_scen]))
                                            {
                                                ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i_scen] = NULL;
                                                cnt--;
                                            }
                                        }
                                    DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]));
                                }
                                DDSIP_node[DDSIP_bb->curnode]->first_sol[j]   = (DDSIP_bb->bestfirst[j]).first_sol;
                                (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] = (DDSIP_bb->bestfirst[j]).cursubsol;
                                (DDSIP_node[DDSIP_bb->curnode]->subbound)[j]  = (DDSIP_bb->bestfirst[j]).subbound;
                                DDSIP_bb->bestfirst[j].first_sol = NULL;
                            }
                        }
                        DDSIP_Free ((void **) &(minfirst));
                        DDSIP_Free ((void **) &(maxfirst));
                        DDSIP_Free ((void **) &(center_point));
                        return 101;
                    }
                    else
                    {
                        if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
                        {
                            fprintf (stderr, "set_new_center_point returned %d\n", status);
                            cb_destruct_problem (&p);
                            DDSIP_Free ((void **) &(minfirst));
                            DDSIP_Free ((void **) &(maxfirst));
                            DDSIP_Free ((void **) &(center_point));
                            return status;
                        }
                        if (DDSIP_param->outlev)
                            fprintf (DDSIP_bb->moreoutfile, "set_new_center_point to bestdual, reevaluate with increased weight\n");
                        // increase weight
                        if (last_weight < 1.e20)
                        {
                            last_weight = next_weight = DDSIP_Dmax(5e-2 * fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 100., 5e+3 * cb_get_last_weight (p));
                            cb_set_next_weight (p, next_weight);
                            if (DDSIP_param->outlev)
                                fprintf (DDSIP_bb->moreoutfile,"########### increased next weight to %g, try again\n", next_weight);
                        }
                        else
                        {
                            last_weight = next_weight = 1e-2;
                            cb_set_next_weight (p, next_weight);
                        }
                        goto NEXT_TRY;
                    }
                }
            }

            if (DDSIP_bb->newTry)
            {
                // increase weight
                if (last_weight < 1.e30)
                {
                    last_weight = next_weight = DDSIP_Dmax(5e-2 * fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 100., 5e+3 * cb_get_last_weight (p));
                    cb_set_next_weight (p, next_weight);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile,"########### increased next weight to %g, try again\n", next_weight);
                    goto NEXT_TRY;
                }
            }
            else
            {
                if ((cb_status = cb_termination_code(p)))
                {
                    fprintf (stderr, "cb_termination_code: %d\n", cb_status);
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "cb_termination_code: %d\n", cb_status);
                }
                /* Get solution information */
                // obj = -cb_get_objval (p);
                // take the best obj value in case the descent step was interrupted by maxsteps or other cause
                obj = DDSIP_bb->dualObjVal;
                next_weight = cb_get_last_weight (p);
                j = (obj <= old_obj);
                // if the step did not increase the bound, increase the weight
                if (((next_weight - last_weight  <= 0.5*last_weight) || (DDSIP_bb->dualdescitcnt == 1 && DDSIP_bb->weight_reset == 1)) &&
                    (!DDSIP_param->cb_inherit || j ||
                    (!j &&
                     ((DDSIP_bb->dualdescitcnt == 1 &&
                       DDSIP_bb->dualitcnt < DDSIP_Imin(DDSIP_param->cb_maxsteps, 11) + init_iters) ||
                     DDSIP_bb->dualdescitcnt > 1))))
                {
                    DDSIP_bb->weight_reset = 0;
                    if (j || obj <= inherited_bound)
                    {
                        if (repeated_increase > -2)
                            repeated_increase = -2;
                        noIncreaseCounter++;
                        many_iters +=2;
/////////
                        if (DDSIP_param->outlev > 10)
                        {
                            fprintf(DDSIP_bb->moreoutfile,"§§§§§§§§§§§§§§§ iters in descent step: %d, up to now: repeated increase= %d, many_iters= %d, weight change by CB: %g, ret-code cb_do_maxsteps: %d §§§§§§§§§§§§§§§§§§\n", DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt, repeated_increase, many_iters, next_weight - last_weight, cb_status);
                        }
/////////
                        if ((DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt >= current_maxsteps*0.6) && DDSIP_param->cb_increaseWeight)
                        {
                            if (noIncreaseCounter == 1)
                            {
                                if (next_weight < 1.1*last_weight)
                                {
                                    last_weight = next_weight;
                                    if (abs(DDSIP_param->riskmod) == 4)
                                        next_weight = 4.0*last_weight;
                                    else if (abs(DDSIP_param->riskmod) == 5)
                                        next_weight = 1.8*last_weight;
                                    else
                                    {
                                        if (many_iters < 4 && last_weight > 1.e-6 * fabs(DDSIP_node[DDSIP_bb->curnode]->bound))
                                            next_weight = 1.3*last_weight;
                                        else
                                            next_weight = 4.0*last_weight;
                                    }
                                }
                                if (DDSIP_bb->dualdescitcnt == 1)
                                {
                                    if (obj + 2.*(fabs(DDSIP_node[DDSIP_bb->curnode]->bound) + 10.) < 0.)
                                        next_weight = 1. + 100.*next_weight;
                                    else
                                        next_weight = 1. + 10.*next_weight;
                                }
                                cb_set_next_weight (p, next_weight);
                                if (DDSIP_param->outlev)
                                    fprintf (DDSIP_bb->moreoutfile,"####### §§§ increased next weight to %g\n", next_weight);
                            }
                            else if (noIncreaseCounter > 1)
                            {
                                repeated_increase--;
                                last_weight = next_weight;
                                if (noIncreaseCounter > 3)
                                {
                                    if (last_weight > 20000.)
                                    {
                                        next_weight = 0.002*last_weight;
                                    }
                                    else
                                    {
                                        next_weight = 100.*last_weight;
                                    }
                                    noIncreaseCounter = 1;
                                    cycleCnt++;
                                }
                                else if (last_weight < 50.)
                                {
                                    if (DDSIP_bb->dualdescitcnt < 3 && DDSIP_bb->dualitcnt > DDSIP_bb->dualdescitcnt*DDSIP_param->cb_maxsteps + 3)
                                    {
                                        if (last_weight < 2.)
                                        {
                                            next_weight = 1. + 100.*last_weight;
                                            if(DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile, "  next_weight=  1. + %g * %g = %g\n", 100., last_weight, next_weight);
                                        }
                                        else if (last_weight < 2. * fabs (DDSIP_bb->bestbound))
                                        {
                                            next_weight = 50.*last_weight;
                                            if(DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile, "  next_weight= %g * %g = %g\n", 50., last_weight, next_weight);
                                        }
                                        else
                                        {
                                            next_weight = 10.*last_weight;
                                            if(DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile, "  next_weight= %g * %g = %g\n", 10., last_weight, next_weight);
                                        }
                                    }
                                    else
                                    {
                                        cpu_secs=(old_obj - obj)/(fabs(old_obj)+1e-6);
                                        if(DDSIP_param->outlev > 10)
                                            fprintf(DDSIP_bb->moreoutfile, " relative decrease = %-16.12g",cpu_secs);
                                        if (cpu_secs < 6e-4)
                                        {
                                            next_weight = (1.3+1e3*cpu_secs)*last_weight;
                                            if(DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile, "  next_weight= %g * %g = %g\n", 1.3+1e3*cpu_secs, last_weight, next_weight);
                                        }
                                        else
                                        {
                                            next_weight = 4.0*last_weight;
                                            if(DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile, "  next_weight= %g * %g = %g\n", 4., last_weight, next_weight);
                                        }
                                    }
                                }
                                else
                                {
                                    if (DDSIP_bb->dualdescitcnt == 2 && DDSIP_bb->dualitcnt > DDSIP_bb->dualdescitcnt*DDSIP_param->cb_maxsteps + 2 + init_iters)
                                    {
                                        next_weight = 10.*last_weight;
                                        if(DDSIP_param->outlev > 10)
                                            fprintf(DDSIP_bb->moreoutfile, "  next_weight= %g * %g = %g\n", 10., last_weight, next_weight);
                                        if (weight_decreases)
                                        {
                                            reduction_factor = 0.8*reduction_factor + 0.18;
///////////     ///////////
                                            if (DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile,"############### repeated increase= %d, reduction factor increased to %g ##################\n",repeated_increase, reduction_factor);
///////////     ///////////
                                        }
                                    }
                                    else if (last_weight < 1.e4)
                                    {
                                        if (many_iters < 4)
                                            next_weight = 4.0*last_weight;
                                        else
                                            next_weight = 8.0*last_weight;
                                    }
                                    else
                                        next_weight = 1.2*last_weight;
                                }
                                cb_set_next_weight (p, next_weight);
                                /////////
                                if (DDSIP_param->outlev)
                                {
                                    fprintf (DDSIP_bb->moreoutfile,"####### §§§ increased next weight to %g\n", next_weight);
                                }
                                /////////
                                if (DDSIP_bb->dualdescitcnt == 1)
                                    old_obj = obj;
                            }
                            if (weight_decreases)
                            {
                                reduction_factor = 0.8*reduction_factor + 0.18;
///////////     ///////////
                                if (DDSIP_param->outlev > 10)
                                    fprintf(DDSIP_bb->moreoutfile,"############### repeated increase= %d, reduction factor increased to %g ##################\n",repeated_increase, reduction_factor);
///////////     ///////////
                            }
                        }
                        if (DDSIP_param->cb_changetol && !limits_reset && (obj < DDSIP_node[DDSIP_bb->curnode]->bound - 1.e-4*fabs(DDSIP_node[DDSIP_bb->curnode]->bound)))
                        {
                            // the branched problems produced a bound worse than already known - reset temporarily cplex tolerance and time limit
                            limits_reset = 2;
                            if (DDSIP_param->outlev)
                            {
                                printf(" reset CPLEX relgap and/or time limit.\n");
                                fprintf(DDSIP_bb->moreoutfile, " reset CPLEX relgap and/or time limit.\n");
                            }
                            for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual; cpu_hrs++)
                            {
                                if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_TILIM)
                                {
                                    old_cpxtimelim = DDSIP_param->cpxdualwhat[cpu_hrs];
                                    DDSIP_param->cpxdualwhat[cpu_hrs] = 1.2*old_cpxtimelim;
                                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                                    if (status)
                                    {
                                        printf("Failed to reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                    }
                                    else
                                    {
                                        if (DDSIP_param->outlev)
                                        {
                                            printf("temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                        }
                                    }
                                }
                                if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_EPGAP)
                                {
                                    old_cpxrelgap = DDSIP_param->cpxdualwhat[cpu_hrs];
                                    DDSIP_param->cpxdualwhat[cpu_hrs] = 0.20*old_cpxrelgap;
                                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                                    if (status)
                                    {
                                        printf("Failed to reset CPLEX relgap to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                    }
                                    else
                                    {
                                        if (DDSIP_param->outlev)
                                        {
                                            printf("temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                        }
                                    }
                                }
                            }
                            for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual2; cpu_hrs++)
                            {
                                if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_TILIM)
                                {
                                    old_cpxtimelim2 = DDSIP_param->cpxdualwhat2[cpu_hrs];
                                    DDSIP_param->cpxdualwhat2[cpu_hrs] = 1.5*old_cpxtimelim2;
                                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                    if (status)
                                    {
                                        printf("Failed to reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                    }
                                    else
                                    {
                                        if (DDSIP_param->outlev)
                                        {
                                            printf("temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                        }
                                    }
                                }
                                if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_EPGAP)
                                {
                                    old_cpxrelgap2 = DDSIP_param->cpxdualwhat2[cpu_hrs];
                                    DDSIP_param->cpxdualwhat2[cpu_hrs] = 0.30*old_cpxrelgap2;
                                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                    if (status)
                                    {
                                        printf("Failed to reset CPLEX relgap to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                    }
                                    else
                                    {
                                        if (DDSIP_param->outlev)
                                        {
                                            printf("temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                            fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                        }
                                    }
                                }
                            } // end for
                        }
                        if (!cb_status && !limits_reset && (noIncreaseCounter > 1 || (obj < DDSIP_node[DDSIP_bb->curnode]->bound - 1.e-4*fabs(DDSIP_node[DDSIP_bb->curnode]->bound))))
                        {
                            limits_reset = 1;
                            if (DDSIP_param->cb_changetol)
                            {
                                // the branched problems produced a worse bound than known - reset temporarily cplex tolerance and time limit
                                if (DDSIP_param->outlev)
                                {
                                    printf("reset CPLEX relgap and/or time limit.\n");
                                    fprintf(DDSIP_bb->moreoutfile, " reset CPLEX relgap and/or time limit.\n");
                                }
                                for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual; cpu_hrs++)
                                {
                                    if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_TILIM)
                                    {
                                        old_cpxtimelim = DDSIP_param->cpxdualwhat[cpu_hrs];
                                        DDSIP_param->cpxdualwhat[cpu_hrs] = 1.2*old_cpxtimelim;
                                        status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                                        if (status)
                                        {
                                            printf("Failed to reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                        }
                                        else
                                        {
                                            if (DDSIP_param->outlev)
                                            {
                                                printf("temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                                fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                            }
                                        }
                                    }
                                    if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_EPGAP)
                                    {
                                        old_cpxrelgap = DDSIP_param->cpxdualwhat[cpu_hrs];
                                        DDSIP_param->cpxdualwhat[cpu_hrs] = 0.20*old_cpxrelgap;
                                        status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                                        if (status)
                                        {
                                            printf("Failed to reset CPLEX relgap to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                        }
                                        else
                                        {
                                            if (DDSIP_param->outlev)
                                            {
                                                printf("temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                                fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                            }
                                        }
                                    }
                                } // end for
                                for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual2; cpu_hrs++)
                                {
                                    if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_TILIM)
                                    {
                                        old_cpxtimelim2 = DDSIP_param->cpxdualwhat2[cpu_hrs];
                                        DDSIP_param->cpxdualwhat2[cpu_hrs] = 1.5*old_cpxtimelim2;
                                        status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                        if (status)
                                        {
                                            printf("Failed to reset CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                        }
                                        else
                                        {
                                            if (DDSIP_param->outlev)
                                            {
                                                printf("temporarily reset CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                                fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                            }
                                        }
                                    }
                                    if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_EPGAP)
                                    {
                                        old_cpxrelgap2 = DDSIP_param->cpxdualwhat2[cpu_hrs];
                                        DDSIP_param->cpxdualwhat2[cpu_hrs] = 0.30*old_cpxrelgap2;
                                        status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                        if (status)
                                        {
                                            printf("Failed to reset CPLEX relgap 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                        }
                                        else
                                        {
                                            if (DDSIP_param->outlev)
                                            {
                                                printf("temporarily reset CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                                fprintf(DDSIP_bb->moreoutfile, "temporarily reset CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                            }
                                        }
                                    }
                                } // end for
                            }
                        }
//
                        if(DDSIP_param->outlev > 10)
                            fprintf(DDSIP_bb->moreoutfile, " obj - old_obj = %-16.12g - %-16.12g = %-10.6g (%g%%), \tnoIncreaseCounter= %d, last_weight= %g, next_weight= %g\n",obj, old_obj, obj - old_obj,100.*(obj - old_obj)/(fabs(old_obj)+1e-10), noIncreaseCounter, last_weight, next_weight);
//
                        // the evaluation in bestdual might have increased the bound
                        obj = DDSIP_bb->dualObjVal;
                    }
                    else
                    {
                        noIncreaseCounter = 0;
                        cycleCnt = 0;
                        next_weight = cb_get_last_weight (p);
                        cur_iters = DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt;
                        cb_get_center (p,center_point);
///////////////
                        // if the center point is not local_bestdual (may occur when maxsteps reached) - set center point to local_bestdual
                        if (cur_iters > 2 && DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 2] &&
                            memcmp(DDSIP_bb->local_bestdual, center_point, sizeof (double) * (DDSIP_bb->dimdual)) &&
                            (int) DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 2] != DDSIP_bb->dualitcnt &&
                            (obj - old_obj)/(fabs(old_obj) + 1.e-10) > 1.e-8 &&
                            DDSIP_bb->dualdescitcnt < DDSIP_bb->current_itlim)
                        {
                            if (DDSIP_param->outlev > 10)
                                fprintf(DDSIP_bb->moreoutfile,"############### set center point to bestdual from iter %g ##################\n",DDSIP_bb->local_bestdual[DDSIP_bb->dimdual + 2]);
                            if ((status = cb_set_new_center_point (p, DDSIP_bb->local_bestdual)))
                            {
                                fprintf (stderr, "set_new_center_point returned %d\n", status);
                                cb_destruct_problem (&p);
                                DDSIP_Free ((void **) &(minfirst));
                                DDSIP_Free ((void **) &(maxfirst));
                                DDSIP_Free ((void **) &(center_point));
                                return status;
                            }
                            obj = DDSIP_bb->dualObjVal;
                            cpu_hrs++;
                        }
                        start_weight = next_weight - last_weight;
                        if (DDSIP_bb->weight_reset == -1)
                        {
                            if (DDSIP_param->outlev > 10)
                                fprintf(DDSIP_bb->moreoutfile,"############### weight_reset 1\n");
                            start_weight = 0.;
                            next_weight = last_weight;
                        }
///////////////
                        if (DDSIP_param->outlev > 10)
                            fprintf(DDSIP_bb->moreoutfile,"############### iters in descent step: %d, up to now: repeated increase= %d, many_iters= %d, weight change by CB: %g ##################\n", cur_iters, repeated_increase, many_iters,  start_weight);
/////////
                        if (start_weight > 0.)
                        {
                            if (DDSIP_bb->dualdescitcnt == 1)
                                reduction_factor = 0.55;
                            else
                                reduction_factor = DDSIP_Dmax (0.65, reduction_factor);
                        }
                        if (DDSIP_param->cb_reduceWeight && last_weight >= DDSIP_Dmin(0.5,reduction_factor)*next_weight && start_weight == 0. &&
                                (cur_iters < 4 || (DDSIP_bb->dualdescitcnt == 1 && cur_iters < 10)))
                        {
                            if (repeated_increase < -1)
                                repeated_increase++;
                            many_iters = 0;
                            repeated_increase += DDSIP_Imax(0, 4 - cur_iters);
                            if (cur_iters > 3)
                            {
                                repeated_increase++;
                                if (DDSIP_bb->dualdescitcnt == 1)
                                {
                                    last_weight = next_weight;
                                    next_weight = last_weight * reduction_factor;
                                    if (DDSIP_bb->curnode < 3)
                                    {
                                        if (cur_iters < 5)
                                        {
                                            next_weight *= 0.5;
                                            weight_decreases++;
                                        }
                                        else if (cur_iters < 7)
                                            next_weight *= 0.65;
                                    }
                                    else if (cur_iters < 7)
                                        next_weight *= 0.8;
                                    cb_set_next_weight (p, next_weight);
                                    weight_decreases++;
///////////     ///////////
                                    if (DDSIP_param->outlev > 10)
                                        fprintf(DDSIP_bb->moreoutfile,"#############0. reduced next weight to %g,  repeated increase= %d ##################\n",next_weight,repeated_increase);
///////////     ///////////
                                }
                            }
                            else if ((abs(DDSIP_param->riskmod) <= 3 && next_weight > 5.e-5) || (abs(DDSIP_param->riskmod) > 3 && next_weight > 1.e-4))
                            {
                                if (repeated_increase > 2 + (next_weight < 0.05?2:0) + (next_weight < 0.005?2:0))
                                {
                                    if (weight_decreases > 1 && (next_weight > 0.049 || reduction_factor > 0.6))
                                    {
                                        reduction_factor = 0.75*reduction_factor + 0.05;
                                        weight_decreases = 1;
///////////     ///////////
                                        if (DDSIP_param->outlev > 10)
                                            fprintf(DDSIP_bb->moreoutfile,"############### repeated increase= %d, reduction factor decreased to %g ##################\n",repeated_increase, reduction_factor);
///////////     ///////////
                                    }
                                    else
                                        weight_decreases++;
                                    last_weight = next_weight;
                                    next_weight = last_weight * reduction_factor;
                                    repeated_increase = 0;
                                    if (DDSIP_bb->dualdescitcnt == 1)
                                    {
                                        next_weight *= 0.65;
                                        weight_decreases++;
                                        repeated_increase = 1;
                                    }
                                    cb_set_next_weight (p, next_weight);
///////////     ///////////
                                    if (DDSIP_param->outlev > 10)
                                        fprintf(DDSIP_bb->moreoutfile,"#############1. reduced next weight to %g,  repeated increase= %d ##################\n",next_weight,repeated_increase);
///////////     ///////////
                                }
                            }
                        }
                        else if (DDSIP_param->cb_increaseWeight && DDSIP_bb->dualdescitcnt > 1 && start_weight == 0. && last_weight*8. >= next_weight && last_weight > 0.99*next_weight)
                        {
                            if (cur_iters > 4)
                            {
                                if (cur_iters > 10 && (obj - old_obj)/(fabs(old_obj)+1e-6) < 1e-13 && last_weight > fabs(old_obj))
                                {
                                    last_weight = next_weight;
                                    next_weight = last_weight * 0.1;
                                    cb_set_next_weight (p, next_weight);
///////////     ///////////
                                    if (DDSIP_param->outlev > 10)
                                        fprintf(DDSIP_bb->moreoutfile,"############10. reduced next weight to %g ##################\n",next_weight);
///////////     ///////////
                                }
                                else if (cur_iters > 5)
                                {
                                    many_iters++;
                                    if (weight_decreases)
                                    {
                                        reduction_factor = 0.7*reduction_factor + 0.27;
///////////     ///////////
                                        if (DDSIP_param->outlev > 10)
                                            fprintf(DDSIP_bb->moreoutfile,"############### repeated increase= %d, reduction factor increased to %g ##################\n",repeated_increase, reduction_factor);
///////////     ///////////
                                    }
                                    if (cur_iters > 6)
                                    {
                                        repeated_increase--;
                                        many_iters++;
                                        if (cur_iters > 9 || (cur_iters > 7 && !weight_decreases))
                                        {
                                            if (next_weight < 1.10*last_weight)
                                            {
                                                last_weight = next_weight;
                                                if (last_weight > 1.)
                                                    next_weight = last_weight * 1.08;
                                                else if (last_weight > 5e-3)
                                                    next_weight = last_weight * 1.25;
                                                else
                                                    next_weight = last_weight * 2.50;
                                                cb_set_next_weight (p, next_weight);
///////////     ///////////
                                                if (DDSIP_param->outlev > 10)
                                                    fprintf(DDSIP_bb->moreoutfile,"#############1. increased next weight to %g,  repeated increase= %d ##################\n",next_weight,repeated_increase);
///////////     ///////////
                                            }
                                            if (cur_iters > 9)
                                            {
                                                many_iters++;
                                                if (repeated_increase > 0 || cur_iters > 11)
                                                    repeated_increase--;
                                            }
                                        }
                                        else if(many_iters && last_weight >= next_weight)
                                        {
                                            last_weight = next_weight;
                                            next_weight = last_weight * 1.04;
                                            cb_set_next_weight (p, next_weight);
///////////     ///////////
                                            if (DDSIP_param->outlev > 10)
                                                fprintf(DDSIP_bb->moreoutfile,"#############2. increased next weight to %g, current iters = %d,  many_iters = %d ##################\n",next_weight,cur_iters,many_iters);
///////////     ///////////
                                        }
                                    }
                                    else if(many_iters)
                                    {
                                        last_weight = next_weight;
                                        next_weight = last_weight * 1.008;
                                        cb_set_next_weight (p, next_weight);
///////////     ///////////
                                        if (DDSIP_param->outlev > 10)
                                            fprintf(DDSIP_bb->moreoutfile,"#############3. increased next weight to %g, current iters = %d,  many_iters = %d ##################\n",next_weight,cur_iters,many_iters);
///////////     ///////////
                                    }
                                    weight_decreases = 0;
                                }
                                else if (next_weight < 1.e-2 || (many_iters > 1 && many_iters < 5))
                                {
                                    last_weight = next_weight;
                                    if (last_weight < 1e5)
                                        next_weight = last_weight * 1.01;
                                    else
                                        next_weight = last_weight * 1.005;
                                    cb_set_next_weight (p, next_weight);
///////////     ///////////
                                    if (DDSIP_param->outlev > 10)
                                        fprintf(DDSIP_bb->moreoutfile,"#############4. increased next weight to %g, current iters = %d,  many_iters = %d ##################\n",next_weight,cur_iters,many_iters);
///////////     ///////////
                                }
                            }
                            else if (next_weight < 1.e-2) // the case cur_iters == 4
                            {
                                last_weight = next_weight;
                                next_weight = last_weight * 1.005;
                                cb_set_next_weight (p, next_weight);
///////////     ///////////
                                if (DDSIP_param->outlev > 10)
                                    fprintf(DDSIP_bb->moreoutfile,"#############4a increased next weight to %g, current iters = %d,  many_iters = %d ##################\n",next_weight,cur_iters,many_iters);
///////////     ///////////
                            }
                            if (many_iters > 4)
                            {
                                last_weight = next_weight;
                                next_weight *= 1.65;
                                cb_set_next_weight (p, next_weight);
                                repeated_increase = -1;
///////////     ///////////
                                if (DDSIP_param->outlev > 10)
                                    fprintf(DDSIP_bb->moreoutfile,"############### many_iters=%d-> too many iters, increased next weight to %g,  repeated increase= %d ##################\n",many_iters,next_weight,repeated_increase);
///////////     ///////////
                                many_iters = 0;
                            }
                        }
                        else if (DDSIP_bb->dualdescitcnt > 1 && DDSIP_param->cb_increaseWeight && (next_weight - last_weight) > 0.2*last_weight && cur_iters > 5)
                        {
                            many_iters++;
                            if (repeated_increase > -3)
                                repeated_increase--;
                            reduction_factor = 0.8*reduction_factor + 0.199;
///////////     ///////////
                            if (DDSIP_param->outlev > 10)
                                fprintf(DDSIP_bb->moreoutfile,"############### repeated increase= %d, reduction factor increased to %g ##################\n",repeated_increase, reduction_factor);
///////////     ///////////
                        }
///////////     ///////////
                        else if (DDSIP_bb->dualdescitcnt == 1 && DDSIP_param->cb_increaseWeight && (next_weight - last_weight) < 0.1*last_weight && cur_iters > DDSIP_param->cb_maxsteps + 4)
                        {
                                next_weight *= 1.1;
                                cb_set_next_weight (p, next_weight);
///////////     ///////////
                                if (DDSIP_param->outlev > 10)
                                    fprintf(DDSIP_bb->moreoutfile,"###############  increased next weight to %g ##################\n",next_weight);
///////////     ///////////
                        }
                        if (DDSIP_param->outlev > 10)
                        {
                            fprintf(DDSIP_bb->moreoutfile,"############### repeated_increase= %d  many_iters= %d    reduction_factor= %g ##################\n",repeated_increase,many_iters, reduction_factor);
                            fprintf(DDSIP_bb->moreoutfile, " obj - old_obj = %-16.12g - %-16.12g = %-10.6g (%g%%), \tnoIncreaseCounter= %d, last_weight= %g, next_weight= %g\n",obj, old_obj, obj - old_obj,100*(obj - old_obj)/(fabs(old_obj)+1e-6), noIncreaseCounter, last_weight, next_weight);
                        }
///////////     ///////////
                        old_obj = obj;

                    }
                }
                else
                {
                    start_weight = next_weight - last_weight;
                    if (DDSIP_bb->weight_reset == -1)
                    {
                        if (DDSIP_param->outlev > 10)
                            fprintf(DDSIP_bb->moreoutfile,"############### weight_reset 2\n");
                        next_weight = last_weight;
                        start_weight = 0.;
                        if (DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt > 7)
                        {
                            last_weight = next_weight = 1.3*last_weight;
                            cb_set_next_weight (p, next_weight);
                        }
                    }
                    if (start_weight > 0.)
                    {
                        if (DDSIP_bb->dualdescitcnt == 1)
                            reduction_factor = 0.55;
                        else
                            reduction_factor = DDSIP_Dmax (0.65, reduction_factor);
                    }
///////////////
                    if (DDSIP_param->outlev > 10)
                    {
                        fprintf(DDSIP_bb->moreoutfile,"############### iters in descent step: %d, up to now: repeated increase= %d, many_iters= %d, weight change by CB: %g ##################\n", DDSIP_bb->dualitcnt - DDSIP_bb->last_dualitcnt, repeated_increase, many_iters,  start_weight);
                        fprintf(DDSIP_bb->moreoutfile,"############### repeated_increase= %d  many_iters= %d    reduction_factor= %g ##################\n",repeated_increase,many_iters, reduction_factor);
                        fprintf(DDSIP_bb->moreoutfile, " obj - old_obj = %-16.12g - %-16.12g = %-10.6g (%g%%), \tnoIncreaseCounter= %d, last_weight= %g, next_weight= %g\n",obj, old_obj, obj - old_obj,100*(obj - old_obj)/(fabs(old_obj)+1e-6), noIncreaseCounter, last_weight, next_weight);
                    }
/////////
                }
                if (DDSIP_bb->dualObjVal > old_obj)
                {
                    old_obj = obj;
                    obj = DDSIP_bb->dualObjVal;
                }
            }

            if (DDSIP_bb->weight_reset != -1)
            {
                last_weight = cb_get_last_weight (p);
            }
            else
            {
                DDSIP_bb->weight_reset = 0;
                start_weight = 0.;
            }


            if (limits_reset && obj > DDSIP_node[DDSIP_bb->curnode]->bound - 1.e-5)
            {
                limits_reset = 0;
                if (DDSIP_param->cb_changetol)
                {
                    // revert to original cplex tolerance and time limit
                    for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual; cpu_hrs++)
                    {
                        if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_TILIM)
                        {
                            DDSIP_param->cpxdualwhat[cpu_hrs] = old_cpxtimelim;
                            status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                            if (status)
                            {
                                printf("Failed to revert CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                            }
                            else
                            {
                                if (DDSIP_param->outlev)
                                {
                                    printf("Reverted CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                    fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                }
                            }
                        }
                        if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_EPGAP)
                        {
                            DDSIP_param->cpxdualwhat[cpu_hrs] = old_cpxrelgap;
                            status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                            if (status)
                            {
                                printf("Failed to revert CPLEX relgap to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                            }
                            else
                            {
                                if (DDSIP_param->outlev)
                                {
                                    printf("Reverted CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                    fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                                }
                            }
                        }
                    } // end for
                    for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual2; cpu_hrs++)
                    {
                        if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_TILIM)
                        {
                            DDSIP_param->cpxdualwhat2[cpu_hrs] = old_cpxtimelim2;
                            status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            if (status)
                            {
                                printf("Failed to revert CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            }
                            else
                            {
                                if (DDSIP_param->outlev)
                                {
                                    printf("Reverted CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                    fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                }
                            }
                        }
                        if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_EPGAP)
                        {
                            DDSIP_param->cpxdualwhat2[cpu_hrs] = old_cpxrelgap2;
                            status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            if (status)
                            {
                                printf("Failed to revert CPLEX relgap 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            }
                            else
                            {
                                if (DDSIP_param->outlev)
                                {
                                    printf("Reverted CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                    fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                                }
                            }
                        }
                    }
                }
            }
            // Print iteration info
            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                time (&DDSIP_bb->cur_time);
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                if (!DDSIP_bb->curnode && DDSIP_bb->bestvalue < DDSIP_infty)
                {
                    if (!DDSIP_Equal (fabs (DDSIP_bb->bestvalue), 0.0))
                        rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / fabs (DDSIP_bb->bestvalue);
                    else
                        rgap = 100. * (DDSIP_bb->bestvalue -DDSIP_node[DDSIP_bb->curnode]->bound) / (fabs (DDSIP_bb->bestvalue) + DDSIP_param->accuracy);
                    rgap = DDSIP_Dmin (rgap, 100.0);
                    printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g               %10.4g%%  %3dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, rgap, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
                else
                {
                    printf ("  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                            DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                    fprintf (DDSIP_outfile, "  | %16d  %7d  %-16.12g %-11.6g                    %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                             DDSIP_bb->dualdescitcnt, DDSIP_bb->dualitcnt, DDSIP_bb->dualObjVal, last_weight, DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
                }
            }
            DDSIP_bb->last_dualitcnt = DDSIP_bb->dualitcnt;
            last_weight = next_weight;
            time (&DDSIP_bb->cur_time);
        }
        // store multipliers from node with highest dual bound up to now - if they are different
        if (DDSIP_bb->dualdescitcnt && DDSIP_bb->local_bestdual[DDSIP_bb->dimdual+1] == DDSIP_bb->curnode &&
                (!DDSIP_bb->bestdual || DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestdual_max ||
                 (obj > DDSIP_bb->bestvalue - 5.*fabs(DDSIP_bb->bestvalue)*(DDSIP_param->relgap)) ||
                  DDSIP_bb->curnode < 6)
           )
        {
            bbest_t * tmp_bestdual;
            cpu_hrs = 1;
            // check whether multipliers are zero
            for (wall_hrs = 0; wall_hrs < DDSIP_bb->dimdual; wall_hrs++)
            {
                if (DDSIP_bb->local_bestdual[wall_hrs] != 0.)
                {
                    break;
                }
            }
            if (!(wall_hrs < DDSIP_bb->dimdual))
                cpu_hrs = 0;
            // check whether we have already stored these multipliers (might happen if CB model cannot be improved and center remains constant)
            tmp_bestdual = DDSIP_bb->bestdual;
            while (tmp_bestdual && cpu_hrs)
            {
                if (!DDSIP_MultEqual (DDSIP_bb->local_bestdual, tmp_bestdual->dual))
                    tmp_bestdual = tmp_bestdual->next;
                else
                {
                    cpu_hrs = 0;
                    break;
                }
            }
            if (cpu_hrs)
            {
                tmp_bestdual = DDSIP_Alloc(sizeof (bbest_t), 1, "bestdual entry (DDSIP_DualOpt)");
                tmp_bestdual->dual = DDSIP_Alloc(sizeof (double), DDSIP_bb->dimdual, "bestdual entry (DDSIP_DualOpt)");
                memcpy (tmp_bestdual->dual, DDSIP_bb->local_bestdual, sizeof (double) * (DDSIP_bb->dimdual));
                tmp_bestdual->node_nr = DDSIP_bb->curnode;
                if (DDSIP_bb->curnode)
                    tmp_bestdual->bound   = DDSIP_Dmin (DDSIP_node[DDSIP_bb->curnode]->bound, 0.2*DDSIP_bb->bestvalue + 0.8*DDSIP_bb->bestbound);
                else
                    tmp_bestdual->bound   = DDSIP_node[DDSIP_bb->curnode]->bound;
                tmp_bestdual->weight  = DDSIP_bb->local_bestdual[DDSIP_bb->dimdual];
                DDSIP_bb->bestdual_cnt++;
                tmp_bestdual->next = DDSIP_bb->bestdual;
                DDSIP_bb->bestdual = tmp_bestdual;
                if (obj > DDSIP_bb->bestvalue - 5.*DDSIP_param->relgap*fabs(DDSIP_bb->bestvalue))
                {
                    DDSIP_bb->bestdual_max = DDSIP_Dmin(DDSIP_bb->bestdual_max, 0.5*DDSIP_bb->bestvalue + 0.5*DDSIP_bb->bestbound);
                }
                else
                {
                    DDSIP_bb->bestdual_max = DDSIP_Dmax(DDSIP_bb->bestdual_max, tmp_bestdual->bound);
                }
                if (DDSIP_param->outlev > 10)
                    fprintf (DDSIP_bb->moreoutfile, " ## added mult. from node %d to bestdual list, #entries: %d, bestdual_max: %17.14g, mult.=%g\n", DDSIP_bb->curnode, DDSIP_bb->bestdual_cnt, DDSIP_bb->bestdual_max, DDSIP_bb->bestdual->weight);
//######################
#ifdef DEBUG
if (DDSIP_param->outlev > 10)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
#endif
//######################
                if (DDSIP_bb->bestdual_cnt > DDSIP_param->cb_bestdualListLength)
                {
                    bbest_t * tmp_previous;
                    // find last entry in list and delete it
                    tmp_bestdual = tmp_previous = DDSIP_bb->bestdual;
                    while (tmp_bestdual->next)
                    {
                        tmp_previous = tmp_bestdual;
                        tmp_bestdual = tmp_bestdual->next;
                    }
                    tmp_previous->next = NULL;
                    if (DDSIP_param->outlev > 10)
                     fprintf (DDSIP_bb->moreoutfile, " ##         set %p->next of entry from node %d of bestdual list to %p\n", tmp_previous, tmp_previous->node_nr, tmp_previous->next);
                    DDSIP_bb->bestdual_cnt--;
                    if (DDSIP_param->outlev > 10)
                     fprintf (DDSIP_bb->moreoutfile, " ## delete last entry from node %d from bestdual list, #entries: %d\n", tmp_bestdual->node_nr, DDSIP_bb->bestdual_cnt);
                    DDSIP_Free ((void **) &(tmp_bestdual->dual));
                    DDSIP_Free ((void **) &(tmp_bestdual));
//######################
#ifdef DEBUG
if (DDSIP_param->outlev > 10)
{
int jj = 1;
bbest_t * tmp1_bestdual = DDSIP_bb->bestdual;
while (tmp1_bestdual)
{
   fprintf (DDSIP_bb->moreoutfile, " ##   %d:  node %2d  bound %16.12g, weight %g, (%p next= %p)\n",jj++,tmp1_bestdual->node_nr,tmp1_bestdual->bound, tmp1_bestdual->weight, tmp1_bestdual, tmp1_bestdual->next);
   tmp1_bestdual = tmp1_bestdual->next;
}
}
#endif
                }
                memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->local_bestdual, sizeof (double) * (DDSIP_bb->dimdual + 3));
//######################
                if (DDSIP_param->outlev > 10)
                    fprintf (DDSIP_bb->moreoutfile," +++-> bb->bestdual updated in node %d, desc. it. %d  (weight %g, bound %g)\n", DDSIP_bb-> curnode, DDSIP_bb->dualdescitcnt, DDSIP_bb->bestdual->weight, DDSIP_bb->bestdual->bound);
//######################
            }
        }
#ifdef DEBUG
        else
        {
            if (DDSIP_param->outlev > 20)
                fprintf (DDSIP_bb->moreoutfile, " ##  mult. of node %2d  bound %16.12g, bestdual_max %16.12g not added to bestdual list after desc.it. %d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestdual_max, DDSIP_bb->dualdescitcnt);
        }
#endif
//######################
        if (DDSIP_param->outlev)
        {
            printf ("   --------- ");
            fprintf (DDSIP_outfile, "   --------- ");
            i_scen = cb_termination_code (p);
            if (i_scen == 32)
            {
                 DDSIP_node[DDSIP_bb->curnode]->cbReturn32 = DDSIP_node[DDSIP_bb->curnode]->depth;
                 if (DDSIP_param->outlev > 20)
                    fprintf (DDSIP_bb->moreoutfile, "########## DDSIP_node[%d]->cbReturn32 = %d\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->depth);
            }
            //if (obj > DDSIP_bb->bestvalue)
            if ((obj - DDSIP_bb->bestvalue) > DDSIP_param->accuracy*(fabs(DDSIP_bb->bestvalue)+1.e-10))
            {
                DDSIP_Print2 ("termination status: cutoff. --------------------------------------------------------------", "\n", 0, 0);
                DDSIP_bb->skip = 2;
                DDSIP_bb->cutoff++;
                DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
            }
            else if (!DDSIP_bb->violations)
            {
                DDSIP_Print2 ("termination status: no violation of nonanticipativity. ------------------------------------", "\n", 0, 0);
                DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
            }
            //else if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - DDSIP_Dmin (0.5*DDSIP_param->relgap, 2.5e-9)*(fabs(DDSIP_bb->bestvalue) + 1e-12))
            else if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - DDSIP_Dmin (0.5*DDSIP_param->relgap, 8e-9)*(fabs(DDSIP_bb->bestvalue) + 1e-14))
            {
                if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue +  2.e-10 * (fabs(DDSIP_bb->bestvalue) + 1e-12))
                {
                    DDSIP_bb->skip = 2;
                    DDSIP_bb->cutoff++;
                    if (DDSIP_param->outlev > 20)
                        fprintf (DDSIP_bb->moreoutfile, "########## skip=+2 for DDSIP_bb->bestsol_in_curnode=%d, DDSIP_node[%d]->bound (%20.15g) - bestvalue = %.8g\n", DDSIP_bb->bestsol_in_curnode, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue);
                }
                else if (DDSIP_bb->found_optimal_node && /* !DDSIP_bb->bestsol_in_curnode && */
                     DDSIP_node[DDSIP_bb->curnode]->violations > 0.4*DDSIP_param->scenarios && (DDSIP_bb->no_reduced_front > 1))
                {
                     DDSIP_bb->skip = -2;
                     if (DDSIP_param->outlev > 20)
                         fprintf (DDSIP_bb->moreoutfile, "########## skip=-2 for DDSIP_bb->bestsol_in_curnode=%d, DDSIP_node[%d]->bound (%20.15g) - bestvalue = %.8g\n", DDSIP_bb->bestsol_in_curnode, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue);
                }
                if (DDSIP_bb->no_reduced_front == 1)
                    DDSIP_Print2 ("termination status: gap reached. --------------------------------------------------------------", "\n", 0, 0);
                else
                {
                    DDSIP_Print2 ("termination status: within relative gap. Number of violations of nonanticipativity: ", " ------------------------------------------\n", 1.*DDSIP_node[DDSIP_bb->curnode]->violations, 1);
                    if (DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue - 1.e-11*fabs(DDSIP_bb->bestvalue))
                        DDSIP_bb->skip = -2;
                    if ((DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bestvalue*nfactor) &&
                        (!(DDSIP_bb->found_optimal_node) || (DDSIP_bb->found_optimal_node && DDSIP_node[DDSIP_bb->curnode]->bound > DDSIP_bb->bound_optimal_node)))
                    {
    //////////////////////////////////////////////////////////////////
                           if (DDSIP_param->outlev > 20)
                               fprintf (DDSIP_bb->moreoutfile, "########## setting found_optimal_node: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bestvalue (%20.15g) = %.8g, - bestvalue*nfactor (%20.15g) = %.8g\n",
                                        DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bestvalue, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue,
                                        DDSIP_bb->bestvalue*nfactor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*nfactor);
    //////////////////////////////////////////////////////////////////
                        DDSIP_bb->found_optimal_node = DDSIP_bb->curnode;
                        DDSIP_bb->bound_optimal_node = DDSIP_node[DDSIP_bb->curnode]->bound;
                    }
    //////////////////////////////////////////////////////////////////
                    else
                    {
                       if (DDSIP_param->outlev > 20)
                           fprintf (DDSIP_bb->moreoutfile, "########## else: found_optimal_node= %d, DDSIP_node[%d]->bound (%20.15g) - bound_optimal_node (%20.15g) = %.8g, - bestvalue*nfactor (%20.15g) = %.8g\n",
                                    DDSIP_bb->found_optimal_node, DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->bound, DDSIP_bb->bound_optimal_node, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bound_optimal_node,
                                    DDSIP_bb->bestvalue*nfactor, DDSIP_node[DDSIP_bb->curnode]->bound - DDSIP_bb->bestvalue*nfactor);
                    }
    //////////////////////////////////////////////////////////////////
                }
                DDSIP_node[DDSIP_bb->curnode]->leaf = 1;
            }
            else if (!i_scen)
            {
                if (DDSIP_bb->dualdescitcnt >= DDSIP_bb->current_itlim)
                    DDSIP_Print2 ("termination status: descent iteration limit exceeded.", " --------------------------------------------------------------\n", 0, 0);
                else if (DDSIP_bb->dualitcnt >= DDSIP_param->cbtotalitlim)
                    DDSIP_Print2 ("termination status: total iteration limit exceeded.", " --------------------------------------------------------------\n", 0, 0);
                else if (cycleCnt >= 2)
                    DDSIP_Print2 ("termination status: no bound increase achieved.", " --------------------------------------------------------------\n", 0, 0);
                else if (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time) > DDSIP_param->timelim)
                    DDSIP_Print2 ("Time limit reached.", "\n", 0, 0);
                else if (DDSIP_killsignal)
                    DDSIP_Print2 ("termination status: User interrupt.", " --------------------------------------------------------------\n", 0, 0);
                else
                    fprintf (DDSIP_outfile, " Unidentified reason for stopping. Termination code of ConicBundle: %d\n", i_scen);
            }
            else
            {
                cb_print_termination_code (p);
                fprintf (DDSIP_outfile, " termination code of ConicBundle: %d ", (i_scen));
                if (i_scen == 1)
                {
                    fprintf (DDSIP_outfile, "Relative precision criterion satisfied.\n");
                    // signal that all scenario solutions should be checked in ub
                    DDSIP_bb->keepSols = -11;
                }
                else if (i_scen == 2)
                    fprintf (DDSIP_outfile, "Time limit exceeded.\n");
                else if (i_scen == 4)
                    fprintf (DDSIP_outfile, "Maximum number of function reevaluations exceeded.\n");
                else if (i_scen == 8)
                    fprintf (DDSIP_outfile, "Maximum number of quadratic subproblem failures exceeded.\n");
                else if (i_scen == 16)
                    fprintf (DDSIP_outfile, "Maximum number of model evaluation failures exceeded.\n");
                else if (i_scen == 32)
                {
                    fprintf (DDSIP_outfile, "Maximum number of failures to increase the augmented model value exceeded.\n");
                    DDSIP_bb->local_bestdual[DDSIP_bb->dimdual] = DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual] = cb_get_last_weight(p);
                    if (DDSIP_param->outlev > 19)
                    {
                        fprintf (DDSIP_bb->moreoutfile, "\nMaximum number of failures to increase the augmented model value exceeded in node %d:   changed weight in bestdual to last weight= %g\n", DDSIP_bb->curnode, DDSIP_bb->local_bestdual[DDSIP_bb->dimdual]);
                    }
                }
                else if (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time) > DDSIP_param->timelim)
                    fprintf (DDSIP_outfile, "Total time limit exceeded.\n");
                else
                    fprintf (DDSIP_outfile, "\n");
            }
        }
////////////////////////////////////////////////////////////////////////////// computation of bounds for premature stop only when cutoff not yet decided
        if (DDSIP_param->prematureStop && !DDSIP_killsignal && (DDSIP_bb->curnode || DDSIP_bb->initial_multiplier) && !DDSIP_node[DDSIP_bb->curnode]->leaf)
        {
            // Initialize multipliers with zero
            memset (DDSIP_bb->startinfo_multipliers, '\0', sizeof (double) * (DDSIP_bb->dimdual));
            DDSIP_bb->keepSols = 0;
            if ((status = cb_set_new_center_point (p, DDSIP_bb->startinfo_multipliers)))
            {
                fprintf (stderr, "set_new_center_point returned %d\n", status);
                if (DDSIP_param->outlev)
                    fprintf (DDSIP_bb->moreoutfile, "set_new_center_point returned %d\n", status);
                cb_destruct_problem (&p);
                DDSIP_Free ((void **) &(minfirst));
                DDSIP_Free ((void **) &(maxfirst));
                DDSIP_Free ((void **) &(center_point));
                return status;
            }
            memcpy (DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag, DDSIP_node[DDSIP_bb->curnode]->subbound, DDSIP_param->scenarios*sizeof(double));
            DDSIP_node[DDSIP_bb->curnode]->BoundNoLag = DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[0] * DDSIP_data->prob[0];
            for (j = 1; j < DDSIP_param->scenarios; j++)
                DDSIP_node[DDSIP_bb->curnode]->BoundNoLag += DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j] * DDSIP_data->prob[j];
        //////////////////////////////
            if (DDSIP_param->outlev > 21)
            {
                for (j = 0; j < DDSIP_param->scenarios; j++)
                {
                    fprintf (DDSIP_bb->moreoutfile, " ## bound zero mult. scen %3d:  %24.15g\n", j+1, DDSIP_node[DDSIP_bb->curnode]->scenBoundsNoLag[j]);
                }
                fprintf (DDSIP_bb->moreoutfile, " ## -----------------------\n ## node %d  BoundNoLag:         %24.15g\n ## -----------------------\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->BoundNoLag);
            }
        //////////////////////////////
            obj = DDSIP_bb->dualObjVal;
            if (DDSIP_param->outlev)
            {
                DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
                time (&DDSIP_bb->cur_time);
                DDSIP_translate_time (difftime(DDSIP_bb->cur_time,DDSIP_bb->start_time),&wall_hrs,&wall_mins,&wall_secs);
                fprintf (DDSIP_outfile, "  |--------------------- %4d  %-16.12g ----------------- zero mult.   %-20.14g %30dh %02d:%02.0f  %3dh %02d:%02.0f\n",
                                             DDSIP_bb->dualitcnt, DDSIP_bb->currentDualObjVal,
                                             DDSIP_node[DDSIP_bb->curnode]->bound, wall_hrs,wall_mins,wall_secs, cpu_hrs,cpu_mins,cpu_secs);
            }
        
            // most probably this never happens...
            if (obj > DDSIP_node[DDSIP_bb->curnode]->bound)
            {
                memcpy (DDSIP_bb->local_bestdual, DDSIP_bb->startinfo_multipliers, sizeof (double) * (DDSIP_bb->dimdual));
    ////////////////////////////
    if (DDSIP_param->outlev && inherited_bound > -DDSIP_infty)
       fprintf (DDSIP_outfile, " ###### zero mult. increased bound %22.14g by %g\n", DDSIP_node[DDSIP_bb->curnode]->bound, obj-DDSIP_node[DDSIP_bb->curnode]->bound);
    ////////////////////////////
    
                DDSIP_node[DDSIP_bb->curnode]->bound = obj;
            }
            if (DDSIP_node[DDSIP_bb->curnode]->depth < DDSIP_Imax(DDSIP_param->cb_depth, 3))
            {
                // do upper bounding for these solutions, too
                if (!DDSIP_bb->curnode && DDSIP_bb->cutAdded && DDSIP_param->outlev)
                {
                    fprintf (DDSIP_outfile, " %6d%101d cuts\n", DDSIP_bb->curnode, DDSIP_bb->cutAdded);
                }
                DDSIP_bb->cutAdded = 0;
                DDSIP_EvaluateScenarioSolutions (&comb);
                if (DDSIP_bb->cutAdded && DDSIP_param->outlev)
                {
                    fprintf (DDSIP_outfile, " %6d%101d cuts\n", DDSIP_bb->curnode, DDSIP_bb->cutAdded);
                }
                DDSIP_bb->skip = 0;
            }
            if (DDSIP_param->outlev)
            {
                fprintf (DDSIP_outfile, "   ------------------------------------------------------------------------------------------------------------------------------\n");
            }
        }
////////////////////////////////////////////////////////////////////////////// computation of bounds for premature stop only when no cutoff
    }
    memcpy (DDSIP_node[DDSIP_bb->curnode]->dual, DDSIP_bb->local_bestdual, sizeof (double) * (DDSIP_bb->dimdual));

    if (DDSIP_bb->dualdescitcnt)
        DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual + 1] = DDSIP_bb->curnode;
    if (DDSIP_param->outlev > 19)
    {
        fprintf (DDSIP_bb->moreoutfile, "\nFinal lambda for node %d:   weight in bestdual= %g, last weight= %g\n MULTIPLIER\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_bb->dimdual], last_weight);
        for (j = 0; j < DDSIP_bb->dimdual; j++)
        {
            fprintf (DDSIP_bb->moreoutfile, " %14.8g", DDSIP_node[DDSIP_bb->curnode]->dual[j]);
            if (!((j+1)%10))
                fprintf (DDSIP_bb->moreoutfile, "\n");
        }
        fprintf (DDSIP_bb->moreoutfile, "\n");
    }

    cb_destruct_problem (&p);
    for (j = 0; j < DDSIP_bb->firstvar; j++)
    {
        minfirst[j] = DDSIP_infty;
        maxfirst[j] = -DDSIP_infty;
    }
    // Reset first stage solutions to the ones that gave the best bound
    for (j = 0; j < DDSIP_param->scenarios; j++)
    {
        if (DDSIP_bb->bestfirst[j].first_sol)
        {
            if (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])
            {
                if ((cnt = (int) ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j])[DDSIP_bb->firstvar] - 0.9)))
                    for (i_scen = j + 1; cnt && i_scen < DDSIP_param->scenarios; i_scen++)
                    {
                        if ((((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]) == (((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i_scen]))
                        {
                            ((DDSIP_node[DDSIP_bb->curnode])->first_sol)[i_scen] = NULL;
                            cnt--;
                        }
                    }
                DDSIP_Free ((void **) &(((DDSIP_node[DDSIP_bb->curnode])->first_sol)[j]));
            }
            DDSIP_node[DDSIP_bb->curnode]->first_sol[j]   = (DDSIP_bb->bestfirst[j]).first_sol;
            (DDSIP_node[DDSIP_bb->curnode]->cursubsol)[j] = (DDSIP_bb->bestfirst[j]).cursubsol;
            (DDSIP_node[DDSIP_bb->curnode]->subbound)[j]  = (DDSIP_bb->bestfirst[j]).subbound;
            DDSIP_bb->bestfirst[j].first_sol = NULL;
            for (cnt=0; cnt<DDSIP_bb->firstvar; cnt++)
            {
                minfirst[cnt]=DDSIP_Dmin (minfirst[cnt],DDSIP_node[DDSIP_bb->curnode]->first_sol[j][cnt]);
                maxfirst[cnt]=DDSIP_Dmax (maxfirst[cnt],DDSIP_node[DDSIP_bb->curnode]->first_sol[j][cnt]);
            }
        }
    }
    // revert changes to  cplex tolerance and time limit
    if (limits_reset)
    {
        // revert to original cplex tolerance and time limit
        limits_reset = 0;
        if (DDSIP_param->cb_changetol)
        {
            for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual; cpu_hrs++)
            {
                if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_TILIM)
                {
                    DDSIP_param->cpxdualwhat[cpu_hrs] = old_cpxtimelim;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to revert CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("Reverted CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX time limit to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                        }
                    }
                }
                if (DDSIP_param->cpxdualwhich[cpu_hrs] == CPX_PARAM_EPGAP)
                {
                    DDSIP_param->cpxdualwhat[cpu_hrs] = old_cpxrelgap;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich[cpu_hrs], DDSIP_param->cpxdualwhat[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to revert CPLEX relgap to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("Reverted CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX relgap tolerance to %g.\n", DDSIP_param->cpxdualwhat[cpu_hrs]);
                        }
                    }
                }
            } // end for
            for (cpu_hrs=0; cpu_hrs < DDSIP_param->cpxnodual2; cpu_hrs++)
            {
                if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_TILIM)
                {
                    DDSIP_param->cpxdualwhat2[cpu_hrs] = old_cpxtimelim2;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to revert CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("Reverted CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX time limit 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                        }
                    }
                }
                if (DDSIP_param->cpxdualwhich2[cpu_hrs] == CPX_PARAM_EPGAP)
                {
                    DDSIP_param->cpxdualwhat2[cpu_hrs] = old_cpxrelgap2;
                    status = CPXsetdblparam (DDSIP_env, DDSIP_param->cpxdualwhich2[cpu_hrs], DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    if (status)
                    {
                        printf("Failed to revert CPLEX relgap 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                    }
                    else
                    {
                        if (DDSIP_param->outlev)
                        {
                            printf("Reverted CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                            fprintf(DDSIP_bb->moreoutfile, "Reverted CPLEX relgap tolerance 2nd opt. to %g.\n", DDSIP_param->cpxdualwhat2[cpu_hrs]);
                        }
                    }
                }
            }
        }
    }
    // Reset cplex parameters
    if (DDSIP_param->cpxnodual)
    {
        status = DDSIP_SetCpxPara (DDSIP_param->cpxno, DDSIP_param->cpxisdbl, DDSIP_param->cpxwhich, DDSIP_param->cpxwhat);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to reset CPLEX parameters (DualOpt) \n");
            DDSIP_Free ((void **) &(minfirst));
            DDSIP_Free ((void **) &(maxfirst));
            DDSIP_Free ((void **) &(center_point));
            return status;
        }
    }
    //determine variable to branch on
    diff = -1.;
    if (!DDSIP_node[DDSIP_bb->curnode]->leaf && !DDSIP_killsignal)
    {
        for (j = 0; j < DDSIP_bb->firstvar; j++)
        {
            maxfirst[j] -= minfirst[j];
            diff = DDSIP_Dmax (diff, maxfirst[j]);
            if (fabs(maxfirst[j])>DDSIP_param->nulldisp)
            {
                if (DDSIP_param->outlev>40)
                    fprintf (DDSIP_bb->moreoutfile," ---- Deviation of variable %d : %g\n",j,maxfirst[j]);
            }
        }
        DDSIP_node[DDSIP_bb->curnode]->dispnorm = diff;
        status = DDSIP_GetBranchIndex (maxfirst);
    }
    DDSIP_Free ((void **) &(minfirst));
    DDSIP_Free ((void **) &(maxfirst));
    DDSIP_Free ((void **) &(center_point));
    return 0;
}

#endif
