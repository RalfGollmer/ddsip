/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
    Last modification: 06.02.2016
	Description:
	The procedures in this file are used to modify the problem according to
	the current scenarios during the branch-and-bound algorithm.

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

#include <DDSIP.h>
#include <DDSIPconst.h>

/*==========================================================================
  Function changes bounds according to current node of b&b tree */
int
DDSIP_ChgBounds (int print)
{
    int i, j, status;

    int *index = (int *) DDSIP_Alloc (sizeof (int), DDSIP_bb->curbdcnt, "index(chgbounds)");

    for (i = 0; i < DDSIP_bb->curbdcnt; i++)
        index[i] = DDSIP_bb->firstindex[DDSIP_bb->curind[i]];
    if (print && DDSIP_param->outlev > 3)
    {
        if (DDSIP_bb->curbdcnt)
            fprintf (DDSIP_bb->moreoutfile, "New bounds\n nr  variable_nr      .    lb                ub   variable_name\n");
        for (i = 0; i < DDSIP_bb->curbdcnt; i++)
        {
            j = index[i];
            status = CPXgetcolname (DDSIP_env, DDSIP_lp, DDSIP_bb->name_buffer, DDSIP_bb->n_buffer, DDSIP_bb->n_buffer_len, &j, j, j);
            if (status)
                fprintf (stderr," Error when querying name of variable %d: %d\n",j,status);
            fprintf (DDSIP_bb->moreoutfile, "%3d  %6d  %16.14g  %16.14g   %s\n", i, index[i], DDSIP_bb->curlb[i], DDSIP_bb->curub[i], DDSIP_bb->n_buffer);
        }
    }

    status = CPXchgbds (DDSIP_env, DDSIP_lp, DDSIP_bb->curbdcnt, index, DDSIP_bb->lbident, DDSIP_bb->curlb);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change lower bounds\n");
        return status;
    }

    status = CPXchgbds (DDSIP_env, DDSIP_lp, DDSIP_bb->curbdcnt, index, DDSIP_bb->ubident, DDSIP_bb->curub);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change upper bounds\n");
        return status;
    }

    DDSIP_Free ((void **) &(index));
    return status;
}

//==========================================================================
// Function is used to change problem parameters for each scenario
// scen is the number of the current scenarios
int
DDSIP_ChgProb (int scen)
{
    int j, i, status = 0;
    int m = DDSIP_Imax (DDSIP_param->stocmat, DDSIP_Imax (DDSIP_param->stocrhs, DDSIP_param->stoccost));
    //int *index = (int *) DDSIP_Alloc (sizeof (int), m, "index(chgprob)");
    char **colname;
    char *colstore;

    double *value = (double *) DDSIP_Alloc (sizeof (double), m, "value(chgprob)");
    double *cost = (double *) DDSIP_Alloc (sizeof (double), DDSIP_bb->firstvar, "cost(chgprob)");

    // Use multipliers passed by CB
    if (DDSIP_param->cb && DDSIP_bb->DDSIP_step != neobj && DDSIP_bb->DDSIP_step != adv && DDSIP_bb->DDSIP_step != eev)
    {
        // Original costs
        for (i = 0; i < DDSIP_bb->firstvar; i++)
            cost[i] = DDSIP_bb->cost[i];

#ifdef CONIC_BUNDLE
        if (DDSIP_bb->DDSIP_step == dual)
            // The additional costs changes for ConicBundle iterations
        {
            if (DDSIP_param->outlev >= DDSIP_current_lambda_outlev)
            {
                if(!scen)
                {
                    fprintf (DDSIP_bb->moreoutfile, "\n Lagrangemult. in OF: in CB from %p ---------------------\n", DDSIP_node[DDSIP_bb->curnode]->dual);
                    fprintf (DDSIP_bb->moreoutfile, "DDSIP_ChgProb: Current lambda for node %d:\n", DDSIP_bb->curnode);
                    for (i = 0; i < DDSIP_bb->dimdual; i++)
                    {
                        fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                        if (!((i+1)%10))
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    fprintf (DDSIP_bb->moreoutfile, "\n");
                }
            }
            if (DDSIP_param->outlev > 74)
            {
                colname = (char **) DDSIP_Alloc (sizeof (char *), (DDSIP_bb->firstvar + DDSIP_bb->secvar), "colname(Change)");
                colstore = (char *) DDSIP_Alloc (sizeof (char), (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, "colstore(Change)");
                status = CPXgetcolname (DDSIP_env, DDSIP_lp, colname, colstore,
                                        (DDSIP_bb->firstvar + DDSIP_bb->secvar) * DDSIP_ln_varname, &j, 0, DDSIP_bb->firstvar + DDSIP_bb->secvar - 1);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to get column names (Change)\n");
                    fprintf (DDSIP_outfile, "ERROR: Failed to get column names (Change)\n");
                    if (DDSIP_param->outlev)
                        fprintf (DDSIP_bb->moreoutfile, "ERROR: Failed to get column names (Change)\n");
                    return 100;
                }
            }
            for (i = 0; i < DDSIP_bb->firstvar; i++)
            {
                for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                        j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                {
                    if (fabs (DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]]) > 4e-16)
                    {
                        if (DDSIP_param->outlev > 74)
                        {
                            fprintf (DDSIP_bb->moreoutfile, " first-stage var: %d ", i);
                            fprintf (DDSIP_bb->moreoutfile, " (%s),", colname[DDSIP_bb->firstindex[i]]);
                            fprintf (DDSIP_bb->moreoutfile, " scen: %d, j=%d, \tnaind[j]= %d: koeff: %g, dual=%g, \tcost add: %g", scen,j, DDSIP_data->naind[j],  DDSIP_data->naval[j], DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]],  DDSIP_data->naval[j] * DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]] / DDSIP_data->prob[scen]);
                            fprintf (DDSIP_bb->moreoutfile, " \t%16.12g changed to", cost[i]);
                        }

                        cost[i] += DDSIP_data->naval[j] * DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]] / DDSIP_data->prob[scen];
                        if (DDSIP_param->outlev > 74)
                        {
                            fprintf (DDSIP_bb->moreoutfile, " %16.12g\n", cost[i]);
                        }
                    }
                }
            }
            if (DDSIP_param->outlev > 74)
            {
                DDSIP_Free ((void**) &(colname));
                DDSIP_Free ((void **) &(colstore));
            }
        }
        else if (DDSIP_bb->DDSIP_step == solve)
        {
            if (DDSIP_param->outlev >= DDSIP_current_lambda_outlev)
            {
                if(!scen)
                {
                    fprintf (DDSIP_bb->moreoutfile,
                             "\n Lagrangemult. in OF: in solve from %p ---------------------\n", DDSIP_node[DDSIP_bb->curnode]->dual);
                    fprintf (DDSIP_bb->moreoutfile, "DDSIP_ChgProb: Current lambda for node %d:\n", DDSIP_bb->curnode);
                    for (i = 0; i < DDSIP_bb->dimdual; i++)
                    {
                        fprintf (DDSIP_bb->moreoutfile, " %14.8g,", DDSIP_node[DDSIP_bb->curnode]->dual[i]);
                        if (!((i+1)%10))
                            fprintf (DDSIP_bb->moreoutfile, "\n");
                    }
                    fprintf (DDSIP_bb->moreoutfile, "\n");
                }
                if (DDSIP_param->outlev > 74)
                {
                    for (i = 0; i < DDSIP_bb->firstvar; i++)
                        for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                                j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                            if (fabs (DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]]) > 1e-16)
                                fprintf (DDSIP_bb->moreoutfile, " first-stage var: %d, scen: %d, j=%d, naind[j]= %d: koeff: %g, dual=%g, cost add: %g\n", i, scen,j, DDSIP_data->naind[j],  DDSIP_data->naval[j], DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]],  DDSIP_data->naval[j] * DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]] / DDSIP_data->prob[scen]);
                    fprintf (DDSIP_bb->moreoutfile, "\n");
                }
            }
            for (i = 0; i < DDSIP_bb->firstvar; i++)
                for (j = DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i];
                        j < DDSIP_data->nabeg[scen * DDSIP_bb->firstvar + i] + DDSIP_data->nacnt[scen * DDSIP_bb->firstvar + i]; j++)
                    cost[i] += DDSIP_data->naval[j] * DDSIP_node[DDSIP_bb->curnode]->dual[DDSIP_data->naind[j]] / DDSIP_data->prob[scen];
        }
#endif

        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_bb->firstvar, DDSIP_bb->firstindex, cost);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to update objective coefficients\n");
            return status;
        }
    }
    // Change rhs
    if (DDSIP_param->stocrhs)
    {
        if (scen == -1)
            // Do this for ExpValProb
        {
            for (j = 0; j < DDSIP_param->stocrhs; j++)
            {
                value[j] = 0.0;
                for (i = 0; i < DDSIP_param->scenarios; i++)
                    value[j] += DDSIP_data->prob[i] * DDSIP_data->rhs[i * DDSIP_param->stocrhs + j];
                // Numerical errors ?
                if (DDSIP_Equal (value[j], 0.0))
                    value[j] = 0.0;
            }
        }
        else
            // Do this for LowerBound, UB
        {
            for (j = 0; j < DDSIP_param->stocrhs; j++)
            {
                value[j] = DDSIP_data->rhs[scen * DDSIP_param->stocrhs + j];
            }
        }

        status = CPXchgrhs (DDSIP_env, DDSIP_lp, DDSIP_param->stocrhs, DDSIP_data->rhsind, value);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change rhs, return code %d\n", status);
            return status;
        }

    }
    // Change costs if required
    if (DDSIP_param->stoccost)
    {
        if (scen == -1)
        {
            // Do this for ExpValProb
            for (j = 0; j < DDSIP_param->stoccost; j++)
            {
                value[j] = 0.0;
                for (i = 0; i < DDSIP_param->scenarios; i++)
                    value[j] += DDSIP_data->prob[i] * DDSIP_data->cost[i * DDSIP_param->stoccost + j];
                // Numerical errors ?
                if (DDSIP_Equal (value[j], 0.0))
                    value[j] = 0.0;
            }
        }
        else
        {
            // Do this for LowerBound
            for (j = 0; j < DDSIP_param->stoccost; j++)
                if (DDSIP_param->riskmod == 3)
                    value[j] = (1 - DDSIP_param->riskweight) * DDSIP_data->cost[scen * DDSIP_param->stoccost + j];
                else
                    value[j] = DDSIP_data->cost[scen * DDSIP_param->stoccost + j];
        }

        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_param->stoccost, DDSIP_data->costind, value);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective, return code %d\n", status);
            return status;
        }
    }
    // Change matrix if required
    if (DDSIP_param->stocmat)
    {
        if (scen == -1)
            // Do this for ExpValProb
            for (j = 0; j < DDSIP_param->stocmat; j++)
            {
                value[j] = 0.0;
                for (i = 0; i < DDSIP_param->scenarios; i++)
                    value[j] += DDSIP_data->prob[i] * DDSIP_data->matval[i * DDSIP_param->stocmat + j];
                // Numerical errors ?
                if (DDSIP_Equal (value[j], 0.0))
                    value[j] = 0.0;
            }
        else
            // Do this for LowerBound
            for (j = 0; j < DDSIP_param->stocmat; j++)
                value[j] = DDSIP_data->matval[scen * DDSIP_param->stocmat + j];

        status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_param->stocmat, DDSIP_data->matrow, DDSIP_data->matcol, value);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change matrix, return code %d\n", status);
            return status;
        }
    }

    //DDSIP_Free ((void **) &(index));
    DDSIP_Free ((void **) &(value));
    DDSIP_Free ((void **) &(cost));
    return status;
}
