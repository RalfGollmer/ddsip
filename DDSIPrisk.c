/*  Authors:           Andreas M"arkert, Ralf Gollmer
	Copyright to:      University of Duisburg-Essen
    Language:          C
    Last modification: 06.02.2016
	Description:
	This file implements the (mean-)risk models supported by DDSIP.

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

int DDSIP_ExpExcess (void);
int DDSIP_ExcessProb (void);
int DDSIP_SemDev (void);
int DDSIP_WorstCase (void);
int DDSIP_TVaR (void);


//==========================================================================
// Risk modelling: expected shortfall below target
int
DDSIP_ExpExcess (void)
{
    int i, status = 0;
    int *rowlist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1),
                                        "rowlist(RiskModel)");
    int *collist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1),
                                        "collist(RiskModel)");
    double *vallist = (double *) DDSIP_Alloc (sizeof (double),
                      (DDSIP_data->novar + 1),
                      "vallist (RiskModel)");
    double *rhs = (double *) DDSIP_Alloc (sizeof (double), 1, "rhs(RiskModel)");
    double *obj = (double *) DDSIP_Alloc (sizeof (double), 1, "obj(RiskModel)");
    double *lb = (double *) DDSIP_Alloc (sizeof (double), 1, "lb(RiskModel)");
    double *ub = (double *) DDSIP_Alloc (sizeof (double), 1, "ub(RiskModel)");
    char *sense = (char *) DDSIP_Alloc (sizeof (char), 1, "sense(RiskModel)");
    char *ctype = (char *) DDSIP_Alloc (sizeof (char), 1, "ctype(RiskModel)");
    char **colname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(RiskModel)");
    char *colstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(RiskModel)");

    // Add variable
    if (DDSIP_param->riskmod > 0)
        obj[0] = DDSIP_param->riskweight;
    else
        obj[0] = 1;

    lb[0] = 0;			//DDSIP_param->risktarget;
    ub[0] = DDSIP_infty;
    ctype[0] = 'C';
    colname[0] = & (colstore[0]);
    sprintf (colname[0], "DDSIP_v_aux_02");

    status = CPXnewcols (DDSIP_env, DDSIP_lp, 1, obj, lb, ub, ctype, colname);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk variables (BbInit) \n");
        return status;
    }

    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colstore));
    DDSIP_Free ((void **) &(ctype));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(obj));

    // Add constraint
    rhs[0] = -DDSIP_param->risktarget;	//0.0;
    sense[0] = 'G';

    status = CPXnewrows (DDSIP_env, DDSIP_lp, 1, rhs, sense, NULL, NULL);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk constraint (BbInit) \n");
        return status;
    }
    DDSIP_Free ((void **) &(rhs));
    DDSIP_Free ((void **) &(sense));

    // New coefficients: v_i >= c x+q y_i - target
    status = CPXgetobj (DDSIP_env, DDSIP_lp, vallist, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    for (i = 0; i < DDSIP_data->novar; i++)
    {
        rowlist[i] = DDSIP_data->nocon;
        collist[i] = i;
        vallist[i] = -vallist[i];
    }

    rowlist[i] = DDSIP_data->nocon;
    collist[i] = DDSIP_data->novar;
    vallist[i] = 1;

    status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_data->novar + 1, rowlist, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change coefficients for risk model (RiskModel) \n");
        return status;
    }
    // Minimize risk term, only
    if (DDSIP_param->riskmod < 0)
    {
        for (i = 0; i < DDSIP_data->novar; i++)
        {
            vallist[i] = 0.0;
            collist[i] = i;
        }

        // Set objective function coefficients to 0
        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_data->novar, collist, vallist);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective coefficients (Risk model 3)\n");
            return status;
        }

    }

    DDSIP_Free ((void **) &(vallist));
    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(rowlist));

    return status;
} // DDSIP_ExpExcess

//==========================================================================
// Risk modelling: counting scenarios above target
// In model 3 we minimize only Q_p (the probability of excess)
// In model 2 we minimize a weighted sum of Q_p and Q_e (expected value)
int
DDSIP_ExcessProb (void)
{
    int i, status = 0;
    int *rowlist;
    int *collist;
    double *vallist;
    double *rhs;
    double *obj = (double *) DDSIP_Alloc (sizeof (double), 1, "obj(RiskModel)");
    double *lb = (double *) DDSIP_Alloc (sizeof (double), 1, "lb(RiskModel)");
    double *ub = (double *) DDSIP_Alloc (sizeof (double), 1, "ub(RiskModel)");
    char *sense;
    char *ctype = (char *) DDSIP_Alloc (sizeof (char), 1, "ctype(RiskModel)");
    char **colname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(RiskModel)");
    char *colstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(RiskModel)");

    // Add #scenarios binary variables to the problem
    if (DDSIP_param->riskmod > 0)
        obj[0] = DDSIP_param->riskweight;
    else
        obj[0] = 1;
    lb[0] = 0;
    ub[0] = 1;
    ctype[0] = 'B';
    colname[0] = & (colstore[0]);
    sprintf (colname[0], "DDSIP_u02");

    status = CPXnewcols (DDSIP_env, DDSIP_lp, 1, obj, lb, ub, ctype, colname);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk variables (BbInit) \n");
        return status;
    }

    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colstore));
    DDSIP_Free ((void **) &(ctype));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(obj));

    // Add constraint to problem
    rhs = (double *) DDSIP_Alloc (sizeof (double), 1, "rhs(RiskModel)");
    sense = (char *) DDSIP_Alloc (sizeof (char), 1, "sense(RiskModel)");

    rhs[0] = DDSIP_param->risktarget;
    sense[0] = 'L';

    status = CPXnewrows (DDSIP_env, DDSIP_lp, 1, rhs, sense, NULL, NULL);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk constraint (BbInit) \n");
        return status;
    }
    DDSIP_Free ((void **) &(rhs));
    DDSIP_Free ((void **) &(sense));

    // New coefficients
    rowlist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1), "rowlist(RiskModel)");
    collist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1), "collist(RiskModel)");
    vallist = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_data->novar + 1), "vallist(RiskModel)");

    // Get objective function coefficients
    status = CPXgetobj (DDSIP_env, DDSIP_lp, vallist, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    for (i = 0; i < DDSIP_data->novar; i++)
    {
        rowlist[i] = DDSIP_data->nocon;
        collist[i] = i;
    }

    rowlist[i] = DDSIP_data->nocon;
    collist[i] = DDSIP_data->novar;
    vallist[i] = -DDSIP_param->riskM;

    status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_data->novar + 1, rowlist, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change coefficients for risk model (BbInit) \n");
        return status;
    }
    // Minimize risk term, only
    if (DDSIP_param->riskmod < 0)
    {
        for (i = 0; i < DDSIP_data->novar; i++)
        {
            vallist[i] = 0.0;
            collist[i] = i;
        }

        // Set objective function coefficients to 0
        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_data->novar, collist, vallist);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective coefficients (Risk model 3)\n");
            return status;
        }

    }


    DDSIP_Free ((void **) &(vallist));
    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(rowlist));

    return status;
} // DDSIP_ExcessProb

//==========================================================================
// Risk modelling: absolute semideviation
int
DDSIP_SemDev (void)
{
    int i, status = 0;
    int *rowlist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1),
                                        "rowlist(SemDev)");
    int *collist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1),
                                        "collist(SemDev)");
    double *vallist = (double *) DDSIP_Alloc (sizeof (double),
                      (DDSIP_data->novar + 1),
                      "vallist (SemDev)");
    double *rhs = (double *) DDSIP_Alloc (sizeof (double), 1, "rhs(SemDev)");
    double *obj = (double *) DDSIP_Alloc (sizeof (double), 1, "obj(SemDev)");
    double *lb = (double *) DDSIP_Alloc (sizeof (double), 1, "lb(SemDev)");
    double *ub = (double *) DDSIP_Alloc (sizeof (double), 1, "ub(SemDev)");
    char *sense = (char *) DDSIP_Alloc (sizeof (char), 1, "sense(SemDev)");
    char *ctype = (char *) DDSIP_Alloc (sizeof (char), 1, "ctype(SemDev)");
    char **colname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(SemDev)");
    char *colstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(SemDev)");

    // Add variable
    obj[0] = DDSIP_param->riskweight;
    lb[0] = DDSIP_param->risktarget;
    ub[0] = DDSIP_infty;
    ctype[0] = 'C';
    colname[0] = & (colstore[0]);
    sprintf (colname[0], "DDSIP_v_aux_02");

    status = CPXnewcols (DDSIP_env, DDSIP_lp, 1, obj, lb, ub, ctype, colname);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk variables (SemDev) \n");
        return status;
    }

    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colstore));
    DDSIP_Free ((void **) &(ctype));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(obj));

    // Add constraint
    rhs[0] = 0.0;
    sense[0] = 'G';

    status = CPXnewrows (DDSIP_env, DDSIP_lp, 1, rhs, sense, NULL, NULL);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk constraint (SemDev) \n");
        return status;
    }
    DDSIP_Free ((void **) &(rhs));
    DDSIP_Free ((void **) &(sense));

    // New coefficients: v_i >= c x+q y_i
    status = CPXgetobj (DDSIP_env, DDSIP_lp, vallist, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    for (i = 0; i < DDSIP_data->novar; i++)
    {
        rowlist[i] = DDSIP_data->nocon;
        collist[i] = i;
        vallist[i] = -vallist[i];
    }

    rowlist[i] = DDSIP_data->nocon;
    collist[i] = DDSIP_data->novar;
    vallist[i] = 1;

    status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_data->novar + 1, rowlist, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change coefficients for risk model (SemDev) \n");
        return status;
    }
    // Change objective (1-riskweight) * EX + riskweight * E max(X,EX)
    // Note the vallist=-vallist above
    for (i = 0; i < DDSIP_data->novar; i++)
        vallist[i] = (DDSIP_param->riskweight - 1) * vallist[i];

    status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_data->novar, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change objective coefficients (SemDev)\n");
        return status;
    }

    DDSIP_Free ((void **) &(vallist));
    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(rowlist));

    return status;
} // DDSIP_SemDev

//==========================================================================
// Risk modelling: semideviation
int
DDSIP_SemDevChgBd (void)
{
    int status;
    int *index = (int *) malloc (sizeof (int));
    double *value = (double *) malloc (sizeof (double));
    char *lu = (char *) malloc (sizeof (char));

    // Get target
    //if (DDSIP_bb->curnode)
    DDSIP_SemDevGetNodeTarget ();

    if (DDSIP_param->outlev > 3)
        fprintf (DDSIP_bb->moreoutfile, "target[%d]=%f\n", DDSIP_bb->curnode, DDSIP_node[DDSIP_bb->curnode]->target);

    // Last variable gets new bound
    index[0] = DDSIP_data->novar;
    value[0] = DDSIP_node[DDSIP_bb->curnode]->target;
    lu[0] = 'L';

    status = CPXchgbds (DDSIP_env, DDSIP_lp, 1, index, lu, value);
    if (status)
        fprintf (stderr, "ERROR: Failed to change bound (SemDevChgBd) \n");

    DDSIP_Free ((void **) &(index));
    DDSIP_Free ((void **) &(value));
    DDSIP_Free ((void **) &(lu));

    return status;
} // DDSIP_SemDevChgBd

//==========================================================================
int
DDSIP_SemDevGetNodeTarget (void)
{
    double objval, bobjval, tmpbestbound = 0.0;

    int status = 0, mipstatus, scen;

    int *index = (int *) malloc (sizeof (int));
    double *value = (double *) malloc (sizeof (double));
    char *lu = (char *) malloc (sizeof (char));

    index[0] = DDSIP_data->novar;
    value[0] = -DDSIP_infty;
    lu[0] = 'L';

    status = CPXchgbds (DDSIP_env, DDSIP_lp, 1, index, lu, value);
    if (status)
        fprintf (stderr, "ERROR: Failed to change bound (SemDevGetNodeTarget) \n");

    DDSIP_Free ((void **) &(index));
    DDSIP_Free ((void **) &(value));
    DDSIP_Free ((void **) &(lu));

    if (DDSIP_param->outlev)
    {
        fprintf (DDSIP_bb->moreoutfile, "\n----------------------\n");
        fprintf (DDSIP_bb->moreoutfile, "Target for ASD...\n");
    }

    if (DDSIP_bb->curbdcnt)
    {
        status = DDSIP_ChgBounds ();
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change problem \n");
            return status;
        }
    }
    // LowerBound problem for each scenario
    //****************************************************************************
    for (scen = 0; scen < DDSIP_param->scenarios; scen++)
    {

        if (DDSIP_killsignal)
            return 1;

        if (DDSIP_param->cpxscr)
        {
            printf ("*****************************************************************************\n");
            printf ("Solving scenario problem %d....\n", scen + 1);
        }

        status = DDSIP_ChgProb (scen);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change problem \n");
            return status;
        }
        // LP-files for debugging are written under sophisticated conditions
        if (DDSIP_param->files > 4)
        {
            char s[DDSIP_ln_fname];
            sprintf (s, "%s/asd_%d_node_%d%s", DDSIP_outdir, scen + 1, DDSIP_bb->curnode, DDSIP_param->coretype);

            status = CPXwriteprob (DDSIP_env, DDSIP_lp, s, NULL);
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to write problem\n");
                return status;
            }
        }

        if (DDSIP_param->cpxscr)
            printf ("--->\n");

//  if (DDSIP_param->cpxtime > DDSIP_param->timelim - GetCpuTime ()) {
//    status = CPXsetdblparam (DDSIP_env, CPX_PARAM_TILIM, DDSIP_Dmax (0, DDSIP_param->timelim - GetCpuTime ()));
//    if (status) {
//      fprintf (stderr, "ERROR: Failed to set cplex parameter (LowerBound) \n");
//      return status;
//    }
//    else{
//      printf ("   time limit nearly reached, set time limit for current CPLEX call to %g.\n",DDSIP_Dmax (0, DDSIP_param->timelim - GetCpuTime ()));
//    }
//  }

        status = CPXmipopt (DDSIP_env, DDSIP_lp);
        // We handle some errors separately (blatant infeasible, error in scenario problem)
        if (DDSIP_Error (status))
        {
            fprintf (stderr, "ERROR: Failed to optimize problem.(LowerBound)%d\n", status);
            return status;
        }
        // Did one of the separately handled errors occure ?
        if (DDSIP_Infeasible (status))
        {
            // If the current node is not the root node --> restore problem to original one
            // Restore types if relaxed
            status = DDSIP_RestoreBoundAndType ();
            if (status)
            {
                fprintf (stderr, "ERROR: Failed to restore bounds \n");
                return status;
            }

            return 1;
        }
        else
        {
            mipstatus = CPXgetstat (DDSIP_env, DDSIP_lp);

            // We did not detect infeasibility
            // If we couldn't find a feasible solution we can at least obtain a lower bound
            if (DDSIP_NoSolution (mipstatus))
            {
                objval = DDSIP_infty;

                status = CPXgetbestobjval (DDSIP_env, DDSIP_lp, &bobjval);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to get value of best remaining node (LowerBound)\n");
                    return status;
                }
            }
            else			// Solution exists
            {
                status = CPXgetobjval (DDSIP_env, DDSIP_lp, &objval);
                if (status)
                {
                    fprintf (stderr, "ERROR: Failed to get best objective value \n");
                    return status;
                }

                if (mipstatus == CPXMIP_OPTIMAL)
                {
                    bobjval = objval;
                }
                else
                {
                    status = CPXgetbestobjval (DDSIP_env, DDSIP_lp, &bobjval);
                    if (status)
                    {
                        fprintf (stderr, "ERROR: Failed to get value of best remaining node\n");
                        return status;
                    }
                }
            }
        }

        // If the tree was exhausted bobjval is huge
        if (fabs (bobjval) > DDSIP_infty)
            bobjval = objval;

        // Debugging information
        if (DDSIP_param->outlev)
        {
            int    cpu_hrs, cpu_mins;
            double cpu_secs;
            DDSIP_translate_time (DDSIP_GetCpuTime(),&cpu_hrs,&cpu_mins,&cpu_secs);
            fprintf (DDSIP_bb->moreoutfile,
                     "LB %3.0d:  Best=%-18.14g\tBound=%-18.14g\tStatus=%3.0d\ttime= %dh %02d:%02.0f\n",
                     scen + 1, objval, bobjval, mipstatus, cpu_hrs,cpu_mins,cpu_secs);
        }

        // Temporary bound is infinity if a scenario problem is infeasible
        if (DDSIP_Equal (DDSIP_Dmin (bobjval, objval), DDSIP_infty))
            tmpbestbound = DDSIP_infty;
        else
            tmpbestbound += DDSIP_data->prob[scen] * DDSIP_Dmin (bobjval, objval);

    }				// end for scen
    //****************************************************************************

    //  if (tmpbestbound<DDSIP_infty)
    DDSIP_node[DDSIP_bb->curnode]->target = DDSIP_Dmax (DDSIP_node[DDSIP_bb->curnode]->target, tmpbestbound);

    if (DDSIP_param->outlev > 4)
    {
        printf ("Target[%d] = %f\n",DDSIP_bb->curnode,DDSIP_node[DDSIP_bb->curnode]->target);
        fprintf (DDSIP_bb->moreoutfile, "Target[%d] = %f\n",DDSIP_bb->curnode,DDSIP_node[DDSIP_bb->curnode]->target);
    }

    // Restore problem to original one
    status = DDSIP_RestoreBoundAndType ();
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to restore bounds \n");
        return status;
    }

    return status;
} // DDSIP_SemDevGetNodeTarget

//==========================================================================
// Risk modelling: worst case costs
int
DDSIP_WorstCase (void)
{
    int i, status = 0;
    int *rowlist;
    int *collist;
    double *vallist;
    double *rhs;
    double *obj = (double *) DDSIP_Alloc (sizeof (double), 1, "obj(RiskModel)");
    double *lb = (double *) DDSIP_Alloc (sizeof (double), 1, "lb(RiskModel)");
    double *ub = (double *) DDSIP_Alloc (sizeof (double), 1, "ub(RiskModel)");
    char *sense;
    char *ctype = (char *) DDSIP_Alloc (sizeof (char), 1, "ctype(RiskModel)");
    char **colname = (char **) DDSIP_Alloc (sizeof (char *), 1, "colname(RiskModel)");
    char *colstore = (char *) DDSIP_Alloc (sizeof (char), DDSIP_ln_varname, "colstore(RiskModel)");

    // Add variable
    if (DDSIP_param->riskmod > 0)
        obj[0] = DDSIP_param->riskweight;
    else
        obj[0] = 1;

    lb[0] = -DDSIP_param->riskM;
    ub[0] = DDSIP_param->riskM;
    ctype[0] = 'C';
    colname[0] = & (colstore[0]);
    if (DDSIP_param->prefix)
    {
        if (!(strlen(DDSIP_param->prefix)))
        {
            fprintf (stderr," *** ERROR: The prefix for the first stage variables has to have a positive length.\n");
            exit (1);
        }
        sprintf (colname[0], "%sDDSIP_n_aux01",DDSIP_param->prefix);
    }
    else
    {
        if (!(strlen(DDSIP_param->postfix)))
        {
            fprintf (stderr," *** ERROR: The postfix for the first stage variables has to have a positive length.\n");
            exit (1);
        }
        sprintf (colname[0], "DDSIP_n_aux01_%s",DDSIP_param->postfix);
    }
    status = CPXnewcols (DDSIP_env, DDSIP_lp, 1, obj, lb, ub, ctype, colname);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk variables (WC) \n");
        return status;
    }

    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colstore));
    DDSIP_Free ((void **) &(ctype));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(obj));

    // Add constraint to problem
    rhs = (double *) DDSIP_Alloc (sizeof (double), 1, "rhs(RiskModel)");
    sense = (char *) DDSIP_Alloc (sizeof (char), 1, "sense(RiskModel)");

    rhs[0] = 0;
    sense[0] = 'G';

    status = CPXnewrows (DDSIP_env, DDSIP_lp, 1, rhs, sense, NULL, NULL);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk constraint (WC) \n");
        return status;
    }
    DDSIP_Free ((void **) &(rhs));
    DDSIP_Free ((void **) &(sense));

    // New coefficients
    rowlist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1), "rowlist(RiskModel)");
    collist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 1), "collist(RiskModel)");
    vallist = (double *) DDSIP_Alloc (sizeof (double), (DDSIP_data->novar + 1), "vallist(RiskModel)");

    // Get objective function coefficients
    status = CPXgetobj (DDSIP_env, DDSIP_lp, vallist, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    for (i = 0; i < DDSIP_data->novar; i++)
    {
        rowlist[i] = DDSIP_data->nocon;
        collist[i] = i;
        vallist[i] = -vallist[i];
    }

    rowlist[i] = DDSIP_data->nocon;
    collist[i] = DDSIP_data->novar;
    vallist[i] = 1;

    status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_data->novar + 1, rowlist, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change coefficients for risk model (BbInit) \n");
        return status;
    }
    // Minimize risk term, only
    if (DDSIP_param->riskmod < 0)
    {
        for (i = 0; i < DDSIP_data->novar; i++)
        {
            vallist[i] = 0.0;
            collist[i] = i;
        }

        // Set objective function coefficients to 0
        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_data->novar, collist, vallist);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective coefficients (Risk model 3)\n");
            return status;
        }
    }

    DDSIP_Free ((void **) &(vallist));
    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(rowlist));

    return status;
} // DDSIP_WorstCase

//==========================================================================
// Risk modelling: Quantil (Value-at-Risk)
int
DDSIP_TVaR (void)
{
    int i, status = 0;
    int *rowlist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 2),
                                        "rowlist(RiskModel)");
    int *collist = (int *) DDSIP_Alloc (sizeof (int), (DDSIP_data->novar + 2),
                                        "collist(RiskModel)");
    double *vallist = (double *) DDSIP_Alloc (sizeof (double),
                      (DDSIP_data->novar + 2),
                      "vallist (RiskModel)");
    double *rhs = (double *) DDSIP_Alloc (sizeof (double), 1, "rhs(RiskModel)");
    double *obj = (double *) DDSIP_Alloc (sizeof (double), 2, "obj(RiskModel)");
    double *lb = (double *) DDSIP_Alloc (sizeof (double), 2, "lb(RiskModel)");
    double *ub = (double *) DDSIP_Alloc (sizeof (double), 2, "ub(RiskModel)");
    char *sense = (char *) DDSIP_Alloc (sizeof (char), 1, "sense(RiskModel)");
    char *ctype = (char *) DDSIP_Alloc (sizeof (char), 2, "ctype(RiskModel)");
    char **colname = (char **) DDSIP_Alloc (sizeof (char *), 2, "colname(RiskModel)");
    char *colstore = (char *) DDSIP_Alloc (sizeof (char), 2 * DDSIP_ln_varname,
                                           "colstore (RiskModel)");

    // Add 2 variables to the problem
    // Lower bound: target
    if (DDSIP_param->riskmod > 0)
        obj[0] = DDSIP_param->riskweight;
    else
        obj[0] = 1;
    lb[0] = -DDSIP_param->riskM;
    ub[0] = DDSIP_param->riskM;
    ctype[0] = 'C';
    colname[0] = & (colstore[0]);
    if (DDSIP_param->prefix)
    {
        if (!(strlen(DDSIP_param->prefix)))
        {
            fprintf (stderr," *** ERROR: The prefix for the first stage variables has to have a positive length.\n");
            exit (1);
        }
        sprintf (colname[0], "%sDDSIP_n_aux01",DDSIP_param->prefix);
    }
    else
    {
        if (!(strlen(DDSIP_param->postfix)))
        {
            fprintf (stderr," *** ERROR: The postfix for the first stage variables has to have a positive length.\n");
            exit (1);
        }
        sprintf (colname[0], "DDSIPn_aux01_%s",DDSIP_param->postfix);
    }

    if (!(1 - DDSIP_param->risklevel > DDSIP_param->accuracy))
        return 129;

    if (DDSIP_param->riskmod > 0)
        obj[1] = DDSIP_param->riskweight / (1 - DDSIP_param->risklevel);
    else
        obj[1] = 1 / (1 - DDSIP_param->risklevel);

    lb[1] = 0;
    ub[1] = DDSIP_infty;
    ctype[1] = 'C';
    colname[1] = & (colstore[DDSIP_ln_varname]);
    sprintf (colname[1], "DDSIP_v_aux_02");

    status = CPXnewcols (DDSIP_env, DDSIP_lp, 2, obj, lb, ub, ctype, colname);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk variables () \n");
        return status;
    }

    DDSIP_Free ((void **) &(colname));
    DDSIP_Free ((void **) &(colstore));
    DDSIP_Free ((void **) &(ctype));
    DDSIP_Free ((void **) &(ub));
    DDSIP_Free ((void **) &(lb));
    DDSIP_Free ((void **) &(obj));

    // Add constraint to problem
    rhs[0] = 0.0;
    sense[0] = 'G';

    status = CPXnewrows (DDSIP_env, DDSIP_lp, 1, rhs, sense, NULL, NULL);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to add risk constraint (BbInit) \n");
        return status;
    }
    DDSIP_Free ((void **) &(rhs));
    DDSIP_Free ((void **) &(sense));

    // New coefficients: n + v_i - c x  - q y_i >= 0
    // Change objective sense and coefficients if problem is a maximization
    status = CPXgetobj (DDSIP_env, DDSIP_lp, vallist, 0, DDSIP_data->novar - 1);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to get objective coefficients\n");
        return status;
    }

    for (i = 0; i < DDSIP_data->novar; i++)
    {
        rowlist[i] = DDSIP_data->nocon;
        collist[i] = i;
        vallist[i] = -vallist[i];
    }

    // Coefficient for n_aux_01
    rowlist[i] = DDSIP_data->nocon;
    collist[i] = DDSIP_data->novar;
    vallist[i++] = 1;

    // Coefficient for v_aux_02
    rowlist[i] = DDSIP_data->nocon;
    collist[i] = DDSIP_data->novar + 1;
    vallist[i] = 1;

    status = CPXchgcoeflist (DDSIP_env, DDSIP_lp, DDSIP_data->novar + 2, rowlist, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change coefficients for risk model (RiskModel) \n");
        return status;
    }
    // Minimize risk term, only
    if (DDSIP_param->riskmod < 0)
    {
        for (i = 0; i < DDSIP_data->novar; i++)
        {
            vallist[i] = 0.0;
            collist[i] = i;
        }

        // Set objective function coefficients to 0
        status = CPXchgobj (DDSIP_env, DDSIP_lp, DDSIP_data->novar, collist, vallist);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective coefficients (Risk model 3)\n");
            return status;
        }
    }


    DDSIP_Free ((void **) &(vallist));
    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(rowlist));

    return status;
} // DDSIP_TVaR

//==========================================================================
// Implementation of the lower bound models for the risk measures
// Monotonuous risk measures allow for another lower bound model:
int
DDSIP_RiskModel (void)
{
    int status = 0;
    char fname[DDSIP_ln_fname];

    printf ("Initializing ");
    if (DDSIP_param->riskmod > 0)
    {
        printf ("mean-");
        fprintf (DDSIP_outfile, "mean-");
    }
    printf ("risk model.\n\t\t");
    fprintf (DDSIP_outfile, "risk model:\n\t\t");
    if (DDSIP_param->riskmod > 0)
    {
        if (DDSIP_param->scalarization)
        {
            printf (" Reference point approach for Exp. Value and");
            fprintf (DDSIP_outfile, " Reference point approach for Exp. Value and");
        }
        else
        {
            printf (" Exp. Value  +  %g", DDSIP_param->riskweight);
            fprintf (DDSIP_outfile, " Exp. Value  +  %g", DDSIP_param->riskweight);
        }
    }

    if (abs (DDSIP_param->riskmod) == 1)
    {
        printf (" Exp. excess of target %g.\n", DDSIP_param->risktarget);
        fprintf (DDSIP_outfile, " Exp. excess of target %g.\n", DDSIP_param->risktarget);
        if (!DDSIP_param->scalarization)
        {
            if (!DDSIP_param->riskalg)
                status = DDSIP_ExpExcess ();
            if (status)
                return status;
        }
    }
    else if (abs (DDSIP_param->riskmod) == 2)
    {
        printf (" Excess prob. w.r.t. target %g.\n", DDSIP_param->risktarget);
        fprintf (DDSIP_outfile, " Excess prob. w.r.t. target %g.\n", DDSIP_param->risktarget);
        if (!DDSIP_param->scalarization)
        {
            if (!DDSIP_param->riskalg)
                status = DDSIP_ExcessProb ();
            if (status)
                return status;
        }
    }
    // Absolute Semideviation
    else if (abs (DDSIP_param->riskmod) == 3)
    {
        printf (" Absolute Semideviation.\n");
        fprintf (DDSIP_outfile, " Absolute Semideviation.\n");
        if (!DDSIP_param->riskalg)
            status = DDSIP_SemDev ();
        if (status)
            return status;
    }
    // Worst case
    else if (abs (DDSIP_param->riskmod) == 4)
    {
        printf (" Worst case costs.\n");
        fprintf (DDSIP_outfile, " Worst case costs.\n");
        status = DDSIP_WorstCase ();
        if (status)
            return status;
    }
    // TVaR
    else if (abs (DDSIP_param->riskmod) == 5)
    {
        printf (" Tail value-at-risk at prob. %g\n", DDSIP_param->risklevel);
        fprintf (DDSIP_outfile, " Tail value-at-risk at prob. %g\n", DDSIP_param->risklevel);
        if (!DDSIP_param->riskalg)
            status = DDSIP_TVaR ();
        if (status)
            return status;
    }
    // Quantile
    else if (abs (DDSIP_param->riskmod) == 6)
    {
        printf (" Value-at-risk at prob. %g\n", DDSIP_param->risklevel);
        fprintf (DDSIP_outfile, " Value-at-risk at prob. %g\n", DDSIP_param->risklevel);
        status = 0;			//VaR (param);
        if (status)
            return status;
    }
    // Standard deviation
    else if (abs (DDSIP_param->riskmod) == 7)
    {
        printf (" Standard deviation.\n");
        fprintf (DDSIP_outfile, " Standard deviation.\n");
        status = 0;			//SD (param);
        if (status)
            return status;
    }

    if (DDSIP_param->riskalg == 1)
    {
        printf ("\t\t Using algorithm for fsd-consistent risk measures.\n");
        fprintf (DDSIP_outfile, "\t\t Using algorithm for fsd-consistent risk measures.\n");
    }
    else if (DDSIP_param->riskalg == 2)
    {
        printf ("\t\t Using algorithm for non-fsd-consistent risk measures.\n");
        fprintf (DDSIP_outfile, "\t\t Using algorithm for non-fsd-consistent risk measures.\n");
    }
    else if (abs (DDSIP_param->riskmod) == 3)
    {
        printf ("\t\t Using algorithm for mean-asd-model.\n");
        fprintf (DDSIP_outfile, "\t\t Using algorithm for mean-asd-model.\n");
        if (DDSIP_param->cb)
        {
            printf ("\t\t Disabling ConicBundle.\n");
            fprintf (DDSIP_outfile, "\t\t Disabling ConicBundle.\n");
            (DDSIP_param->cb = 0);
        }
    }
    else
    {
        printf ("\t\t Using scenario decomposition.\n");
        fprintf (DDSIP_outfile, "\t\t Using scenario decomposition.\n");
    }

    if (DDSIP_param->files > 1)
    {
        sprintf (fname, "%s/risk%s", DDSIP_outdir, DDSIP_param->coretype);
        status = CPXwriteprob (DDSIP_env, DDSIP_lp, fname, NULL);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to write problem\n");
            return status;
        }
    }

    return status;
} // DDSIP_RiskModel

//==========================================================================
int
DDSIP_DeleteRiskObj (void)
{
    int cnt = 0, i = 0, status = 0;
    int *collist = NULL;
    double *vallist = NULL;

    // Set coefficients for variables of risk models to 0
    cnt = 1;
    if (abs (DDSIP_param->riskmod) == 5)
        cnt = 2;

    collist = (int *) malloc (sizeof (int) * cnt);
    vallist = (double *) malloc (sizeof (double) * cnt);

    collist[0] = DDSIP_bb->novar - 1;
    vallist[0] = 0.0;

    if (abs (DDSIP_param->riskmod) == 5)
    {
        collist[1] = DDSIP_bb->novar - 2;
        vallist[1] = 0.0;
    }

    status = CPXchgobj (DDSIP_env, DDSIP_lp, cnt, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change objective function (DeleteRiskObj(1)) \n");
        return status;
    }

    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(vallist));

    // For absolute semideviation and pure risk models
    // Set original coefficients
    // ?????  The same is done again IN ANY CASE in UpperBound (sipub.c) after leaving this function -- WHY ???
    //        this omits Lagrange parameters proposed by the dual from the objective.
    //        thus:  the call to DeleteRiskObj is needed just for affecting the additional variables
    //        which means the following lines could be deleted altogether....
    //
    if (DDSIP_param->riskmod < 0 || DDSIP_param->riskmod == 3)
    {
        cnt = DDSIP_data->novar;
        collist = (int *) malloc (sizeof (int) * cnt);

        for (i = 0; i < cnt; i++)
        {
            collist[i] = i;
        }

        status = CPXchgobj (DDSIP_env, DDSIP_lp, cnt, collist, DDSIP_data->cost + DDSIP_param->scenarios * DDSIP_param->stoccost);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective function (DeleteRiskObj(2)) \n");
            return status;
        }

        DDSIP_Free ((void **) &(collist));
    }

    return status;
} // DDSIP_DeleteRiskObj

//==========================================================================
int
DDSIP_UndeleteRiskObj (void)
{
    int cnt = 0, i = 0, status = 0;
    int *collist = NULL;
    double *vallist = NULL;

    // Add coefficients for risk variables
    cnt = 1;
    //TVaR
    if (abs (DDSIP_param->riskmod) == 5)
        cnt = 2;

    collist = (int *) malloc (sizeof (int) * cnt);
    vallist = (double *) malloc (sizeof (double) * cnt);

    collist[0] = DDSIP_bb->novar - 1;
    if (DDSIP_param->riskmod > 0)
        vallist[0] = DDSIP_param->riskweight;
    else
        vallist[0] = 1;

    //TVaR
    if (abs (DDSIP_param->riskmod) == 5)
    {
        collist[0] = DDSIP_bb->novar - 2;
        collist[1] = DDSIP_bb->novar - 1;
        vallist[1] = DDSIP_param->riskweight / (1 - DDSIP_param->risklevel);
    }

    status = CPXchgobj (DDSIP_env, DDSIP_lp, cnt, collist, vallist);
    if (status)
    {
        fprintf (stderr, "ERROR: Failed to change objective function (UndeleteRiskObj(1)) \n");
        return status;
    }

    DDSIP_Free ((void **) &(collist));
    DDSIP_Free ((void **) &(vallist));

    // Set coefficient for original objective to 0 for pure risk model (min RX)
    // Change objective (1-riskweight) * EX + riskweight * E max(X,EX) for absolute semideviation
    if (DDSIP_param->riskmod < 0 || DDSIP_param->riskmod == 3)
    {
        cnt = DDSIP_data->novar;
        collist = (int *) malloc (sizeof (int) * cnt);
        vallist = (double *) malloc (sizeof (double) * cnt);

        for (i = 0; i < cnt; i++)
        {
            collist[i] = i;
            if (DDSIP_param->riskmod < 0)
                vallist[i] = 0.0;
            else
                vallist[i] = (1 - DDSIP_param->riskweight) * DDSIP_data->cost[i];
        }

        status = CPXchgobj (DDSIP_env, DDSIP_lp, cnt, collist, vallist);
        if (status)
        {
            fprintf (stderr, "ERROR: Failed to change objective function (UndeleteRiskObj(2)) \n");
            return status;
        }

        DDSIP_Free ((void **) &(collist));
        DDSIP_Free ((void **) &(vallist));
    }

    return status;
} // DDSIP_UndeleteRiskObj

//==========================================================================
// Compute risk value for a distribution according to scenbd
double
DDSIP_RiskLb (double *scensol)
{
    int i;
    double risk = 0.0L, exp = 0.0L, we, wr, d;

    for (i = 0; i < DDSIP_param->scenarios; i++)
    {
        exp += DDSIP_data->prob[i] * scensol[i];
    }
    DDSIP_bb->curexp = exp;

    if (DDSIP_param->outlev)
        fprintf (DDSIP_bb->moreoutfile, "\tQ_E(LB)= %.10g\t", exp);

    if (DDSIP_param->scalarization)
    {
        we = (DDSIP_bb->curexp - DDSIP_param->ref_point[0]) / DDSIP_param->ref_scale[0];
        //reference point scalarization
        risk = 0.0;
        if (!DDSIP_param->cb || !DDSIP_bb->violations)
        {
            // Expected excess
            if (DDSIP_param->riskmod == 1)
                for (i = 0; i < DDSIP_param->scenarios; i++)
                    risk += DDSIP_data->prob[i] * DDSIP_Dmax (scensol[i] - DDSIP_param->risktarget, 0);

            // Excess Probabilities
            else if (DDSIP_param->riskmod == 2)
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    if (scensol[i] > DDSIP_param->risktarget)
                        risk += DDSIP_data->prob[i];
                }
            else
            {
                fprintf (stderr, "ERROR: Illegal risk model for reference scalarization.\n");
                exit (1);
            }
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "Q_R(LB)= %.10g\n", risk);
        }
        else
            fprintf (DDSIP_bb->moreoutfile, "\n");
        wr = (risk - DDSIP_param->ref_point[1]) / DDSIP_param->ref_scale[1];
        //d = DDSIP_Dmax (we, wr);
        d = we;
        if (d < wr)
        {
            d = wr;
            DDSIP_bb->ref_max = 1;
        }
        else
            DDSIP_bb->ref_max = 0;
        risk = d + DDSIP_param->ref_eps * (we + wr);
// if (DDSIP_bb->violations){
//  DDSIP_bb->currisk = 0.;
//  risk = (1.+DDSIP_param->ref_eps) * we;
// }

    }
    else
    {
        //weighted sum scalarization
        if (DDSIP_param->riskmod > 0)
            risk = exp;
        if (DDSIP_param->riskalg == 1)
        {
            // Expected excess
            if (abs (DDSIP_param->riskmod) == 1)
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    risk += DDSIP_param->riskweight * DDSIP_data->prob[i] * DDSIP_Dmax (scensol[i] - DDSIP_param->risktarget, 0);
                }

            // Excess Probabilities
            else if (abs (DDSIP_param->riskmod) == 2)
                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    if (scensol[i] > DDSIP_param->risktarget)
                        risk += DDSIP_param->riskweight * DDSIP_data->prob[i];
                }
            // Semideviation
            else if (abs (DDSIP_param->riskmod) == 3)
                for (i = 0; i < DDSIP_param->scenarios; i++)
                    risk += DDSIP_param->riskweight * DDSIP_data->prob[i] * DDSIP_Dmax (scensol[i] - exp, 0);

            // Tail Value-at-Risk=5, VaR=6
            else if (abs (DDSIP_param->riskmod) == 5 || abs (DDSIP_param->riskmod) == 6)
            {
                int j;
                double sumprob;

                int *ordind = (int *) malloc (sizeof (int) * DDSIP_param->scenarios);

                for (i = 0; i < DDSIP_param->scenarios; i++)
                {
                    ordind[i] = i;
                }

                DDSIP_qsort_ins_A (scensol, ordind, 0, DDSIP_param->scenarios-1);

                sumprob = DDSIP_data->prob[ordind[0]];
                i = 1;
                while (sumprob < DDSIP_param->risklevel && i < DDSIP_param->scenarios)
                {
                    sumprob += DDSIP_data->prob[ordind[i]];
                    i++;
                }

                i--;

                if (abs (DDSIP_param->riskmod) == 6)
                    risk += DDSIP_param->riskweight * scensol[ordind[i]];
                else
                {
                    sumprob = 1 - DDSIP_param->risklevel;
                    for (j = i + 1; j < DDSIP_param->scenarios; j++)
                    {
                        risk += DDSIP_param->riskweight * DDSIP_data->prob[ordind[j]] / (1 - DDSIP_param->risklevel) * scensol[ordind[j]];
                        sumprob -= DDSIP_data->prob[ordind[j]];
                    }

                    risk += DDSIP_param->riskweight * sumprob / (1 - DDSIP_param->risklevel) * scensol[ordind[i]];
                }

                DDSIP_Free ((void **) & (ordind));
            }
            DDSIP_bb->currisk = (risk - exp) / DDSIP_param->riskweight;
            if (DDSIP_param->outlev)
                fprintf (DDSIP_bb->moreoutfile, "Q_R(LB)= %.10g\n", DDSIP_bb->currisk);
        }
        else if (DDSIP_param->outlev)
            fprintf (DDSIP_bb->moreoutfile, "\n");

    }
    return risk;
} // DDSIP_RiskLb

//==========================================================================
// Compute risk objective: DDSIP_bb->currisk
// For comparison we always compute the values for all risk measures: DDSIP_bb->curriskval[]
int
DDSIP_RiskObjective (double *scensol)
{
    int i, j, status = 0;

    double sumprob;

    int *ordind = (int *) DDSIP_Alloc(sizeof (int), DDSIP_param->scenarios, "ordind(RiskObjective)");
//double *tmpscensol = (double *) malloc (sizeof (double) * DDSIP_param->scenarios);

    // Initialize
    DDSIP_bb->currisk = 0.0;
    for (i = 0; i < DDSIP_maxrisk; i++)
        DDSIP_bb->curriskval[i] = 0.0;

    // Expected excess
    for (i = 0; i < DDSIP_param->scenarios; i++)
    {
        DDSIP_bb->curriskval[0] += DDSIP_data->prob[i] * DDSIP_Dmax (scensol[i] - DDSIP_param->risktarget, 0);
    }

    // Excess Probabilities
    for (i = 0; i < DDSIP_param->scenarios; i++)
        if (scensol[i] > DDSIP_param->risktarget)
            DDSIP_bb->curriskval[1] += DDSIP_data->prob[i];

    // Semideviation
    for (i = 0; i < DDSIP_param->scenarios; i++)
        DDSIP_bb->curriskval[2] += DDSIP_data->prob[i] * DDSIP_Dmax (scensol[i] - DDSIP_bb->curexp, 0);

    // Worst Case Costs
    DDSIP_bb->curriskval[3] = scensol[0];
    for (i = 1; i < DDSIP_param->scenarios; i++)
        DDSIP_bb->curriskval[3] = DDSIP_Dmax (DDSIP_bb->curriskval[3], scensol[i]);

    // (Tail) Value-at-Risk
    for (i = 0; i < DDSIP_param->scenarios; i++)
    {
        ordind[i] = i;
    }

    DDSIP_qsort_ins_A (scensol, ordind, 0, DDSIP_param->scenarios-1);

    sumprob = DDSIP_data->prob[ordind[0]];
    i = 1;
    while (sumprob < DDSIP_param->risklevel && i < DDSIP_param->scenarios)
    {
        sumprob += DDSIP_data->prob[ordind[i]];
        i++;
    }
    i--;

    // VaR
    DDSIP_bb->curriskval[5] = scensol[ordind[i]];

    //TVaR
    sumprob = 1 - DDSIP_param->risklevel;
    for (j = i + 1; j < DDSIP_param->scenarios; j++)
    {
        DDSIP_bb->curriskval[4] += DDSIP_data->prob[ordind[j]] / (1 - DDSIP_param->risklevel) * scensol[ordind[j]];
        sumprob -= DDSIP_data->prob[ordind[j]];
    }

    DDSIP_bb->curriskval[4] += sumprob / (1 - DDSIP_param->risklevel) * scensol[ordind[i]];

    DDSIP_Free ((void **) &(ordind));

    // Variance
    for (i = 0; i < DDSIP_param->scenarios; i++)
        DDSIP_bb->curriskval[6] += DDSIP_data->prob[i] * pow (scensol[i] - DDSIP_bb->curexp, 2);
    // Stand. dev.
    DDSIP_bb->curriskval[6] = sqrt (DDSIP_bb->curriskval[6]);

    switch (abs (DDSIP_param->riskmod))
    {
        // Expected excess
    case 1:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[0];
        break;
        // Excess Probabilities
    case 2:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[1];
        break;
        // Semideviation
    case 3:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[2];
        break;
        // Worst Case Costs
    case 4:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[3];
        break;
        // Tail Value-at-Risk
    case 5:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[4];
// Last first-stage variable is the value-at-risk
// but why should the heuristics value enter the sol. from lower bounds?? And this variable is not fixed...so the action doesn't do anything
//    if (!DDSIP_param->riskalg){
//      for (i = 0; i < DDSIP_param->scenarios; i++)
//        (DDSIP_node[DDSIP_bb->curnode]->first_sol)[i][DDSIP_bb->firstvar - 1] = DDSIP_bb->curriskval[5];
//    }
        break;
        // VaR
    case 6:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[5];
        break;
        // Standard deviation
    case 7:
        DDSIP_bb->currisk = DDSIP_bb->curriskval[6];
        break;
        // Error
    default:
        status = 1;
    }

    return status;
} // DDSIP_RiskObjective
